"""End-to-end public-source reproduction and offline conditional-result checks."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np

from io_data import HERE,verify_sources,digest
from baselines import sodium_baseline,force_baseline,macroscopicity_reproduction
from combinations import combinations,xenon_audit,radiation_response,white_scan_diagnostics
from response import layered_force_coefficient
from author_oracle import oracle
from inference import calibration_demo
from spectral_bounds import affine_certificate,finite_bound
from report import generate


def canonical(value):
    if isinstance(value,dict): return {k:canonical(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [canonical(v) for v in value]
    if isinstance(value,(float,np.floating)):
        if not np.isfinite(value): raise ValueError("Nonfinite result")
        return float(f"{value:.12g}")
    if isinstance(value,np.integer): return int(value)
    if isinstance(value,np.bool_): return bool(value)
    return value


def serial(value):
    return json.dumps(canonical(value),indent=2,sort_keys=True,allow_nan=False)+"\n"


def metadata():
    code = hashlib.sha256()
    for path in sorted(HERE.glob("*.py")):
        code.update(path.name.encode()+b"\0"+path.read_bytes())
    code.update((HERE/"requirements.txt").read_bytes())
    return {"code_sha256":code.hexdigest(),"protocol_sha256":digest((HERE/"protocol.json").read_bytes()),
            "sufficiency_sha256":digest((HERE/"sufficiency.json").read_bytes()),
            "manifest_sha256":digest((HERE/"manifest.json").read_bytes()),
            "source_hashes":{x["name"]:x["sha256"] for x in json.loads((HERE/"manifest.json").read_text())["files"]},
            "revision_association":"content hashes; no circular output/commit hash",
            "units":"SI except explicitly named keV, mK, mW, atomic mass units",
            "empirical_joint_confidence_enabled":False}


def fits_probe(path):
    raw = Path(path).read_bytes(); pos=0; summaries=[]
    while pos<len(raw):
        fields={}
        while True:
            block=raw[pos:pos+2880]; pos+=2880
            if len(block)!=2880: raise ValueError("Truncated FITS header")
            stop=False
            for j in range(0,2880,80):
                card=block[j:j+80].decode("ascii")
                if card.startswith("END "): stop=True; break
                if card[8:10]=="= ":
                    value=card[10:].strip()
                    value=value[1:].split("'")[0] if value.startswith("'") else value.split("/")[0].strip()
                    fields[card[:8].strip()]=value.strip()
            if stop: break
        summary={"extension":fields.get("EXTNAME","PRIMARY"),"rows":int(fields.get("NAXIS2","0")),
                 "columns":int(fields.get("TFIELDS","0"))}
        summaries.append(summary)
        n=int(fields.get("NAXIS","0")); size=0 if n==0 else abs(int(fields["BITPIX"]))//8
        if n:
            for j in range(1,n+1): size*=int(fields[f"NAXIS{j}"])
            size=(size+int(fields.get("PCOUNT","0")))*int(fields.get("GCOUNT","1"))
        pos+=((size+2879)//2880)*2880
    return {"evidence":"observed","sha256":digest(raw),"bytes":len(raw),"extensions":summaries,
            "scope":"downloaded access/metadata probe only; no LPF rotational science inference"}


def radiation_summary(directory,protocol):
    audit,e,eff=xenon_audit(directory)
    rates={}
    for rc in protocol["rc_grid_m"]:
        for tau in [0.,1e-12,1e-8,1e-4,1.,1e4,1e8,1e12]:
            for alpha in [1.,1.5]:
                rates[f"{rc}:{tau}:{alpha}"]=radiation_response(e,eff,rc,tau,alpha)
    return {"audit":audit,"counts_per_lambda":rates,"evidence":"observed",
            "interpretation":"derived forward coefficients, identity energy-response assumption"}


def certificate_result():
    args={"knots":[0,1,2],"kernels":[[1,1,0],[0,1,1]],"macro":[1,2,1],
          "tail_slopes":[0,0],"macro_tail_slope":0,"weights":[1,1],"bounds":[2,3]}
    return {"evidence":"formal_conditional","input":args,"validation":affine_certificate(**args),
            "finite_primal_dual":finite_bound(args["kernels"],args["bounds"],args["macro"]),
            "escape_example":finite_bound([[1,0,0],[0,1,0]],[2,3],[1,1,1]),
            "empirical_exclusion_enabled":False}


def compare(expected,actual,where="root"):
    if isinstance(expected,dict):
        if expected.keys()!=actual.keys(): raise AssertionError("Key mismatch "+where)
        for k in expected: compare(expected[k],actual[k],where+"."+k)
    elif isinstance(expected,list):
        if len(expected)!=len(actual): raise AssertionError("Length mismatch "+where)
        for j,(a,b) in enumerate(zip(expected,actual)): compare(a,b,where+f"[{j}]")
    elif isinstance(expected,(int,float)) and not isinstance(expected,bool):
        if not np.isclose(expected,actual,rtol=5e-5,atol=1e-290): raise AssertionError(f"Numeric mismatch {where}: {expected} != {actual}")
    elif expected!=actual: raise AssertionError("Mismatch "+where)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir",type=Path)
    parser.add_argument("--output-dir",type=Path,default=HERE/"results")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    protocol=json.loads((HERE/"protocol.json").read_text())
    sufficiency=json.loads((HERE/"sufficiency.json").read_text())
    if protocol["empirical_joint_exclusion_enabled"] or any(sufficiency[k]["enabled"] for k in
        ["white_joint_confidence_region", "colored_joint_physics", "arbitrary_spectrum_exclusion", "radiation_parameter_exclusion", "novel_physics_result"]):
        raise ValueError("This implementation does not support enabling empirical exclusions or novelty")
    meta=metadata()
    names=["baselines","compatibility","certificate","validation"]
    expected={n:json.loads((HERE/"results"/(n+".json")).read_text()) for n in names} if args.check or args.data_dir is None else None
    if args.data_dir:
        verify_sources(args.data_dir)
        baseline={"sodium":sodium_baseline(args.data_dir),"force":force_baseline(args.data_dir),
                  "macroscopicity":macroscopicity_reproduction(args.data_dir),
                  "radiation":radiation_summary(args.data_dir,protocol),
                  "lpf_access_probe":fits_probe(args.data_dir/"lpf_drs_sample.fits"),
                  "white_scan_diagnostic":white_scan_diagnostics(args.data_dir,protocol)}
        author=oracle(args.data_dir)
    else:
        if not args.check: parser.error("Use --data-dir to generate observed-source results")
        baseline=expected["baselines"]["data"]
        author=expected["validation"]["data"]["author_oracle"]
    compatibility=combinations(None,protocol,radiation_input=baseline["radiation"],
                               original_diagnostic=baseline["white_scan_diagnostic"])
    compatibility["layer_thickness_sensitivity"]=[
        {"rc_m":rc,"thickness_m":d,"two_sided_coefficient":layered_force_coefficient(rc,d)}
        for rc in protocol["rc_grid_m"] for d in [366e-9,370e-9,374e-9]]
    validation={"evidence":"simulated","simulation":calibration_demo(**protocol["simulation"]),
                "author_oracle":author,"scope":"simulation calibration independent of observational likelihood"}
    observed={"baselines":baseline,"compatibility":compatibility,"certificate":certificate_result(),"validation":validation}
    results={name:{"metadata":meta,"data":canonical(value)} for name,value in observed.items()}
    if args.check:
        for name in names:
            if expected[name]["metadata"]!=meta: raise AssertionError("Stale input/code association: "+name)
            compare(expected[name]["data"],results[name]["data"],name)
        hashes=json.loads((HERE/"results/artifact-hashes.json").read_text())
        for name,sha in hashes.items():
            if digest((HERE/"results"/name).read_bytes())!=sha: raise AssertionError("Artifact hash mismatch: "+name)
        with tempfile.TemporaryDirectory() as tmp:
            generate(expected,Path(tmp))
            for name in ["report.md","diagnostics.svg"]:
                if (Path(tmp)/name).read_bytes()!=(HERE/"results"/name).read_bytes():
                    raise AssertionError("Report regeneration differs: "+name)
        if args.data_dir: verify_sources(args.data_dir)
        print("Source-associated calculations agree; report/SVG and artifact hashes reproduce"+(" (full sources)" if args.data_dir else " (offline derived inputs)"))
    else:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        for name,value in results.items(): (args.output_dir/(name+".json")).write_text(serial(value))
        # Render from the serialized values: regeneration has identical rounding.
        results={n:json.loads((args.output_dir/(n+".json")).read_text()) for n in names}
        generate(results,args.output_dir)
        hashes={p.name:digest(p.read_bytes()) for p in sorted(args.output_dir.iterdir()) if p.is_file() and p.name!="artifact-hashes.json"}
        (args.output_dir/"artifact-hashes.json").write_text(serial(hashes))
        verify_sources(args.data_dir)
        print("Regenerated observed baselines, two candidate comparisons, certificates and report")


if __name__=="__main__": main()
