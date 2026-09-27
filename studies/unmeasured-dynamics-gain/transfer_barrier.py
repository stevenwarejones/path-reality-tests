"""Exact limitation of the transfer-inequality relaxation, NOT a physical witness."""
import argparse,json
from fractions import Fraction as F
from dynamics_common import *
from dynamics_verify import TARGET,record
B=F(21,100)

def assigned(row):return B+(1-2*B)*F(row['k'],row['n'])

def scalar_probability(word,known,endpoint):
    if endpoint not in [0,1]:raise ValueError('endpoint must be zero or one')
    if word==TARGET:return F(endpoint)
    if word in known:return assigned(known[word])
    return F(2+endpoint,5)

def algebra_checks(b,d,q):
    # Non-target inputs: h >= 2 sqrt(b(1-b)) >= 1-b.
    if not F(1,5)<=b<=F(1,2):raise ValueError('range argument fails')
    for endpoint in [0,1]:
        c=F(2+endpoint,5)
        # Ordinary proper extensions of target; the other output may be target.
        for gap in [abs(c-endpoint),abs(c-b),abs(c-(1-b))]:
            if gap*gap>b:raise ValueError('ordinary extension bound fails')
        # Exceptional observed extension TY, against a non-target output.
        for gap in [abs(d-b),abs(d-(1-b))]:
            if gap*gap>b:raise ValueError('exceptional extension range fails')
        # TY compared with target itself: other input must be IIY.
        if (d-endpoint)**2>(q if endpoint==0 else 1-q):raise ValueError('exceptional cancellation bound fails')
    return True

def run(cache=None):
    rows=json.loads((BASE/'results/observations.json').read_text());rows=rows['rows'] if isinstance(rows,dict) else rows
    if cache:
        raw=sources(cache)
        if len(rows)!=len(raw) or any(any(a[k]!=b[k] for k in ['family','row','word_sha256','k','n','length']) for a,b in zip(rows,raw)):raise ValueError('source binding mismatch')
        extensions=[i for i,r in enumerate(raw) if r['family']=='GST' and r['word'][:len(TARGET)]==TARGET]
        if extensions!=[121] or raw[121]['word']!=TARGET+(2,):raise ValueError('exceptional extension list changed')
        if raw[119]['word']!=TARGET[:-1]:raise ValueError('cancellation prefix changed')
    if (rows[121]['word_sha256'],rows[119]['word_sha256'])!=(target_hash(TARGET+(2,)),target_hash(TARGET[:-1])):raise ValueError('exception word changed')
    if target_hash(TARGET) in {r['word_sha256'] for r in rows}:raise ValueError('target became measured')
    gst_hashes=[r['word_sha256'] for r in rows if r['family']=='GST']
    if len(set(gst_hashes))!=len(gst_hashes):raise ValueError('GST scalar assignment requires unique recorded words')
    R=region(rows);slacks=[];gst=[i for i,r in enumerate(rows) if r['family']=='GST'];p={i:assigned(rows[i]) for i in gst}
    for i in gst:slacks.extend([p[i]-F(R['lo_num'][i],GRID),F(R['hi_num'][i],GRID)-p[i]])
    for g in R['groups']:
        if g['family']!='GST':continue
        ids=[i for i in gst if g['min_length']<=rows[i]['length']<=g['max_length']]
        mean=sum(rows[i]['n']*p[i] for i in ids)/g['n']
        slacks.extend([mean-F(g['lower_numerator'],g['lower_denominator']),F(g['upper_numerator'],g['upper_denominator'])-mean])
    if min(slacks)<=0:raise ValueError('GST constraint violated')
    d=p[121];q=p[119];algebra_checks(B,d,q)
    return {'claim':'Complete-GST transfer-inequality relaxation admits the entire target interval [0,1]',
        'physical_model':False,'does_not_prove':'A physical GST-only crossing witness, a full-class GST bound, or impossibility of combination gain',
        'target':list(TARGET),'target_sha256':target_hash(TARGET),'b':record(B),
        'recorded_GST_assignment':'b+(1-2b)*k/n',
        'endpoint_assignments':[{'target':0,'other_unrecorded_words':[2,5]},{'target':1,'other_unrecorded_words':[3,5]}],
        'exceptional_recorded_proper_extensions':[{'index':121,'word':[0,0,2,2,2],'probability':record(d)}],
        'cancellation_prefix':{'index':119,'word':[0,0,2],'probability':record(q)},
        'checked_two_sided_GST_constraints':len(slacks)//2,'minimum_GST_slack':record(min(slacks)),
        'all_word_proof':'algebra_checks plus the exhaustive proper-extension fact, verified against pinned words by --cache; proof in followup.md',
        'interval_completion':'Geometric-mean concavity makes the transfer-inequality feasible set convex; interpolate the two endpoint assignments',
        'confidence':'Original region unchanged; this is a relaxation assignment, not a fitted quantum model'}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--cache');a.add_argument('--check',action='store_true');args=a.parse_args();out=run(args.cache);p=HERE/'results/transfer-barrier.json'
    if args.check:
        if json.loads(p.read_text())!=out:raise ValueError('barrier artifact mismatch')
    else:save(p,out)
    print(out['claim'],out['minimum_GST_slack']['decimal'])
