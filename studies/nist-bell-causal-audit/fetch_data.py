"""Fetch pinned NIST inputs outside the repository; verify before use."""
import argparse
import json
from pathlib import Path
import urllib.request
from reconstruct import sha256

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def fetch(directory, keys, verify_only=False):
    directory=Path(directory).resolve()
    if directory.is_relative_to(ROOT):
        raise ValueError('large external inputs must be stored outside the repository')
    manifest=json.loads((HERE/'manifest.json').read_text())
    directory.mkdir(parents=True,exist_ok=True)
    for key in keys:
        spec=manifest['files'][key]
        target=directory/(key+('.hdf5' if key=='hdf5' else '.zip'))
        if target.exists():
            if target.stat().st_size != spec['bytes'] or sha256(target)!=spec['sha256']:
                raise ValueError(f'{target}: checksum mismatch; original retained')
        elif verify_only:
            raise FileNotFoundError(target)
        else:
            temporary=target.with_suffix(target.suffix+'.partial')
            try:
                with urllib.request.urlopen(spec['url'],timeout=90) as response, temporary.open('wb') as out:
                    size=0
                    while block:=response.read(1024*1024):
                        size+=len(block)
                        if size>spec['bytes']:
                            raise ValueError('response exceeds pinned input size')
                        out.write(block)
                if size!=spec['bytes'] or sha256(temporary)!=spec['sha256']:
                    raise ValueError('download does not match pinned input')
                temporary.replace(target)
            finally:
                temporary.unlink(missing_ok=True)
        print(f'{key}: verified {spec["bytes"]} bytes')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--keys',nargs='+',choices=['hdf5','alice','bob','code'],default=['hdf5','alice','bob'])
    p.add_argument('--verify-only',action='store_true')
    a=p.parse_args()
    fetch(a.directory,a.keys,a.verify_only)
