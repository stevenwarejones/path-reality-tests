"""Compare a complete local source reconstruction to the reviewed artifacts."""
import argparse
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def raw_summary(raw):
    result={}
    for side,value in raw.items():
        result[side]={k:v for k,v in value.items() if k not in ('patches','fixtures')}
        result[side]['changed_words']={k:len(v) for k,v in value['patches'].items()}
        result[side]['nominal_corrections']=value['patches']['nominal']
    return result


def verify(directory):
    directory=Path(directory)
    for name in ('counts','analysis'):
        if json.loads((directory/(name+'.json')).read_text())!=json.loads((HERE/'results'/f'{name}.json').read_text()):
            raise ValueError(f'{name}: source regeneration differs')
    raw=json.loads((directory/'reconstruction.json').read_text())
    if raw_summary(raw)!=json.loads((HERE/'results/reconstruction.json').read_text()):
        raise ValueError('raw reconstruction summary differs')
    print('Full source reconstruction, count tables and conditional analysis match')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    a=p.parse_args();verify(a.directory)
