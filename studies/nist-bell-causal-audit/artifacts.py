"""Keep generated numeric arrays compact and one record per line."""
import json


def write_json(path, value):
    lines=[]
    for key,item in value.items():
        prefix='  '+json.dumps(key)+': '
        if isinstance(item,list) and item:
            text='[\n'+',\n'.join('    '+json.dumps(row,separators=(',',':'),allow_nan=False) for row in item)+'\n  ]'
        else:
            text=json.dumps(item,separators=(',',':'),allow_nan=False)
        lines.append(prefix+text)
    path.write_text('{\n'+',\n'.join(lines)+'\n}\n')
