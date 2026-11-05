"""Revalidate saved responses with current parsers; never request missing batches."""
from pathlib import Path
import argparse,json,hashlib
import pandas as pd
from scripts.acquisition.mois_diel_environment import population_batch,viirs_batch

def run(root,queue,out):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    viirs=[];population=[];ncache=0
    for path in sorted(root.glob('*/cache/*.json')):
        record=json.loads(path.read_text());p=record['payload']
        signature=hashlib.sha256(json.dumps({'url':record['url'],'payload':p},sort_keys=True).encode()).hexdigest()
        if signature!=record['request_sha256']:raise ValueError('Cache request digest differs')
        if path.name.startswith('worldpop2_'):
            cells=[f['id'] for f in p['features']['features']]
            population.extend(population_batch(cells,p['years'],path.parent))
        elif path.name.startswith('viirs_'):
            # Infer the exact cached cells from the regular nine-point geometry.
            points=json.loads(p['geometry'])['points'];cells=[]
            for lon,lat in points[::9]:
                import math
                cells.append(f'{math.floor((lon+180)/.25)}_{math.floor((lat+90)/.25)}')
            year=int(json.loads(p['mosaicRule'])['where'].split('=')[1])
            # Reconstructed batch identity must match the existing cache filename.
            key=hashlib.sha256(('|'.join(cells)).encode()).hexdigest()[:16]
            if path.name!=f'viirs_{year}_{key}.json':raise ValueError('VIIRS cache geometry identity differs')
            viirs.extend(viirs_batch(cells,year,path.parent))
        else:continue
        ncache+=1
    v=pd.DataFrame(viirs);p=pd.DataFrame(population);q=pd.read_csv(queue);keys=['cell025','year']
    if not len(v) or not len(p):raise ValueError('Environmental response caches missing')
    if v.duplicated(keys).any() or p.duplicated(keys).any():raise ValueError('Duplicate environment keys')
    expected=set(map(tuple,q[keys].to_numpy()))
    if set(map(tuple,v[keys].to_numpy()))!=expected:raise ValueError('Incomplete VIIRS requests; cannot fit partial global result')
    if not expected.issubset(set(map(tuple,p[keys].to_numpy()))):raise ValueError('Incomplete population requests')
    e=q.merge(v,on=keys,validate='one_to_one').merge(p,on=keys,validate='one_to_one')
    e.to_csv(out/'sameyear_environment.csv',index=False)
    r={'requested_cell_years':len(q),'rows':len(e),'response_caches_revalidated':ncache,'valid_radiance':int(e.radiance.notna().sum()),'valid_population':int(e.population_density.notna().sum()),'sameyear_verified':True}
    (out/'environment_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--queue',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.root,a.queue,a.out)
