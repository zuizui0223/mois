"""Same-year, response-verified VIIRS/WorldPop acquisition with durable cache."""
from pathlib import Path
import argparse,datetime,hashlib,json,time
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import pandas as pd
import requests

VIIRS='https://di-lawimagery1.img.arcgis.com/arcgis/rest/services/NighttimeLightsMDNB/ImageServer'
POP='https://worldpop.arcgis.com/arcgis/rest/services/WorldPop_Population_Density_1km/ImageServer'
POP3='https://api.worldpop.org/v3'

def year_ms(year):return int(datetime.datetime(year,1,1,tzinfo=datetime.timezone.utc).timestamp()*1000)

def parse_sample(s,year,month,kind,lon,lat):
    loc=s.get('location',{})
    if not np.allclose([loc.get('x',np.nan),loc.get('y',np.nan)],[lon,lat],atol=1e-7,rtol=0):raise ValueError('Returned point differs')
    a=s.get('attributes') or {}
    if kind=='viirs':
        if a.get('year')!=year or a.get('month')!=month:raise ValueError('VIIRS year/month differs')
    elif a.get('StdTime')!=year_ms(year):raise ValueError('Population year differs')
    try:v=list(map(float,str(s.get('value','')).split()))
    except ValueError:return {'value':None,'cloud_free':None,'status':'NO_DATA'}
    if not v or not np.isfinite(v).all() or abs(v[0])>1e30:return {'value':None,'cloud_free':None,'status':'NO_DATA'}
    if kind=='viirs':
        if len(v)!=2:return {'value':None,'cloud_free':None,'status':'BAD_BANDS'}
        if v[1]<2:return {'value':None,'cloud_free':v[1],'status':'LOW_CLOUD_FREE'}
        return {'value':v[0],'cloud_free':v[1],'status':'OK'}
    return {'value':v[0] if v[0]>=0 else None,'cloud_free':None,'status':'OK' if v[0]>=0 else 'NO_DATA'}

def cached_request(url,payload,path,post=False,timeout=120):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    signature=hashlib.sha256(json.dumps({'url':url,'payload':payload},sort_keys=True).encode()).hexdigest()
    if path.exists():
        x=json.loads(path.read_text());assert x['request_sha256']==signature;return x['response']
    for retry in range(4):
        try:
            if post=='json':r=requests.post(url,json=payload,timeout=timeout)
            elif post:r=requests.post(url,data=payload,timeout=timeout)
            else:r=requests.get(url,params=payload,timeout=timeout)
            if r.status_code==429:raise RuntimeError('Provider rate limit; preserve checkpoint and resume later: '+r.text[:150])
            r.raise_for_status();d=r.json()
            if d.get('error'):raise ValueError(str(d['error']))
            break
        except (requests.RequestException,ValueError):
            if retry==3:raise
            time.sleep(2*(retry+1))
    record={'url':url,'payload':payload,'request_sha256':signature,'response':d,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'remaining':r.headers.get('X-RateLimit-Remaining')}
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(record,separators=(',',':')),encoding='utf-8');temp.replace(path)
    return d

def cell_points(cell):
    gx,gy=map(int,cell.split('_'));west=gx*.25-180;south=gy*.25-90
    return [(west+.25*(i+.5)/3,south+.25*(j+.5)/3) for j in range(3) for i in range(3)]

def viirs_batch(cells,year,cache):
    points=[point for cell in cells for point in cell_points(cell)]
    key=hashlib.sha256(('|'.join(cells)).encode()).hexdigest()[:16]
    p={'f':'json','geometry':json.dumps({'points':points,'spatialReference':{'wkid':4326}}),'geometryType':'esriGeometryMultipoint','returnFirstValueOnly':'false','outFields':'year,month,Name','interpolation':'RSP_NearestNeighbor','mosaicRule':json.dumps({'where':f'year={year}'}),'renderingRule':json.dumps({'rasterFunction':'None'})}
    d=cached_request(VIIRS+'/getSamples',p,Path(cache)/f'viirs_{year}_{key}.json',post=True)
    samples=[];seen=set()
    for s in d.get('samples',[]):
        loc=int(s['locationId']);assert 0<=loc<len(points)
        month=int(s['attributes']['month']);assert 1<=month<=12
        assert (loc,month) not in seen;seen.add((loc,month))
        v=parse_sample(s,year,month,'viirs',*points[loc])
        samples.append({'cell025':cells[loc//9],'point':loc%9,'year':year,'month':month,'radiance':v['value'],'cloud_free':v['cloud_free'],'status':v['status']})
    raw=pd.DataFrame(samples)
    rows=[]
    for cell in cells:
        r=raw.loc[raw.cell025.eq(cell)] if len(raw) else raw
        ok=r.loc[r.status.eq('OK')] if len(r) else r
        medians=[];n_months=[]
        for point in range(9):
            vals=ok.loc[ok.point.eq(point),'radiance'] if len(ok) else pd.Series(dtype=float)
            n_months.append(len(vals))
            if len(vals)>=8:medians.append(float(vals.median()))
        rows.append({'cell025':cell,'year':year,'radiance':float(np.mean(medians)) if len(medians)>=6 else np.nan,'n_valid_points':len(medians),'n_valid_point_months':sum(n_months),'n_returned_point_months':len(r),'radiance_spatial_sd':float(np.std(medians)) if medians else np.nan})
    return rows

def polygon(cell):
    gx,gy=map(int,cell.split('_'));x=gx*.25-180;y=gy*.25-90
    return {'type':'Feature','id':cell,'properties':{},'geometry':{'type':'Polygon','coordinates':[[[x,y],[x+.25,y],[x+.25,y+.25],[x,y+.25],[x,y]]]}}

def population_batch(cells,years,cache):
    key=hashlib.sha256(('|'.join(cells)+'|'+','.join(map(str,years))).encode()).hexdigest()[:16]
    payload={'features':{'type':'FeatureCollection','features':[polygon(c) for c in cells]},'years':years,'include':'population'}
    d=cached_request(POP3+'/batch',payload,Path(cache)/f'worldpop2_{key}.json',post='json',timeout=240)
    rows=[];seen=set()
    for feature in d.get('features',[]):
        cell=feature['id'];assert cell in cells and cell not in seen;seen.add(cell)
        for value in feature.get('years',[]):
            year=value['data_year'];assert year in years
            value_ok=not feature.get('skipped',False) and np.isfinite(value['population_density']) and value['population_density']>=0
            rows.append({'cell025':cell,'year':year,'population_density':value['population_density'] if value_ok else np.nan,'population_area_km2':feature['area_km2']})
    assert seen==set(cells) and len(rows)==len(cells)*len(years)
    return rows

def run(queue,out,cache,kind,workers=3,shard=0,shards=1):
    q=pd.read_csv(queue,dtype={'cell025':str});assert not q.duplicated(['cell025','year']).any()
    out=Path(out);out.mkdir(parents=True,exist_ok=True);jobs=[]
    if kind=='viirs':
        for year,frame in q.groupby('year',sort=True):
            cells=sorted(frame.cell025.unique())
            for start in range(0,len(cells),20):jobs.append((cells[start:start+20],int(year)))
    else:
        # <=80 full 0.25deg cells remain under the 100,000 km2 anonymous area limit.
        cells=sorted(q.cell025.unique());years=sorted(map(int,q.year.unique()))
        for start in range(0,len(cells),80):jobs.append((cells[start:start+80],years))
        if len(jobs)>900:raise ValueError('Population requests exceed conservative daily budget; use bulk rasters')
    jobs=[j for i,j in enumerate(jobs) if i%shards==shard]
    allrows=[];done=0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures=[ex.submit(viirs_batch if kind=='viirs' else population_batch,*job,cache) for job in jobs]
        for f in as_completed(futures):
            allrows.extend(f.result());done+=1
            if done%10==0 or done==len(jobs):
                pd.DataFrame(allrows).to_csv(out/f'{kind}_shard{shard:03d}.csv',index=False)
                print(json.dumps({'kind':kind,'batches_done':done,'batches_total':len(jobs),'rows':len(allrows),'shard':shard}),flush=True)
    return allrows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--queue',required=True);p.add_argument('--out',required=True);p.add_argument('--cache',required=True);p.add_argument('--kind',choices=['viirs','population'],required=True);p.add_argument('--workers',type=int,default=3);p.add_argument('--shard',type=int,default=0);p.add_argument('--shards',type=int,default=1);a=p.parse_args();run(a.queue,a.out,a.cache,a.kind,a.workers,a.shard,a.shards)
