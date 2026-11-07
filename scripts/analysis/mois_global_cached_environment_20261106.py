"""Worldwide mainland+island acquisition, without area or latitude thresholds.
Numerical environment only: no photographs or ROI inference. Four monthly VIIRS
values are not an official annual product. WorldPop request year is 2020; its
response does not echo a year. Raw responses and missing outcomes are preserved.
"""
from __future__ import annotations
import argparse,calendar,concurrent.futures,csv,datetime as dt,hashlib,json,math,statistics,threading,time
from pathlib import Path
import requests

ENDPOINTS={
 'viirs':'https://di-lawimagery1.img.arcgis.com/arcgis/rest/services/NighttimeLightsMDNB/ImageServer/identify',
 'population_global1':'https://worldpop.arcgis.com/arcgis/rest/services/WorldPop_Population_Density_1km/ImageServer/identify'}
MONTHS=(1,4,7,10)
STOP={'AUTH_OR_RATE_LIMIT','NETWORK_UNAVAILABLE'}
REUSABLE={'OK','NO_DATA','LOW_CLOUD_COUNT'}

def digest(b):return hashlib.sha256(b).hexdigest()
def numeric(x):
    if x is None or isinstance(x,bool):return None
    try:v=float(x)
    except (ValueError,TypeError):return None
    return v if math.isfinite(v) else None

def spec(product,year,month,lon,lat):
    if product not in ENDPOINTS:raise ValueError('Unknown product')
    if product=='population_global1' and not 2000<=int(year)<=2020:raise ValueError('Population year unsupported; no year substitution')
    if product=='viirs' and not 1<=int(month)<=12:raise ValueError('Invalid month')
    if not(-180<=float(lon)<=180 and -90<=float(lat)<=90):raise ValueError('Invalid coordinates')
    return {'product':product,'endpoint':ENDPOINTS[product],'year':int(year),
            'month':int(month) if product=='viirs' else 0,'longitude':float(lon),'latitude':float(lat),
            'rendering':'Average Monthly Radiance (Raw Values)' if product=='viirs' else 'native_population_density'}

def key(s):return digest(json.dumps(s,sort_keys=True,separators=(',',':')).encode())
def times(year,month):
    a=dt.datetime(year,month,1,tzinfo=dt.timezone.utc)
    b=dt.datetime(year,month,calendar.monthrange(year,month)[1],23,59,59,tzinfo=dt.timezone.utc)
    return f'{int(a.timestamp()*1000)},{int(b.timestamp()*1000)}'

def parameters(s):
    q={'f':'json','geometry':json.dumps({'x':s['longitude'],'y':s['latitude'],'spatialReference':{'wkid':4326}}),
       'geometryType':'esriGeometryPoint','returnGeometry':'false','returnCatalogItems':'true'}
    if s['product']=='viirs':
        q['time']=times(s['year'],s['month']);q['renderingRule']=json.dumps({'rasterFunction':s['rendering']})
    else:q['time']=times(s['year'],1).split(',')[0]+','+times(s['year'],12).split(',')[1]
    return q

def parse(obj,product,year,month):
    if obj.get('error'):
        err=obj['error'];code=err.get('code') if isinstance(err,dict) else None
        return None,'AUTH_OR_RATE_LIMIT' if code in (401,403,429,498,499) else 'API_ERROR'
    v=numeric(obj.get('value'))
    if v is None:return None,'NO_DATA'
    if product=='population_global1':return (v,'OK') if v>=0 else (None,'NO_DATA')
    features=(obj.get('catalogItems') or {}).get('features') or []
    if len(features)!=1:return None,'SOURCE_MONTH_NOT_SINGLE'
    at=features[0].get('attributes') or {}
    if numeric(at.get('year'))!=year or numeric(at.get('month'))!=month:return None,'WRONG_YEAR_MONTH'
    vals=(obj.get('properties') or {}).get('Values') or []
    if not vals or len(str(vals[0]).split())<2:return None,'NO_CLOUD_COUNT'
    rad,n=map(numeric,str(vals[0]).split()[:2])
    if rad is None or n is None:return None,'INVALID_CLOUD_COUNT'
    if abs(rad-v)>0.001:return None,'VALUE_MISMATCH'
    if n<2:return None,'LOW_CLOUD_COUNT'
    return v,'OK'

def save_response(folder,s,raw,http,origin='new_request'):
    folder.mkdir(parents=True,exist_ok=True);k=key(s)
    (folder/(k+'.response')).write_bytes(raw)
    if http in (401,403,429):val,status=None,'AUTH_OR_RATE_LIMIT'
    elif http!=200:val,status=None,'HTTP_ERROR'
    else:
        try:val,status=parse(json.loads(raw),s['product'],s['year'],s['month'])
        except (ValueError,TypeError,AttributeError):val,status=None,'INVALID_JSON'
    r={'request_key':k,'request':s,'params':parameters(s),'http_status':http,
       'response_sha256':digest(raw),'value':val,'status':status,'origin':origin,
       'retrieved_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
       'population_response_year_echo':False if s['product']=='population_global1' else None}
    (folder/(k+'.json')).write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    return r

def read_cache(folder,s):
    p=folder/(key(s)+'.json')
    if not p.exists():return None
    r=json.loads(p.read_text())
    if r.get('request')!=s:raise ValueError('Request identity changed')
    raw=folder/(key(s)+'.response')
    if not raw.is_file() or digest(raw.read_bytes())!=r.get('response_sha256'):raise ValueError('Cached raw source SHA mismatch')
    if r['status'] not in REUSABLE:return None
    return r

def summary(records):
    vals=[r['value'] for r in records if r['product']=='viirs' and r['status']=='OK']
    pop=[r['value'] for r in records if r['product']=='population_global1' and r['status']=='OK']
    return {'viirs_valid_months':len(vals),'viirs_median':statistics.median(vals) if len(vals)>=3 else None,
            'population_density':pop[0] if pop else None,'complete':len(vals)>=3 and bool(pop)}

def import_preflight(folder,seed):
    if seed is None:return 0
    root=Path(seed);r=json.loads((root/'preflight_receipt.json').read_text());count=0
    for row in r['requests']:
        s=spec(row['product'],row['year'],row['month'],r['longitude'],r['latitude'])
        raw=(root/row['response_file']).read_bytes()
        if digest(raw)!=row['response_sha256']:raise ValueError('Preflight bytes changed')
        if not (folder/(key(s)+'.json')).exists():
            save_response(folder,s,raw,row['http_status'],origin='preflight_run_42190860854');count+=1
    return count

def run(targets,out,seed=None,workers=4,offline=False):
    if not 1<=workers<=4:raise ValueError('Max four concurrent provider requests')
    out=Path(out);out.mkdir(parents=True,exist_ok=True);cache=out/'cache';cache.mkdir(exist_ok=True)
    rows=list(csv.DictReader(Path(targets).open(newline='')))
    if not rows or len(rows)>256 or len({r['target_id'] for r in rows})!=len(rows):raise ValueError('Invalid bounded target manifest')
    keys=[]
    for row in rows:
        y=int(row['target_id'].split('|')[0]);lon=float(row['query_longitude']);lat=float(row['query_latitude'])
        for p,m in [('viirs',m) for m in MONTHS]+[('population_global1',0)]:keys.append(key(spec(p,y,m,lon,lat)))
    if len(keys)!=len(set(keys)):raise ValueError('Duplicate source requests in batch')
    imported=import_preflight(cache,seed);stop=threading.Event();lock=threading.Lock()
    ledger=out/'request_ledger.jsonl';request_rows=[];results=[]
    def collect(row):
        if stop.is_set():return None
        y=int(row['target_id'].split('|')[0]);observed=[]
        for product,month in [('viirs',m) for m in MONTHS]+[('population_global1',0)]:
            if stop.is_set():break
            s=spec(product,y,month,row['query_longitude'],row['query_latitude'])
            record=read_cache(cache,s);reused=record is not None
            if record is None:
                if offline:raise ValueError('Offline replay lacks source response: '+key(s))
                try:
                    response=requests.get(s['endpoint'],params=parameters(s),timeout=25)
                    record=save_response(cache,s,response.content,response.status_code)
                except requests.ConnectionError as exc:
                    record={'request_key':key(s),'request':s,'status':'NETWORK_UNAVAILABLE','value':None,'error':str(exc)[:300]}
                except requests.RequestException as exc:
                    record={'request_key':key(s),'request':s,'status':'REQUEST_FAILED','value':None,'error':str(exc)[:300]}
            rec={**record,'target_id':row['target_id'],'product':product,'month':month,'reused':reused}
            if rec['status'] in STOP:stop.set()
            observed.append(rec)
            with lock:
                request_rows.append(rec)
                with ledger.open('a') as f:f.write(json.dumps(rec,sort_keys=True)+'\n')
            if not reused:time.sleep(.08)
        if not observed:return None
        result={**row,**summary(observed)}
        for r in observed:result[f"{r['product']}_{r['month']}_status"]=r['status'];result[f"{r['product']}_{r['month']}_value"]=r['value']
        return result
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(collect,rows):
            if result is not None:results.append(result)
    fields=list(dict.fromkeys(k for r in results for k in r))
    values=out/('replay_values.csv' if offline else 'values.csv')
    with values.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(results)
    statuses={k:sum(r['status']==k for r in request_rows) for k in sorted({r['status'] for r in request_rows})}
    receipt={'schema':'mois_global_no_area_real_collection_v2','target_sha256':digest(Path(targets).read_bytes()),
      'targets_planned':len(rows),'targets_attempted':len(results),'targets_complete':sum(r['complete'] for r in results),
      'successful_numerical_values_total':sum(r['status']=='OK' for r in request_rows),
      'new_successful_values':sum(r['status']=='OK' and not r['reused'] for r in request_rows),
      'requests_new':sum(not r['reused'] for r in request_rows),'requests_reused':sum(r['reused'] for r in request_rows),
      'seed_records_imported':imported,'statuses':statuses,'stopped_for_provider_or_network':stop.is_set(),
      'target_classes':{k:sum(r['land_class']==k for r in results) for k in ('mainland','island')},
      'complete_classes':{k:sum(r['land_class']==k and r['complete'] for r in results) for k in ('mainland','island')},
      'new_photo_requests':0,'ROI_reinference':0,'area_filter':None,'latitude_filter':None,
      'population_year_requested':2020,'population_year_echoed_in_response':False,
      'four_month_summary_is_not_official_annual_product':True,'offline_replay':offline,
      'values_sha256':digest(values.read_bytes())}
    (out/('replay_receipt.json' if offline else 'receipt.json')).write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print('MOIS_GLOBAL_ACTUAL_ENVIRONMENT_VALUES',json.dumps(receipt,sort_keys=True),flush=True)
    if stop.is_set():raise RuntimeError('Stopped after provider/network failure; checkpoints retained')
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--targets',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path);p.add_argument('--seed',type=Path)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--offline',action='store_true')
    a=p.parse_args();run(a.targets,a.out,a.seed,a.workers,a.offline)
