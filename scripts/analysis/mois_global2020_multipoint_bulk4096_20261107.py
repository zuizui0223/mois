"""Source-audited worldwide 2020 VIIRS/WorldPop multi-location extraction.

For each 4096 target block:
- 512 VIIRS sites per API query, four months (cloudfree band verified),
- 128 WorldPop Global1 sites per API query, matching 2020 only.
No island-size/latitude cutoffs; no photo/ROI operations or zero imputation.
Original JSON response bytes and their SHA are retained, and can be replayed
offline to produce a byte-identical values.csv without further API requests.
"""
from __future__ import annotations
import argparse,calendar,concurrent.futures,csv,datetime as dt,hashlib,io,json,math,statistics
from pathlib import Path
import requests
from scripts.analysis.mois_global_multipoint_numeric_20261107 import (
    extract_samples,VIIRS,POP,MONTHS)

MANIFEST_SHA="fa06170ce59e41b6d7e486c236acf7d004091a035f7cbddbb3632bcd788ffa65"
SERVICE={"viirs":VIIRS,"population":POP}

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_4096_GIS: "+msg)

def digest(b):return hashlib.sha256(b).hexdigest()

def daterange(y,m):
    a=dt.datetime(y,m,1,tzinfo=dt.timezone.utc)
    b=dt.datetime(y,m,calendar.monthrange(y,m)[1],23,59,59,tzinfo=dt.timezone.utc)
    return f"{int(a.timestamp()*1000)},{int(b.timestamp()*1000)}"

def validate_source_targets(rows):
    require(1<=len(rows)<=40000,"out-of-range source target cohort")
    require(len({x['target_id'] for x in rows})==len(rows),"duplicate target_id")
    require(all(str(x['target_id']).startswith("2020|") for x in rows),
            "only source-year 2020 is supported with Global1 population")
    require(all(x['land_class'] in ("island","mainland") for x in rows),
            "unresolved shoreline point is not a supported known-land target")
    require(all(-180<=float(x['query_longitude'])<=180 and
                -90<=float(x['query_latitude'])<=90 for x in rows),
            "invalid geographic coordinates")
    require(len({(float(x['query_longitude']),float(x['query_latitude'])) for x in rows})==len(rows),
            "duplicate raster request coordinates")
    return rows

def plan_requests(rows):
    validate_source_targets(rows)
    plans=[]
    for product,month,block in [("viirs",m,512) for m in MONTHS]+[("population",0,128)]:
        for start in range(0,len(rows),block):
            n=min(block,len(rows)-start)
            plans.append({"product":product,"month":month,"start":start,"n":n,
                         "cache_name":f"{product}_2020_{month:02d}_{start:05d}_{n:03d}.response"})
    return plans

def provider_params(product,month,rows,year=2020):
    require(year==2020,"2020 only: WorldPop Global1 2024 extrapolation prohibited")
    require(product in SERVICE,"unknown raster source")
    require(1<=len(rows)<= (512 if product=="viirs" else 128),
            "unsafe source multipoint count (URL 414 for 512 population)")
    require(product!="viirs" or month in MONTHS,"invalid satellite month")
    points=[[float(r['query_longitude']),float(r['query_latitude'])] for r in rows]
    geom={"points":points,"spatialReference":{"wkid":4326}}
    return {"f":"json","geometryType":"esriGeometryMultipoint",
            "geometry":json.dumps(geom,separators=(',',':')),
            "returnFirstValueOnly":"true","outFields":"*",
            "interpolation":"RSP_NearestNeighbor",
            "time":daterange(year,month) if product=="viirs" else (
                   daterange(year,1).split(',')[0]+','+daterange(year,12).split(',')[1])}

def summarize_one(months,pop):
    vals=[float(v) for status,v in months.values() if status=="OK"]
    density=float(pop[1]) if pop[0]=="OK" and pop[1] is not None else None
    require(density is None or density>=0,"population source below zero")
    return {"viirs_valid_months":len(vals),
            "viirs_median":statistics.median(vals) if len(vals)>=3 else None,
            "population_density":density,"complete":len(vals)>=3 and density is not None}

def _one(plan,rows,out,offline):
    first,last=plan['start'],plan['start']+plan['n']
    points=rows[first:last]
    p=out/'raw'/plan['cache_name']
    if offline or p.is_file():
        require(p.is_file(),"offline original source response missing: "+str(p))
        raw=p.read_bytes()
        reused=True
    else:
        reused=False
        params=provider_params(plan['product'],plan['month'],points,year=2020)
        try:
            response=requests.get(SERVICE[plan['product']],params=params,timeout=85,
                headers={"User-Agent":"MOIS source-audited global all-land numeric sampling 2026-11-07"})
        except requests.RequestException as exc:
            raise RuntimeError("source provider network error "+str(plan)) from exc
        # Store raw bytes also on non-200 responses, for audit.
        raw=response.content
        p.parent.mkdir(parents=True,exist_ok=True)
        if response.status_code!=200:
            (out/'failed').mkdir(parents=True,exist_ok=True)
            (out/'failed'/plan['cache_name']).write_bytes(raw)
            raise ValueError("provider HTTP "+str(response.status_code)+" "+str(plan))
        p.write_bytes(raw)
    try:obj=json.loads(raw)
    except ValueError as exc:raise ValueError("malformed raster response: "+str(plan)) from exc
    require(not obj.get('error'),"provider source error: "+str(obj.get('error')))
    received=obj.get('samples')
    require(isinstance(received,list),"missing provider multipoint samples")
    parsed=extract_samples(plan['product'],plan['month'],points,obj,year=2020)
    require(len(parsed)==len(points),"provider missed target ID positions")
    rec={"product":plan['product'],"month":plan['month'],
         "start":first,"n":len(points),
         "original_response_sha256":digest(raw),
         "response_bytes":len(raw),
         "returned_samples":len(received),
         "value_statuses":{z:sum(x['status']==z for x in parsed)
                          for z in sorted({x['status'] for x in parsed})},
         "source_replayed_not_refetched":reused,
         "cache_name":plan['cache_name']}
    return plan,parsed,rec

def run(targets,out,offline=False,workers=4,manifest_sha=MANIFEST_SHA,expected_rows=4096):
    require(1<=workers<=4,"maximum 4 concurrent provider queries")
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    source=Path(targets).read_bytes()
    require(digest(source)==manifest_sha,"frozen global target selection SHA differs")
    rows=list(csv.DictReader(io.StringIO(source.decode())))
    validate_source_targets(rows)
    plans=plan_requests(rows)
    require(len(rows)==expected_rows,"source cohort count changed")
    require(len(plans)==4*math.ceil(len(rows)/512)+math.ceil(len(rows)/128),
            "multipoint API plan count changed")
    # Each network request is for a disjoint product/month/position group.
    result={x['target_id']:dict(x) for x in rows}
    recs=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for plan,vals,rec in pool.map(lambda p:_one(p,rows,out,offline),plans):
            for i,v in enumerate(vals,start=plan['start']):
                r=result[rows[i]['target_id']]
                k=f"viirs_{plan['month']}" if plan['product']=="viirs" else "population_global1_0"
                r[k+'_status']=v['status']
                r[k+'_value']=v['value']
                if plan['product']=="viirs":r[k+'_cloudfree']=v['cloud']
            recs.append(rec)
    for row in rows:
        r=result[row['target_id']]
        status={m:(r[f"viirs_{m}_status"],r[f"viirs_{m}_value"]) for m in MONTHS}
        pop=(r["population_global1_0_status"],r["population_global1_0_value"])
        r.update(summarize_one(status,pop))
    final=list(result.values())
    output=io.StringIO(newline='')
    w=csv.DictWriter(output,fieldnames=list(dict.fromkeys(k for r in final for k in r)),
                     lineterminator='\n')
    w.writeheader();w.writerows(final)
    data=output.getvalue().encode()
    (out/('replay_values.csv' if offline else 'values.csv')).write_bytes(data)
    state_count={str(i):sum(r['viirs_valid_months']==i for r in final) for i in range(5)}
    receipt={"schema":"MOIS_2020_global_4096_multipoint_original_API_v1",
        "input_manifest_SHA256":digest(source),
        "targets_attempted":len(rows),
        "targets_complete":sum(r['complete'] for r in final),
        "targets_quality_incomplete":sum(not r['complete'] for r in final),
        "attempted_by_land":{"island":sum(r['land_class']=="island" for r in final),
                             "mainland":sum(r['land_class']=="mainland" for r in final)},
        "complete_by_land":{"island":sum(r['complete'] and r['land_class']=="island" for r in final),
                            "mainland":sum(r['complete'] and r['land_class']=="mainland" for r in final)},
        "viirs_valid_months_histogram":state_count,
        "viirs_valid_month_values":sum(
            r[f'viirs_{m}_status']=="OK" for r in final for m in MONTHS),
        "population_valid_values":sum(r["population_global1_0_status"]=="OK" for r in final),
        "planned_API_requests":len(plans),
        "API_requests_new":sum(not x["source_replayed_not_refetched"] for x in recs),
        "API_requests_reused_from_saved_raw":sum(x["source_replayed_not_refetched"] for x in recs),
        "responses_source_SHA256":{p['cache_name']:x['original_response_sha256']
           for p,x in zip(plans,recs)},
        "original_raw_responses_retained":True,
        "data_SHA256":digest(data),
        "offline":offline,
        "island_area_minimum":None,
        "latitude_cutoff":None,
        "population_year_requested_2020_but_source_year_not_echoed":True,
        "four_month_raw_VIIRS_not_official_annual_composite":True,
        "new_photo_downloads_or_ROI_inferences":0,
        "flower_unknowns_not_imputed":True,
        "new_ecological_coefficients_fitted":False}
    (out/('replay_receipt.json' if offline else 'receipt.json')).write_text(
        json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print('MOIS_ORIGINAL_4096_GLOBAL_MULTIPOINT_SOURCE_RESULT',
        json.dumps({k:receipt[k] for k in ('targets_attempted','targets_complete',
          'targets_quality_incomplete','complete_by_land','viirs_valid_month_values',
          'population_valid_values','API_requests_new','API_requests_reused_from_saved_raw',
          'offline')},sort_keys=True),flush=True)
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--targets',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--offline',action='store_true')
    p.add_argument('--expected-manifest-sha',default=MANIFEST_SHA)
    p.add_argument('--expected-rows',type=int,default=4096)
    a=p.parse_args()
    run(a.targets,a.out,a.offline,a.workers,a.expected_manifest_sha,a.expected_rows)
