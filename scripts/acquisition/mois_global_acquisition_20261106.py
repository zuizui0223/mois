"""MOIS global acquisition: mainland + islands of every mapped size; no colour imputation.

Photo intake is a separate recent-creation API arm, NOT the frozen 93M census.
Environmental intake operates on all existing source designs and preserves their
identity. Downloading environmental covariates together does not pool outcomes.
"""
from __future__ import annotations
import argparse,calendar,concurrent.futures,datetime as dt,gzip,hashlib,io,json,math,os,re,time
from pathlib import Path
from urllib.parse import urlsplit,urlunsplit,urljoin
import numpy as np
import pandas as pd
import requests

CONFIG={'include_mainland':True,'minimum_island_area_km2':None,
 'retain_unresolved_geography':True,'latitude_restriction':None,
 'photo_get_attempts':1,'new_strict9_automatically_claimed':0,
 'photo_batch_limit':512,'photo_total_bytes_limit':1073741824,
 'photo_per_file_limit':33554432,'photo_download_workers':4,
 'source_window_created_after':'2026-10-28','source_window_created_before':'2026-11-06',
 'source_snapshot_arm':'live_recent_creation_20261106_separate_from_frozen93m',
 'native_species_only_filter':False,'land_area_used_for_acquisition':False}
CANONICAL_SHA='f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c'
BASELINE_SHA='3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835'
VIIRS='https://di-lawimagery1.img.arcgis.com/arcgis/rest/services/NighttimeLightsMDNB/ImageServer'
LICENSES={'cc0','cc-by','cc-by-sa','cc-by-nc','cc-by-nc-sa'}

def require(ok,msg):
 if not bool(ok):raise ValueError(msg)
def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save_json(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(obj,sort_keys=True,indent=2,default=str,allow_nan=False)+'\n')
def land_scope(area):
 if not np.isfinite(area) or area<=0:return 'unresolved'
 return 'mainland' if area>7000000 else 'island'
def global_cells(frame,year):
 require(frame.photo_id.is_unique and frame.observation_id.is_unique,'Duplicate original IDs')
 x=frame.loc[frame.observed_on.astype(str).str[:4].eq(str(year))].copy()
 valid=x.latitude.between(-90,90)&x.longitude.between(-180,180)
 n_bad=int((~valid).sum());x=x.loc[valid]
 x['gx']=np.floor((x.longitude.clip(upper=np.nextafter(180.,-np.inf))+180)/.05).clip(0,7199).astype(int)
 x['gy']=np.floor((x.latitude.clip(upper=np.nextafter(90.,-np.inf))+90)/.05).clip(0,3599).astype(int)
 x['cell_id']=x.gx.astype(str)+'_'+x.gy.astype(str)
 t=x.groupby('cell_id',as_index=False).agg(n_attempts=('photo_id','size'),gx=('gx','first'),gy=('gy','first'))
 t['longitude']=(t.gx+.5)*.05-180;t['latitude']=(t.gy+.5)*.05-90;t['year']=year
 t['work_key']=t.cell_id.map(lambda k:hashlib.sha256(('MOIS-global-env-20261106|'+k).encode()).hexdigest())
 t=t.sort_values('work_key').reset_index(drop=True)
 return t,{'source_year':year,'all_coordinate_cells':len(t),'source_photo_ids_valid_coordinates':len(x),
           'invalid_coordinates_retained':n_bad,'no_island_area_filter':True,'no_colour_success_filter':True,
           'mainland_and_unresolved_geography_retained':True}

def parse_viirs(j,y,m):
 if j.get('error'):return None,None,'api_error'
 try:v=float(j['value'])
 except (KeyError,ValueError,TypeError):return None,None,'nodata'
 fs=(j.get('catalogItems') or {}).get('features') or []
 if len(fs)!=1:return None,None,'scene_not_unique'
 a=fs[0].get('attributes') or {}
 if int(a.get('year',-1))!=y or int(a.get('month',-1))!=m:return None,None,'wrong_scene_date'
 try:p,c=map(float,str(j['properties']['Values'][0]).split()[:2])
 except (KeyError,TypeError,ValueError,IndexError):return None,None,'no_quality_count'
 if not all(np.isfinite([v,p,c])) or abs(v-p)>.001:return None,None,'invalid_value'
 if c<2:return None,c,'insufficient_cloudfree'
 return v,c,'ok'

def unseen_photos(rows,prior_ids,prior_obs,already):
 good=[];seen=set(already)
 for r in rows:
  pid=str(r['photo_id']);obs=r['observation_uuid'];url=r['original_url']
  if pid in prior_ids or obs in prior_obs or pid in seen:continue
  if r.get('license_code') not in LICENSES:continue
  u=urlsplit(url)
  if u.scheme!='https' or u.hostname not in {'inaturalist-open-data.s3.amazonaws.com','static.inaturalist.org'}:continue
  if not re.fullmatch(r'/photos/'+re.escape(pid)+r'/original\.(?:jpg|jpeg|png)',u.path,re.I):continue
  seen.add(pid);good.append(r)
 return good

def photos(args):
 from PIL import Image,ImageOps
 root=args.out;root.mkdir(parents=True,exist_ok=True)
 if (root/'receipt.json').exists():
  print('EXISTING_COMPLETED_PHOTO_BATCH_NO_REFETCH');return
 require(digest(args.baseline)==BASELINE_SHA,'1181146 frozen source SHA drift')
 old=pd.read_parquet(args.baseline,columns=['photo_id','observation_id'])
 require(len(old)==1181146 and old.photo_id.is_unique,'Prior source identity count mismatch')
 prior=set(old.photo_id.astype(str));prior_obs=set(old.observation_id.dropna().astype(str))
 ledger=root/'attempts.jsonl';done=set()
 if ledger.exists():
  done={str(json.loads(l)['photo_id']) for l in ledger.read_text().splitlines() if l.strip()}
 pages=root/'api_pages';pages.mkdir(exist_ok=True);media=root/'originals';media.mkdir(exist_ok=True)
 session=requests.Session();session.headers.update({'User-Agent':'MOIS-academic-global-source-research/20261106'})
 records=[];cursor=0;api_count=0;api_stop=None
 for page in range(10):
  f=pages/f'{page:03d}.json'
  params={'taxon_id':47125,'photos':'true','photo_license':'any','per_page':200,
          'order_by':'id','order':'asc','created_d1':CONFIG['source_window_created_after'],
          'created_d2':CONFIG['source_window_created_before'],'id_above':cursor}
  if f.exists():j=json.loads(f.read_text())
  else:
   time.sleep(1.1)
   try:
    resp=session.get('https://api.inaturalist.org/v1/observations',params=params,timeout=35);api_count+=1
    if resp.status_code!=200:
     api_stop=f'HTTP_{resp.status_code}';break
    j=resp.json();save_json(f,j)
   except (requests.RequestException,ValueError) as e:api_stop=type(e).__name__;break
  observations=j.get('results',[])
  if not observations:break
  require(all(int(z['id'])>cursor for z in observations),'iNat cursor did not advance')
  cursor=max(int(z['id']) for z in observations)
  for ob in observations:
   created=str(ob.get('created_at',''))
   if created[:10]<CONFIG['source_window_created_after']:continue
   uuid=ob.get('uuid'); tax=ob.get('taxon') or {};coords=(ob.get('geojson') or {}).get('coordinates') or [None,None]
   if not uuid:continue
   for ph in ob.get('photos',[]):
    pid=str(ph['id']);u=urlsplit(ph.get('url') or '')
    path=re.sub(r'/(?:square|small|medium|large)\.', '/original.',u.path)
    rows={'photo_id':pid,'observation_uuid':uuid,'observation_numeric_id':str(ob['id']),
      'taxon_id':str(tax.get('id','')),'taxon_name':tax.get('name'),'created_at':created,
      'observed_on':ob.get('observed_on'),'longitude':coords[0],'latitude':coords[1],
      'positional_accuracy':ob.get('positional_accuracy'),'geoprivacy':ob.get('geoprivacy'),
      'coordinates_obscured':ob.get('obscured'),'license_code':ph.get('license_code'),
      'attribution':ph.get('attribution'),'original_url':urlunsplit((u.scheme,u.netloc,path,'','')),
      'photograph_original_dimensions':ph.get('original_dimensions'),
      'source_arm':CONFIG['source_snapshot_arm'],'flower_presence_verified':False}
    records.append(rows)
  if len(observations)<200:break
 save_json(root/'source_intake.json',{'request_frame':params,'last_observation_cursor':cursor,
       'complete_source_window':len(observations)<200 if api_stop is None else False,
       'api_stop':api_stop,'no_geographic_filter':True,'candidates':len(records)})
 unique=unseen_photos(records,prior,prior_obs,done)
 unique.sort(key=lambda r:hashlib.sha256(('MOIS-new-original-20261106|'+r['photo_id']).encode()).hexdigest())
 selected=unique[:CONFIG['photo_batch_limit']]
 pd.DataFrame(unique).to_json(root/'pending_photo_id_queue.jsonl',orient='records',lines=True)
 save_json(root/'selected_source_ids.json',selected)
 (root/'licensing.txt').write_text('Photo licenses and creators are recorded per row in selected_source_ids.json. Original bytes are unmodified. No transfer of ownership. CC-BY/SA/NC terms retain their full effect.\n')
 spent=sum(p.stat().st_size for p in media.iterdir())
 def fetch(row):
  pid=row['photo_id'];rr={'photo_id':pid,'observation_uuid':row['observation_uuid'],'status':'not_acquired'}
  target=media/(pid+Path(urlsplit(row['original_url']).path).suffix.lower())
  if target.exists():return {**rr,'status':'already_cached_no_refetch','bytes':target.stat().st_size,'sha256':digest(target)}
  try:
   with requests.get(row['original_url'],stream=True,timeout=(10,40),headers={'User-Agent':'MOIS-academic-source-intake/20261106'}) as r:
    rr['http_status']=r.status_code
    if r.status_code!=200:return {**rr,'status':f'HTTP_{r.status_code}'}
    if int(r.headers.get('Content-Length','0'))>CONFIG['photo_per_file_limit']:return {**rr,'status':'over_file_budget'}
    data=bytearray()
    for part in r.iter_content(262144):
     data.extend(part)
     if len(data)>CONFIG['photo_per_file_limit']:return {**rr,'status':'over_file_budget'}
   with Image.open(io.BytesIO(data)) as im:
    im.load();rgb=im.convert('RGB');rgb.load();rr.update(width=im.width,height=im.height,source_mode=im.mode)
   target.write_bytes(data)
   return {**rr,'status':'downloaded_decoded','sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)}
  except Exception as e:return {**rr,'status':'download_or_decode_error','error':str(e)[:160]}
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for start in range(0,len(selected),4):
   if spent>=CONFIG['photo_total_bytes_limit']:break
   res=list(pool.map(fetch,selected[start:start+4]))
   with ledger.open('a') as f:
    for rr in res:f.write(json.dumps(rr,sort_keys=True)+'\n');spent+=rr.get('bytes',0)
   if any(r.get('http_status') in (403,429) for r in res):api_stop='image_provider_limit';break
 results=[json.loads(l) for l in ledger.read_text().splitlines()] if ledger.exists() else []
 successful=[r for r in results if r['status'] in ('downloaded_decoded','already_cached_no_refetch')]
 success_ids={r['photo_id'] for r in successful}
 pd.DataFrame(records).drop_duplicates('photo_id').to_csv(root/'discovered_photo_metadata.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 require(not success_ids & prior,'New photo overlaps source baseline')
 receipt={'schema':'mois_global_fresh_original_download_only_v1','scope':CONFIG,
   'prior_unique_photo_ids':len(prior),'source_created_window_is_separate_from_frozen93m':True,
   'api_requests':api_count,'observations_cursor':cursor,'metadata_photo_rows':len(records),
   'metadata_distinct_photo_ids':len({r['photo_id'] for r in records}),'eligible_unseen_photo_ids':len(unique),
   'batch_selected_ids':len(selected),'actual_attempted_photo_ids':len({r['photo_id'] for r in results}),
   'actual_downloaded_and_decoded_photo_ids':len(success_ids),'new_strict9_success':0,
   'photo_id_overlap_baseline':0,'baseline_sha256':BASELINE_SHA,'provider_stop':api_stop,
   'whole_project_crossarm_cumulative_NOT_updated':True,'not_a_global_flower_frequency_sample':True,
   'no_MNAR_sensitivity':True,'new_photo_bytes':sum(r.get('bytes',0) for r in successful),
   'attempt_status_counts':pd.Series([r['status'] for r in results],dtype=str).value_counts().to_dict()}
 save_json(root/'receipt.json',receipt);print('MOIS_NEW_GLOBAL_ORIGINALS',json.dumps(receipt,sort_keys=True),flush=True)


def env_get(row,cache):
 result={'cell_id':row['cell_id'],'year':int(row['year']),'longitude':row['longitude'],'latitude':row['latitude']}
 vals=[]
 for month in (1,4,7,10):
  path=cache/f'{row["cell_id"]}_{row["year"]}_{month:02d}.json'
  if path.exists():j=json.loads(path.read_text());state='reused'
  else:
   y=int(row['year']);a=dt.datetime(y,month,1,tzinfo=dt.timezone.utc);b=dt.datetime(y,month,calendar.monthrange(y,month)[1],23,59,59,tzinfo=dt.timezone.utc)
   params={'f':'json','geometry':json.dumps({'x':row['longitude'],'y':row['latitude'],'spatialReference':{'wkid':4326}}),
    'geometryType':'esriGeometryPoint','returnGeometry':'false','returnCatalogItems':'true',
    'time':f'{int(a.timestamp()*1000)},{int(b.timestamp()*1000)}',
    'renderingRule':json.dumps({'rasterFunction':'Average Monthly Radiance (Raw Values)'})}
   try:
    r=requests.get(VIIRS+'/identify',params=params,timeout=(10,25))
    j=r.json() if r.ok else {'error':{'http':r.status_code}}
   except (requests.RequestException,ValueError) as e:j={'error':{'type':type(e).__name__}}
   save_json(path,j);state='new_query'
  v,c,status=parse_viirs(j,int(row['year']),month)
  result[f'radiance_{month:02d}']=v;result[f'cloudfree_{month:02d}']=c
  result[f'status_{month:02d}']=status;result[f'request_{month:02d}']=state
  result[f'raw_sha256_{month:02d}']=digest(path)
  if status=='ok':vals.append(v)
 result['valid_months']=len(vals);result['radiance_median']=float(np.median(vals)) if len(vals)>=3 else None
 return result

def worldpop_2024(out,cells):
 """Find an official global R2025A 2024 raster, download once, then sample counts/km2."""
 import rasterio
 from pyproj import Geod
 base='https://data.worldpop.org/GIS/Population/Global_2015_2030/R2025A/2024/0_Mosaicked/v1/'
 dirs=[base+'1km/',base+'1km/constrained/',base]
 session=requests.Session();urls=[];checks=[]
 for u in dirs:
  try:
   r=session.get(u,timeout=20);checks.append({'url':u,'status':r.status_code})
   if r.ok:
    urls += [urljoin(u,h) for h in re.findall(r'href=[\"\']([^\"\']+\.tif)[\"\']',r.text,re.I)
      if '2024' in h and 'pop' in h.lower() and '1km' in h.lower()]
  except requests.RequestException as e:checks.append({'url':u,'error':type(e).__name__})
  if urls:break
 if not urls:
  save_json(out/'population_source_status.json',{'status':'no_verified_2024_raster_url','checks':checks,'not_replaced_by_2020':True});return None
 url=sorted(set(urls))[0];tif=out/'WorldPop2024_global_source.tif'
 if not tif.exists():
  try:
   with session.get(url,stream=True,timeout=(15,120)) as r:
    r.raise_for_status();size=int(r.headers.get('Content-Length',0));require(size<=1500000000,'population raster >1.5GB')
    count=0
    with tif.open('wb') as f:
     for chunk in r.iter_content(1048576):
      count+=len(chunk);require(count<=1500000000,'population byte budget');f.write(chunk)
  except Exception as e:
   tif.unlink(missing_ok=True);save_json(out/'population_source_status.json',{'status':'raster_fetch_failed','url':url,'error':str(e)[:150]});return None
 geod=Geod(ellps='WGS84');data=[]
 with rasterio.open(tif) as src:
  require(src.crs.to_epsg()==4326,'Need explicit population raster CRS transform')
  for rr in cells.to_dict('records'):
   lon,lat=float(rr['longitude']),float(rr['latitude']);v=next(src.sample([(lon,lat)],masked=True))[0]
   if np.ma.is_masked(v) or not np.isfinite(v) or v<0:pop=area=None
   else:
    row,col=src.index(lon,lat);x0,y0=src.transform*(col,row);x1,y1=src.transform*(col+1,row+1)
    area=abs(geod.polygon_area_perimeter([x0,x1,x1,x0],[y0,y0,y1,y1])[0])/1e6
    pop=float(v)
   data.append({'cell_id':rr['cell_id'],'population_2024_people_per_pixel':pop,
     'population_pixel_area_km2':area,'population_2024_density_per_km2':pop/area if area else None})
  meta={'status':'sampled_original_raster','year':2024,'url':url,'raster_sha256':digest(tif),
   'raster_bytes':tif.stat().st_size,'crs':str(src.crs),'transform':list(src.transform),'count_units_not_density_until_area_conversion':True,
   'population_is_modelled_not_census_observation':True,'raster_nodata':src.nodata}
 pd.DataFrame(data).to_csv(out/'population_2024_point_samples.csv',index=False);save_json(out/'population_source_status.json',meta)
 return data

def environment(args):
 out=args.out;out.mkdir(parents=True,exist_ok=True)
 if (out/'receipt.json').exists():print('EXISTING_ENV_BATCH_NO_REFETCH');return
 require(digest(args.canonical)==CANONICAL_SHA,'Canonical SHA drift')
 d=pd.read_csv(args.canonical,compression='gzip',usecols=['photo_id','observation_id','observed_on','latitude','longitude'],dtype={'photo_id':str,'observation_id':str},low_memory=False)
 cells,meta=global_cells(d,args.year)
 cells.to_csv(out/'global_all_coordinate_cell_queue.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 prev=pd.read_csv(args.previous,dtype={'cell_id':str}) if args.previous else pd.DataFrame()
 exclude=set(prev.cell_id) if 'cell_id' in prev else set()
 pending=cells.loc[~cells.cell_id.isin(exclude)].copy();selected=pending.head(args.limit).copy()
 selected.to_csv(out/'this_global_environment_batch.csv',index=False)
 cache=out/'raw_viirs';cache.mkdir(exist_ok=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
  rows=list(pool.map(lambda r:env_get(r,cache),selected.to_dict('records')))
 measured=pd.DataFrame(rows);measured.to_csv(out/'global_new_nightlight_values.csv',index=False)
 print('MOIS_GLOBAL_LIGHT_DOWNLOADED',len(measured),int(measured.valid_months.ge(3).sum()),flush=True)
 pop=worldpop_2024(out,selected) if args.year==2024 else None
 receipt={'schema':'mois_global_environment_no_area_cutoff_v1','scope':CONFIG,'source_sha256':CANONICAL_SHA,
  **meta,'prior_2024_cell_ids_reused_not_retrieved':len(set(cells.cell_id)&exclude),
  'selected_new_cell_years':len(selected),'source_responses_attempted':len(rows)*4,
  'measured_complete_nightlight_cells':int(measured.valid_months.ge(3).sum()) if len(rows) else 0,
  'source_valid_nightlight_months':int(measured.valid_months.sum()) if len(rows) else 0,
  'new_population_values':sum(r['population_2024_density_per_km2'] is not None for r in pop) if pop else 0,
  'remaining_queued_cell_years':len(pending)-len(selected),'no_new_ROI9':True,
  'no_automatic_ecological_coefficient':True,'source_designs_not_pooled_for_inference':True}
 save_json(out/'receipt.json',receipt);print('MOIS_GLOBAL_ENVIRONMENT',json.dumps(receipt,sort_keys=True),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
 a=s.add_parser('photos');a.add_argument('--baseline',type=Path,required=True);a.add_argument('--out',type=Path,required=True)
 b=s.add_parser('environment');b.add_argument('--canonical',type=Path,required=True);b.add_argument('--previous',type=Path);b.add_argument('--year',type=int,default=2024);b.add_argument('--limit',type=int,default=512);b.add_argument('--out',type=Path,required=True)
 args=p.parse_args();photos(args) if args.command=='photos' else environment(args)
