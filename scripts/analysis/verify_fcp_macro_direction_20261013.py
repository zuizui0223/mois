#!/usr/bin/env python3
"""Independent source re-open + direct row-level weighted least squares verifier.

Does not import the estimator or use its sufficient statistics/bootstrap implementation.
Confirms point coefficients and source membership, NOT independently resampled CIs.
"""
from __future__ import annotations
import argparse,io,json,zipfile,hashlib
from pathlib import Path
import numpy as np,pandas as pd
PHOTO_SHA='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
ENV_SHA='6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca'
HUES=('white','yellow_orange','red_pink','blue_purple')
C26=['wc_elevation_m',*[f'wc_bio{i}' for i in range(1,20)],
     'wc_srad_annual_mean','wc_srad_monthly_cv','wc_wind_annual_mean','wc_wind_monthly_cv',
     'wc_vapr_annual_mean','wc_vapr_monthly_cv']
SOIL=['soil_pH','soil_SOC','soil_N','soil_clay','soil_available_water_proxy',
      'soil_cec_0_30cm_source_raw','soil_sand_0_30cm_source_raw','soil_silt_0_30cm_source_raw',
      'soil_bdod_0_30cm_source_raw','soil_cfvo_0_30cm_source_raw']
AXES=['wc_bio5','wc_bio12','wc_bio14','wc_bio15','wc_srad_annual_mean','wc_vapr_annual_mean']
GEO=['abs_lat','lon_sin','lon_cos']

def read(p,name,cols=None):
 with zipfile.ZipFile(p) as z:return pd.read_csv(io.BytesIO(z.read(name)),compression='gzip',usecols=cols,low_memory=False)

def independent_input(p,e):
 if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=PHOTO_SHA:raise ValueError('photo source differs')
 if hashlib.sha256(Path(e).read_bytes()).hexdigest()!=ENV_SHA:raise ValueError('env source differs')
 c=['inat_taxon_id','observation_id','photo_id','species','morph','measurement_status','flower_effective_pixels','nuisance_pixel_fraction']
 a=read(p,'fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz',c)
 b=read(p,'fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz',c)
 t=pd.concat([a,b],ignore_index=True).drop_duplicates('photo_id')
 q=read(e,'fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz',
   ['inat_taxon_id','observation_id','photo_id','latitude','longitude',*C26,*SOIL])
 u=t.merge(q,on=['inat_taxon_id','observation_id','photo_id'],validate='one_to_one')
 u['abs_lat']=u.latitude.abs()
 u['lon_sin']=np.sin(np.deg2rad(u.longitude));u['lon_cos']=np.cos(np.deg2rad(u.longitude))
 u['genus']=u.species.str.split().str[0]
 cls=u.measurement_status.eq('classified_four_state_morph') & u.morph.isin(HUES)
 if len(u)!=100543 or cls.sum()!=45953:raise ValueError('Photo denominator drift')
 return u,cls

def source_species(src,columns):
 f=src.copy()
 for k in HUES:f['hue_'+k]=f.morph.eq(k).astype(float)
 z=f.groupby('inat_taxon_id',sort=True).agg({**{j:'mean' for j in dict.fromkeys([*columns,'latitude','longitude',*('hue_'+k for k in HUES)])},'genus':'first','photo_id':'size'}).reset_index()
 lat=z.latitude.to_numpy(float);lon=z.longitude.to_numpy(float)
 gy=np.clip(np.floor((np.sin(np.deg2rad(lat))+1)*4.5),0,8).astype(int)
 gx=np.clip(np.floor((lon+180)/20),0,17).astype(int)
 z['cell']=gy*18+gx
 return z

def independent_lstsq(sp,photo_means,photo_sd,anchor,weight,bycell):
 sx=sp.copy();sx['cluster']=sx.genus.astype(str)+'__'+sx.cell.astype(str) if bycell else sx.genus.astype(str)
 n=sx.cluster.value_counts();sx=sx[sx.cluster.isin(n[n>=3].index)].copy()
 cols=[anchor] if bycell else [*GEO,anchor]
 xx=(sx[cols].to_numpy(float)-photo_means[cols].to_numpy(float))/photo_sd[cols].to_numpy(float)
 y=sx[['hue_'+k for k in HUES]].to_numpy(float)
 g,ind=np.unique(sx.cluster.to_numpy(str),return_inverse=True)
 count=np.bincount(ind)
 for j in range(xx.shape[1]):
  av=np.bincount(ind,weights=xx[:,j])/count
  xx[:,j]-=av[ind]
 for j in range(4):
  av=np.bincount(ind,weights=y[:,j])/count
  y[:,j]-=av[ind]
 if weight=='genus_equal':
  ct=sx.genus.value_counts();w=1/ct.loc[sx.genus].to_numpy(float)
 else:w=np.ones(len(sx),float)
 beta=np.linalg.lstsq(xx*np.sqrt(w)[:,None],y*np.sqrt(w)[:,None],rcond=None)[0][-1]
 return beta,int(len(sx)),int(sx.genus.nunique())

def run(photo,env,saved):
 doc=json.loads(Path(saved).read_text());u,cls=independent_input(photo,env)
 jobs=[('climate','all','world','genus_geo'),('climate','all','world','genus_geocell'),
       ('climate','all','E90_E180','genus_geo'),('climate','all','E90_E180','genus_geocell'),
       ('climate','strict','world','genus_geo'),('soil','all','world','genus_geo'),
       ('soil','strict','world','genus_geo')]
 checks=0;maximum=0.;readouts=[]
 for mode,quality,region,adj in jobs:
  feat=C26+(SOIL if mode=='soil' else [])
  allmask=cls & u[['latitude','longitude',*feat]].notna().all(axis=1)
  full=u.loc[allmask]
  means=full[[*GEO,*feat]].mean();sd=full[[*GEO,*feat]].std(ddof=0).where(lambda x:x>=1e-12,1.0)
  mask=allmask
  if quality=='strict':mask=mask & u.flower_effective_pixels.ge(23826) & u.nuisance_pixel_fraction.le(.1375717742195135)
  subset=u.loc[mask]
  if region=='E90_E180':subset=subset.loc[subset.longitude.ge(90)&subset.longitude.le(180)]
  sp=source_species(subset,[*GEO,*feat])
  entry=doc['subcohorts'][f'{mode}__{quality}__{region}']['estimation' if adj=='genus_geo' else 'cell_conditioned']
  for w in ['species_equal','genus_equal']:
   estimates=entry[w]
   assert estimates['status']=='OK'
   for a in AXES+( ['soil_available_water_proxy'] if mode=='soil' else []):
    out,n,ng=independent_lstsq(sp,means,sd,a,w,adj=='genus_geocell')
    assert estimates['n_species']==n and estimates['n_named_genus']==ng,(mode,quality,region,adj,w,a,'denominator',n,ng,estimates['n_species'],estimates['n_named_genus'])
    for j,h in enumerate(HUES):
     row=next(x for x in estimates['estimates'] if x['anchor']==a and x['hue']==h)
     err=abs(row['point_pp_per_original_photo_env_sd']-100*out[j]);maximum=max(err,maximum)
     if err>1e-8:raise RuntimeError(f'photo-source direct WLS mismatch: {mode} {region} {quality} {adj} {w} {a} {h}: {err}')
     checks+=1
   readouts.append({'mode':mode,'quality':quality,'region':region,'adjustment':adj,'weight':w,'source_species':n,'genera':ng})
 return {'schema':'fcp_macro_abiotic_hue_independent_source_direct_wls_v1','status':'PASS','source_SHA256':{'photo':PHOTO_SHA,'environment':ENV_SHA},
   'direct_source_genus_WLS_colour_coefficients_verified':checks,'max_abs_difference_pp':maximum,
   'study_cells':readouts,'bootstrap_ci_resampled_independently':False,
   'notes':['An independently authored input reader and direct row-level sqrt-weight least squares agrees with all saved genus-sufficient-statistics point estimates.','Confidence intervals are reproduced by pinned production source replay, not re-bootstrapped in this verifier.']}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--photo-archive',type=Path,required=True);p.add_argument('--environment-archive',type=Path,required=True)
 p.add_argument('--result-json',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 q=run(a.photo_archive,a.environment_archive,a.result_json);a.output.write_text(json.dumps(q,indent=2)+'\n');print('INDEPENDENT VERIFIED',q['direct_source_genus_WLS_colour_coefficients_verified'],q['max_abs_difference_pp'])
