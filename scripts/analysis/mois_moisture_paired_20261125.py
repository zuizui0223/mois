"""Paired local white-photo-labelled species fractions, not verified local petals.
Run: python mois_moisture_paired_20261125.py --inputs /mnt/data --out /mnt/data/mois_paired_results
Requires numpy, pandas, statsmodels. Sources are prior signed Github artifacts.
"""
import argparse, hashlib, io, json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import statsmodels.api as sm

SOURCES={
 'paired':(14389295373,'3efaf6ee359aca4a49713d55f5bbadd0c86467ccd4d27c300f1e06936245759a','output_sha256','WCVP_native_introduced_matched_100km_white_genus_components.csv.gz'),
 'env':(14397855889,'4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6','numeric_outputs_sha256','all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz'),
 'geo':(14400177816,'cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd','output_sha256','global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz')}
def sha(b): return hashlib.sha256(b).hexdigest()
def read(root,key):
 ident,digest,manifest,member=SOURCES[key]
 p=root/f'github-actions-artifact-{ident}.zip'
 assert sha(p.read_bytes())==digest, key
 with ZipFile(p) as z:
  r=json.loads(z.read('receipt.json'))
  for name,h in r[manifest].items(): assert sha(z.read(name))==h,(key,name)
  d=pd.read_csv(io.BytesIO(z.read(member)),compression='gzip')
  old=json.loads(z.read('full_models.json')) if key=='paired' else None
 return d,old

def raw_design(d):
 return pd.DataFrame({'NPP':d.NPP_2001_2023_mean_kgC_m2_yr,'log_current_precip':np.log1p(d.present_precip_mm),'current_temp':d.present_temperature_C,'LGM_warming':d.threeGCM_delta_temperature_C_median,'LGM_rain':d.threeGCM_delta_log_precip_median,'NDVI':d.MODIS_NDVI_2001_2021_mean,'log_woody_AGB':np.log1p(d.ESA_AGB_woody_Mg_ha_2020),'abs_lat':abs(d.lat),'log_native_species':np.log1p(d.colour_species_n__native),'log_introduced_species':np.log1p(d.colour_species_n__introduced)})

def fit(d,km,weight):
 x=raw_design(d); assert np.isfinite(x).all().all()
 x=(x-x.mean())/x.std(ddof=1); b=d[f'block_{km}km'].astype(str)
 w=pd.Series(np.ones(len(d)) if weight=='equal' else np.minimum(d.colour_species_n__native,d.colour_species_n__introduced),index=d.index)
 assert (w>0).all()
 def demean(a):
  if isinstance(a,pd.Series): return a-(a*w).groupby(b).transform('sum')/w.groupby(b).transform('sum')
  return a-a.mul(w,axis=0).groupby(b).transform('sum').div(w.groupby(b).transform('sum'),axis=0)
 xx=demean(x); keep=b.map(b.value_counts())>1; xx=xx.loc[keep]
 assert np.linalg.matrix_rank(xx)==xx.shape[1]
 ys={'native':d.observed_white__native,'introduced':d.observed_white__introduced,'difference':d.observed_native_minus_introduced,'all_recorded':d.white_species_fraction}
 # exact contrast under common design and weights; reference prior coefficient checked separately.
 models={}
 for name,y in ys.items():
  assert np.isfinite(y).all() and (y.between(0,1).all() if name!='difference' else True)
  f=sm.WLS(demean(y).loc[keep],xx,weights=w.loc[keep]).fit(cov_type='cluster',cov_kwds={'groups':b.loc[keep]},use_t=False)
  models[name]=f
 assert np.max(np.abs(models['native'].params-models['introduced'].params-models['difference'].params))<1e-10
 rows=[]
 for name,f in models.items():
  for term in x.columns:
   ci=f.conf_int().loc[term]*100
   rows.append(dict(response=name,block_km=km,weight=weight,n_cells=int(keep.sum()),n_blocks=int(b.loc[keep].nunique()),term=term,beta_pp=float(f.params[term]*100),ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),p=float(f.pvalues[term])))
 return rows

def run(root,out):
 p,old=read(root,'paired');env,_=read(root,'env');geo,_=read(root,'geo')
 assert len(env)==12295 and len(geo)==18577
 assert p.groupby('WCVP_geo_mode').size().to_dict()=={'moderate':3602,'strict':3106}
 keys=['grid_x','grid_y']; rows=[]
 for mode,sub in p.groupby('WCVP_geo_mode'):
  g=geo[keys+[f'{mode}_native_white_species_fraction',f'{mode}_introduced_white_species_fraction',f'{mode}_native_n_photo9_species',f'{mode}_introduced_n_photo9_species']]
  d=sub.merge(env[keys+['white_species_fraction']],on=keys,validate='one_to_one').merge(g,on=keys,validate='one_to_one').reset_index(drop=True)
  assert len(d)==len(sub)
  for status in ['native','introduced']:
   assert np.allclose(d[f'observed_white__{status}'],d[f'{mode}_{status}_white_species_fraction'],atol=1e-10,rtol=0)
   assert (d[f'colour_species_n__{status}']==d[f'{mode}_{status}_n_photo9_species']).all()
  assert np.allclose(d.observed_white__native-d.observed_white__introduced,d.observed_native_minus_introduced,atol=1e-10,rtol=0)
  for km in [500,1000]:
   for weight in ['equal','minority_richness']:
    for r in fit(d,km,weight): rows.append({'mode':mode,**r})
 results=pd.DataFrame(rows)
 # Reproduce old unweighted 500 km contrast for all ten terms, not just rainfall.
 for mode in ['strict','moderate']:
  v=old[mode]['primary_500km_equal_cell']['components']['observed_native_minus_introduced']['coefficient_vector']
  r=results.query("mode==@mode and block_km==500 and weight=='equal' and response=='difference'").set_index('term')
  assert max(abs(r.loc[k,'beta_pp']-100*x) for k,x in v.items())<1e-8
 out.mkdir(parents=True,exist_ok=False)
 dest=out/'paired_white_moisture_coefficients.csv';results.to_csv(dest,index=False)
 receipt={'source_artifact_sha256':{k:v[1] for k,v in SOURCES.items()},'old_10_coefficients_reproduced_both_modes':True,'native_minus_intro_coefficient_identity':True,'observed_native_and_intro_fractions_crosschecked':True,'new_estimates_local_grid_white_label_species_fraction_not_local_flower_color_or_organ_truth':True,'no_direct_GBIF_record_effort_adjustment':True,'all_recorded_same_paired_cells':True,'minority_species_count_weight_is_sensitivity_not_inverse_variance':True,'scope':'matched cells; both groups >=30 species; strict and moderate overlap; no independent replication; current precip primary; other terms exploratory; no multiplicity correction or causal claim','output_sha256':{dest.name:sha(dest.read_bytes())}}
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(results.query("term=='log_current_precip'").to_string(index=False))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--inputs',type=Path,required=True);a.add_argument('--out',type=Path,required=True);s=a.parse_args();run(s.inputs,s.out)
