"""Broader strict/moderate native grid sensitivity, reuses paired audit module."""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from mois_moisture_paired_20261125 import read,sha

def run(root,out):
 env,_=read(root,'env');geo,_=read(root,'geo');k=['grid_x','grid_y'];d=env.merge(geo,on=k,validate='one_to_one')
 assert len(d)==12295
 rows=[]
 for mode,count in [('strict',7201),('moderate',8119)]:
  q=d.loc[d[f'{mode}_native_n_photo9_species'].ge(30)].reset_index(drop=True)
  assert len(q)==count
  for cohort in ['all_native30','native100','native30_intro_lt30']:
   s=q.copy() if cohort=='all_native30' else q.loc[q[f'{mode}_native_n_photo9_species'].ge(100) if cohort=='native100' else q[f'{mode}_introduced_n_photo9_species'].lt(30)].reset_index(drop=True)
   x=pd.DataFrame({'log_current_precip':np.log1p(s.present_precip_mm),'current_temp':s.present_temperature_C,'LGM_warming':s.threeGCM_delta_temperature_C_median,'LGM_rain':s.threeGCM_delta_log_precip_median,'NPP':s.NPP_2001_2023_mean_kgC_m2_yr,'NDVI':s.MODIS_NDVI_2001_2021_mean,'log_woody_AGB':np.log1p(s.ESA_AGB_woody_Mg_ha_2020),'abs_lat':s.lat.abs(),'log_native_species':np.log1p(s[f'{mode}_native_n_photo9_species']),'log_all_species':np.log1p(s.n_accepted_species)})
   x=(x-x.mean())/x.std(ddof=1);assert np.isfinite(x).all().all()
   y=s[f'{mode}_native_white_species_fraction']
   for km in [500,1000]:
    b=s[f'block_{km}km'].astype(str);keep=b.map(b.value_counts()).gt(1)
    for weighting in ['equal','native_richness']:
     w=pd.Series(np.ones(len(s)) if weighting=='equal' else s[f'{mode}_native_n_photo9_species'],index=s.index)
     den=w.groupby(b).transform('sum');xx=x-x.mul(w,axis=0).groupby(b).transform('sum').div(den,axis=0);yy=y-(y*w).groupby(b).transform('sum')/den
     assert np.linalg.matrix_rank(xx.loc[keep])==xx.shape[1]
     f=sm.WLS(yy.loc[keep],xx.loc[keep],weights=w.loc[keep]).fit(cov_type='cluster',cov_kwds={'groups':b.loc[keep]},use_t=False)
     v='log_current_precip';ci=f.conf_int().loc[v]*100
     rows.append(dict(mode=mode,cohort=cohort,n_source=len(s),block_km=km,weight=weighting,n_cells=int(keep.sum()),n_blocks=b.loc[keep].nunique(),beta_pp=f.params[v]*100,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=f.pvalues[v]))
 r=pd.DataFrame(rows);out.mkdir(parents=True,exist_ok=False);p=out/'native_grid_rainfall_sensitivity.csv';r.to_csv(p,index=False)
 (out/'receipt.json').write_text(json.dumps({'source':'same verified env and geo artifact as paired analysis','not_independent_of_paired_cohort':True,'no_direct_GBIF_record_effort_adjustment':True,'native_richness_weight_not_inverse_variance':True,'no_multiple_test_correction':True,'output_sha256':{p.name:sha(p.read_bytes())}},indent=2)+'\n');print(r.to_string(index=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.inputs,a.out)
