#!/usr/bin/env python3
"""Independent reanalysis of signed MOIS 2026-11-17 Hmsc held-out predictions.

Inputs are connector-downloaded, SHA-verified source ZIPs. No model refitting:
this is an exploratory stratified prediction audit, not a fourth-corner test.
Requires pandas and numpy. Set MOIS_DATA_DIR to the directory with the four ZIPs.
"""
import hashlib, io, json, zipfile, os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(os.environ.get('MOIS_DATA_DIR', '.')).resolve()
PREFLIGHT=(ROOT/'mois_Hmsc_preflight_20261117.zip','e9136c8c79a55986ebb6034273027554c02aaf391c6cf60e182942a792e8edb6')
SOURCES=[
 ('Europe_temperate',ROOT/'mois_Hmsc_14293759536.zip','feac5c7e6f5be700e98fe4251e24c9f73721b1051560ebce2addeb2d71da495c'),
 ('North_America_temperate',ROOT/'mois_Hmsc_14294502608.zip','ff3fb04ea6c410d60dbcc0f8d4dfd80aa473e68cc4bc1f96d2a436f434e4a1bc'),
 ('East_Asia_temperate',ROOT/'mois_Hmsc_14294487737.zip','596e8573bd2a60f034519d4b04f869dc17c4c2ca6c2d1ac9e82a5db676a371ea'),
]
REPS=20000
rng=np.random.default_rng(20261117)

def verify(path,expected):
    found=hashlib.sha256(path.read_bytes()).hexdigest()
    if found!=expected: raise ValueError(f'SHA mismatch {path}: {found}')

def csv(z,name):
    return pd.read_csv(io.BytesIO(z.read(name)))

def block_ci(df,group_col,level,delta_col):
    sel=df[df[group_col]==level]
    grouped=sel.groupby('block_500km')[delta_col].agg(['sum','count'])
    if len(grouped)<3: raise ValueError('too few spatial blocks')
    ix=rng.integers(0,len(grouped),size=(REPS,len(grouped)))
    sums=grouped['sum'].to_numpy()[ix].sum(axis=1)
    ns=grouped['count'].to_numpy()[ix].sum(axis=1)
    return [float(np.percentile(sums/ns,2.5)),float(np.percentile(sums/ns,97.5))],len(grouped)

def contrast_ci(df,delta_col):
    g=df.groupby(['block_500km','white_proxy'])[delta_col].agg(['sum','count'])
    blocks=sorted(df.block_500km.unique())
    ss=np.zeros((len(blocks),2));nn=np.zeros((len(blocks),2))
    for i,b in enumerate(blocks):
        for j in (0,1):
            if (b,j) in g.index: ss[i,j]=g.loc[(b,j),'sum'];nn[i,j]=g.loc[(b,j),'count']
    ix=rng.integers(0,len(blocks),size=(REPS,len(blocks)))
    sums=ss[ix].sum(axis=1);ns=nn[ix].sum(axis=1)
    ok=(ns>0).all(axis=1)
    dif=sums[ok,1]/ns[ok,1]-sums[ok,0]/ns[ok,0]
    return [float(np.percentile(dif,2.5)),float(np.percentile(dif,97.5))],len(blocks)

def species_equal_sensitivity(df,delta_col):
    species=sorted(df.species_key.unique())
    blocks=sorted(df.block_500km.unique())
    si={k:i for i,k in enumerate(species)}
    bi={k:i for i,k in enumerate(blocks)}
    ones=np.zeros((len(blocks),len(species)))
    sums=np.zeros_like(ones)
    for (block,species_key),g in df.groupby(['block_500km','species_key']):
        i=bi[block];j=si[species_key]
        sums[i,j]=g[delta_col].sum()
        ones[i,j]=len(g)
    color=dict(zip(df.species_key,df.white_proxy))
    mask=np.array([color[k] for k in species],dtype=int)
    def est_for_index(ix):
        ss=sums[ix].sum(axis=0);nn=ones[ix].sum(axis=0)
        means=np.divide(ss,nn,where=nn>0,out=np.full_like(ss,np.nan))
        return [np.nanmean(means[mask==c]) for c in (0,1)]
    point=est_for_index(np.arange(len(blocks)))
    draws=rng.integers(0,len(blocks),size=(REPS,len(blocks)))
    diffs=np.zeros((REPS,3))
    for k,idx in enumerate(draws):
        a,b=est_for_index(idx)
        diffs[k,:]=a,b,b-a
    loo=[]
    for i in range(len(blocks)):
        a,b=est_for_index(np.array([j for j in range(len(blocks)) if i!=j]))
        loo.append([a,b,b-a])
    loo=np.asarray(loo)
    rows=[]
    for ind,grp in enumerate(('nonwhite','white','white_minus_nonwhite')):
        value=point[ind] if ind<2 else point[1]-point[0]
        rows.append(dict(group=grp,comparison=delta_col,n=len(species) if ind==2 else int(sum(mask==ind)),n_500km_blocks=len(blocks),mean_diff_logloss=float(value),CI95_lower=float(np.percentile(diffs[:,ind],2.5)),CI95_upper=float(np.percentile(diffs[:,ind],97.5)),loo_min=float(loo[:,ind].min()),loo_max=float(loo[:,ind].max()),p_value='not_computed',interpretation='negative_is_better' if ind<2 else 'negative_means_larger_gain_for_white',weighting='equal_species'))
    return rows

def compute():
    verify(*PREFLIGHT)
    for _,path,expected in SOURCES: verify(path,expected)
    allrows=[];diagnostics=[]
    with zipfile.ZipFile(PREFLIGHT[0]) as pref:
      for region,path,_ in SOURCES:
        with zipfile.ZipFile(path) as z:
            tr=csv(pref,f'{region}/TR_photo9_species_traits.csv')
            env=csv(pref,f'{region}/X_true_GBIF_100km_environment.csv')
            Ys=csv(pref,f'{region}/Y_apparent_recorded_not_absence.csv')
            preds=[csv(z,n+'_spatial_holdout_predictions.csv') for n in ('current_climate_only','plus_NPP_AGB_NDVI','plus_LGM_history')]
            names=('climate','vegetation','LGM')
            base=preds[0].drop(columns='variant').rename(columns={'p_fixed_only':'p_climate'})
            for name,p in zip(names[1:],preds[1:]):
                assert not p[['site_row','species_key']].duplicated().any()
                assert len(p)==len(base)
                p=p.drop(columns='variant').rename(columns={'p_fixed_only':'p_'+name,'recorded':'Y_check'})
                base=base.merge(p,on=['site_row','species_key'],validate='one_to_one')
                same=base['recorded'].eq(base.Y_check)|(base['recorded'].isna()&base.Y_check.isna())
                assert same.all(),f'Y discrepancy {region}'
                base=base.drop(columns='Y_check')
            assert not base[['site_row','species_key']].duplicated().any()
            base['preflight_site_id']=base.site_row.astype(int)-1  # R is 1-based; source X/Y 0-based
            assert set(base['preflight_site_id']).issubset(set(env.site_id))
            assert set(base.species_key).issubset(set(tr.key))
            base['white_proxy']=base.species_key.map(dict(zip(tr.key,tr.white_proxy)))
            assert not base.white_proxy.isna().any()
            base['block_500km']=base.preflight_site_id.map(dict(zip(env.site_id,env.block_500km)))
            assert base.block_500km.notna().all()
            orig=np.array([Ys.at[int(i),s] for i,s in zip(base.preflight_site_id,base.species_key)],dtype=float)
            curr=base.recorded.to_numpy(dtype=float)
            assert np.all((np.isnan(orig)&np.isnan(curr))|(orig==curr)),f'Original Y wrong {region}'
            valid=base.dropna(subset=['recorded']).copy()
            assert set(valid.recorded.unique())=={0,1}
            for name in names:
                pp=np.clip(valid[f'p_{name}'].to_numpy(dtype=float),1e-12,1-1e-12)
                yy=valid.recorded.to_numpy(dtype=float)
                valid[f'll_{name}']=-(yy*np.log(pp)+(1-yy)*np.log1p(-pp))
            valid['delta_veg']=valid.ll_vegetation-valid.ll_climate
            valid['delta_LGM']=valid.ll_LGM-valid.ll_vegetation
            comparison=csv(z,'500km_unseen_block_model_comparison.csv').set_index('variant')
            for var,code in [('current_climate_only','climate'),('plus_NPP_AGB_NDVI','vegetation'),('plus_LGM_history','LGM')]:
                known=float(comparison.loc[var,'log_loss'])
                observed=float(valid[f'll_{code}'].mean())
                assert abs(observed-known)<1e-6,(region,var,known,observed)
            diagnostics.append(dict(region=region,number_species=base.species_key.nunique(),heldout_sites=base.site_row.nunique(),valid_species_cells=len(valid),n500kmblocks=valid.block_500km.nunique(),source_Y_disagreements=0))
            for colour in (0,1):
                group=valid[valid.white_proxy==colour]
                for v in ('delta_veg','delta_LGM'):
                    ci,blocks=block_ci(valid,'white_proxy',colour,v)
                    allrows.append(dict(region=region,group='white' if colour else 'nonwhite',comparison=v,n=int(len(group)),n_500km_blocks=blocks,mean_diff_logloss=float(group[v].mean()),CI95_lower=ci[0],CI95_upper=ci[1],p_value='not_computed',interpretation='negative_is_better'))
            for v in ('delta_veg','delta_LGM'):
                ci,blocks=contrast_ci(valid,v)
                diff=valid.loc[valid.white_proxy==1,v].mean()-valid.loc[valid.white_proxy==0,v].mean()
                allrows.append(dict(region=region,group='white_minus_nonwhite',comparison=v,n=len(valid),n_500km_blocks=blocks,mean_diff_logloss=float(diff),CI95_lower=ci[0],CI95_upper=ci[1],p_value='not_computed',interpretation='negative_means_larger_gain_for_white'))
            for v in ('delta_veg','delta_LGM'):
                for rr in species_equal_sensitivity(valid,v):
                    rr['region']=region
                    allrows.append(rr)
            for y in (0,1):
                group=valid[valid.recorded==y]
                allrows.append(dict(region=region,group='recorded_'+str(y),comparison='delta_LGM',n=len(group),n_500km_blocks=group.block_500km.nunique(),mean_diff_logloss=float(group.delta_LGM.mean()),CI95_lower=np.nan,CI95_upper=np.nan,p_value='not_computed',interpretation='negative_is_better'))
    out=pd.DataFrame(allrows)
    out['weighting']=out['weighting'].fillna('equal_species_cells')
    dest=ROOT/'mois_Hmsc_LGM_colour_stratified_audit_20261117.csv'
    out.to_csv(dest,index=False,float_format='%.10f')
    receipt=dict(method='exploratory_reanalysis_not_independent_holdout',preflight_sha=PREFLIGHT[1],sources={r:s for r,_,s in SOURCES},n_block_bootstrap=REPS,rng_seed=20261117,verification='source Y incl NA, paired heldout keys; full LL rederived within 1e-6 due exported probability precision',original_panels_selected_using_all_region_GBIF_recording=True,photo_white_is_unverified_proxy=True,uncertainty='single_holdout_split_conditions_only',regions=diagnostics)
    with open(ROOT/'mois_Hmsc_LGM_colour_stratified_audit_20261117_receipt.json','w') as f: json.dump(receipt,f,ensure_ascii=False,indent=2)
    print(out.to_string(index=False,float_format=lambda x:f'{x:+.5f}'))
    print('SHA csv',hashlib.sha256(dest.read_bytes()).hexdigest())
    return out

if __name__=='__main__':compute()
