#!/usr/bin/env python3
"""Source-locked, no-human moisture association refit on existing photographs.

Selection weights describe technical ROI availability under conditional MAR.
They do NOT correct volunteer recruitment, validate petals, or recover native
flora abundance. Climate predictors are annual normals, not observation weather.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA='mois_existing_observation_moisture_MAR_refit_20261030_v1'
HUES=('white','yellow_orange','red_pink','blue_purple')
CLIMATE=('wc2_bio05_10m','wc2_bio12_10m','wc2_bio14_10m','wc2_annual_vapr_10m','annual_solar_Wm2','wc2_elev_2p5m')
CONTROLS=('wc2_bio05_10m','annual_solar_Wm2','wc2_elev_2p5m')
OUTCOMES=('white_photo_label','white_ROI_fraction')
SOURCE_SHA='faf6951505f1002a1327f64a2ced13431b562cf13893707125fafff918160f62'
CLIMATE_SHA='ec39b49649acc2ca4a468c00c86d05800e33ea2db1f91a3b119076be38432735'

def require(condition, message):
    if not bool(condition): raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda:fh.read(4<<20),b''):h.update(block)
    return h.hexdigest()

def validate_frame(d):
    require(d.photo_id.notna().all() and d.observation_id.notna().all(),'Missing source ID')
    require(d.photo_id.is_unique and d.observation_id.is_unique,'Duplicate source IDs')
    require(pd.api.types.is_bool_dtype(d.strict9.dtype) and d.strict9.notna().all(),'ROI must be boolean')
    ok=d.strict9.to_numpy()
    require(d.loc[ok,'original_roi_four_hue'].isin(HUES).all(),'Invalid successful hue')
    require(d.loc[ok,'strict9_white'].between(0,1).all(),'Invalid successful white fraction')
    require(d.loc[~ok,'strict9_white'].isna().all(),'Technical unknown recoloured')

def strata_weights(frame):
    """Return observed counts and N/K weights on successes ONLY; no imputation."""
    validate_frame(frame)
    z=frame[['stratum','block','q','strict9','original_roi_four_hue','strict9_white']].copy()
    z['W']=z.strict9 & z.original_roi_four_hue.eq('white')
    z['S']=z.strict9_white.where(z.strict9,0.)
    t=z.groupby('stratum',sort=True,observed=True).agg(
        N=('strict9','size'),K=('strict9','sum'),W=('W','sum'),S=('S','sum'),
        block=('block','first'),q=('q','first'))
    t['weight']=t.N/t.K.replace(0,np.nan)
    w=frame.stratum.map(t.weight).to_numpy(float)
    w[~frame.strict9.to_numpy()]=np.nan
    require(int(t.N.sum())==len(frame),'Wrong stratum denominator')
    return t,w

def within_block_sums(x,y,w,groups,blocks,nblocks=None):
    """Weighted within-group least squares, aggregated to spatial blocks."""
    x=np.asarray(x,float);y=np.asarray(y,float);w=np.asarray(w,float)
    groups=np.asarray(groups,int);blocks=np.asarray(blocks,int)
    require(x.ndim==2 and y.ndim==2 and len(x)==len(y)==len(w),'Invalid arrays')
    require(np.isfinite(x).all() and np.isfinite(y).all() and np.isfinite(w).all() and (w>0).all(),'Nonfinite model inputs')
    den=np.bincount(groups,weights=w)
    xc=x.copy();yc=y.copy()
    for j in range(x.shape[1]):
        mean=np.divide(np.bincount(groups,weights=w*x[:,j]),den,out=np.zeros_like(den),where=den>0)
        xc[:,j]-=mean[groups]
    for j in range(y.shape[1]):
        mean=np.divide(np.bincount(groups,weights=w*y[:,j]),den,out=np.zeros_like(den),where=den>0)
        yc[:,j]-=mean[groups]
    nb=int(blocks.max())+1 if nblocks is None else int(nblocks)
    p=x.shape[1];n_y=y.shape[1]
    A=np.zeros((nb,p,p));b=np.zeros((nb,p,n_y))
    for j in range(p):
        for k in range(p): A[:,j,k]=np.bincount(blocks,weights=w*xc[:,j]*xc[:,k],minlength=nb)
        for k in range(n_y): b[:,j,k]=np.bincount(blocks,weights=w*xc[:,j]*yc[:,k],minlength=nb)
    return A,b

def solve_bootstrap(A,b,draws):
    At=np.einsum('rb,bjk->rjk',draws,A,optimize=True)
    bt=np.einsum('rb,bjk->rjk',draws,b,optimize=True)
    ans=np.full_like(bt,np.nan)
    for i in range(len(draws)):
        try: ans[i]=np.linalg.solve(At[i],bt[i])
        except np.linalg.LinAlgError: pass
    return ans

def interval(x):
    require(np.isfinite(x).sum()>=.95*len(x),'Too many degenerate spatial resamples')
    return [float(z) for z in np.nanquantile(x,[.025,.975])]

def run(joined,receipt,out,*,replicates=999):
    joined,receipt,out=map(Path,(joined,receipt,out))
    require(not out.exists(),'Refuse to overwrite an analysis result')
    rec=json.loads(receipt.read_text())
    require(rec['source_photo_sha256']==SOURCE_SHA and rec['climate_sha256']==CLIMATE_SHA,'Wrong frozen source')
    require(sha(joined)==rec['output_sha256'][joined.name],'Joined input SHA changed')
    d=pd.read_csv(joined,dtype={'photo_id':str,'observation_id':str,'inat_taxon_id':str,'source_subarm':str},low_memory=False)
    validate_frame(d)
    require(len(d)==1032882 and int(d.strict9.sum())==401703,'Canonical input census changed')
    geo=d.latitude.between(-90,90)&d.longitude.between(-180,180)
    valid=d.climate_complete & geo & np.isfinite(d[list(CLIMATE)].to_numpy(float)).all(axis=1)
    valid &= d.wc2_bio12_10m.ge(0) & d.wc2_annual_vapr_10m.ge(0)
    f=d.loc[valid].copy().reset_index(drop=True)
    require(len(f)>0 and f.climate_join.eq('both').all(),'No genuine climate support')
    require(int(f.strict9.sum())>100000,'Unexpected loss of original successful images')
    x5=np.clip(np.floor((f.longitude.to_numpy(float)+180)/5).astype(int),0,71)
    y5=np.clip(np.floor((f.latitude.to_numpy(float)+90)/5).astype(int),0,35)
    f['cell5']=pd.Series(y5).astype(str)+'_'+pd.Series(x5).astype(str)
    f['block']=pd.Series(y5//2).astype(str)+'_'+pd.Series(x5//2).astype(str)
    f['group']=f.source_design.astype(str)+'|'+f.source_subarm.fillna('unknown')+'|'+f.cell5+'|'+f.local_season.fillna('unknown')
    block_codes,block_labels=pd.factorize(f.block,sort=True);nb=len(block_labels)
    require(nb>=30,'Insufficient spatial blocks')
    rng=np.random.default_rng(20261030)
    draws=rng.multinomial(nb,np.full(nb,1/nb),size=replicates)
    ok=f.strict9.to_numpy()
    Y=np.c_[f.loc[ok,'original_roi_four_hue'].eq('white').to_numpy(float),f.loc[ok,'strict9_white'].to_numpy(float)]
    groups,_=pd.factorize(f.loc[ok,'group'],sort=True)
    control=f[list(CONTROLS)].to_numpy(float)
    ctrmean=control.mean(axis=0);ctrsd=control.std(axis=0)
    require((ctrsd>0).all(),'Constant control variable')
    control=(control-ctrmean)/ctrsd
    out.mkdir(parents=True)
    coeff_rows=[];profile_rows=[];contrast_rows=[];diagnostics={}
    for exposure,col,transform in [('rainfall','wc2_bio12_10m','log1p'),('vapour','wc2_annual_vapr_10m','identity')]:
        raw=f[col].to_numpy(float)
        cuts=np.quantile(raw,[.2,.4,.6,.8])
        require(len(np.unique(cuts))==4,'Exposure bins have tied boundaries')
        f['q']=np.searchsorted(cuts,raw,side='right')
        f['stratum']=f['group']+'|Q'+f.q.astype(str)
        t,ipw=strata_weights(f)
        t.to_csv(out/(exposure+'_response_strata.csv.gz'),compression={'method':'gzip','mtime':0})
        supported=t.loc[t.K.gt(0)]
        z=np.log1p(raw) if transform=='log1p' else raw.copy()
        mean=float(z.mean());sd=float(z.std());require(sd>0,'Constant exposure')
        X=np.c_[(z-mean)/sd,control][ok]
        q=f.q.to_numpy()[ok];bc=block_codes[ok]
        n_unsupported=int(t.loc[t.K.eq(0),'N'].sum())
        diagnostics[exposure]={'raw_quintile_boundaries':cuts.tolist(),'transform':transform,'transformed_mean':mean,'transformed_sd_all_eligible':sd,'response_strata':len(t),'unsupported_strata':int(t.K.eq(0).sum()),'unsupported_observations':n_unsupported,'supported_observations':int(supported.N.sum()),'strata_with_1_to_4_success':int(t.K.between(1,4).sum()),'observations_in_strata_with_1_to_4_success':int(t.loc[t.K.between(1,4),'N'].sum()),'max_technical_response_weight':float(np.nanmax(ipw)),'effective_success_N':float(np.nansum(ipw)**2/np.nansum(ipw**2))}
        boots={};contrastboots={}
        for mode,weights in [('unweighted',np.ones(ok.sum())),('MAR_IPW',ipw[ok])]:
            require(np.isfinite(weights).all() and (weights>=1).all(),'Invalid response weight')
            den=np.bincount(bc*5+q,weights=weights,minlength=nb*5).reshape(nb,5)
            nums=np.stack([np.bincount(bc*5+q,weights=weights*Y[:,j],minlength=nb*5).reshape(nb,5) for j in range(2)],axis=2)
            for qi in range(5):
                totalN=int(t.loc[t.q.eq(qi),'N'].sum());totalK=int(t.loc[t.q.eq(qi),'K'].sum())
                for j,outcome in enumerate(OUTCOMES):
                    profile_rows.append({'exposure':exposure,'outcome':outcome,'method':mode,'quintile':qi+1,'all_observations':totalN,'strict9_success':totalK,'classifiability':totalK/totalN,'effective_denominator':float(den[:,qi].sum()),'fraction':float(nums[:,qi,j].sum()/den[:,qi].sum())})
            bootden=draws@den
            bootnum=np.einsum('rb,bqk->rqk',draws,nums,optimize=True)
            diffs=100*(bootnum[:,4,:]/bootden[:,4,None]-bootnum[:,0,:]/bootden[:,0,None])
            contrastboots[mode]=diffs
            for j,outcome in enumerate(OUTCOMES):
                contrast=100*(nums[:,4,j].sum()/den[:,4].sum()-nums[:,0,j].sum()/den[:,0].sum())
                contrast_rows.append({'exposure':exposure,'outcome':outcome,'method':mode,'wettest_minus_driest_pp':float(contrast),'bootstrap95_pp':interval(diffs[:,j]),'spatial_blocks':nb})
            A,b=within_block_sums(X,Y,weights,groups,bc,nblocks=nb)
            totalA=A.sum(0);totalb=b.sum(0)
            require(np.linalg.matrix_rank(totalA)==4,'Within model not identified')
            est=np.linalg.solve(totalA,totalb)
            boot=solve_bootstrap(A,b,draws)
            boots[mode]=boot
            for j,outcome in enumerate(OUTCOMES):
                coeff_rows.append({'exposure':exposure,'outcome':outcome,'method':mode,'slope_pp_per_sd':float(100*est[0,j]),'bootstrap95_pp':interval(100*boot[:,0,j]),'control_coefficients':{k:float(100*est[i+1,j]) for i,k in enumerate(CONTROLS)},'condition_number':float(np.linalg.cond(totalA)),'successful_observations':int(ok.sum()),'fixed_effect_groups':int(groups.max()+1),'spatial_blocks':nb,'bootstrap_finite':int(np.isfinite(boot[:,0,j]).sum())})
        diagnostics[exposure]['paired_slope_change_IPW_minus_raw_95_pp']={o:interval(100*(boots['MAR_IPW'][:,0,j]-boots['unweighted'][:,0,j])) for j,o in enumerate(OUTCOMES)}
        diagnostics[exposure]['paired_wet_dry_change_IPW_minus_raw_95_pp']={o:interval(contrastboots['MAR_IPW'][:,j]-contrastboots['unweighted'][:,j]) for j,o in enumerate(OUTCOMES)}
    pd.DataFrame(profile_rows).to_csv(out/'quintile_profiles.csv',index=False)
    summary={'schema':SCHEMA,'joined_input_sha256':sha(joined),'join_receipt_sha256':sha(receipt),'source_photo_sha256':SOURCE_SHA,'climate_sha256':CLIMATE_SHA,'canonical_observations':len(d),'canonical_roi9_success':int(d.strict9.sum()),'climate_id_matched':int(d.climate_join.eq('both').sum()),'analysis_eligible_observations':len(f),'analysis_roi9_success':int(ok.sum()),'analysis_roi9_missing':int((~ok).sum()),'not_in_climate_analysis':len(d)-len(f),'analysis_source_design_counts':{str(k):int(v) for k,v in f.source_design.value_counts().items()},'climate_source':rec['climate_source'],'bootstrap_replicates':replicates,'bootstrap_seed':20261030,'bootstrap_unit':'10-degree latitude/longitude blocks; not equal-area; approximate spatial replication','fixed_effects':'original source design/subarm x 5-degree cell x hemisphere-local season','control_columns':list(CONTROLS),'control_means':ctrmean.tolist(),'control_sds':ctrsd.tolist(),'diagnostics':diagnostics,'within_environment_slopes':coeff_rows,'source_wet_dry_contrasts':contrast_rows,'new_images':0,'human_labels':0,'unknown_colours_imputed':False,'MAR_assumption_not_tested':True,'native_world_flora_estimand_identified':False,'independent_rain_vapour_causal_effects_identified':False,'multiple_comparisons':'primary rainfall; vapour and continuous ROI fraction supporting descriptive checks; nominal intervals not multiplicity-adjusted','uncertainty_limit':'Spatial intervals do not include MNAR, classification error or volunteer recruitment bias.'}
    summary['output_sha256']={p.name:sha(p) for p in sorted(out.iterdir())}
    (out/'moisture_bias_refit_receipt.json').write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('MOIS_MOISTURE_REFIT_COMPLETED',json.dumps({'eligible':len(f),'success':int(ok.sum()),'models':len(coeff_rows)},sort_keys=True),flush=True)
    for row in coeff_rows:print(json.dumps(row,sort_keys=True),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--joined',required=True,type=Path);p.add_argument('--receipt',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();run(a.joined,a.receipt,a.output)
if __name__=='__main__':main()
