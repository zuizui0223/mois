#!/usr/bin/env python3
"""Source-locked taxon composition and environmental associations; numeric data only.

No new images, ROI inference, MNAR scenarios or unknown-label imputation.
Palette composition is NOT botanical inconspicuousness or calibrated colour.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

PALETTE=('white','yellow','orange','red','pink','magenta','purple','blue','bronze')
FRACTIONS=tuple('strict9_'+c for c in PALETTE)
METRICS=FRACTIONS+('palette_entropy','palette_dominance')
MAIN='historical_K1_FCP_source_heterogeneous'
PINS={
 'canonical':'f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c',
 'palette':'faf6951505f1002a1327f64a2ced13431b562cf13893707125fafff918160f62',
 'geography':'b60d643fb255abeec4f88ef4cfb38a040fdd24eb635baf26c16fbccb5711c562'}
CLIMATE=('wc2_elev_2p5m','wc2_bio05_10m','wc2_bio12_10m','annual_solar_Wm2')

def require(ok,msg):
    if not bool(ok):raise ValueError('MOIS_COMPOSITION: '+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()

def join_palette(canonical,palette):
    keys=['photo_id','observation_id','inat_taxon_id']
    a=canonical.copy();b=palette.copy()
    for table in (a,b):
        for col in keys:table[col]=table[col].fillna('').astype(str)
        require(table.photo_id.is_unique,'Duplicate photo ID')
        require(table[keys].ne('').all().all(),'Missing exact source identity')
        require(pd.api.types.is_bool_dtype(table.strict9.dtype),'Nonboolean strict9')
    require(a.observation_id.is_unique,'Canonical observation UUID repeated')
    b=b.set_index('photo_id')
    require(a.photo_id.isin(b.index).all(),'Canonical photo absent from palette')
    aligned=b.loc[a.photo_id]
    for col in ('observation_id','inat_taxon_id','strict9'):
        require(np.array_equal(a[col].to_numpy(),aligned[col].to_numpy()),'ID or strict9 drift: '+col)
    require(np.allclose(a.strict9_white.to_numpy(float),aligned.strict9_white.to_numpy(float),atol=1e-7,equal_nan=True),'White fraction drift')
    for col in FRACTIONS:a[col]=aligned[col].to_numpy(float)
    return palette_traits(a)

def palette_traits(d):
    d=d.copy();ok=d.strict9.to_numpy(bool);v=d[list(FRACTIONS)].to_numpy(float)
    require(np.isnan(v[~ok]).all(),'Technical unknown contains colours')
    require(np.isfinite(v[ok]).all() and ((v[ok]>=0)&(v[ok]<=1)).all(),'Invalid strict9 components')
    require(np.allclose(v[ok].sum(1),1,atol=1e-5),'Nine palette does not close to one')
    entropy=np.full(len(d),np.nan);dominance=entropy.copy()
    entropy[ok]=-np.sum(v[ok]*np.log(np.maximum(v[ok],1e-30)),axis=1)/np.log(9)
    dominance[ok]=v[ok].max(axis=1)
    d['palette_entropy']=entropy;d['palette_dominance']=dominance
    return d

def region_masks(d):
    x=pd.to_numeric(d.latitude,errors='coerce')
    return {'arctic':x.between(66.5622,90),'north_subpolar':x.ge(55)&x.lt(66.5622),
            'north_temperate':x.ge(23.5)&x.lt(55),'south_temperate':x.gt(-55)&x.le(-23.5),
            'island_cached':d.land_type.eq('island'),'mainland_cached':d.land_type.eq('mainland')}

def decompose(a,b,metric):
    """Exact equal-taxon four-term identity; not a causal or phylogenetic partition."""
    require(a.index.is_unique and b.index.is_unique and len(a)>0 and len(b)>0,'Invalid taxon summaries')
    require(np.isfinite(a[metric]).all() and np.isfinite(b[metric]).all(),'Missing measured taxon trait')
    common=a.index.intersection(b.index);only_a=a.index.difference(b.index);only_b=b.index.difference(a.index)
    wa=1/len(a);wb=1/len(b)
    within=float((.5*(wa+wb)*(a.loc[common,metric]-b.loc[common,metric])).sum())
    weights=float((.5*(wa-wb)*(a.loc[common,metric]+b.loc[common,metric])).sum())
    exclusive_a=float(wa*a.loc[only_a,metric].sum());exclusive_b=float(-wb*b.loc[only_b,metric].sum())
    delta=float(a[metric].mean()-b[metric].mean())
    require(abs(delta-within-weights-exclusive_a-exclusive_b)<1e-10,'Decomposition does not reconcile')
    return {'focal_taxa':len(a),'reference_taxa':len(b),'shared_taxa':len(common),
            'focal_taxa_shared_fraction':len(common)/len(a),'reference_taxa_shared_fraction':len(common)/len(b),
            'focal_observations':int(a.n.sum()),'reference_observations':int(b.n.sum()),
            'observation_weighted_difference_pp':100*float(np.average(a[metric],weights=a.n)-np.average(b[metric],weights=b.n)),
            'taxon_equal_focal':float(a[metric].mean()),'taxon_equal_reference':float(b[metric].mean()),
            'taxon_equal_difference_pp':100*delta,'shared_taxon_outcome_pp':100*within,
            'shared_taxon_weight_pp':100*weights,'focal_only_taxa_pp':100*exclusive_a,
            'reference_only_taxa_pp':100*exclusive_b,
            'common_taxa_equal_difference_pp':100*float((a.loc[common,metric]-b.loc[common,metric]).mean()) if len(common) else None}

def within_fit(x,y,w,groups,blocks):
    """Weighted FE regression and nominal spatial-cluster CR1 covariance."""
    x=np.asarray(x,float);y=np.asarray(y,float);w=np.asarray(w,float)
    groups=np.asarray(groups,int);blocks=np.asarray(blocks,int)
    require(x.ndim==2 and y.ndim==2 and len(x)==len(y)==len(w),'Bad model shape')
    require(np.isfinite(x).all() and np.isfinite(y).all() and np.isfinite(w).all() and (w>0).all(),'Invalid model values')
    ng=int(groups.max())+1;nb=int(blocks.max())+1;n,p=x.shape
    require(nb>=2 and n>ng+p,'Insufficient residual support')
    den=np.bincount(groups,weights=w,minlength=ng)
    require((den>0).all(),'Noncontiguous empty FE codes')
    xc=x.copy();yc=y.copy()
    for j in range(p):xc[:,j]-=(np.bincount(groups,weights=w*x[:,j])/den)[groups]
    for j in range(y.shape[1]):yc[:,j]-=(np.bincount(groups,weights=w*y[:,j])/den)[groups]
    gram=xc.T@(w[:,None]*xc)
    require(np.linalg.matrix_rank(gram)==p,'Environmental model not identified')
    inv=np.linalg.inv(gram);beta=inv@(xc.T@(w[:,None]*yc));res=yc-xc@beta
    factor=nb/(nb-1)*(n-1)/(n-ng-p);se=np.empty_like(beta)
    for j in range(y.shape[1]):
        scores=np.stack([np.bincount(blocks,weights=w*xc[:,k]*res[:,j],minlength=nb) for k in range(p)],axis=1)
        cov=factor*inv@(scores.T@scores)@inv
        se[:,j]=np.sqrt(np.maximum(np.diag(cov),0))
    crit=student_t.ppf(.975,nb-1)
    return {'beta':beta,'se':se,'low':beta-crit*se,'high':beta+crit*se,
            'within_vif':np.diag(gram)*np.diag(inv),'n':n,'groups':ng,'blocks':nb,
            'condition_number':float(np.linalg.cond(gram)),
            'elevation_informative_rows':int((np.abs(xc[:,0])>1e-10).sum())}

def regional_tables(d,out):
    masks=region_masks(d);coverage=[];profiles=[];tables={}
    for sector,mask in masks.items():
        f=d.loc[mask]
        for (arm,subarm),g in f.groupby(['source_design','source_subarm'],dropna=False,observed=True):
            coverage.append({'sector':sector,'source_design':arm,'source_subarm':str(subarm),
                             'observations':len(g),'strict9_success':int(g.strict9.sum()),
                             'unknown':int((~g.strict9).sum()),'attempted_taxa':g.inat_taxon_id.nunique(),
                             'successful_taxa':g.loc[g.strict9,'inat_taxon_id'].nunique()})
        good=f.loc[f.source_design.eq(MAIN)&f.strict9]
        t=good.groupby('inat_taxon_id',observed=True)[list(METRICS)].mean()
        t['n']=good.groupby('inat_taxon_id',observed=True).size();tables[sector]=t
        if len(t):
            q=t.reset_index();q.insert(0,'sector',sector);profiles.append(q)
    pd.DataFrame(coverage).to_csv(out/'region_source_coverage.csv',index=False)
    pd.concat(profiles,ignore_index=True).to_csv(out/'taxon_region_palette.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    contrasts=[]
    for a,b in [('arctic','north_subpolar'),('south_temperate','north_temperate'),('island_cached','mainland_cached')]:
        if not len(tables[a]) or not len(tables[b]):continue
        for metric in METRICS:
            row={'focal':a,'reference':b,'metric':metric,**decompose(tables[a],tables[b],metric)}
            contrasts.append(row)
    pd.DataFrame(contrasts).to_csv(out/'taxon_composition_contrasts.csv',index=False)
    return contrasts

def elevation_models(d,out):
    climate=d.loc[d.source_design.eq(MAIN)&d.climate_complete&d.latitude.between(-90,90)&d.longitude.between(-180,180)].copy()
    climate=climate.loc[np.isfinite(climate[list(CLIMATE)].to_numpy(float)).all(1)&climate.wc2_bio12_10m.ge(0)]
    raw=climate[list(CLIMATE)].to_numpy(float);raw[:,2]=np.log1p(raw[:,2])
    means=raw.mean(0);sds=raw.std(0);require((sds>0).all(),'Constant predictor')
    standardized=(raw-means)/sds
    ok=climate.strict9.to_numpy(bool);f=climate.loc[ok].copy()
    X=standardized[ok];Y=f[list(METRICS)].to_numpy(float)
    taxon_n=f.groupby('inat_taxon_id',observed=True).photo_id.transform('size').to_numpy(float);w=1/taxon_n
    x5=np.clip(np.floor((f.longitude.to_numpy()+180)/5).astype(int),0,71)
    y5=np.clip(np.floor((f.latitude.to_numpy()+90)/5).astype(int),0,35)
    cells=pd.Series(y5).astype(str)+'_'+pd.Series(x5).astype(str)
    group=f.source_subarm.fillna('unknown').reset_index(drop=True)+'|'+cells+'|'+f.local_season.fillna('unknown').reset_index(drop=True)
    blocks=pd.Series(y5//2).astype(str)+'_'+pd.Series(x5//2).astype(str)
    g,_=pd.factorize(group,sort=True);b,_=pd.factorize(blocks,sort=True)
    rows=[];diagnostics={}
    for name,cols in [('regional_elevation',[0]),('regional_environment',[0,1,2,3])]:
        result=within_fit(X[:,cols],Y,w,g,b)
        require(np.allclose(result['beta'][:,:9].sum(1),0,atol=1e-7),'Compositional slopes do not sum to zero')
        diagnostics[name]={k:result[k] for k in ('n','groups','blocks','condition_number','elevation_informative_rows')}
        for k,i in enumerate(cols):
            for j,metric in enumerate(METRICS):
                rows.append({'model':name,'predictor':CLIMATE[i],'outcome':metric,
                             'effect_pp_per_source_sd':100*float(result['beta'][k,j]),
                             'nominal_spatial95_low_pp':100*float(result['low'][k,j]),
                             'nominal_spatial95_high_pp':100*float(result['high'][k,j]),
                             'within_vif':float(result['within_vif'][k]),'successful_observations':len(f),
                             'taxa':f.inat_taxon_id.nunique(),'fixed_effect_groups':result['groups'],'spatial_blocks':result['blocks']})
    pd.DataFrame(rows).to_csv(out/'elevation_environment_coefficients.csv',index=False)
    return {'attempted_climate_complete_observations':len(climate),'strict9_success':len(f),'nominal_taxa':int(f.inat_taxon_id.nunique()),
            'predictor_columns':list(CLIMATE),'predictor_means':means.tolist(),'predictor_sds':sds.tolist(),
            'predictor_correlations':np.corrcoef(raw,rowvar=False).tolist(),'weights':'equal nominal taxon total weight among successful observations',
            'models':diagnostics,'same_observations_weights_and_FE_in_both_models':True,'causal_mediation_identified':False},rows

def run(canonical,palette,geography,output):
    paths={'canonical':Path(canonical),'palette':Path(palette),'geography':Path(geography)}
    out=Path(output);require(not out.exists(),'Refuse to overwrite results')
    hashes={k:sha(p) for k,p in paths.items()};require(hashes==PINS,'Frozen source SHA changed')
    a=pd.read_csv(paths['canonical'],dtype={'photo_id':str,'observation_id':str,'inat_taxon_id':str},low_memory=False)
    p=pd.read_parquet(paths['palette'])
    require(len(a)==1032882 and int(a.strict9.sum())==401703,'Canonical source census changed')
    require(len(p)==1181146 and int(p.strict9.sum())==451643,'Original palette census changed')
    d=join_palette(a,p);del p,a
    geo=pd.read_csv(paths['geography'],dtype={'photo_id':str},low_memory=False)
    require(geo.photo_id.is_unique and len(geo)==313535,'Geography cache drift')
    d=d.merge(geo,on='photo_id',how='left',validate='one_to_one',indicator='_geo')
    d.land_type=d.land_type.fillna('not_in_old_cache')
    out.mkdir(parents=True)
    contrasts=regional_tables(d,out);models,coefficients=elevation_models(d,out)
    capabilities={'measured':list(METRICS),'not_measured':['green_as_independent_component','botanically_verified_brown_petals','calibrated_Lab_chroma','background_contrast','physical_display_size','perianth_reduction'],
                  'inconspicuousness_score_created':False,'bronze_not_verified_brown':True,'white_not_low_conspicuousness':True,
                  'taxonomic_level':'original nominal inat_taxon_id, no genus/phylogeny crosswalk yet',
                  'plant_abundance_or_native_flora_census':False}
    (out/'trait_availability.json').write_text(json.dumps(capabilities,indent=2)+'\n')
    receipt={'schema':'mois_original_canonical_taxon_palette_environment_20261030_v1','input_sha256':hashes,
             'canonical_observations':len(d),'strict9_success':int(d.strict9.sum()),'unknown':int((~d.strict9).sum()),
             'main_source_observations':int(d.source_design.eq(MAIN).sum()),'main_source_success':int((d.source_design.eq(MAIN)&d.strict9).sum()),
             'exact_cached_geography_matches':int(d._geo.eq('both').sum()),'island_analysis_not_expanded_beyond_old_cache':True,
             'new_photos':0,'new_inference':0,'unknown_colours_imputed':False,'MNAR_scenarios_run':False,'adaptive_arms_pooled_into_primary':False,
             'primary_target':'classifiable photographed nominal taxa in the historical source, NOT worldwide flora',
             'environment':models,'regional_metrics':len(contrasts),'environmental_coefficients':len(coefficients),
             'interpretation':'Exploratory descriptive composition and conditional associations; no botanical dullness, phylogenetic or causal attribution',
             'output_sha256':{f.name:sha(f) for f in sorted(out.iterdir())}}
    (out/'composition_environment_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('MOIS_COMPOSITION_ENVIRONMENT_COMPLETED',json.dumps({k:receipt[k] for k in ('canonical_observations','strict9_success','exact_cached_geography_matches','regional_metrics','environmental_coefficients')},sort_keys=True))
    for row in contrasts:
        if row['metric'] in ('strict9_white','strict9_bronze','palette_entropy'):print('REGIONAL_RESULT',json.dumps(row,sort_keys=True))
    for row in coefficients:
        if row['predictor']=='wc2_elev_2p5m' and row['outcome'] in ('strict9_white','palette_entropy'):print('ELEVATION_RESULT',json.dumps(row,sort_keys=True))
    return receipt

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('canonical','palette','geography','output'):parser.add_argument('--'+name,required=True,type=Path)
    args=parser.parse_args();run(args.canonical,args.palette,args.geography,args.output)
