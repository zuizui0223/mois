#!/usr/bin/env python3
"""FCP macroecological ABIO-HUE direction atlas (2026-10-13, source exposed).

Tests DIRECTION of photographed four-colour proportions along a *frozen* suite
of environmental axes. The population is original species means, with equal
species or equal-genus weights, controlling nominal genus and geographic trend.
Global-to-sector sign portability, coarse same-cell conditioning, strict-photo
quality and complete-climate/soil opportunity are mandatory negative controls.
These are retrospective observational PHOTO traits, NOT selection/fitness or
frequency of genetic colour morphs. Do not claim geographic causal independence.
"""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_PHOTO_SHA='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
SOURCE_ENV_SHA='6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca'
CLASSES=('white','yellow_orange','red_pink','blue_purple')
BLOCKS={
 'elevation':('wc_elevation_m',),
 'temperature':tuple(f'wc_bio{i}' for i in range(1,12)),
 'precipitation':tuple(f'wc_bio{i}' for i in range(12,20)),
 'solar_radiation':('wc_srad_annual_mean','wc_srad_monthly_cv'),
 'wind':('wc_wind_annual_mean','wc_wind_monthly_cv'),
 'vapor_pressure':('wc_vapr_annual_mean','wc_vapr_monthly_cv'),
 'soil':('soil_pH','soil_SOC','soil_N','soil_clay','soil_available_water_proxy',
         'soil_cec_0_30cm_source_raw','soil_sand_0_30cm_source_raw',
         'soil_silt_0_30cm_source_raw','soil_bdod_0_30cm_source_raw',
         'soil_cfvo_0_30cm_source_raw'),
}
ALL_FEATURES=tuple(v for block in BLOCKS.values() for v in block)
assert len(ALL_FEATURES)==36
# Ecologically interpretable directional anchors chosen from variable semantics,
# NOT from outcome-dependent significant-feature scanning.
ANCHOR_CLIMATE=('wc_bio5','wc_bio12','wc_bio14','wc_bio15',
                'wc_srad_annual_mean','wc_vapr_annual_mean')
ANCHOR_SOIL=ANCHOR_CLIMATE+('soil_available_water_proxy',)
GEO=('abs_lat','lon_sin','lon_cos')
SECTORS=('W180_W90','W90_E0','E0_E90','E90_E180')
BOOT_SEED=2026101348
VERSION='fcp_macroecology_colour_environment_direction_genus_geosector_v1'


def sha_ok(path:Path,expected:str)->None:
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if h!=expected:raise RuntimeError(f'original archive SHA mismatch {path}: {h}')


def read_archive(path:Path,member:str,usecols=None)->pd.DataFrame:
    with zipfile.ZipFile(path) as z:
        return pd.read_csv(io.BytesIO(z.read(member)),compression='gzip',usecols=usecols,low_memory=False)


def original_source(photo_zip:Path,env_zip:Path):
    sha_ok(photo_zip,SOURCE_PHOTO_SHA);sha_ok(env_zip,SOURCE_ENV_SHA)
    cols=['inat_taxon_id','observation_id','photo_id','species','morph',
          'measurement_status','flower_effective_pixels','nuisance_pixel_fraction']
    a=read_archive(photo_zip,'fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz',cols)
    b=read_archive(photo_zip,'fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz',cols)
    d=pd.concat((a,b),ignore_index=True)
    if len(a)!=42111 or len(b)!=85337 or d.photo_id.duplicated().sum()!=26905:
        raise RuntimeError('photo source denominator drift')
    if d[d.photo_id.duplicated(keep=False)].groupby('photo_id').inat_taxon_id.nunique().max()!=1:
        raise RuntimeError('different source taxon on identical photo')
    d=d.drop_duplicates('photo_id').copy()
    if len(d)!=100543 or d.inat_taxon_id.nunique()!=42111:raise RuntimeError('source taxon opportunity drift')
    ok=d.measurement_status.eq('classified_four_state_morph') & d.morph.isin(CLASSES)
    if ok.sum()!=45953:raise RuntimeError('source photo-colour classifier drift')
    ext=read_archive(env_zip,'fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz',
      ['inat_taxon_id','observation_id','photo_id','latitude','longitude',*ALL_FEATURES])
    if len(ext)!=100543 or ext.photo_id.nunique()!=100543:raise RuntimeError('environment source denominator drift')
    merged=d.merge(ext,on=['inat_taxon_id','observation_id','photo_id'],validate='one_to_one',indicator=True)
    if not merged._merge.eq('both').all():raise RuntimeError('photo-environment source identity mismatch')
    merged=merged.drop(columns=['_merge'])
    merged['genus']=merged.species.str.split().str[0]
    merged['abs_lat']=merged.latitude.abs()
    merged['lon_sin']=np.sin(np.deg2rad(merged.longitude))
    merged['lon_cos']=np.cos(np.deg2rad(merged.longitude))
    merged['sector_index']=np.minimum(3,np.maximum(0,np.floor((merged.longitude+180)/90))).astype('Int64')
    merged['photo_label_ok']=merged.measurement_status.eq('classified_four_state_morph') & merged.morph.isin(CLASSES)
    good_geo=merged.latitude.between(-90,90)&merged.longitude.between(-180,180)
    climate=tuple(x for g,b in BLOCKS.items() if g!='soil' for x in b)
    assert len(climate)==26
    c=merged.photo_label_ok & good_geo & merged[[*GEO,*climate]].notna().all(axis=1)
    s=c & merged[list(BLOCKS['soil'])].notna().all(axis=1)
    if (int(c.sum()),int(s.sum()),merged.loc[c,'inat_taxon_id'].nunique(),merged.loc[s,'inat_taxon_id'].nunique())!=(45836,31925,25028,20786):
        raise RuntimeError('source classifiable climate / soil population drift')
    quality=merged.flower_effective_pixels.ge(23826)&merged.nuisance_pixel_fraction.le(.1375717742195135)
    return merged,{'climate':c,'soil':s},quality


def geo_cell(lat,lon):
    # Original source's reproducible coarse approx equal-area cell index,
    # a spatial sensitivity only, NOT true fine-scale site matched null.
    y=np.clip(np.floor((np.sin(np.deg2rad(lat))+1)*4.5),0,8).astype(int)
    x=np.clip(np.floor((lon+180)/20),0,17).astype(int)
    return 18*y+x


def species_frame(d:pd.DataFrame,features:tuple[str,...]):
    src=d.copy()
    for h in CLASSES:src['hue_'+h]=src.morph.eq(h).astype(float)
    groupcols=[*GEO,*features,*['hue_'+h for h in CLASSES], 'latitude','longitude']
    agg={c:'mean' for c in dict.fromkeys(groupcols)}
    agg.update({'genus':'first','photo_id':'size'})
    s=src.groupby('inat_taxon_id',sort=True).agg(agg).rename(columns={'photo_id':'n_photos'}).reset_index()
    if not s.inat_taxon_id.is_unique:raise RuntimeError('source species repeated')
    s['cell']=geo_cell(s.latitude.to_numpy(float),s.longitude.to_numpy(float))
    return s


def cluster_suff(species:pd.DataFrame,anchors,means,scales,geos=True,cell_conditioned=False):
    xnames=(*GEO,*anchors) if geos else anchors
    src=species.copy()
    if cell_conditioned:
        src['group']=src.genus.astype(str)+'__'+src.cell.astype(str)
    else:src['group']=src.genus.astype(str)
    n=src.group.value_counts()
    src=src[src.group.isin(n[n>=3].index)].copy()
    if len(src)<120 or src.genus.nunique()<30:
        return None,{'n_species':len(src),'n_genera':int(src.genus.nunique()),'n_groups':int(src.group.nunique())}
    xx=(src[list(xnames)].to_numpy(float)-means[list(xnames)].to_numpy(float))/scales[list(xnames)].to_numpy(float)
    yy=src[['hue_'+h for h in CLASSES]].to_numpy(float)
    names,gi=np.unique(src.group.to_numpy(str),return_inverse=True)
    weights=src.groupby('group').size().loc[names].to_numpy(int)
    xbar=np.stack([np.bincount(gi,weights=xx[:,i])/weights for i in range(xx.shape[1])],axis=1)
    ybar=np.stack([np.bincount(gi,weights=yy[:,i])/weights for i in range(4)],axis=1)
    cx=xx-xbar[gi];cy=yy-ybar[gi]
    groups={}
    for k,name in enumerate(names):
        take=gi==k
        z=cx[take];u=cy[take]
        groups[name]=(z.T@z,z.T@u,len(z),str(src.loc[src.group.eq(name),'genus'].iloc[0]))
    return groups,{'n_species':int(len(src)),'n_genera':int(src.genus.nunique()),'n_groups':int(len(names)),
                   'n_photo_rows':int(src.n_photos.sum()),'n_geo_cells':int(src.cell.nunique()),
                   'min_qualified_species_per_cluster':3}


def regress_boot(groups,anchors,nboot,seed,equal_genus=False,cell_conditioned=False):
    """Return all signed hue×anchor coefficients with same-genus multinomial draws.

    All anchors in a scope share genus draws. The conservative maxT family is
    4 photos colours x all nominated anchors WITHIN THAT scope/weight/population;
    it does not select-correct scopes, modalities or original hypothesis.
    """
    gen=set(v[3] for v in groups.values());gnames=sorted(gen);ng=len(gnames)
    if ng<30:return {'status':'HOLD_INSUFFICIENT_NAMED_GENUS','n_genera':ng}
    arrays=sorted(groups.items())
    gr=[v[3] for _,v in arrays]
    nm=np.array([v[2] for _,v in arrays],int)
    wg=np.ones(len(arrays)) if not equal_genus else np.array([1/max(1,int(sum(nm[j] for j,g in enumerate(gr) if g==v))) for v in gr])
    # Weights within *named genus* not within genus×cell group so each genus has one vote.
    mats=np.stack([v[0] for _,v in arrays])*wg[:,None,None]
    rhs=np.stack([v[1] for _,v in arrays])*wg[:,None,None]
    ngroup=len(arrays);ngen=len(gnames);gid=np.array([gnames.index(g) for g in gr],int)
    gm=np.stack([mats[gid==i].sum(axis=0) for i in range(ngen)])
    gy=np.stack([rhs[gid==i].sum(axis=0) for i in range(ngen)])
    rng=np.random.default_rng(seed)
    draws=rng.multinomial(ngen,np.repeat(1/ngen,ngen),size=nboot)
    # Stage arrays for every anchor x 4 hues, aligned on common genus resamples.
    vals=np.empty((nboot,len(anchors),4),float)
    points=np.empty((len(anchors),4),float)
    if cell_conditioned:
        # after identical genus×spatial-cell conditioning, each anchor is a
        # one-variable slope. No imposed directional geographic regressors.
        totalM=gm.sum(axis=0);totalY=gy.sum(axis=0)
        bm=np.einsum('bg,gij->bij',draws,gm,optimize=True)
        by=np.einsum('bg,gik->bik',draws,gy,optimize=True)
        for j in range(len(anchors)):
            den=float(totalM[j,j])
            if den<1e-10:return {'status':'HOLD_NO_WITHIN_CELL_ENVIRONMENT_VARIATION','anchor':anchors[j]}
            points[j]=totalY[j]/den
            good=bm[:,j,j]>1e-10
            if not good.all():return {'status':'HOLD_DEGENERATE_CLUSTER_BOOTSTRAP','anchor':anchors[j]}
            vals[:,j,:]=by[:,j,:]/bm[:,j,j,None]
    else:
        totalM=gm.sum(axis=0);totalY=gy.sum(axis=0)
        bm=np.einsum('bg,gij->bij',draws,gm,optimize=True)
        by=np.einsum('bg,gik->bik',draws,gy,optimize=True)
        for j in range(len(anchors)):
            cols=np.array([0,1,2,3+j],int)
            mat=totalM[np.ix_(cols,cols)]
            if np.linalg.cond(mat)>1e10:return {'status':'HOLD_ILL_CONDITIONED','anchor':anchors[j]}
            points[j]=np.linalg.solve(mat,totalY[cols])[-1,:]
            sub=bm[:,cols,:][:,:,cols]
            try:v=np.linalg.solve(sub,by[:,cols,:])
            except np.linalg.LinAlgError:return {'status':'HOLD_SINGULAR_BOOTSTRAP','anchor':anchors[j]}
            vals[:,j]=v[:,-1,:]
    deviations=vals-points[None,:,:]
    sd=deviations.std(axis=0,ddof=1)
    valid=sd>1e-12
    if not valid.any():return {'status':'HOLD_DEGENERATE_FAMILY_VARIANCE'}
    maxstat=np.max(np.abs(deviations[:,valid]/sd[valid]),axis=1)
    crit=float(np.quantile(maxstat,.95))
    out=[]
    for j,key in enumerate(anchors):
        for k,hue in enumerate(CLASSES):
            lo,hi=np.quantile(vals[:,j,k],[.025,.975])
            simul=[float(points[j,k]-crit*sd[j,k]),float(points[j,k]+crit*sd[j,k])]
            out.append({'anchor':key,'hue':hue,'point_pp_per_original_photo_env_sd':float(points[j,k]*100),
                'cluster_ci95_pp':[float(100*lo),float(100*hi)],
                'scope_anchor_hue_maxT_ci95_pp':[float(100*simul[0]),float(100*simul[1])],
                'maxT_excludes_zero':bool(simul[0]>0 or simul[1]<0),
                'nominal_ci_excludes_zero':bool(lo>0 or hi<0),
                'bootstrap_positive_fraction':float((vals[:,j,k]>0).mean())})
    # The 4 mutually exclusive photos class-probability coefficients sum to 0.
    if not np.allclose(points.sum(axis=1),0,atol=1e-8):raise RuntimeError('Four-hue compositional conservation drift')
    return {'status':'OK','n_named_genus':ngen,'n_genus_subgroups':ngroup,
            'n_bootstrap':nboot,'anchors':list(anchors),'n_colour_anchor_family':len(out),
            'scope_colour_anchor_family_maxT95':crit,'estimates':out,
            'maxT_discoveries':int(sum(x['maxT_excludes_zero'] for x in out))}


def conditional_full(species,features,anchors,means,scales,equal_genus=False):
    src=species.copy();n=src.genus.value_counts();src=src[src.genus.isin(n[n>=3].index)]
    if len(src)<120 or src.genus.nunique()<30:return {'status':'HOLD'}
    names=(*GEO,*features)
    xx=(src[list(names)].to_numpy(float)-means[list(names)].to_numpy(float))/scales[list(names)].to_numpy(float)
    yy=src[['hue_'+x for x in CLASSES]].to_numpy(float)
    g,idx=np.unique(src.genus.to_numpy(str),return_inverse=True)
    ct=np.bincount(idx);xm=np.stack([np.bincount(idx,weights=xx[:,j])/ct for j in range(xx.shape[1])],axis=1)
    ym=np.stack([np.bincount(idx,weights=yy[:,j])/ct for j in range(4)],axis=1)
    z=xx-xm[idx];u=yy-ym[idx]
    w=np.ones(len(src)) if not equal_genus else 1/ct[idx]
    a=z.T@(w[:,None]*z);b=z.T@(w[:,None]*u)
    # A fixed regularization fraction rather than data-driven selected ridge.
    # Same penalty-equivalent species scale as FCP ridge40 / ~16000 taxa.
    alpha=.0025*w.sum()
    beta=np.linalg.solve(a+alpha*np.eye(len(names)),b)
    return {'status':'OK','ridge_penalty_fraction_per_weight_unit':.0025,
       'target_coefficients_pp':{key:{h:float(100*beta[names.index(key),i]) for i,h in enumerate(CLASSES)} for key in anchors},
       'note':'Descriptive all26/all36 partial coefficient with correlated environmental features, not causal or inferential CI'}


def run(photo_zip,env_zip,nboot=399,sector_scopes=('world','W180_W90','W90_E0','E0_E90','E90_E180'),quality_scopes=('all','strict')):
    d,pools,technical_quality=original_source(photo_zip,env_zip)
    r={'schema':VERSION,'created_jst':'2026-10-13','status':'RETROSPECTIVE_MACRO_PHOTOCOLOUR_ABIOTIC_DIRECTION',
       'source_photo_sha256':SOURCE_PHOTO_SHA,'source_environment_sha256':SOURCE_ENV_SHA,
       'source_photo_ids':100543,'source_taxon_opportunities':42111,'classified_source_photo_ids':45953,
       'full_source_environment_features':len(ALL_FEATURES),'original_macro_opportunities':{
        'climate26_classified_photos':45836,'all36_classified_photos':31925,
        'climate26_taxa':25028,'all36_taxa':20786},
       'original_colour_classes':list(CLASSES),'original_environment_block_inventory':{k:list(v) for k,v in BLOCKS.items()},
       'frozen_primary_directional_anchors':list(ANCHOR_CLIMATE),'soil_additional_anchor':'soil_available_water_proxy',
       'source_original_photo_quality':{'flower_effective_pixels_atleast':23826,'nuisance_pixel_fraction_atmost':.1375717742195135},
       'subcohorts':{},'scientific_nonclaims':[
        'This measures photographic species-level colour allocation across environment, not within-species genetically identified morph frequencies.',
        'Genus conditioning is not direct phylogenetic branch-length adjustment.',
        'Broad coordinate and genus-by-original-cell comparisons are not a full fine-spatial random-field correction.',
        'Feature coefficients under multicollinearity and source MNAR image classifiability are not causal environment-selection effects.',
        'No gradient is promoted as global selection if sign or familywise/photo-quality/independent-botanical tests fail.',
        'No pigment cost, UV physiology, pollinator visitation, plant genotypes or relative reproductive success was measured.',
       ]}
    for mode,mask in pools.items():
        features=tuple(x for name,block in BLOCKS.items() if mode=='soil' or name!='soil' for x in block)
        anchors=ANCHOR_SOIL if mode=='soil' else ANCHOR_CLIMATE
        available=d.loc[mask].copy()
        norm_columns=(*GEO,*features)
        means=available[list(norm_columns)].mean();scales=available[list(norm_columns)].std(ddof=0)
        scales=scales.where(scales>=1e-12,1.0)
        for quality in quality_scopes:
            src=available if quality=='all' else available.loc[technical_quality.loc[available.index]].copy()
            for region in sector_scopes:
                sub=src if region=='world' else src.loc[src.sector_index.eq(SECTORS.index(region))]
                spp=species_frame(sub,features)
                key=f'{mode}__{quality}__{region}'
                info={'source_classified_photos':int(len(sub)),'n_photo_species':int(len(spp)),
                      'n_source_genera':int(spp.genus.nunique()),'source_species_photo_count_median':float(spp.n_photos.median()) if len(spp) else 0,
                      'source_normalization':'same photo source/global environmental mean and SD within climate26 vs soil36, frozen across region/quality',
                      'source_species_grain':'one named original taxon per region/photography opportunity',
                      'estimation':{},'cell_conditioned':{}}
                for weight in ['species_equal','genus_equal']:
                    g,qc=cluster_suff(spp,anchors,means,scales,geos=True)
                    if g is None:est={'status':'HOLD_INSUFFICIENT_GENUS_OR_SPECIES',**qc}
                    else:
                        est=regress_boot(g,anchors,nboot,BOOT_SEED+len(region)+len(mode)*11+(25 if quality=='strict' else 0),equal_genus=weight=='genus_equal')
                        est.update(qc)
                    if est['status']=='OK':
                        est['all_features_conditional_ridge_direction']=conditional_full(spp,features,anchors,means,scales,equal_genus=weight=='genus_equal')
                    info['estimation'][weight]=est
                    # Spatially matched genus×coarse equal-area cell, interpreted as stricter sensitivity only.
                    if quality=='all' and mode=='climate':
                        gg,meta=cluster_suff(spp,anchors,means,scales,geos=False,cell_conditioned=True)
                        if gg is None:col={'status':'HOLD_INSUFFICIENT_LOCAL_CONGENERIC_GROUPS',**meta}
                        else:
                            col=regress_boot(gg,anchors,nboot,BOOT_SEED+len(region)+900,equal_genus=weight=='genus_equal',cell_conditioned=True)
                            col.update(meta)
                        info['cell_conditioned'][weight]=col
                r['subcohorts'][key]=info
                print('COHORT',key,'species',len(spp),'genera',spp.genus.nunique(),
                      'signed',info['estimation']['species_equal']['status'],flush=True)
    return r


def flatten(result):
    rows=[]
    for key,case in result['subcohorts'].items():
        mode,quality,region=key.split('__')
        for kind in ('estimation','cell_conditioned'):
            for w,dat in case[kind].items():
                if dat['status']!='OK':continue
                full=dat.get('all_features_conditional_ridge_direction',{}).get('target_coefficients_pp',{})
                for x in dat['estimates']:
                    rows.append({'mode':mode,'quality':quality,'region':region,'adjustment':'genus_geo' if kind=='estimation' else 'genus_geocell',
                                 'weight':w,'n_species':dat.get('n_species'), 'n_genera':dat.get('n_named_genus'),
                                 'anchor':x['anchor'],'hue':x['hue'],'direction_pp_per_env_global_sd':x['point_pp_per_original_photo_env_sd'],
                                 'cluster_ci_lo_pp':x['cluster_ci95_pp'][0],'cluster_ci_hi_pp':x['cluster_ci95_pp'][1],
                                 'maxT_ci_lo_pp':x['scope_anchor_hue_maxT_ci95_pp'][0],
                                 'maxT_ci_hi_pp':x['scope_anchor_hue_maxT_ci95_pp'][1],
                                 'maxT_excludes_zero':x['maxT_excludes_zero'],
                                 'full26_36_ridge_partial_pp':full.get(x['anchor'],{}).get(x['hue'],None),
                                 'bootstrap_positive_fraction':x['bootstrap_positive_fraction']})
    return pd.DataFrame(rows)


def main():
    a=argparse.ArgumentParser()
    a.add_argument('--photo-archive',required=True,type=Path)
    a.add_argument('--environment-archive',required=True,type=Path)
    a.add_argument('--bootstrap',type=int,default=399)
    a.add_argument('--outdir',required=True,type=Path)
    args=a.parse_args()
    r=run(args.photo_archive,args.environment_archive,nboot=args.bootstrap)
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/'result.json').write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    tab=flatten(r);tab.to_csv(args.outdir/'direction_ledger.csv',index=False)
    print('DIRECTION RESULTS',json.dumps({'cohorts':len(r['subcohorts']),'tested_direction_rows':len(tab),
       'scope_MAXT_positive':int(tab.maxT_excludes_zero.sum()),
       'macro_photo_only':True},ensure_ascii=False),flush=True)

if __name__=='__main__':main()
