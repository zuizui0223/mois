#!/usr/bin/env python3
"""Independent FCP raw-photo reader and direct row WLS for matched spatial support.

Does NOT import the primary matched-support analysis (or its sufficient-stat
bootstrap) and does NOT independently re-bootstrap confidence limits.
"""
from __future__ import annotations
import argparse,hashlib,io,zipfile,json
from pathlib import Path
import numpy as np,pandas as pd
PHOTO_SHA='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
ENV_SHA='6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca'
ANCHORS=('wc_bio5','wc_bio12','wc_bio14','wc_bio15','wc_srad_annual_mean','wc_vapr_annual_mean')
COLOURS=('white','yellow_orange','red_pink','blue_purple')
GEO=('abs_lat','lon_sin','lon_cos')
ALLFEATURES=['wc_elevation_m',*[f'wc_bio{j}' for j in range(1,20)],'wc_srad_annual_mean','wc_srad_monthly_cv','wc_wind_annual_mean','wc_wind_monthly_cv','wc_vapr_annual_mean','wc_vapr_monthly_cv']


def read(path,member,cols):
    with zipfile.ZipFile(path) as z:
        return pd.read_csv(io.BytesIO(z.read(member)),compression='gzip',usecols=cols,low_memory=False)


def source(photo,environment):
    if hashlib.sha256(Path(photo).read_bytes()).hexdigest()!=PHOTO_SHA:raise ValueError('photo SHA mismatch')
    if hashlib.sha256(Path(environment).read_bytes()).hexdigest()!=ENV_SHA:raise ValueError('env SHA mismatch')
    ids=['photo_id','inat_taxon_id','observation_id']
    cc=[*ids,'species','morph','measurement_status']
    a=read(photo,'fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz',cc)
    b=read(photo,'fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz',cc)
    assert len(a)==42111 and len(b)==85337
    t=pd.concat([a,b],ignore_index=True).drop_duplicates('photo_id')
    assert len(t)==100543
    ee=read(environment,'fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz',[*ids,'latitude','longitude',*ALLFEATURES])
    assert len(ee)==100543
    t=t.merge(ee,on=ids,validate='one_to_one')
    valid=t.measurement_status.eq('classified_four_state_morph') & t.morph.isin(COLOURS)
    valid&=t[['latitude','longitude',*ALLFEATURES]].notna().all(axis=1)
    valid&=t.latitude.between(-90,90) & t.longitude.between(-180,180)
    t=t.loc[valid].copy()
    assert len(t)==45836 and t.inat_taxon_id.nunique()==25028
    t['abs_lat']=t.latitude.abs();t['lon_sin']=np.sin(np.deg2rad(t.longitude));t['lon_cos']=np.cos(np.deg2rad(t.longitude))
    t['genus']=t.species.astype(str).str.split().str[0]
    norm=t[[*GEO,*ALLFEATURES]]
    mu=norm.mean();sd=norm.std(ddof=0).where(lambda a:a>=1e-12,1.0)
    return t,mu,sd


def one_species_frame(rows):
    d=rows.copy()
    for colour in COLOURS:d[colour]=d.morph.eq(colour).astype(float)
    vars=[*GEO,*ALLFEATURES,*COLOURS,'latitude','longitude']
    agg={x:'mean' for x in vars};agg['genus']='first'
    s=d.groupby('inat_taxon_id',sort=True).agg(agg).reset_index()
    sy=np.sin(np.deg2rad(s.latitude.to_numpy(float)))
    iy=np.clip(np.floor((sy+1)*4.5),0,8).astype(int)
    ix=np.clip(np.floor((s.longitude.to_numpy(float)+180)/20),0,17).astype(int)
    s['cell']=iy*18+ix
    return s


def stage_species(s):
    group=s.genus.astype(str)+'|'+s.cell.astype(str)
    nc=group.value_counts()
    b=s.loc[group.map(nc).ge(3)].copy()
    g=set(b.genus)
    a=s[s.genus.isin(g)].copy()
    assert a.genus.nunique()==b.genus.nunique() and b.inat_taxon_id.isin(a.inat_taxon_id).all()
    return [a,b,b]


def wls(src:pd.DataFrame,mu:pd.Series,sd:pd.Series,axis:str,colour:str,weight:str,stage:int):
    df=src.copy()
    cols=[axis] if stage==2 else [*GEO,axis]
    label=df.genus if stage!=2 else df.genus+'__'+df.cell.astype(str)
    X=(df[cols]-mu[cols])/sd[cols]
    Y=df[colour]
    xx=X-X.groupby(label).transform('mean')
    yy=Y-Y.groupby(label).transform('mean')
    gn=df.groupby('genus').genus.transform('size')
    ww=np.ones(len(df)) if weight=='species_equal' else 1.0/gn.to_numpy(float)
    w=np.sqrt(ww)
    fit=np.linalg.lstsq(xx.to_numpy(float)*w[:,None],yy.to_numpy(float)*w,rcond=None)
    if fit[2]<len(cols):raise ValueError('rank deficient source WLS')
    return float(fit[0][-1])


def run(photo,env,receipt):
    z=json.loads(Path(receipt).read_text())
    assert z['schema']=='fcp_macro_hue_same_species_spatial_conditioning_20261014_v1'
    t,mu,sd=source(photo,env)
    N=0;maxdiff=0;skips=[]
    for region in ('world','E90_E180'):
        subset=t if region=='world' else t[t.longitude.between(90.,180.)]
        s=one_species_frame(subset)
        A,B,C=stage_species(s)
        gold=z['regions'][region]
        assert len(B)==gold['support']['n_matched_species']
        assert len(A)==gold['support']['n_original_species_all_genera']
        for weight in ('species_equal','genus_equal'):
            entries=gold['fit'][weight]['estimate']
            for stage,key in enumerate(('full_photo_in_matched_genera','same_species_with_broad_geography','same_species_with_genus_cell')):
                for j,axis in enumerate(ANCHORS):
                    for k,col in enumerate(COLOURS):
                        new=wls([A,B,C][stage],mu,sd,axis,col,weight,stage)
                        old=entries[key]['point'][j][k]
                        delta=abs(new-old)
                        if delta>maxdiff:maxdiff=delta
                        if delta>1e-10:raise AssertionError((region,weight,stage,axis,col,new,old,delta))
                        N+=1
            dif=gold['fit'][weight]['stage_differences']['spatial_adjustment_on_same_species']['point']
            for j in range(len(ANCHORS)):
                for k in range(len(COLOURS)):
                    a=entries['same_species_with_broad_geography']['point'][j][k]
                    b=entries['same_species_with_genus_cell']['point'][j][k]
                    assert abs(a-b-dif[j][k])<1e-12
    return {'schema':'fcp_spatial_support_independent_original_row_wls_v1','status':'PASS',
            'source_photo_sha256':PHOTO_SHA,'source_abiotic_sha256':ENV_SHA,
            'independent_direct_WLS_slope_reconstructions':N,'max_abs_coeff_difference_probability_units':maxdiff,
            'paired_bootstrap_CI_resampled_independently':False,
            'matched_B_C_species_equal_and_rows_unchanged':True,
            'interpretation':'Purely computational replication, not independent biological data'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--photo-archive',type=Path,required=True)
    p.add_argument('--environment-archive',type=Path,required=True)
    p.add_argument('--result',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    x=run(a.photo_archive,a.environment_archive,a.result)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2)+'\n')
    print('INDEPENDENT_SOURCE_WLS',json.dumps(x))
if __name__=='__main__':main()
