#!/usr/bin/env python3
"""Within-genus white/photo moisture slopes across genus environmental contexts.

Owner-only original FCP photo/environment ZIPs; outputs aggregate numerical
results only. These retrospective slopes are NOT plant fitness, pigment cost,
causal colour adaptation, or proof of invariance across all environments.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

PHOTO_SHA = '2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
ENV_SHA = '6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca'
CLIMATE = [*(f'wc_bio{i}' for i in range(1,20)), 'wc_srad_annual_mean',
           'wc_srad_monthly_cv', 'wc_wind_annual_mean', 'wc_wind_monthly_cv',
           'wc_vapr_annual_mean', 'wc_vapr_monthly_cv', 'wc_elevation_m']
GEO = ('abs_lat', 'lon_sin', 'lon_cos')
AXES = ('wc_bio12', 'wc_vapr_annual_mean')
MODERATORS = ('wc_srad_annual_mean', 'wc_bio5', 'wc_bio12')
BOOTSTRAPS = 1999
SEED = 2026101608


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(1 << 20), b''):
            h.update(part)
    return h.hexdigest()


def read(zpath: Path, member: str, cols: list[str]) -> pd.DataFrame:
    with zipfile.ZipFile(zpath) as zipped:
        with zipped.open(member) as content:
            return pd.read_csv(content, compression='gzip', usecols=cols, low_memory=False)


def load(photo_zip: Path, env_zip: Path):
    for zip_path, expected in ((photo_zip, PHOTO_SHA), (env_zip, ENV_SHA)):
        if digest(zip_path) != expected:
            raise ValueError('Source archive SHA mismatch')
        with zipfile.ZipFile(zip_path) as zipped:
            if zipped.testzip() is not None:
                raise ValueError('Original archive ZIP member integrity failed')
    f = ['inat_taxon_id', 'observation_id', 'photo_id', 'species', 'morph', 'measurement_status']
    breadth = read(photo_zip, 'fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz', f)
    cell = read(photo_zip, 'fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz', f)
    raw = pd.concat((breadth, cell), ignore_index=True)
    if (len(breadth), len(cell), int(raw.photo_id.duplicated().sum())) != (42111, 85337, 26905):
        raise AssertionError('Original photo opportunity drift')
    if raw.groupby('photo_id').inat_taxon_id.nunique().max() != 1:
        raise AssertionError('Original photo ID reused for multiple taxa')
    raw = raw.drop_duplicates('photo_id')
    original = read(env_zip,
        'fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz',
        ['inat_taxon_id','observation_id','photo_id','latitude','longitude',*CLIMATE])
    df = raw.merge(original, on=['inat_taxon_id','observation_id','photo_id'], validate='one_to_one')
    if len(df) != 100543:
        raise AssertionError('Original all-photo cohort drift')
    eligibility = (
        df.measurement_status.eq('classified_four_state_morph') &
        df.morph.isin(('white','yellow_orange','red_pink','blue_purple')) &
        df.latitude.between(-90, 90) & df.longitude.between(-180,180) &
        df[CLIMATE].notna().all(axis=1)
    )
    df = df.loc[eligibility].copy()
    if len(df) != 45836:
        raise AssertionError('Original climate26 complete source drift')
    df['genus'] = df.species.str.split().str[0]
    df['abs_lat'] = df.latitude.abs()
    df['lon_sin'] = np.sin(np.deg2rad(df.longitude))
    df['lon_cos'] = np.cos(np.deg2rad(df.longitude))
    df['white'] = df.morph.eq('white').astype(float)
    features = list(dict.fromkeys([*GEO,*AXES,*MODERATORS]))
    norm = {k: {'mean':float(df[k].mean()), 'sd':float(df[k].std(ddof=0))} for k in features}
    fields = {k:'mean' for k in features}
    fields.update({'genus':'first','white':'mean','photo_id':'size'})
    species = (df.groupby('inat_taxon_id',sort=True).agg(fields)
               .rename(columns={'photo_id':'n_photos'}).reset_index())
    counts = species.genus.value_counts()
    species = species.loc[species.genus.isin(counts.index[counts>=3])].reset_index(drop=True)
    if (len(species),species.genus.nunique(),int(species.n_photos.sum())) != (20694,1934,38498):
        raise AssertionError('Qualified world genus-centred source support drift')
    return species,norm


def moments(species: pd.DataFrame, norm: dict, axis: str):
    genera, gi = np.unique(species.genus.to_numpy(), return_inverse=True)
    n = np.bincount(gi).astype(float)
    w = 1/n[gi]
    names = (*GEO, axis)
    xx = np.stack([(species[k].to_numpy(float)-norm[k]['mean'])/norm[k]['sd'] for k in names],axis=1)
    y = species.white.to_numpy(float)
    mean_x = np.stack([np.bincount(gi,weights=xx[:,j])/n for j in range(xx.shape[1])],axis=1)
    mean_y = np.bincount(gi,weights=y)/n
    centered_x = xx - mean_x[gi]
    centered_y = y - mean_y[gi]
    aa = np.zeros((len(genera),4,4))
    bb = np.zeros((len(genera),4))
    np.add.at(aa,gi,centered_x[:,:,None]*centered_x[:,None,:]*w[:,None,None])
    np.add.at(bb,gi,centered_x*centered_y[:,None]*w[:,None])
    return genera,aa,bb


def fit(aa: np.ndarray,bb: np.ndarray, seed: int, nboot: int):
    n = len(aa)
    point = float(100*np.linalg.solve(aa.sum(axis=0),bb.sum(axis=0))[-1])
    rng=np.random.default_rng(seed)
    draws=rng.multinomial(n,np.ones(n)/n,size=nboot)
    draw_a=np.einsum('bg,gij->bij',draws,aa,optimize=True)
    draw_b=np.einsum('bg,gi->bi',draws,bb,optimize=True)
    sampled=100*np.linalg.solve(draw_a,draw_b[...,None])[:,-1,0]
    return point,sampled


def get_moderation(species: pd.DataFrame,norm: dict, moderator: str, seed: int,nboot: int):
    genus_mean=species.groupby('genus',sort=True)[moderator].mean()
    median=float(genus_mean.median())
    high=(genus_mean.to_numpy()>=median)
    if high.sum()!=967:
        raise AssertionError('Nonbalanced environmental genus split')
    report={'median_genus_mean_in_archived_source_unit':round(median,6),
            'n_low_mean_genera':int((~high).sum()),'n_high_mean_genera':int(high.sum()),
            'axes':{}}
    for axis in AXES:
        genus,aa,bb=moments(species,norm,axis)
        if not np.array_equal(genus,genus_mean.index.to_numpy()):
            raise AssertionError('Genus group order drift')
        lower,lower_boot=fit(aa[~high],bb[~high],seed+50,nboot)
        upper,upper_boot=fit(aa[high],bb[high],seed+150,nboot)
        delta=upper_boot-lower_boot
        report['axes'][axis]={
            'low_genus_pp_per_original_photo_environment_sd':round(lower,6),
            'low_percentile_genus_bootstrap95_pp':np.round(np.quantile(lower_boot,[.025,.975]),6).tolist(),
            'high_genus_pp_per_original_photo_environment_sd':round(upper,6),
            'high_percentile_genus_bootstrap95_pp':np.round(np.quantile(upper_boot,[.025,.975]),6).tolist(),
            'high_minus_low_pp':round(float(upper-lower),6),
            'high_minus_low_percentile_bootstrap95_pp':np.round(np.quantile(delta,[.025,.975]),6).tolist(),
            'bootstrap_fraction_high_greater_than_low':round(float((delta>0).mean()),6),
        }
    return report


def build(photo_zip: Path,env_zip: Path, nboot: int = BOOTSTRAPS):
    species,norm=load(photo_zip,env_zip)
    for axis,original in (('wc_bio12',2.254545235123587),
                          ('wc_vapr_annual_mean',3.951850325152891)):
        genus,aa,bb=moments(species,norm,axis)
        observed=100*np.linalg.solve(aa.sum(axis=0),bb.sum(axis=0))[-1]
        if abs(observed-original)>1e-9:
            raise AssertionError('Did not reproduce original genus-equal world source coefficient')
    per={}
    for index,key in enumerate(MODERATORS):
        per[key]=get_moderation(species,norm,key,SEED+index*10,nboot)
    return {
        'schema':'mois_moisture_white_genus_climate_regime_v1',
        'analysis_date':'2026-10-15',
        'source_photo_sha256':PHOTO_SHA,'source_environment_sha256':ENV_SHA,
        'n_original_photo_opportunities':100543,'n_original_climate26_classifiable_photos':45836,
        'n_primary_qualified_species':20694,'n_primary_qualified_genera':1934,
        'n_primary_qualified_classifiable_photos':38498,
        'source_original_photo_moisture_sd':{k:round(norm[k]['sd'],9) for k in AXES},
        'source_original_genus_equal_within_coefficients_pp':{
            'wc_bio12':2.254545235,'wc_vapr_annual_mean':3.951850325},
        'frozen_original_species_and_genus_support':True,
        'source_photo_colours_relabelled':False,
        'bootstrap_seed_base':SEED,'bootstrap_replicates':nboot,
        'moderator_definitions':{
            'wc_srad_annual_mean':'annual mean surface solar radiation (archived source units)',
            'wc_bio5':'warmest month maximum temperature (archived source units)',
            'wc_bio12':'annual precipitation (archived source units)',
        },
        'design':'Split the 1934 source-eligible genera into equal halves at the median of equally weighted eligible species environmental averages per genus. For each half, fit the original genus-centred photographic-white/moisture slope with genus-equal species weighting, globally fixed original photo-site moisture/geography SD and broad within-genus geography covariates. Draw 1999 multinomial genus bootstrap replicates independently within halves and compare high minus low. All three moderator splits are retrospective source-exposed, not prospective ecological treatments.',
        'results_by_genus_mean_environment':per,
        'limits':['Within each split, 95% percentile intervals are ordinary genus bootstraps, not original 24-hue-axis maxT nor study-wide post-selection adjustment.',
                  'The observed environmental-median cutpoints are fixed during bootstrap; all environmental regimes are nonrandom samples of genera.',
                  'Group coefficients concern different sampled taxa and different within-group covariate supports; no causal moderator or relative plant fitness is estimated.',
                  'The dataset lacks calibrated floral/UV pigment, photosynthetic energy budgets, pollen or seed success, genotypes and frequency-dependent fitness.',
                  'No original photo IDs, species names, genus names, geographic observations, original photo media or site coordinates are exported in this report.']
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--photo-zip',type=Path,required=True)
    parser.add_argument('--env-zip',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write',action='store_true')
    modes.add_argument('--check',action='store_true')
    args=parser.parse_args()
    out=build(args.photo_zip,args.env_zip)
    if args.write:
        args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print('WROTE aggregate only climate-regime moderation result')
    else:
        frozen=json.loads(args.output.read_text(encoding='utf-8'))
        if out!=frozen:
            raise AssertionError('Original source climate-regime result differs from committed aggregate')
        print('PASS source-exact climate-regime replay and percentiles')


if __name__=='__main__':
    main()
