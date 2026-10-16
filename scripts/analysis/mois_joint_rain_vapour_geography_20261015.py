#!/usr/bin/env python3
"""Re-fit rainfall and atmospheric vapour on the SAME original source species.

The two source ZIP files are SHA256-verified in the imported prior study loader.
Outputs aggregate slopes and genus bootstrap uncertainty only; no sensitive
photo IDs, taxon/genus names, coordinates or per-photo data are exported.
Exploratory covariate comparisons do not identify causal climate or fitness.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np

ORIGINAL=("abs_lat","lon_sin","lon_cos","wc_elevation_m","wc_bio12","wc_vapr_annual_mean")
MODELS={
    "geo3_rain_only": [0,1,2,4],
    "geo3_vapour_only": [0,1,2,5],
    "geo3_joint": [0,1,2,4,5],
    "geo4_rain_only": [0,1,2,3,4],
    "geo4_vapour_only": [0,1,2,3,5],
    "geo4_joint": [0,1,2,3,4,5],
}
SEED=202610151330
DRAWS=1999
EXPECTED_SHA={
    "photo":"2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0",
    "environment":"6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca",
}


def load_owner_source(photo_zip:Path,env_zip:Path):
    source_path=Path(__file__).with_name("mois_geography_vapour_ovb_20261015.py")
    if not source_path.exists():
        raise FileNotFoundError("Source-verified archived FCP original-photo loader missing: "+str(source_path))
    spec=importlib.util.spec_from_file_location("mois_prior_source_geography",source_path)
    source=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    assert tuple(source.VARS)==ORIGINAL
    assert (source.PHOTO_SHA,source.ENV_SHA)==(EXPECTED_SHA["photo"],EXPECTED_SHA["environment"])
    taxa,norms=source.source_species(photo_zip,env_zip)
    ga,gy=source.genus_moments(taxa,norms)
    assert (len(taxa),int(taxa.genus.nunique()),int(taxa.n_photo.sum()),len(ga))==(20694,1934,38498,1934)
    return norms,ga,gy


def point_and_boot(ga:np.ndarray,gy:np.ndarray,nboot:int,seed:int):
    n=len(ga)
    A=ga.sum(axis=0);B=gy.sum(axis=0)
    rng=np.random.default_rng(seed)
    samples=rng.multinomial(n,np.full(n,1/n),size=nboot)
    bm=np.einsum('bg,gij->bij',samples,ga,optimize=True)
    by=np.einsum('bg,gi->bi',samples,gy,optimize=True)
    point,boot={},{}
    for name,indices in MODELS.items():
        idx=np.array(indices,int)
        m=A[np.ix_(idx,idx)]
        p=100*np.linalg.solve(m,B[idx])
        sim=100*np.linalg.solve(bm[:,idx,:][:,:,idx],by[:,idx,None])[:,:,0]
        point[name]=p;boot[name]=sim
    return A,point,boot


def two_contrast_summary(point:np.ndarray,draws:np.ndarray):
    assert len(point)==2 and draws.shape[1]==2
    sd=np.std(draws-point[None,:],axis=0,ddof=1)
    if np.any(sd<=0):raise ValueError('Source genus-bootstrap slope uncertainty degenerate')
    z=np.max(np.abs((draws-point[None,:])/sd),axis=1)
    critical=float(np.quantile(z,.95))
    names=('first','second')
    return {
        'maxT_two_contrast_95_critical':critical,
        'results':{
            name: {
                'point_pp':float(point[j]),
                'ordinary_bootstrap95_pp':np.quantile(draws[:,j],[.025,.975]).tolist(),
                'two_contrast_simultaneous95_pp':[float(point[j]-critical*sd[j]),float(point[j]+critical*sd[j])],
                'bootstrap_fraction_positive':float((draws[:,j]>0).mean())}
            for j,name in enumerate(names)
        }
    }


def partial_correlation_after_geography(A:np.ndarray,geo:list[int]):
    moisture=[4,5]
    G=np.asarray(geo)
    W=np.asarray(moisture)
    sub=A[np.ix_(W,W)]-A[np.ix_(W,G)]@np.linalg.solve(A[np.ix_(G,G)],A[np.ix_(G,W)])
    rho=float(sub[0,1]/np.sqrt(sub[0,0]*sub[1,1]))
    if not -1<rho<1:raise ValueError('Invalid partial correlation')
    return rho


def compute(photo_zip:Path,env_zip:Path,nboot:int=DRAWS,seed:int=SEED):
    norms,ga,gy=load_owner_source(photo_zip,env_zip)
    A,p,boot=point_and_boot(ga,gy,nboot,seed)
    if not np.isclose(p['geo3_rain_only'][-1],2.254545235123587,atol=1e-10):
        raise AssertionError('Original white annual-precipitation source coefficient drift')
    if not np.isclose(p['geo3_vapour_only'][-1],3.951850325152891,atol=1e-10):
        raise AssertionError('Original white atmospheric-vapour source coefficient drift')
    result={}
    rotate=np.array([[1,1],[1,-1]],dtype=float)/np.sqrt(2)
    for geo in ('geo3','geo4'):
        rain_only=p[geo+'_rain_only'][-1]
        vap_only=p[geo+'_vapour_only'][-1]
        pt=p[geo+'_joint'][-2:]
        draws=boot[geo+'_joint'][:,-2:]
        water=two_contrast_summary(pt,draws)
        rotated=two_contrast_summary(rotate@pt,draws@rotate.T)
        group={
            'n_species':20694,'n_genera':1934,
            'nuisance_covariates':['abs_lat','lon_sin','lon_cos']+(['wc_elevation_m'] if geo=='geo4' else []),
            'single_axis_pp_per_original_photo_site_sd':{
                'annual_precipitation':float(rain_only),'atmospheric_vapour':float(vap_only)},
            'both_axes_each_adjusted_for_other':{
                'annual_precipitation':water['results']['first'],
                'atmospheric_vapour':water['results']['second'],
                'two_water_axis_maxT95_critical':water['maxT_two_contrast_95_critical']},
            'rotated_equal_standardized_moisture_axes':{
                'common_equal_increase':rotated['results']['first'],
                'rain_minus_vapour_contrast':rotated['results']['second'],
                'two_rotated_axis_maxT95_critical':rotated['maxT_two_contrast_95_critical'],
                'formula': 'common=(source_z_BIO12+source_z_vapour)/sqrt(2); difference=(source_z_BIO12-source_z_vapour)/sqrt(2); both are in same full model'},
            'water_axis_partial_correlation_given_geography':partial_correlation_after_geography(A,[0,1,2]+([3] if geo=='geo4' else [])),
            'normal_equation_condition_number':float(np.linalg.cond(A[np.ix_(MODELS[geo+'_joint'],MODELS[geo+'_joint'])])),
        }
        assert np.isclose(group['rotated_equal_standardized_moisture_axes']['common_equal_increase']['point_pp'],
                          (pt[0]+pt[1])/np.sqrt(2),atol=1e-12)
        result[geo]=group
    return {
        'schema':'mois_original_20694_joint_moisture_and_geography_v1',
        'date_jst':'2026-10-15',
        'original_source_sha256':EXPECTED_SHA,
        'source_support':{'unique_original_photos':100543,'eligible_climate26_photos':45836,
                          'main_eligible_photos':38498,'main_eligible_species':20694,'main_eligible_genera':1934},
        'source_environment_sd':{v:float(norms[v]['sd']) for v in ORIGINAL},
        'bootstrap':{'draws':nboot,'seed':seed,'genus_cluster_resampling':'same multinomial genus draws across all six model specifications',
                     'interval_definition':'ordinary percentile 95 and centred-error standardized within-specification two-contrast maxT95',
                     'programme_wide_multiplicity_correction':False},
        'models':result,
        'limits':[
            'The two standardised water axes are derived from annual rainfall and mean atmospheric vapour. Their equal-direction rotated axis is an exploratory climate indicator, NOT soil moisture, water potential or a plant physiological measurement.',
            'Each 2-contrast maxT refers to that selected model and selected paired contrasts, NOT to all historic climate/colour comparisons; the original 24-contrast maxT estimates remain primary.',
            'The fitted single- and two-moisture models retain exactly the same 20694 eligible species, genus-equal weights and source SD; coefficients are observational, not biological causal components.',
            'An interval containing zero does not prove an independent contribution is absent; correlated predictors and nonrandom photographic sampling limit mechanistic identification.',
            'No original photo IDs, taxon/genus identities, site coordinates, photo media, rights/attribution or individual source rows were exported.',
        ],
        'biological_interpretation':'Original-world photographed-white composition covaries positively with a joint annual precipitation plus atmospheric vapour climate direction after genus and broad geography/elevation adjustment; unique partial attribution to either moisture axis is not stable under the two-contrast model. No pigment cost-benefit trade-off, natural selection, or evolutionary 50:50 equilibrium established.'
    }


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--photo-zip',required=True,type=Path)
    cli.add_argument('--env-zip',required=True,type=Path)
    cli.add_argument('--output',required=True,type=Path)
    action=cli.add_mutually_exclusive_group(required=True)
    action.add_argument('--write',action='store_true')
    action.add_argument('--check',action='store_true')
    args=cli.parse_args()
    report=compute(args.photo_zip,args.env_zip)
    if args.write:
        args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
        print('WROTE SHA256-replayed original-source joint-moisture coefficient audit')
    else:
        expected=json.loads(args.output.read_text(encoding='utf8'))
        if expected!=report:raise AssertionError('Original-source joint moisture receipts differ from replay')
        print('PASS original-source 20,694-species joint rainfall/vapour models and rotated moisture axes')

if __name__=='__main__':main()
