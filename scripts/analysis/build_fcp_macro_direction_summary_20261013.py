#!/usr/bin/env python3
"""Produce a small decision receipt from the SHA-locked 20-cohort result.
No outcome-driven candidate selection: only frozen six climate anchors, all modes.
"""
import argparse,json
from pathlib import Path
import pandas as pd

ORDER=[('climate','all','world','genus_geo'),('climate','strict','world','genus_geo'),
       ('climate','all','world','genus_geocell'),('climate','all','E90_E180','genus_geo'),
       ('climate','all','E90_E180','genus_geocell'),('climate','strict','E90_E180','genus_geo'),
       ('soil','all','world','genus_geo'),('soil','strict','world','genus_geo'),
       ('soil','all','E90_E180','genus_geo')]
KEY_AXES=['wc_bio12','wc_bio14','wc_bio5','wc_vapr_annual_mean','wc_srad_annual_mean','soil_available_water_proxy']

def build(json_source,ledger):
    a=json.loads(Path(json_source).read_text());d=pd.read_csv(ledger)
    assert a['schema']=='fcp_macroecology_colour_environment_direction_genus_geosector_v1'
    assert len(a['subcohorts'])==20 and len(d)==1280
    assert a['source_photo_ids']==100543 and a['full_source_environment_features']==36
    key=[]
    for mode,quality,region,adjustment in ORDER:
        for axis in KEY_AXES:
            if axis=='soil_available_water_proxy' and mode!='soil':continue
            z=d[(d['mode']==mode)&(d['quality']==quality)&(d['region']==region)&
                (d['adjustment']==adjustment)&(d['weight']=='genus_equal') &
                (d['anchor']==axis)&(d['hue']=='white')]
            assert len(z)==1,(mode,quality,region,adjustment,axis)
            r=z.iloc[0]
            key.append({'mode':mode,'quality':quality,'region':region,'adjustment':adjustment,
                'axis':axis,'hue':'white','n_photo_species':int(r.n_species),'n_genera':int(r.n_genera),
                'point_pp':float(r.direction_pp_per_env_global_sd),
                'family_95_interval_pp':[float(r.maxT_ci_lo_pp),float(r.maxT_ci_hi_pp)],
                'family_interval_excludes_zero':bool(r.maxT_excludes_zero),
                'conditional_all26_or36_ridge_slope_pp':None if pd.isna(r.full26_36_ridge_partial_pp) else float(r.full26_36_ridge_partial_pp)})
    info=d.groupby(['mode','quality','adjustment','weight']).maxT_excludes_zero.agg(['sum','count'])
    families=[]
    for idx,q in info.iterrows():
        families.append({'mode':idx[0],'quality':idx[1],'adjustment':idx[2],
             'weight':idx[3],'familywise_interval_zero_exclusions':int(q['sum']),
             'comparisons':int(q['count'])})
    out={'schema':'fcp_directional_macroecology_20261013_decision_receipt_v1',
         'source_source_SHA256':{'original_photo':a['source_photo_sha256'],'all36_abiotic':a['source_environment_sha256']},
         'original_photo_opportunities':a['source_photo_ids'],
         'source_photo_classified':a['classified_source_photo_ids'],
         'environmental_anchors':a['frozen_primary_directional_anchors'],
         'n_categorical_hues':4,'n_geo_sectors':4,'n_macro_scope_cohorts':20,
         'directional_estimates_and_sensitivities':len(d),
         'scope_familywise_maxT_zero_exclusions':int(d.maxT_excludes_zero.sum()),
         'key_white_directional_observations':key,
         'maxT_summary_by_source_design':families,
         'scientific_verdict':{
             'broad_photographic_colour_environment_association':'SUPPORTED_AS_SMALL_CONTEXT_DEPENDENT_PREDICTION',
             'specific_marginal_rain_and_vapour_white_slope':'SUPPORTED_IN_UNFILTERED_PHOTO_ONLY',
             'same_direction_under_genus_geocell_space_photo_quality_and_external_source':'NOT_ESTABLISHED',
             'universal_global_abiotic_hue_direction':'NOT_SUPPORTED_BY_ALL_REQUIRED_SENSITIVITIES',
             'world_common_hotspot_map':'NOT_ESTABLISHED_PREVIOUS_FCP_RGFCA_G1_P_0_070',
             'climate_driven_universal_selection':'UNIDENTIFIABLE_NO_FITNESS_GENOTYPES_OR_CAUSAL_CONTROL',
             'genetic_maintenance_of_morphs':'UNIDENTIFIABLE_NO_POPULATION_GENOTYPE_OR_RELATIVE_FITNESS',
         },
         'independent_botanical_source_previous_pr191':{'east_white_blue_direction':'NOT_INDEPENDENTLY_CONFIRMED',
            'n_original_source_bothcolours_taxa':2476,'n_fixed_east_vs_other_matched_genera':58,
            'genus_equal_east_pp':3.56,'bootstrap_ci95_pp':[-15.34,22.60]},
         'disclosure':['All candidates are post-outcome source-exposed macroecological PHOTO descriptions. No p-values are claimed for a universal genetic climate selection mechanism.',
         'Each maxT family is restricted to anchors times four hues WITHIN a cohort/weight/scope; it does not family-correct the research programme across all regions and model choices.',
         'The strict photo and shared coarse-geocell source filters change their included taxa and environment support; absent inference is not evidence the true response is absent.',
         'The prior NHM 191 photographic-site comparison reused FCP photography locations; its null does not independently resolve genetics or habitat history.']}
    return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);p.add_argument('--ledger',type=Path,required=True);p.add_argument('--out',type=Path,required=True);x=p.parse_args()
 z=build(x.result,x.ledger);x.out.write_text(json.dumps(z,indent=2,ensure_ascii=False)+'\n')
 print(json.dumps({'cohorts':z['n_macro_scope_cohorts'],'estimates':z['directional_estimates_and_sensitivities'],
   'family_within_cohort_zero_exclusions':z['scope_familywise_maxT_zero_exclusions'],
   'claim':z['scientific_verdict']['universal_global_abiotic_hue_direction']}))
