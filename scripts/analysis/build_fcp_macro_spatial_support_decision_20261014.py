#!/usr/bin/env python3
"""Generate strict, no-rescue ecological interpretation from already-frozen source result."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA='fcp_spatial_support_macroecological_decision_20261014_v1'

def build(jpath,lpath,vpath):
    r=json.loads(Path(jpath).read_text())
    v=json.loads(Path(vpath).read_text())
    rows=pd.read_csv(lpath)
    assert r['schema']=='fcp_macro_hue_same_species_spatial_conditioning_20261014_v1'
    assert v['status']=='PASS' and v['independent_direct_WLS_slope_reconstructions']==288
    assert len(rows)==576
    assert (rows['type'].eq('stage').sum(),rows['type'].eq('paired_difference').sum())==(288,288)
    assert r['climate26_eligible_source_photos']==45836
    result={'schema':SCHEMA,'outcome_exposure':'retrospective_same_original_FCP_photo_source',
        'n_original_photo_opportunities':100543,'n_classified_climate26_photo_observations':45836,
        'n_environmental_anchors':6,'n_photographic_hue_outcomes':4,
        'n_source_test_regions':2,'n_weightings':2,
        'n_environmental_direction_rows':576,
        'independent_direct_source_WLS_checks':288,
        'max_abs_wls_discrepancy_probability':v['max_abs_coeff_difference_probability_units'],
        'source_sha256_photo':r['photo_source_sha256'],
        'source_sha256_environment':r['abiotic_source_sha256'],
        'base_pr193_conditions_comparable':True,
        'claims':{'wide_scale_env_colour_direction':'OBSERVATIONAL_ASSOCIATION_IN_PHOTOGRAPHIC_COLOUR_COMPOSITION',
            'same_species_spatial_adjustment_difference':'NOT_FAMILYWISE_ROBUST_IN_EITHER_GENUS_EQUAL_GLOBAL_OR_EAST',
            'environment_only_photo_slope':'UNIDENTIFIABLE_FROM_OBSERVATIONAL_SOURCE',
            'within_species_colour_adaptation':'NOT_TESTED',
            'common_world_atlas_hotspot':'NOT_SUPPORTED_BY_PREVIOUS_RGFCA_P_HOT',
            'universal_climate_selection':'NOT_TESTED_NO_W_OR_CAUSAL_COLOUR_ENV_INTERVENTION',
            'genetic_colour_polymorphism_maintenance':'NOT_TESTED_NO_GENOTYPES_FREQUENCY_FITNESS'},
        'BIO12_white_by_region_and_weight':{},'negative_gates':[]}
    for region in ('world','E90_E180'):
        c=r['regions'][region]
        assert c['status']=='OK' and c['support']['only_cell_adjustment_changes_between_last_two_stages']
        j=r['frozen_anchors'].index('wc_bio12');k=r['frozen_four_hues'].index('white')
        for weight in ('species_equal','genus_equal'):
            f=c['fit'][weight]
            key=region+'|'+weight
            a=[float(f['estimate'][stage]['point'][j][k])*100 for stage in ('full_photo_in_matched_genera','same_species_with_broad_geography','same_species_with_genus_cell')]
            b=f['stage_differences']['spatial_adjustment_on_same_species']
            ci=np.array(b['paired_genus_bootstrap_CI95'][j][k])*100
            maxci=np.array(b['family_maxT_95_CI'][j][k])*100
            assert abs((a[1]-a[2])-b['point'][j][k]*100)<1e-10
            result['BIO12_white_by_region_and_weight'][key]={
                'n_matched_genera':c['support']['n_matched_named_genera'],
                'n_matched_species':c['support']['n_matched_species'],
                'matched_photographic_individuals':c['support']['n_matched_original_photos'],
                'cohort_same_matched_genera_all_species_pp':a[0],
                'same_species_with_broad_coordinate_pp':a[1],
                'same_species_within_genus_cell_pp':a[2],
                'paired_adjustment_difference_pp':float(b['point'][j][k])*100,
                'paired_nominal_95_interval_pp':ci.tolist(),
                'paired_24_hue_axis_maxT_95_interval_pp':maxci.tolist(),
                'difference_excludes_zero_after_24_maxT':bool(maxci[0]>0 or maxci[1]<0)}
            if weight=='genus_equal' and not(maxci[0]<0<maxci[1]):
                raise RuntimeError('All main genus-equal spatial adjustment differences are supposed to be HOLD under 24-axis family')
        result['negative_gates'].append({'region':region,'strictly_matched_n_species':c['support']['n_matched_species'],
            'unmeasured_real_local_soil_or_individual_fitness':True,
            'genus_equal_spatial_adjustment_CI_includes_zero':True})
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--result-json',required=True,type=Path)
    p.add_argument('--result-csv',required=True,type=Path)
    p.add_argument('--independent',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    a=p.parse_args()
    obj=build(a.result_json,a.result_csv,a.independent)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
    print('FCP M1 MATCHED-GENUS DECISION',json.dumps(obj['claims'],ensure_ascii=False))
