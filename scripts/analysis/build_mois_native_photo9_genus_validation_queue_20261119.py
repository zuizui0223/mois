#!/usr/bin/env python3
"""Prioritize genuine-photo botanical validation for currently one-sided genus contrasts.
Unverified photo ROI colour != a verified petal colour. No photos are fetched.
"""
import io,zipfile,hashlib,json,argparse
from pathlib import Path
import pandas as pd
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--source',required=True)
p.add_argument('--out',required=True)
args=p.parse_args()
source=Path(args.source)
assert hashlib.sha256(source.read_bytes()).hexdigest()=='955f4c0fdf2ba1be49f3195e065db52b048e28400b4a7938dc8c97222c4e280d'
with zipfile.ZipFile(source) as z:
    metadata=json.loads(z.read('receipt.json'))
    x=pd.read_csv(io.BytesIO(z.read('WCVP_strict_native_species_LGM_3GCM_climatic_analogue_all_eligible.csv.gz')),compression='gzip')
assert len(x)==10882 and x.source_specieskey.is_unique
x=x[x.source_hue9.ne('bronze')].copy()
x['white']=x.source_hue9.eq('white')
x['good']=x.hue_status.eq('sampling_consensus_single')
counts=x.groupby('genus').white.agg(['min','max'])
genera=counts.index[counts['min'].ne(counts['max'])]
x=x[x.genus.isin(genera)]
good=x[x.good].groupby(['genus','white']).size().unstack(fill_value=0)
good=good.reindex(index=genera,columns=[False,True],fill_value=0)
# One side has consensus; other has none. These genera require ONE additional genuine verified species
one=good.index[((good[True]>0)^(good[False]>0))]
assert len(genera)==688 and len(one)==199
rows=[]
for genus in one:
    pool=x[x.genus.eq(genus)]
    has_consensus_white=good.loc[genus,True]>0
    wanted_white=not has_consensus_white
    candidates=pool[(pool.white==wanted_white)&(~pool.good)].copy()
    assert len(candidates)>0
    # Existing source confidence / ecological support, not LGM response or apparent effect size.
    candidates['quality_rank']=candidates.hue_status.map({
        'provisional_single':0,'sampled_multicolour':1,'multicolour_candidate':2,
    }).fillna(99).astype(int)
    candidates['native_rich']=candidates.n_modern_native_occupied_cells.ge(100).astype(int)
    candidates['lift_good']=candidates.modern_heldout_native_to_region_background_lift.ge(1.2).astype(int)
    candidates=candidates.sort_values(['quality_rank','native_rich','lift_good','n_modern_native_occupied_cells','source_specieskey'],ascending=[True,False,False,False,True])
    best=candidates.iloc[0]
    prior_good=good.loc[genus,not wanted_white]
    rows.append({'genus':genus,
        'missing_confirmed_class':'white' if wanted_white else 'nonwhite',
        'candidate_source_specieskey':best.source_specieskey,
        'candidate_species':best.species,
        'source_photo9_label':best.source_hue9,
        'source_hue_status':best.hue_status,
        'n_native_current_100km_cells':int(best.n_modern_native_occupied_cells),
        'n_native_500km_blocks':int(best.n_modern_native_500km_blocks),
        'n_current_TDWG3_native_regions':int(best.n_WCVP_TDWG3_native_regions),
        'existing_opposite_photo_consensus_species_in_genus':int(prior_good),
        'n_alternative_unverified_species_in_needed_class':len(candidates),
        'modern_heldout_native_to_region_background_lift':float(best.modern_heldout_native_to_region_background_lift),
        'is_provisional_single_photo_status':bool(best.hue_status=='provisional_single'),
        'not_confirmed_petal_colour':True,
        'desired_action':'Independent botanical petal ROI validation from original photos, 3 observations per species minimum; record organ, polymorphism and observation IDs',
    })
result=pd.DataFrame(rows)
result['provisional_priority']=result.is_provisional_single_photo_status.astype(int)
result['native100_priority']=result.n_native_current_100km_cells.ge(100).astype(int)
result=result.sort_values(['provisional_priority','native100_priority','n_native_current_100km_cells','genus'],ascending=[False,False,False,True]).reset_index(drop=True)
result.insert(0,'rank',np.arange(1,len(result)+1))
result.insert(1,'batch',np.where(result['rank']<=75,'NEXT_75_GENERA','LATER'))
out=Path(args.out)
assert not out.exists(),'Do not overwrite previous source-backed queue'
out.mkdir(parents=True)
a=out/'prioritized_missing_genus_colour_photo9_QA_targets.csv'
result.to_csv(a,index=False)
summary={
    'source_zip_SHA256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'source_accepted_photo9_species':10882,
    'excluded_bronze_from_primary_7_nonwhite_hues':True,
    'genera_white_and_nonwhite_photo9_nonbronze':len(genera),
    'genera_both_sides_consensus_photo':int(((good[True]>0)&(good[False]>0)).sum()),
    'genera_one_side_consensus_only':len(one),
    'genera_neither_side_consensus_photo':int(((good[True]==0)&(good[False]==0)).sum()),
    'priority_batch_unique_genera':75,
    'priority_species_provisional_single':int(result[result['rank']<=75].is_provisional_single_photo_status.sum()),
    'priority_hue_counts':result[result['rank']<=75].source_photo9_label.value_counts().to_dict(),
    'priority_missing_class':result[result['rank']<=75].missing_confirmed_class.value_counts().to_dict(),
    'ranking_has_not_used_any_LGM_response':True,
    'no_original_images_fetched_or_reannotated':True,
    'requires_3_independent_original_observations_and_manual_petal_organ_QA':True,
    'estimated_additional_comparable_genera_if_every_selected_species_confirmed':75,
    'not_guaranteed_additional_genera_due_to_petal_colour_misclassification_or_polymorphism':True,
    'not_independent_fossil_refugia_evidence':True,
    'csv_sha256':hashlib.sha256(a.read_bytes()).hexdigest(),
}
(out/'receipt.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print('\nTOP10:')
print(result.head(10)[['rank','genus','missing_confirmed_class','candidate_species','source_photo9_label','source_hue_status']].to_string(index=False))