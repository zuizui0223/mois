#!/usr/bin/env python3
"""Incrementally match only newly hue-labelled MOIS taxa to GBIF COL XR.

Source: signed 42,111 source taxon atlas after third genuine photo wave,
and completed all-hue GBIF COL XR matches for the PREVIOUS atlas.
Never requery or re-count already matched species. Preserve synonyms,
non-exact matches, missing petals and source photographic uncertainty.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_gbif_colxr_batch_match import (
    match_batch, parse_row, source_sha
)
from scripts.analysis.mois_gbif_colxr_nine_palette_name_pilot import (
    PALETTE, NAME_RULE, COL_XR, require
)
from scripts.analysis.mois_gbif_colxr_nine_colour_traits_bridge import build_traits
from scripts.analysis.mois_gbif_colxr_sql_isea3h_requests import make_requests

SCHEMA = "mois_incremental_actual_wave3_photo9_gbif_colxr_match_v1"
ATLAS_AFTER2304_SHA = "4369818640dd557530274077610bd91637e8e86f3a6094dbeb12b5fd857aa697"
EXPECTED_SOURCE_TAXA = 42111
PRIOR_NAMED = 28309
PRIOR_EXACT = 26852
PRIOR_SOURCE_HUES = 28312
UPDATED_SOURCE_HUES = 28591

def exact_name_source(atlas):
    require({"inat_taxon_id","species","hue_status","dominant_hue9"} <= set(atlas.columns),
            "Updated real photo hue atlas incomplete")
    d = atlas.copy()
    for c in ("inat_taxon_id","species","hue_status","dominant_hue9"):
        d[c] = d[c].fillna("").astype(str).str.strip()
    require(d.inat_taxon_id.is_unique and len(d)==EXPECTED_SOURCE_TAXA,
            "Updated immutable 42,111 source taxa corrupted or duplicated")
    require(int(d.dominant_hue9.isin(PALETTE).sum()) == UPDATED_SOURCE_HUES,
            "Authentic new 279 photo hue labels absent")
    return d.loc[d.dominant_hue9.isin(PALETTE) &
                 d.species.str.match(NAME_RULE)].copy()

def select_new_named(updated, prior):
    """Conserve old identity mappings, then match ONLY taxon IDs absent before."""
    require({"inat_taxon_id","species","source_hue9","hue_status","status","key",
             "COL_XR_checklistKey"} <= set(prior.columns),
            "Frozen full GBIF prior matching ledger lacks required evidence columns")
    require(len(prior)==PRIOR_NAMED and prior.inat_taxon_id.is_unique and
            int(prior.status.eq("colxr_exact_accepted_species").sum()) == PRIOR_EXACT,
            "Earlier source-complete GBIF naming receipt drifted")
    require(prior.COL_XR_checklistKey.eq(COL_XR).all(),
            "GBIF reference checklist dataset changed")
    accepted = prior.loc[prior.status.eq("colxr_exact_accepted_species")]
    require(accepted.key.ne("").all() and accepted.key.is_unique,
            "Earlier accepted COL XR reference key collision")
    known = updated.set_index("inat_taxon_id").loc[prior.inat_taxon_id]
    require((known.species.to_numpy() == prior.species.to_numpy()).all() and
            (known.dominant_hue9.to_numpy() == prior.source_hue9.to_numpy()).all() and
            (known.hue_status.to_numpy() == prior.hue_status.to_numpy()).all(),
            "Known name/hue states changed between signed photo atlases")
    new = updated.loc[~updated.inat_taxon_id.isin(prior.inat_taxon_id)].copy()
    require(len(new)<=UPDATED_SOURCE_HUES-PRIOR_SOURCE_HUES,
            "New named-photo hue source exceeds true incremental source opportunity")
    return new.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)

def combine_names(prior,new_results):
    allnames = pd.concat([prior,new_results],ignore_index=True,sort=False)
    require(allnames.inat_taxon_id.is_unique,"Source taxon names repeated on incremental append")
    accepted = allnames.loc[allnames.status.eq("colxr_exact_accepted_species")].copy()
    require(accepted.key.ne("").all(),"Accepted GBIF key missing from actual result")
    collision = accepted.loc[accepted.key.duplicated(keep=False),"key"].unique()
    if len(collision):
        allnames.loc[allnames.key.isin(collision) &
                     allnames.status.eq("colxr_exact_accepted_species"),
                     "status"]="duplicate_accepted_species_key_review"
    return allnames, allnames.loc[allnames.status.eq("colxr_exact_accepted_species")].copy()

def classify_new_rows(added,*,batch_size=200,pause=1.0):
    require(1<=batch_size<=250 and 0.5<=pause<=10,"Unsupported GBIF request size or pace")
    frames=[]
    for index in range(0,len(added),batch_size):
        batch=added.iloc[index:index+batch_size]
        try:
            actual=match_batch(batch.species.tolist())
            rows=[parse_row(row,result) for row,result in zip(batch.itertuples(index=False),actual)]
        except Exception as err:
            rows=[{"inat_taxon_id":row.inat_taxon_id,
                   "species":row.species,
                   "source_hue9":row.dominant_hue9,
                   "hue_status":row.hue_status,
                   "COL_XR_checklistKey":COL_XR,
                   "key":"",
                   "status":"api_error_review",
                   "error_note":str(err)[:400]}
                  for row in batch.itertuples(index=False)]
        frames.append(pd.DataFrame(rows))
        print("MOIS_GBIF_NEW_ONLY_BATCH",index+len(batch),len(added),
              "accepted",sum(x["status"]=="colxr_exact_accepted_species" for x in rows),
              flush=True)
        time.sleep(pause)
    if frames:return pd.concat(frames,ignore_index=True,sort=False)
    return pd.DataFrame(columns=["inat_taxon_id","species","source_hue9",
                                 "hue_status","COL_XR_checklistKey","status","key"])

def run(atlas_path,prior_match_path,output,*,batch_size=200,pause=1.0):
    output=Path(output)
    require(not output.exists(),"Cannot overwrite previous frozen incremental GBIF reconciliation")
    require(source_sha(atlas_path)==ATLAS_AFTER2304_SHA,
            "Source-authenticated 42,111 photo atlas after third wave SHA drifted")
    atlas=pd.read_csv(atlas_path,compression="infer",dtype={"inat_taxon_id":str})
    earlier=pd.read_csv(prior_match_path,compression="infer",dtype=str,keep_default_na=False)
    updated=exact_name_source(atlas)
    new=select_new_named(updated,earlier)
    names_new=classify_new_rows(new,batch_size=batch_size,pause=pause)
    full,accepted=combine_names(earlier,names_new)
    # Validated nine original fractions stay joined to the same source taxon and
    # current accepted COL XR key. Old 26,852 names require no new API calls.
    traits=build_traits(atlas,accepted,verify_counts=False)
    require(len(traits)==len(accepted) and len(traits)>=PRIOR_EXACT,
            "Original numeric nine-fraction species bridge lost known accepted taxa")
    output.mkdir(parents=True)
    paths={}
    output_tables=(
        ("mois_incremental_new_colour_taxa_colxr_name_results.csv.gz",names_new),
        ("mois_all_current_named_hue_taxa_colxr_matches.csv.gz",full),
        ("mois_current_exact_accepted_colxr_species.csv.gz",accepted),
        ("mois_current_colxr_species_real9_fraction_traits.csv.gz",traits),
    )
    for filename,table in output_tables:
        path=output/filename
        table.to_csv(path,index=False,compression="gzip",float_format="%.7f")
        paths[filename]=source_sha(path)
    # This creates requests for external authenticated GBIF SQL downloads but
    # does NOT submit them, as that requires the user's GBIF credentials.
    requests=make_requests(full,chunk_size=250)
    sql_dir=output/"gbif_SQL_unsubmitted"
    sql_dir.mkdir()
    for i,request in enumerate(requests):
        path=sql_dir/f"colxr_isea3h_species_wave3_{i:03d}.json"
        path.write_text(json.dumps(request["GBIF_download_request"],indent=2)+"\n")
        paths[str(path.relative_to(output))]=source_sha(path)
    note={"schema":SCHEMA,
          "source_original_42111_taxa":EXPECTED_SOURCE_TAXA,
          "source_photo_dominant_hue_atlas_28591":UPDATED_SOURCE_HUES,
          "previous_named_hue_taxa":len(earlier),
          "new_named_hue_taxa_matched":len(new),
          "new_exact_accepted_colxr_species":int(
              names_new.status.eq("colxr_exact_accepted_species").sum()),
          "cumulative_named_hue_taxa":len(full),
          "cumulative_exact_accepted_colxr_species":len(accepted),
          "cumulative_preserved_nine_fraction_traits":len(traits),
          "source_names_without_automatic_accepted_match":len(full)-len(accepted),
          "new_name_match_status_counts":{
              str(k):int(v) for k,v in names_new.status.value_counts().items()},
          "GBIF_SQL_ISEA3H_request_files_prepared":len(requests),
          "GBIF_authenticated_download_submitted":False,
          "GBIF_occurrence_records_retrieved":0,
          "species_wide_hue_fixation_proven":False,
          "original_roi9_petal_truth_verified":False,
          "previous_GBIF_matches_byte_sha256":source_sha(prior_match_path),
          "updated_species_photo_atlas_byte_sha256":source_sha(atlas_path),
          "outputs_SHA256":paths}
    (output/"receipt.json").write_text(json.dumps(note,sort_keys=True,indent=2)+"\n")
    print("MOIS_GBIF_INCREMENTAL_NEW_9HUE_COLXR_REAL",json.dumps({
        k:note[k] for k in ("new_named_hue_taxa_matched",
          "new_exact_accepted_colxr_species","cumulative_exact_accepted_colxr_species",
          "GBIF_SQL_ISEA3H_request_files_prepared")},sort_keys=True),flush=True)
    return note


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--atlas",required=True)
    p.add_argument("--prior-matches",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--batch-size",type=int,default=200)
    p.add_argument("--pause",type=float,default=1.0)
    z=p.parse_args()
    run(z.atlas,z.prior_matches,z.output,batch_size=z.batch_size,pause=z.pause)


if __name__=="__main__":
    main()
