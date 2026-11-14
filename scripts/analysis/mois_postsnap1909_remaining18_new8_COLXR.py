#!/usr/bin/env python3
"""GBIF COL XR incremental exact name match of up to 8 new clear photographed hues after the second 18 independent photos.

The previous 37,095 COL XR accepted species identities are reused from a signed
photographic nine-fraction bridge. Eight new clear original photo hues are
potentially name-matchable, but species-wide hue fixation is NOT proven.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
from scripts.analysis.mois_gbif_colxr_new_hue_incremental import (
    classify_new_rows,combine_names
)
from scripts.analysis.mois_gbif_colxr_batch_match import source_sha
from scripts.analysis.mois_gbif_colxr_nine_palette_name_pilot import (
    PALETTE,NAME_RULE,COL_XR,require
)
from scripts.analysis.mois_gbif_colxr_nine_colour_traits_bridge import build_traits
from scripts.analysis.mois_gbif_colxr_sql_isea3h_requests import make_requests

SCHEMA="mois_postsnapshot_new1909_remaining18_real_original_COLXR_new8_v1"
APPEND_SCHEMA="mois_remaining1909_remaining18_original_photo_ROI9_pipeline_v2"
PRIOR_MATCH_SCHEMA="mois_postsnapshot_new1909_first64_real_original_COLXR_new17_v1"
PRIOR_ACCEPTED=37095
PRIOR_HUES=39164
AFTER_HUES=39172
TAXA=42111


def select_new(atlas,old):
    require(len(atlas)==TAXA and atlas.inat_taxon_id.is_unique,
            "42111 source taxon atlas missing or duplicate")
    for field in ("inat_taxon_id","species","dominant_hue9","hue_status"):
        atlas[field]=atlas[field].fillna("").astype(str).str.strip()
    require(int(atlas.dominant_hue9.isin(PALETTE).sum())==AFTER_HUES,
            "Eight newly photographed dominant hue taxa not in current atlas")
    for field in ("inat_taxon_id","species","source_hue9","hue_status",
                  "COL_XR_checklistKey","status","key"):
        old[field]=old[field].fillna("").astype(str).str.strip()
    require(old.inat_taxon_id.is_unique and
            int(old.status.eq("colxr_exact_accepted_species").sum())==PRIOR_ACCEPTED
            and old.COL_XR_checklistKey.eq(COL_XR).all()
            and old.loc[old.status.eq("colxr_exact_accepted_species"),"key"].is_unique,
            "Previous 37095 exact unique COL XR accepted species key source changed")
    indexed=atlas.set_index("inat_taxon_id").reindex(old.inat_taxon_id)
    require(indexed.species.notna().all() and
            (indexed.species.to_numpy()==old.species.to_numpy()).all() and
            (indexed.dominant_hue9.to_numpy()==old.source_hue9.to_numpy()).all() and
            (indexed.hue_status.to_numpy()==old.hue_status.to_numpy()).all(),
            "Previously COL XR named species/hue statuses drifted")
    new=atlas.loc[atlas.dominant_hue9.isin(PALETTE)&
                  atlas.species.str.match(NAME_RULE)&
                  ~atlas.inat_taxon_id.isin(old.inat_taxon_id)].copy()
    require(len(new)<=AFTER_HUES-PRIOR_HUES and new.inat_taxon_id.is_unique,
            "Unexpected number of newly clear original-photo hue names")
    return new.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)


def run(atlas_path,append_receipt,prior_names,prior_receipt,out,*,pause=1.0):
    out=Path(out);require(not out.exists(),"Never overwrite signed new COL XR name match ledger")
    app=json.loads(Path(append_receipt).read_text())
    require(app.get("schema")==APPEND_SCHEMA and
            app.get("stage")=="full_42111_taxon_atlas_new1909_true_photo9_appended" and
            app.get("nominal_taxa_with_clear_sampled_photo_hue")==AFTER_HUES
            and app.get("nominal_taxa_with_genuine_nineROI")==39911
            and app.get("new_strict9_numeric_original_photos")==8
            and app.get("new_clear_photo_hue_taxa")==8
            and app.get("prior_original_photo_and_observation_UUID_collisions")==0
            and app.get("output_SHA256",{}).get(Path(atlas_path).name)==source_sha(atlas_path)
            and source_sha(atlas_path)=="588f74124ff2982af5ada2fe8c0258b0e28f9d91010e0451b0b3023e7c67fbab",
            "Signed remainder18 original image9 atlas/receipt provenance mismatch")
    prior=json.loads(Path(prior_receipt).read_text())
    require(prior.get("schema")==PRIOR_MATCH_SCHEMA and
            prior.get("cumulative_exact_accepted_COLXR_species")==PRIOR_ACCEPTED and
            prior.get("signed_outputs_SHA256",{}).get(Path(prior_names).name)==source_sha(prior_names),
            "Previous 37095 accepted-source COL XR names SHA mismatch")
    atlas=pd.read_csv(atlas_path,compression="infer",dtype={"inat_taxon_id":str})
    prev=pd.read_csv(prior_names,compression="infer",dtype=str,keep_default_na=False)
    added=select_new(atlas,prev)
    matched=classify_new_rows(added,batch_size=50,pause=pause)
    full,accepted=combine_names(prev,matched)
    traits=build_traits(atlas,accepted,verify_counts=False)
    require(len(accepted)==len(traits)>=PRIOR_ACCEPTED and
            len(full)==len(prev)+len(added),
            "GBIF exact accepted species had lost source-measured 9 colour fractions")
    out.mkdir(parents=True)
    files={}
    for name,df in (
        ("new_postsnapshot_only_COLXR_name_matching.csv.gz",matched),
        ("current_all_photo9_COLXR_names.csv.gz",full),
        ("current_exact_accepted_COLXR_species.csv.gz",accepted),
        ("current_exact_accepted_COLXR_numeric9_traits.csv.gz",traits),
    ):
        p=out/name
        df.to_csv(p,index=False,compression="gzip",float_format="%.7f")
        files[name]=source_sha(p)
    requests=make_requests(full,chunk_size=250)
    sql=out/"GBIF_ISEA3H8_SQL_unsubmitted";sql.mkdir()
    for i,r in enumerate(requests):
        p=sql/f"colxr_photo9_r8_{i:03d}.json"
        p.write_text(json.dumps(r["GBIF_download_request"],indent=2)+"\n")
        files[str(p.relative_to(out))]=source_sha(p)
    report={
        "schema":SCHEMA,
        "source_original_MOIS_nominal_taxa":TAXA,
        "source_genuine_numeric9_taxa":39911,
        "source_clear_sampled_photo_hue_taxa":AFTER_HUES,
        "previous_exact_accepted_COLXR_species":PRIOR_ACCEPTED,
        "new_clear_photo_hue_taxa_name_eligible":len(added),
        "new_exact_accepted_COLXR_species":int(
            matched.status.eq("colxr_exact_accepted_species").sum()),
        "cumulative_exact_accepted_COLXR_species":len(accepted),
        "current_exact_accepted_numeric9_trait_species":len(traits),
        "new_taxonomic_match_status_counts":{
            str(k):int(v) for k,v in matched.status.value_counts().items()},
        "GBIF_ISEA3H8_SQL_unsubmitted_files":len(requests),
        "GBIF_occurrence_authenticated_download_submitted":False,
        "GBIF_occurrence_records_downloaded":0,
        "species_wide_photo_hue_fixed_proven":False,
        "source_photo_atlas_SHA256":source_sha(atlas_path),
        "prior_GBIF_colxr_name_SHA256":source_sha(prior_names),
        "signed_outputs_SHA256":files,
    }
    (out/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_POSTSNAPSHOT_REMAINING18_NEW8_OFFICIAL_COLXR_NAME_MATCH",
          json.dumps(report,sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("atlas","append-receipt","prior-names","prior-receipt","output"):
        p.add_argument("--"+k,required=True)
    p.add_argument("--pause",type=float,default=1.0)
    a=p.parse_args()
    run(a.atlas,a.append_receipt,a.prior_names,a.prior_receipt,a.output,pause=a.pause)
