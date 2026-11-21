#!/usr/bin/env python3
"""Append independently audited post-snapshot original inaturalist ROI9 to 42111 taxa.

Five 2026-10-28+ created-at CC-original cohorts were independently source-
signed, FCP-v4 measured, and photo/observation-ID matched to the old complete
1,231,244 original photo registry. This final increment must not imply that
image SHA overlap against historical 1.2m physical pixels was fully checked.
"""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path
from urllib.parse import urlsplit

import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_palette9 import PALETTE,class_labels,sha256
from scripts.analysis.mois_species_first_wave3_incremental import add_to_catalogue
from scripts.acquisition.mois_zero9_recent_inat_new_created_photo_source import CC,ATLAS_SHA

SCHEMA="mois_post20261028_five_true_new_original_photo_ROI9_42111_append_v1"
AUDIT_SCHEMA="mois_post20261028_fullID_gated_actual_RGB9_independent_audit_v1"
WAVES={0:(2,1,1),64:(1,0,0),128:(3,1,0),192:(1,0,0),256:(4,1,1)}
PRIOR_ATTEMPTED=1231244
PRIOR_NUMERIC=467536
PRIOR_SUPPORTED=39882
PRIOR_LABELLED=39145


def require(ok,reason):
    if not bool(ok):raise ValueError("MOIS_POSTSNAPSHOT_APPEND: "+reason)


def load_five_audits(roots):
    require(set(roots)==set(WAVES),"Require exactly five independent audited hash strata")
    frames=[];proof={}
    for offset,(n,strict_count,clear_count) in WAVES.items():
        root=Path(roots[offset])
        receipt=root/"receipt.json"
        csv=root/"mois_new_inat_post20261028_actual_original_nine_roi.csv.gz"
        r=json.loads(receipt.read_text())
        require(r.get("schema")==AUDIT_SCHEMA
                and r.get("original_licensed_photo_attempts")==n
                and r.get("actual_genuine_strict9_numeric_success")==strict_count
                and r.get("actual_clear_photo_nine_palette_labels")==clear_count
                and r.get("actual_technical_unknowns")==n-strict_count
                and r.get("photo_observation_taxon_collision_count")==0
                and r.get("new_first_strict9_nominal_taxa")==strict_count
                and r.get("full_historical_raw_image_SHA_collision_check_possible") is False
                and r.get("official_iNat_42111_original_frame_auto_appended") is False
                and sha256(csv)==r.get("audited_derivative_sha256"),
                "Actual original-photo wave SHA or numeric count drift: "+str(offset))
        d=pd.read_csv(csv,compression="gzip",dtype=str,keep_default_na=False)
        require(len(d)==n and d.photo_id.is_unique and d.observation_id.is_unique
                and d.inat_taxon_id.is_unique and
                d.photo_id.str.fullmatch(r"\d+").all()
                and d.observation_id.str.fullmatch(
                    r"[0-9a-fA-F]{8}(?:-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}").all(),
                "New original photo/observation/taxon IDs duplicated or malformed")
        require(d.photo_license.isin(CC).all() and
                d.source_observation_created_at.ne("").all() and
                d.photo_url_large.str.startswith("https://inaturalist-open-data.s3.amazonaws.com/photos/").all(),
                "Post-snapshot original iNaturalist URL/license/date gate changed")
        for t in d.itertuples(index=False):
            from datetime import datetime,timezone
            dt=datetime.fromisoformat(t.source_observation_created_at.replace("Z","+00:00"))
            url=urlsplit(t.photo_url_large)
            require(dt.tzinfo is not None and
                    dt.astimezone(timezone.utc)>=datetime(2026,10,28,tzinfo=timezone.utc)
                    and url.hostname=="inaturalist-open-data.s3.amazonaws.com"
                    and url.path.startswith("/photos/"+t.photo_id+"/"),
                    "Photo URL or new observation source date inconsistent")
        cols=["flower_fraction_"+h for h in PALETTE]
        vec=d[cols].apply(pd.to_numeric,errors="coerce").to_numpy(float)
        numeric=(d.morph.isin(("white","yellow_orange","red_pink","blue_purple")).to_numpy(bool)
                 &np.isfinite(vec).all(axis=1)&(vec>=0).all(axis=1)
                 &(vec<=1.00001).all(axis=1)
                 &(abs(np.nan_to_num(vec).sum(axis=1)-1)<=.001))
        recorded=d.strict9_original_photo_numeric.map({"True":True,"False":False})
        require(recorded.notna().all() and np.array_equal(numeric,recorded.to_numpy(bool))
                and int(numeric.sum())==strict_count,
                "Attested strict9 status does not match true nine photo fractions")
        label=np.full(n,"",dtype=object)
        if numeric.any():label[numeric]=class_labels(vec[numeric])
        require(np.array_equal(label,d.dominant_hue9_label.to_numpy())
                and int((label!="").sum())==clear_count,
                "Attested dominant 9 hue was not measured from same photo")
        delta=d[["inat_taxon_id","photo_id","observation_id","image_sha256"]].copy()
        delta["strict9"]=numeric
        delta["dominant_label9"]=label
        for i,c in enumerate(PALETTE):
            delta["strict9_"+c]=np.where(numeric,vec[:,i],np.nan)
        frames.append(delta)
        proof[str(offset)]={"audited_receipt_SHA256":sha256(receipt),
                            "audited_true_numeric_photo_CSV_SHA256":sha256(csv),
                            "source_photo_manifest_SHA256":r["source_photo_manifest_sha256"],
                            "new_source_image_sha":r["new_photo_image_SHA256"]}
    result=pd.concat(frames,ignore_index=True)
    require(len(result)==11 and
            result.photo_id.is_unique and result.observation_id.is_unique and
            result.inat_taxon_id.is_unique and
            result.image_sha256.ne("").all() and
            result.image_sha256.is_unique,
            "Repeated source original photo, observation, taxon, or returned raw image SHA")
    require(int(result.strict9.sum())==3 and int(result.dominant_label9.ne("").sum())==2,
            "Combined five actual original RGB9 results changed")
    return result,proof


def append(atlas_path,historical_registry,roots,out):
    out=Path(out)
    require(not out.exists(),"Never overwrite signed post-snapshot original photo append")
    require(sha256(atlas_path)==ATLAS_SHA,"Signed base after-wave56 42111 taxon atlas SHA changed")
    atlas=pd.read_csv(atlas_path,compression="gzip",dtype={"inat_taxon_id":str})
    require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique
            and int(atlas.after_strict9.gt(0).sum())==PRIOR_SUPPORTED
            and int(atlas.dominant_hue9.fillna("").ne("").sum())==PRIOR_LABELLED,
            "Prior taxon-level numeric nine fractions/source hue ledger drifted")
    new,proof=load_five_audits(roots)
    history=pd.read_csv(historical_registry,compression="gzip",
                        usecols=["photo_id","observation_id"],dtype=str,keep_default_na=False)
    require(len(history)==PRIOR_ATTEMPTED and history.photo_id.is_unique,
            "Original frozen complete source photo ID registry absent")
    require(not set(history.photo_id).intersection(new.photo_id)
            and not (set(history.observation_id)-{""}).intersection(new.observation_id),
            "Post-snapshot new source photo/observation ID was previously attempted")
    atlas_id=atlas.set_index("inat_taxon_id")
    require(new.inat_taxon_id.isin(atlas_id.index).all()
            and atlas_id.loc[new.inat_taxon_id,"after_strict9"].eq(0).all(),
            "New post-snapshot photo belongs to a previously measured taxon")
    updated=add_to_catalogue(atlas,new,expected_photos=len(new))
    strict=int(new.strict9.sum())
    require(int(updated.after_strict9.gt(0).sum())==PRIOR_SUPPORTED+strict
            and int(updated.dominant_hue9.fillna("").ne("").sum())==PRIOR_LABELLED+2,
            "New post-snapshot photo source did not add exactly 3 numeric taxa and 2 clear hues")
    out.mkdir(parents=True)
    paths={
      "updated_atlas":out/"mois_42111_taxa_after_post_snapshot_5_signed_waves.csv.gz",
      "real_original_photo_delta":out/"mois_new11_iNat_original_photo9_audited_delta.parquet",
      "new_gbif_colxr_name_candidates":out/"mois_new2_clear9_source_taxa_for_GBIF_COLXR.csv",
    }
    updated.to_csv(paths["updated_atlas"],index=False,compression="gzip")
    new.to_parquet(paths["real_original_photo_delta"],index=False)
    two=updated.loc[updated.inat_taxon_id.isin(
        new.loc[new.dominant_label9.ne(""),"inat_taxon_id"])]
    two[["inat_taxon_id","species","hue_status","dominant_hue9"]].to_csv(
        paths["new_gbif_colxr_name_candidates"],index=False)
    report={
       "schema":SCHEMA,
       "cumulative_iNaturalist_original_photo_attempts":PRIOR_ATTEMPTED+len(new),
       "cumulative_iNaturalist_original_numeric_ROI9":PRIOR_NUMERIC+strict,
       "cumulative_nominal_taxa_with_numeric_nineROI":PRIOR_SUPPORTED+strict,
       "cumulative_nominal_taxa_with_clear_photo_hue":PRIOR_LABELLED+2,
       "nominal_taxa_without_numeric_nineROI":42111-PRIOR_SUPPORTED-strict,
       "new_genuine_original_photo_attempts":len(new),
       "new_nine_fraction_numeric_photo_ROI":strict,
       "new_clear_nine_palette_photo_labels":2,
       "new_technical_unknowns":len(new)-strict,
       "new_taxa_with_first_numeric_nineROI":strict,
       "source_iNaturalist_global_93m_index_untouched":True,
       "historical_original_photo_and_observation_ID_overlap":0,
       "additional_duplicate_photo_image_SHA_among_new11":0,
       "raw_image_SHA_overlap_with_historical_1231244_not_fully_verifiable":True,
       "species_wide_colour_fixed_proven":False,
       "GBIF_worldwide_occurrence_download_performed":False,
       "old_full_history_original_registry_SHA256":sha256(historical_registry),
       "prior_original_taxon_atlas_SHA256":sha256(atlas_path),
       "five_independent_actual_RGB9_audits":proof,
       "outputs_SHA256":{k:sha256(p) for k,p in paths.items()}
    }
    (out/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_POSTSNAPSHOT_SIGNED_11_ORIGINAL_ROI9_APPEND_42111",
          json.dumps(report,sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--atlas",required=True)
    p.add_argument("--history",required=True)
    p.add_argument("--audit",action="append",required=True,nargs=2,metavar=("OFFSET","DIRECTORY"))
    p.add_argument("--output",required=True)
    a=p.parse_args()
    roots={int(k):v for k,v in a.audit}
    append(a.atlas,a.history,roots,a.output)
