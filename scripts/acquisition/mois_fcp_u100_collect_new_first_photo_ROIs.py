#!/usr/bin/env python3
"""Collect never-before-measured FCP U100 depth1 photos after K1 resource gate.

FCP's first coverage decision is one original observation and one original
photo per each of 20,242 eligible nominal taxa from all 93,209,351 candidate
photo IDs. Previous source-verified ROIs for 15,483 exact photo IDs are reused;
1,618 identical photo IDs are already reserved in the independent active
80,035-observation K1 run; only 3,141 NEW photo IDs need inference. This
collector validates those 3,141 only, never counts 20,242 metadata photo IDs
as if their flower hues were known.

Every photo's inclusion probability is exactly 1/N_taxon/M_observation, with
N and M from the immutable first-prefix signed global 93m source photo roster.
The workflow downloads raw photo pixels once in separate model workers. This
collector is data-only: no new S3 fetch, YOLO+EfficientSAM model launch, or
historical 0.9m-row Parquet rewrite.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import collect_wave, require, sha, HUES
from scripts.acquisition.mois_fcp_u100_depth1_prior_k1_reuse_queue import (
    PREFIX_ARM, FULL_TARGET,
)

SCHEMA="mois_fcp_u100_depth1_3141_new_original_ROI_photo_real_completion_v1"
PREFIX_PHOTOS_SHA="0403905f04d2def07d3f2e0caee28728a54923ef187f1fc2c63318857012b713"
SOURCE_SHA="fe49344241fbc051eebbbcbc8eb627a040198dc2b3ee56c89401876b58d28a6a"
TOTAL_SOURCE=93_209_351
NEW=3141
PRIOR_REUSE=15483
RESERVED_K1=1618


def check_prefix_probabilities(prefix_csv, *, prefix_sha=PREFIX_PHOTOS_SHA,
                               expected_eligible=FULL_TARGET):
    require(sha(prefix_csv)==prefix_sha,
            "Globally original 93m U100 FCP per-taxon prefix photo roster SHA changed")
    source={}
    taxa=set();observations=set()
    with Path(prefix_csv).open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        required={"photo_id","observation_id","inat_taxon_id",
                  "photos_in_original_observation",
                  "full_93m_distinct_original_observations_of_taxon",
                  "depth1_joint_original_photo_inclusion_probability"}
        require(required<=set(r.fieldnames or ()),
                "Original FCP source photo selection probability fields absent")
        for row in r:
            pid=row["photo_id"];obs=row["observation_id"];taxon=row["inat_taxon_id"]
            require(pid.isdecimal() and taxon.isdecimal() and bool(obs)
                    and pid not in source and taxon not in taxa and obs not in observations,
                    "Global source original photo/taxon/observation repeated")
            n=int(row["full_93m_distinct_original_observations_of_taxon"])
            m=int(row["photos_in_original_observation"])
            p=float(row["depth1_joint_original_photo_inclusion_probability"])
            require(n>=100 and m>=1 and 0<p<=1 and math.isfinite(p)
                    and math.isclose(p,1/n/m,rel_tol=1e-12,abs_tol=0.0),
                    "Original FCP prefix1 photo probability 1/N/M changed")
            source[pid]=(obs,taxon,p)
            taxa.add(taxon);observations.add(obs)
    require(len(source)==expected_eligible,"Not exactly the original first-photo source taxa")
    return source


def collect(source, workers, prefix_csv, output, *,
            expected_manifest_sha=SOURCE_SHA, prefix_sha=PREFIX_PHOTOS_SHA,
            expected_new=NEW, expected_reuse=PRIOR_REUSE,
            expected_inflight=RESERVED_K1, expected_eligible=FULL_TARGET):
    import pandas as pd
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
    source,workers,output=map(Path,(source,workers,output))
    require(not output.exists(),"Immutable FCP original-photo ROI output must not be overwritten")
    mf=source/"manifest.json"
    require(sha(mf)==expected_manifest_sha,
            "Precolour FCP first-photo actual work manifest SHA changed")
    d=json.loads(mf.read_text(encoding="utf-8"))
    require(d.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and d.get("full_frame_photo_ids")==TOTAL_SOURCE
            and d.get("source_arm")==PREFIX_ARM
            and d.get("selected_FCP_U100_eligible_source_taxa")==expected_eligible
            and d.get("selected_original_photo_observations")==expected_eligible
            and d.get("previously_measured_original_photo_ids_reused")==expected_reuse
            and d.get("inflight_K1_photo_IDs_reserved_without_redownload")==expected_inflight
            and d.get("new_photo_attempts_planned")==expected_new
            and d.get("technical_pilot_photo_IDs_not_refetched")==0
            and d.get("unsupported_original_photo_format_selected_as_unknown")==0
            and d.get("colour_measured") is False
            and d.get("photo_already_reserved_by_K1_never_fetched_twice") is True
            and d.get("FCP_top100_photo_prefix_inclusion_probability_corrected_to_1_over_N_M") is True,
            "FCP first-photo source not the pinned 20242-taxon original probability design")
    require(expected_new+expected_inflight+expected_reuse==expected_eligible,
            "FCP 20,242 original first-photo decisions were not all accounted for")
    prefix=check_prefix_probabilities(prefix_csv,prefix_sha=prefix_sha,
                                      expected_eligible=expected_eligible)
    data,audit=collect_wave(source,workers,METHOD)
    require(audit["wave_complete"] and len(data)==expected_new
            and audit["attempted"]==expected_new
            and audit["planned_photo_attempts"]==expected_new
            and len(audit["completed_units"])==d["work_units"],
            "Original 3,141 new FCP depth1 photo work units have not all been sealed")
    seen_taxa=set();seen_obs=set();seen_photo=set()
    unknown=Counter()
    for row in data:
        pid=row["photo_id"];taxon=row["inat_taxon_id"];obs=row["observation_id"]
        require(pid in prefix
                and pid not in seen_photo and taxon not in seen_taxa
                and obs not in seen_obs
                and prefix[pid][:2]==(obs,taxon)
                and row["source_arm"]==PREFIX_ARM,
                "Original first photo/observation/taxon/source arm altered in ROI result")
        p=float(row["selection_probability"])
        require(math.isfinite(p) and 0<p<=1 and
                math.isclose(p,prefix[pid][2],rel_tol=1e-12,abs_tol=0.0),
                "Wrong real source photo conditional FCP inclusion probability")
        seen_photo.add(pid);seen_taxa.add(taxon);seen_obs.add(obs)
        known=row["morph"] in HUES
        failed=row["measurement_status"]=="image_acquisition_failed"
        require(row["morph"] in (*HUES,"mixed_uncertain")
                and known==(row["measurement_status"]=="classified_four_state_morph")
                and not(known and failed),
                "Technical nonresponse was incorrectly assigned petal hue")
        unknown[row["morph"]]+=1
        unknown["classified"]+=known
        unknown["acquisition_failed"]+=failed
        unknown["unclassifiable_acquired"]+=not known and not failed
    require(unknown["classified"]==audit["classifiable"]
            and unknown["acquisition_failed"]==audit["acquisition_failed"]
            and unknown["unclassifiable_acquired"]==audit["unclassifiable_acquired"],
            "Original FCP frozen ROI-v4 true-photo counts were changed")
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-fcp-U100-depth1-RGB-",dir=output.parent))
    try:
        file=stage/"new_FCP_depth1_original_photo_ROIs.parquet"
        df=pd.DataFrame(data).sort_values("photo_id",key=lambda s:s.map(int)).reset_index(drop=True)
        df.to_parquet(file,compression="zstd",index=False)
        stats={
            "schema":SCHEMA,
            "original_93m_photo_ID_source_population":TOTAL_SOURCE,
            "eligible_global_first_FCP_nominal_taxon_IDs":expected_eligible,
            "first_original_photo_candidates_retained":expected_eligible,
            "previously_measured_original_image_ROIs_reused_elsewhere":expected_reuse,
            "inflight_K1_original_photo_IDs_reserved_elsewhere":expected_inflight,
            "newly_true_image_measured_original_photo_attempts":expected_new,
            "newly_true_image_FCP_four_hue_classified":unknown["classified"],
            "newly_true_image_acquisition_failed":unknown["acquisition_failed"],
            "newly_true_image_unclassifiable_acquired":unknown["unclassifiable_acquired"],
            "new_original_photo_hue_counts":{h:unknown[h] for h in (*HUES,"mixed_uncertain")},
            "new_distinct_source_taxa_and_observations":len(seen_obs),
            "frozen_global_FCP_U100_first_photo_SHA256":prefix_sha,
            "frozen_FCP_new_original_photo_source_SHA256":expected_manifest_sha,
            "actual_new_image_parquet_SHA256":sha(file),
            "signed_original_image_source_photo_SHA_audits":audit["completed_unit_audit"],
            "individual_photo_design_probabilities_1_over_N_over_M_preserved":True,
            "photo_colour_inferred_by_unchanged_frozen_FCP_ROI_v4":True,
            "new_images_downloaded_OR_ROI_inferred_by_this_data_only_collector":0,
            "previous_real_photo_pixel_Parquets_rewritten":False,
            "technical_unknowns_kept_unclassified":True,
            "individual_botanical_native_wild_anthesis_verified":False,
            "globally_unbiased_white_flora_fraction_estimated":False,
        }
        (stage/"receipt.json").write_text(json.dumps(stats,indent=2,sort_keys=True)+"\n")
        os.replace(stage,output)
        print("MOIS_FCP_U100_3141_NEVERBEFOREMEASURED_ORIGINAL_ROI_COMPLETE",
              json.dumps({k:stats[k] for k in (
                  "eligible_global_first_FCP_nominal_taxon_IDs",
                  "previously_measured_original_image_ROIs_reused_elsewhere",
                  "inflight_K1_original_photo_IDs_reserved_elsewhere",
                  "newly_true_image_measured_original_photo_attempts",
                  "newly_true_image_FCP_four_hue_classified")},sort_keys=True),flush=True)
        return stats
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("source","workers","prefix","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    collect(a.source,a.workers,a.prefix,a.output)


if __name__=="__main__":
    main()
