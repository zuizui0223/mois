#!/usr/bin/env python3
"""FCP-style nested source depth 5 from SHAs of already-selected top-100 photos.

For each of 20,242 eligible nominal taxa, retain one original photo from each
of the five highest-priority INDEPENDENT original iNaturalist observations.
The rank-1 photo is the EXACT prior FCP-U100 depth-1 source member. Ranks 2-5
are incremental candidate metadata ONLY, not automatically new image tasks:
all existing FCP/MOIS original-photo IDs MUST be re-used first.

Select by the same globally pinned original-observation hash priority as the
full 93,209,351-photo source. Preserve true nested inclusion probability
P(photo in top-5 taxon prefix)=5/N_taxon/M_photo_per_observation. The
incremental 2..5 candidate's unconditional inclusion probability is
4/N_taxon/M; these are two different estimands, both explicitly recorded.
No photos downloaded, no ROI model inference or colour-dependent filtering.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_fcp_u100_progressive_depth1_photo_source import (
    SOURCE_FIELDS, ARM as FIRST_ARM, GLOBAL_ELIGIBLE_TAXA, GROUPS,
    make_first_row,
)
from scripts.acquisition.mois_fcp_u100_full93m_source_onepass_group import (
    SCHEMA as PARENT_SCHEMA, SOURCE_ARM as PARENT_ARM,
    FCP_TARGET, original_observation_priority,
)
from scripts.acquisition.mois_fcp_u100_full93m_global_eligibility_ledger import (
    SCHEMA as GLOBAL_SCHEMA,
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    EXPECTED_PHOTOS, EXPECTED_TAXA, require, sha256,
)

SCHEMA = "mois_full93m_FCP_U100_nested_5_original_observation_photos_per_taxon_v1"
ARM = "mois_full93m_FCP_U100_nested_prefix5_source_v1"
DEPTH = 5
GLOBAL_RECEIPT_SHA = "152e77c2ac9c3b414a2d30d1e4e5cc74ee93bbf60ffc2993234551cd28fc39d3"
FIELDS = (
    "photo_id", "observation_id", "inat_taxon_id", "species",
    "observer_id", "observed_on", "recorded_month", "latitude_band",
    "latitude", "longitude", "positional_accuracy", "photo_license",
    "original_image_extension", "photos_in_original_observation",
    "full_93m_distinct_original_observations_of_taxon",
    "U100_source_original_observation_rank",
    "source_rank1_original_photo_id",
    "cumulative_prefix5_observation_probability",
    "cumulative_prefix5_joint_original_photo_probability",
    "rank2to5_incremental_observation_probability",
    "rank2to5_incremental_joint_original_photo_probability",
    "depth1_joint_original_photo_inclusion_probability",
    "parent_U100_joint_original_photo_inclusion_probability",
    "prefix_depth", "incremental_since_depth1",
    "source_arm", "parent_original_photo_U100_source_arm",
    "photo_source_disposition", "photo_url_large",
    "independent_accepted_species_and_native_stage_verified",
)


def source_row(row, rank, first_photo):
    """Validate upstream full-U100 row exactly and retain all probability axes."""
    proof = make_first_row(row)
    nobs = int(row["full_93m_distinct_original_observations_of_taxon"])
    m = int(row["photos_in_original_observation"])
    require(1 <= rank <= DEPTH and nobs >= FCP_TARGET and m >= 1,
            "Original nested stage-five observation/photo multiplicity changed")
    first = float(proof["depth1_joint_original_photo_inclusion_probability"])
    prefix = DEPTH / nobs / m
    delta = (DEPTH - 1) / nobs / m
    parent = float(proof["parent_U100_joint_original_photo_inclusion_probability"])
    require(
        math.isfinite(prefix) and 0 < first <= prefix <= parent <= 1
        and math.isclose(prefix, DEPTH * first, rel_tol=1e-12)
        and math.isclose(delta, (DEPTH - 1) * first, rel_tol=1e-12),
        "Frozen nested FCP source original-photo inclusion probabilities disagree")
    out = {k: row[k] for k in (
        "photo_id", "observation_id", "inat_taxon_id", "species", "observer_id",
        "observed_on", "recorded_month", "latitude_band", "latitude", "longitude",
        "positional_accuracy", "photo_license", "original_image_extension",
        "photos_in_original_observation",
        "full_93m_distinct_original_observations_of_taxon",
        "photo_source_disposition", "photo_url_large",
    )}
    out.update(
        U100_source_original_observation_rank=str(rank),
        source_rank1_original_photo_id=str(first_photo),
        cumulative_prefix5_observation_probability=format(DEPTH/nobs, ".17g"),
        cumulative_prefix5_joint_original_photo_probability=format(prefix, ".17g"),
        rank2to5_incremental_observation_probability=format((DEPTH-1)/nobs, ".17g"),
        rank2to5_incremental_joint_original_photo_probability=format(delta, ".17g"),
        depth1_joint_original_photo_inclusion_probability=format(first, ".17g"),
        parent_U100_joint_original_photo_inclusion_probability=format(parent, ".17g"),
        prefix_depth=str(DEPTH),
        incremental_since_depth1=str(rank > 1).lower(),
        source_arm=ARM,
        parent_original_photo_U100_source_arm=PARENT_ARM,
        independent_accepted_species_and_native_stage_verified="false",
    )
    require(set(out)==set(FIELDS), "Nested full93m FCP source photo field schema drift")
    return {key:out[key] for key in FIELDS}


def freeze_group(source_csv, source_receipt, global_receipt, output, group, *,
                 expected_global_sha=GLOBAL_RECEIPT_SHA,
                 expected_global_taxa=GLOBAL_ELIGIBLE_TAXA,
                 expected_all_taxa=EXPECTED_TAXA, expected_photos=EXPECTED_PHOTOS):
    output=Path(output)
    require(not output.exists(), "Never mutate frozen nested stage5 candidate source")
    require(type(group) is int and 0<=group<GROUPS,
            "Original SHA-authenticated FCP U100 taxon group outside 0..15")
    require(sha256(global_receipt)==expected_global_sha,
            "Worldwide 42111-taxon U100 source global receipt SHA changed")
    root=json.loads(Path(global_receipt).read_text(encoding="utf-8"))
    receipt=json.loads(Path(source_receipt).read_text(encoding="utf-8"))
    witnesses=root.get("group_original_photo_SHA_proof")
    require(
        root.get("schema")==GLOBAL_SCHEMA
        and root.get("source_all_original_photo_metadata_IDs")==expected_photos
        and root.get("all_source_nominal_taxon_IDs")==expected_all_taxa
        and root.get("source_FCP_style_exact100_photo_observation_taxa")==expected_global_taxa
        and root.get("source_FCP_style_exact100_observation_photo_targets")
            ==expected_global_taxa * FCP_TARGET
        and root.get("photo_ROI_pixels_opened_or_reinferred")==0
        and isinstance(witnesses, list) and len(witnesses)==GROUPS,
        "Complete global original FCP U100 93m-source proof invalid")
    refs={x["group"]:x for x in witnesses}
    require(set(refs)==set(range(GROUPS))
            and sha256(source_receipt)==refs[group]["group_receipt_sha256"],
            "FCP original top100 per-group receipt not in the global signed audit")
    require(
        receipt.get("schema")==PARENT_SCHEMA
        and receipt.get("original_taxon_modulo_source_group")==group
        and receipt.get("full_original_photo_source_IDs")==expected_photos
        and receipt.get("full_original_named_plant_taxa")==expected_all_taxa
        and receipt.get("FCP_original_photo_observation_target_per_taxon")==FCP_TARGET
        and receipt.get("no_original_image_downloads") is True
        and receipt.get("no_original_ROI_or_hue_inferences") is True
        and receipt.get("all_old_FCP_MOIS_original_photo_IDs_must_be_reused_before_any_new_ROI") is True
        and receipt.get("source_group_original_photo_ID_SHA",{}).get(
            "original_full_93m_source_registry_SHA256")
            ==root["full93m_original_source_registry_sha256"]
        and sha256(source_csv)==receipt["U100_complete_original_photo_candidate_SHA256"],
        "Original ranked 100-photo taxon source or SHA changed")
    expected_taxa=receipt["FCP_U100_source_eligible_nominal_taxa"]
    source_expected=receipt["FCP_U100_source_original_photo_opportunities"]
    require(source_expected==expected_taxa*FCP_TARGET
            and expected_taxa==refs[group]["eligible_taxa"],
            "Original full-U100 source taxon count contradicts signed census")

    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix=f"mois-fcp-nested5-g{group:02d}-",
                                dir=output.parent))
    prefix=stage/"original_FCP_U100_prefix5_rank1to5_candidates.csv"
    delta=stage/"original_FCP_U100_rank2to5_INCREMENTAL_CANDIDATES_not_new_ROIs.csv"
    visited_taxa=set()
    original_ids=set()
    original_obs=set()
    selected=set()
    selected_obs=set()
    incremental=set()
    first_photo_ids=[]
    count=0
    supported=unsupported=0
    try:
        with (Path(source_csv).open("r",encoding="utf-8",newline="") as stream,
              prefix.open("w",encoding="utf-8",newline="") as all_out,
              delta.open("w",encoding="utf-8",newline="") as new_out):
            reader=csv.DictReader(stream)
            require(tuple(reader.fieldnames or ())==SOURCE_FIELDS,
                    "Original source top100 metadata schema changed")
            a=csv.DictWriter(all_out,fieldnames=FIELDS,lineterminator="\n")
            b=csv.DictWriter(new_out,fieldnames=FIELDS,lineterminator="\n")
            a.writeheader();b.writeheader()
            current=None
            rank=0
            prev_rank=("", "")
            first_id=None
            nobs_for_taxon=None
            for row in reader:
                pid=row["photo_id"]
                obs=row["observation_id"]
                taxon=row["inat_taxon_id"]
                require(pid.isdecimal() and taxon.isdecimal() and bool(obs)
                        and int(taxon)%128//8==group
                        and pid not in original_ids and obs not in original_obs
                        and row["source_arm"]==PARENT_ARM,
                        "Original 100-observation source photo ID duplicated or group changed")
                original_ids.add(pid);original_obs.add(obs)
                if taxon!=current:
                    if current is not None:
                        require(rank==FCP_TARGET,
                                "Prior original FCP taxon did not retain exactly 100 observations")
                    require(taxon not in visited_taxa,
                            "Original FCP-U100 candidate source taxon not contiguous")
                    visited_taxa.add(taxon)
                    current=taxon
                    rank=0
                    prev_rank=("", "")
                    first_id=None
                    nobs_for_taxon=None
                rank+=1
                require(rank<=FCP_TARGET,
                        "Original taxon selected more than 100 observations")
                source_rank=(original_observation_priority(obs),obs)
                require(source_rank>prev_rank,
                        "Original observation priority changed from global hash")
                prev_rank=source_rank
                nobs=int(row["full_93m_distinct_original_observations_of_taxon"])
                if nobs_for_taxon is None:nobs_for_taxon=nobs
                require(nobs_for_taxon==nobs,"Original observation denominator changes within taxon")
                if rank==1:
                    first_id=pid
                    first_photo_ids.append(pid)
                if rank<=DEPTH:
                    out=source_row(row,rank,first_id)
                    require(pid not in selected and obs not in selected_obs,
                            "Nested original selected photo ID/observation repeated")
                    selected.add(pid);selected_obs.add(obs)
                    a.writerow(out)
                    if rank>1:
                        require(pid not in incremental,"Nested delta original photo ID duplicated")
                        incremental.add(pid)
                        b.writerow(out)
                    supported+=out["photo_source_disposition"]=="supported"
                    unsupported+=out["photo_source_disposition"]!="supported"
                count+=1
            if current is not None:
                require(rank==FCP_TARGET,
                        "Final original FCP source taxon missing one of 100 observations")
        require(count==source_expected
                and len(visited_taxa)==expected_taxa
                and len(selected)==DEPTH*expected_taxa
                and len(incremental)==(DEPTH-1)*expected_taxa
                and len(first_photo_ids)==expected_taxa
                and supported+unsupported==len(selected)
                and len(original_ids)==source_expected==len(original_obs),
                "FCP stage5 candidate photo/observation source conservation violated")
        stats={
            "schema":SCHEMA,
            "global_source_original_photo_IDs_metadata_only":expected_photos,
            "global_nominal_source_taxa":expected_all_taxa,
            "global_U100_eligible_taxa":expected_global_taxa,
            "original_taxon_modulo_source_group":group,
            "source_group_original_photo_rank100_receipt_SHA256":sha256(source_receipt),
            "source_group_original_photo_rank100_CSV_SHA256":sha256(source_csv),
            "global_U100_complete_42111_source_receipt_SHA256":sha256(global_receipt),
            "original_photo_rank1_depth1_source_ARM_reused":FIRST_ARM,
            "nested_original_photo_depth5_source_arm":ARM,
            "source_group_eligible_nominal_taxa":expected_taxa,
            "full_100_observation_photo_rank_candidates_audited":count,
            "first5_observation_independent_original_photo_candidates":len(selected),
            "first1_original_photo_IDs_retained":len(first_photo_ids),
            "new_candidate_ranks2to5_NOT_yet_measured":len(incremental),
            "first5_supported_photo_extensions":supported,
            "first5_unsupported_photo_extensions_technical_unknown":unsupported,
            "nested_first5_photo_inclusion_probability":
                "5 / N_original_distinct_observations_of_taxon / M_photos_of_selected_observation",
            "rank2to5_only_inclusion_probability":
                "4 / N_original_distinct_observations_of_taxon / M_photos_of_selected_observation",
            "conditional_on_original_source_taxon_at_least100_observations":True,
            "source_group_original_photo_IDs_distinct":True,
            "source_group_original_observations_distinct":True,
            "source_group_candidate_photo_rank_unchanged":True,
            "first5_csv_SHA256":sha256(prefix),
            "ranks2to5_csv_SHA256":sha256(delta),
            "original_photo_IDs_already_measured_by_old_FCP_MOIS_and_K1_excluded":False,
            "new_image_download_or_frozen_ROI_inference_this_stage":0,
            "wild_native_flowering_photo_or_accepted_taxonomy_independently_verified":False,
            "globally_unbiased_white_flora_frequency_estimated":False,
        }
        (stage/"nested5_source_receipt.json").write_text(
            json.dumps(stats,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_FCP_NESTED5_ORIGINAL_PRECOLOUR_SOURCE",json.dumps({
        "group":group,
        "eligible_taxa":expected_taxa,
        "first5_photo_candidates":len(selected),
        "new_rank2to5_source_candidates_before_prior_ROI_reuse":len(incremental),
        "original_photo_bytes_downloaded":0,
        "colour_inference":0,
    },sort_keys=True),flush=True)
    return stats


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",required=True,type=Path)
    ap.add_argument("--receipt",required=True,type=Path)
    ap.add_argument("--global-receipt",required=True,type=Path)
    ap.add_argument("--group",required=True,type=int)
    ap.add_argument("--output",required=True,type=Path)
    a=ap.parse_args()
    freeze_group(a.source,a.receipt,a.global_receipt,a.output,a.group)


if __name__=="__main__":
    main()
