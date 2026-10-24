#!/usr/bin/env python3
"""FCP-style 93m progressive original-photo coverage: depth 1 before 5/20/100.

Read an ALREADY source-SHA-verified global U100 photo-observation roster for
one of 16 taxon-modulo groups. Original U100 ranks 100 distinct observations
per eligible taxon by a globally seeded observation priority, with one
independently uniformly photo-hashed original image per observation. Instead
of starting two MILLION expensive ROI evaluations at once, freeze the FIRST
hash-ranked original photo for each of the 20,242 eligible taxa.

This first-prefix sample is a random original observation among the N_i
observations of species i, NOT a uniform row among its photographs. For this
selected photo j of observation o, m_o = number of original photo IDs attached
to o, so the exact conditional design probability is:

         pi_1(photo j | eligible taxon i) = 1/N_i * 1/m_o.

This is 1/100th the probability written in its archived top-100 U100 parent
source (100/N_i/m_o). It is a NESTED DESIGN, so future stage depth 5, 20
and 100 can add only previously unmeasured top-ranked observations, reusing
both existing FCP and MOIS original photo ROIs. No colour, flowering, native
status, or unsupported photo format is used to select source candidates.

FCP-U100 eligibility excludes taxa with fewer than 100 original observations:
these remain in the full 42,111-taxon source census and complementary
season/latitude probability sampling arm, NOT declared absent or monochrome.

Zero photo downloads, zero model inferencing and zero photo exclusions based
on green/flower colour occur here.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile

from scripts.acquisition.mois_fcp_u100_full93m_source_onepass_group import (
    FCP_TARGET, SCHEMA as PARENT_SCHEMA,
    SOURCE_ARM as PARENT_ARM, original_observation_priority,
)
from scripts.acquisition.mois_fcp_u100_full93m_global_eligibility_ledger import (
    SCHEMA as GLOBAL_SCHEMA,
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    EXPECTED_PHOTOS, EXPECTED_TAXA, normalized_extension, require, sha256,
)

SCHEMA = "mois_full93m_fcp_u100_progressive_onephoto_per_eligible_taxon_v1"
ARM = "mois_full93m_fcp_u100_prefix1_one_observation_per_eligible_taxon_v1"
GLOBAL_ELIGIBLE_TAXA = 20_242
GLOBAL_FIRST100_PHOTOS = 2_024_200
GROUPS = 16
SOURCE_FIELDS = (
    "photo_id", "observation_id", "inat_taxon_id", "species",
    "observer_id", "observed_on", "recorded_month", "latitude_band",
    "latitude", "longitude", "positional_accuracy", "photo_license",
    "original_image_extension", "photos_in_original_observation",
    "full_93m_distinct_original_observations_of_taxon",
    "conditional_species_observation_probability",
    "conditional_one_photo_probability",
    "conditional_full_93m_photo_selection_probability",
    "photo_source_disposition", "photo_url_large", "source_arm",
)
OUTPUT_FIELDS = (
    "photo_id", "observation_id", "inat_taxon_id", "species",
    "observer_id", "observed_on", "recorded_month", "latitude_band",
    "latitude", "longitude", "positional_accuracy", "photo_license",
    "original_image_extension", "photos_in_original_observation",
    "full_93m_distinct_original_observations_of_taxon",
    "U100_source_original_observation_rank",
    "depth1_observation_inclusion_probability",
    "depth1_one_photo_inclusion_probability",
    "depth1_joint_original_photo_inclusion_probability",
    "parent_U100_joint_original_photo_inclusion_probability",
    "parent_original_photo_U100_source_arm",
    "progressive_source_arm",
    "photo_source_disposition", "photo_url_large",
    "independent_accepted_species_and_native_stage_verified",
)


def float_close(value, expected, label):
    try:
        observed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("Invalid " + label) from exc
    require(math.isfinite(observed) and 0 < observed <= 1
            and abs(observed-expected) <= 1e-12 * max(expected, observed),
            "Source-conditional inclusion probability mismatch: " + label)


def make_first_row(row):
    pid = row["photo_id"]
    obs = row["observation_id"]
    taxon = row["inat_taxon_id"]
    require(pid.isdecimal() and bool(obs) and taxon.isdecimal()
            and row["source_arm"] == PARENT_ARM,
            "Original U100 source photo, observation or source arm invalid")
    nobs = int(row["full_93m_distinct_original_observations_of_taxon"])
    m = int(row["photos_in_original_observation"])
    require(nobs >= FCP_TARGET and 1 <= m <= EXPECTED_PHOTOS,
            "Original full-photo U100 taxon/observation count is impossible")
    old_conditional = FCP_TARGET/nobs
    one_photo = 1/m
    parent_joint = old_conditional * one_photo
    first_conditional = 1/nobs
    joint = first_conditional * one_photo
    float_close(row["conditional_species_observation_probability"],
                old_conditional, "U100 observation")
    float_close(row["conditional_one_photo_probability"],one_photo,
                "uniform original photo within observation")
    float_close(row["conditional_full_93m_photo_selection_probability"],
                parent_joint,"U100 original photo")
    ext = normalized_extension(row["original_image_extension"])
    require(row["photo_source_disposition"] ==
            ("supported" if ext is not None else "unsupported_original_extension"),
            "Original photo extension/source technical status was silently filtered")
    if ext is None:
        require(not row["photo_url_large"],
                "Unsupported photo unexpectedly supplied a valid original-image URL")
    else:
        require(row["photo_url_large"] ==
                f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}",
                "Original-photo URL/ID differs from source")
    output = {key:row[key] for key in (
        "photo_id","observation_id","inat_taxon_id","species","observer_id",
        "observed_on","recorded_month","latitude_band","latitude","longitude",
        "positional_accuracy","photo_license","original_image_extension",
        "photos_in_original_observation",
        "full_93m_distinct_original_observations_of_taxon",
        "photo_source_disposition","photo_url_large",
    )}
    output.update(
        U100_source_original_observation_rank="1",
        depth1_observation_inclusion_probability=format(first_conditional,".17g"),
        depth1_one_photo_inclusion_probability=format(one_photo,".17g"),
        depth1_joint_original_photo_inclusion_probability=format(joint,".17g"),
        parent_U100_joint_original_photo_inclusion_probability=format(parent_joint,".17g"),
        parent_original_photo_U100_source_arm=PARENT_ARM,
        progressive_source_arm=ARM,
        independent_accepted_species_and_native_stage_verified="false",
    )
    return output


def freeze_group(source_csv, group_receipt, global_receipt, output, group, *,
                 expected_taxa=EXPECTED_TAXA,
                 expected_photos=EXPECTED_PHOTOS,
                 expected_global_eligible=GLOBAL_ELIGIBLE_TAXA):
    """One bounded pass over signed U100 metadata; never load image bytes."""
    source_csv,group_receipt,global_receipt,output = map(
        Path,(source_csv,group_receipt,global_receipt,output))
    require(not output.exists(),"Previously frozen full93m FCP stage1 source is immutable")
    require(type(group) is int and 0 <= group < GROUPS,
            "Source taxon group index is outside the frozen 16-part population")
    root=json.loads(global_receipt.read_text(encoding="utf-8"))
    r=json.loads(group_receipt.read_text(encoding="utf-8"))
    witnesses=root.get("group_original_photo_SHA_proof")
    require(root.get("schema")==GLOBAL_SCHEMA
            and root.get("source_all_original_photo_metadata_IDs")==expected_photos
            and root.get("all_source_nominal_taxon_IDs")==expected_taxa
            and root.get("source_FCP_style_exact100_photo_observation_taxa")==expected_global_eligible
            and root.get("source_FCP_style_exact100_observation_photo_targets")==
                    expected_global_eligible*FCP_TARGET
            and root.get("photo_ROI_pixels_opened_or_reinferred")==0
            and type(witnesses) is list and len(witnesses)==GROUPS,
            "Global full93m FCP-U100 42,111-taxon census receipt missing or modified")
    group_receipts={item["group"]:item for item in witnesses}
    require(set(group_receipts)==set(range(GROUPS))
            and sha256(group_receipt)==group_receipts[group]["group_receipt_sha256"],
            "Group source checksum contradicts globally certified full93m U100 source")
    require(r.get("schema")==PARENT_SCHEMA
            and r.get("original_taxon_modulo_source_group")==group
            and r.get("full_original_photo_source_IDs")==expected_photos
            and r.get("full_original_named_plant_taxa")==expected_taxa
            and r.get("FCP_original_photo_observation_target_per_taxon")==FCP_TARGET
            and r.get("no_original_image_downloads") is True
            and r.get("no_original_ROI_or_hue_inferences") is True
            and r.get("all_old_FCP_MOIS_original_photo_IDs_must_be_reused_before_any_new_ROI") is True
            and r.get("source_group_original_photo_ID_SHA",{}).get(
                 "original_full_93m_source_registry_SHA256") ==
                 root["full93m_original_source_registry_sha256"]
            and sha256(source_csv)==r.get("U100_complete_original_photo_candidate_SHA256"),
            "Not the original source-SHA-pinned FCP U100 taxon/photo roster")
    source_expected = r["FCP_U100_source_original_photo_opportunities"]
    eligible_expected = r["FCP_U100_source_eligible_nominal_taxa"]
    require(source_expected == FCP_TARGET*eligible_expected
            and 0 <= eligible_expected <= expected_taxa
            and eligible_expected == group_receipts[group]["eligible_taxa"],
            "FCP source conditional original observation count changed")

    output.parent.mkdir(parents=True,exist_ok=True)
    staged=Path(tempfile.mkdtemp(prefix="mois-fcp-depth1-source-",dir=output.parent))
    selected=staged/"FCP_U100_prefix1_original_photo_per_taxon.csv"
    try:
        total=0
        eligible=0
        supported=0
        unsupported=0
        selected_ids=set()
        used_taxa=set()
        current_taxon=None
        count_for_taxon=0
        old_priority=("", "")
        with (source_csv.open("r",encoding="utf-8",newline="") as input_stream,
              selected.open("w",encoding="utf-8",newline="") as out):
            reader=csv.DictReader(input_stream)
            require(tuple(reader.fieldnames or ())==SOURCE_FIELDS,
                    "FCP U100 source photo columns changed without source receipt migration")
            writer=csv.DictWriter(out,fieldnames=OUTPUT_FIELDS,lineterminator="\n")
            writer.writeheader()
            for row in reader:
                taxon=row["inat_taxon_id"]
                require(taxon.isdecimal() and int(taxon)%128//8==group,
                        "Photo-source taxon not owned by its full93m source part")
                if taxon != current_taxon:
                    if current_taxon is not None:
                        require(count_for_taxon == FCP_TARGET,
                                "Previous original U100 taxon did not retain exactly 100 observations")
                    require(taxon not in used_taxa, "Same taxon source roster is noncontiguous")
                    used_taxa.add(taxon)
                    current_taxon=taxon
                    count_for_taxon=0
                    old_priority=("", "")
                    eligible+=1
                count_for_taxon+=1
                require(count_for_taxon<=FCP_TARGET,
                        "One taxon had more than its frozen 100 original source observations")
                rank=(original_observation_priority(row["observation_id"]),
                      row["observation_id"])
                require(rank > old_priority,
                        "FCP source observation priority order changed inside eligible taxon")
                old_priority=rank
                if count_for_taxon == 1:
                    first=make_first_row(row)
                    require(first["photo_id"] not in selected_ids,
                            "Duplicate selected first-prefix original photo ID")
                    selected_ids.add(first["photo_id"])
                    writer.writerow(first)
                    supported+=first["photo_source_disposition"]=="supported"
                    unsupported+=first["photo_source_disposition"]!="supported"
                total+=1
            if current_taxon is not None:
                require(count_for_taxon==FCP_TARGET,
                        "Final source taxon lost its last original photo/observation")
        require(total==source_expected
                and eligible==eligible_expected
                and len(selected_ids)==eligible
                and supported+unsupported==eligible,
                "Incomplete full93m FCP U100 prefix sample photo source")
        source_stats={
            "schema":SCHEMA,
            "parent_original_FCP_U100_full93m_source_schema":PARENT_SCHEMA,
            "global_full93m_FCP_source_receipt_SHA256":sha256(global_receipt),
            "global_full93m_photo_ID_opportunities":expected_photos,
            "global_full93m_nominal_taxon_IDs":expected_taxa,
            "global_FCP_U100_source_eligible_taxa":expected_global_eligible,
            "original_taxon_modulo_source_group":group,
            "group_original_FCP_U100_source_receipt_SHA256":sha256(group_receipt),
            "group_original_FCP_U100_candidate_SHA256":sha256(source_csv),
            "U100_original_source_candidates_accounted":total,
            "U100_source_eligible_taxa_in_group":eligible,
            "depth1_original_photo_observations_selected":eligible,
            "depth1_supported_original_image_formats":supported,
            "depth1_unsupported_original_image_formats_unknown":unsupported,
            "source_taxa_below100_observations_retained_in_42111_census":
                r["source_nominal_taxa"] - eligible,
            "depth1_original_photo_candidate_SHA256":sha256(selected),
            "original_FCP_prefix_depths_prospectively_nested":[1,5,20,100],
            "original_source_probability_formula":"1/N_taxon / M_original_observation",
            "source_probability_true_for_eligible_taxa_only":True,
            "source_ineligible_taxa_are_not_missing_or_nonflowering":True,
            "randomization_unconditional_on_photo_hue_or_image_quality":True,
            "original_photo_source_prior_ROI_reuse_not_yet_applied":True,
            "original_images_downloaded":0,
            "original_colour_inferences_executed":0,
            "native_wild_or_individual_flowering_photo_verified":False,
            "not_a_global_unbiased_flora_white_flower_census":True,
        }
        (staged/"FCP_U100_prefix1_group_source_receipt.json").write_text(
            json.dumps(source_stats,indent=2,sort_keys=True)+"\n")
        os.replace(staged,output)
        print("MOIS_FCP_FULL93M_U100_PROGRESSIVE_DEPTH1_PER_TAXON",
              json.dumps({k:source_stats[k] for k in (
                  "original_taxon_modulo_source_group",
                  "U100_source_eligible_taxa_in_group",
                  "depth1_original_photo_observations_selected",
                  "depth1_supported_original_image_formats",
                  "depth1_unsupported_original_image_formats_unknown")},sort_keys=True),flush=True)
        return source_stats
    except BaseException:
        shutil.rmtree(staged,ignore_errors=True)
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("source-csv","group-receipt","global-receipt","output"):
        p.add_argument("--"+flag,type=Path,required=True)
    p.add_argument("--group",type=int,required=True)
    x=p.parse_args()
    freeze_group(x.source_csv,x.group_receipt,x.global_receipt,x.output,x.group)


if __name__=="__main__":
    main()
