#!/usr/bin/env python3
"""Assemble independently SHA-signed 16-group FCP U100 nested first5 original source.

100% metadata-only. Preserve exact photo, original observation and taxon keys
and original hash-rank first1→first5 nesting. Distinguish 80,968 incremental
rank2to5 SOURCE candidates from never-before-inferred photos, which require
historical FCP/MOIS/K1 photo-ID dedup before any new image download.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_fcp_u100_depth5_source import (
    SCHEMA as GROUP_SCHEMA, ARM, FIELDS, DEPTH, GROUPS, GLOBAL_RECEIPT_SHA,
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    EXPECTED_PHOTOS, EXPECTED_TAXA, require, sha256,
)

SCHEMA = "mois_full93m_FCP_U100_nested_first5_all20242_original_photo_source_v1"
GLOBAL_ELIGIBLE=20_242
STAGE_FIRST5=101_210
STAGE_INCREMENTAL=80_968
FIRST1_SHA="0403905f04d2def07d3f2e0caee28728a54923ef187f1fc2c63318857012b713"
PREFIX_FILENAME="original_FCP_U100_prefix5_rank1to5_candidates.csv"
DELTA_FILENAME="original_FCP_U100_rank2to5_INCREMENTAL_CANDIDATES_not_new_ROIs.csv"


def read_csv(path):
    with Path(path).open("r",encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        require(tuple(reader.fieldnames or ())==FIELDS,
                "Frozen U100 nested5 source photo field schema changed")
        yield from reader


def global_source(folders, first1, output, *,
                  expected_groups=GROUPS,
                  expected_eligible=GLOBAL_ELIGIBLE,
                  expected_prefix5=STAGE_FIRST5,
                  expected_incremental=STAGE_INCREMENTAL,
                  expected_first1_sha=FIRST1_SHA):
    output=Path(output)
    require(not output.exists(),"Previously certified world U100 nested5 source is immutable")
    require(len(folders)==expected_groups and expected_prefix5==DEPTH*expected_eligible
            and expected_incremental==(DEPTH-1)*expected_eligible,
            "Missing original 16 group five-observation source or population denominator")
    require(sha256(first1)==expected_first1_sha,
            "Exact old FCP prefix1 original photo membership SHA changed")
    first={}
    with Path(first1).open(encoding="utf-8",newline="") as s:
        reader=csv.DictReader(s)
        require({"photo_id","inat_taxon_id","observation_id",
                 "depth1_joint_original_photo_inclusion_probability"}
                <=set(reader.fieldnames or ()),
                "Previously signed global U100 first original photo missing fields")
        for row in reader:
            taxon=row["inat_taxon_id"]
            require(taxon.isdecimal() and taxon not in first
                    and row["photo_id"].isdecimal() and bool(row["observation_id"]),
                    "Nested depth1 original photo/taxon not unique")
            first[taxon]=row
    require(len(first)==expected_eligible,
            "Previously sampled source U100 first photo does not cover all eligible taxa")
    ranks_by_taxon={}
    ids=set()
    observations=set()
    delta_ids=set()
    group_audits={}
    old_taxa_groups=set()
    technical=0
    for folder in folders:
        folder=Path(folder)
        receipt_file=folder/"nested5_source_receipt.json"
        receipt=json.loads(receipt_file.read_text(encoding="utf-8"))
        group=receipt.get("original_taxon_modulo_source_group")
        # Give an exact contract failure rather than a generic source mismatch:
        # the immutable published group metadata may differ from a synthetic
        # unit-test fixture. Any mismatch still FAILS CLOSED before source data.
        checks = {
            "group_valid_unique":type(group) is int and 0<=group<expected_groups
                and group not in group_audits,
            "schema":receipt.get("schema")==GROUP_SCHEMA,
            "global_U100_receipt_SHA":receipt.get(
                "global_U100_complete_42111_source_receipt_SHA256")==GLOBAL_RECEIPT_SHA,
            "source_arm":receipt.get("nested_original_photo_depth5_source_arm")==ARM,
            "source_total_photo_IDs":receipt.get(
                "global_source_original_photo_IDs_metadata_only")==EXPECTED_PHOTOS,
            "source_total_nominal_taxa":receipt.get(
                "global_nominal_source_taxa")==EXPECTED_TAXA,
            "U100_eligible_taxa":receipt.get(
                "global_U100_eligible_taxa")==expected_eligible,
            "source_group_photo_ID_uniqueness":receipt.get(
                "source_group_original_photo_IDs_distinct") is True,
            "source_group_observation_uniqueness":receipt.get(
                "source_group_original_observations_distinct") is True,
            "source_group_photo_rank_unchanged":receipt.get(
                "source_group_candidate_photo_rank_unchanged") is True,
            "source_photo_reuse_not_yet_done":receipt.get(
                "original_photo_IDs_already_measured_by_old_FCP_MOIS_and_K1_excluded") is False,
            "no_actual_new_colour_pixels":receipt.get(
                "new_image_download_or_frozen_ROI_inference_this_stage")==0,
        }
        violated=[name for name, ok in checks.items() if not ok]
        require(not violated,
                "Unsigned/modified or repeated full93m original five-photo source taxon group: "
                + ",".join(violated))

        a=folder/PREFIX_FILENAME
        b=folder/DELTA_FILENAME
        require(sha256(a)==receipt["first5_csv_SHA256"]
                and sha256(b)==receipt["ranks2to5_csv_SHA256"],
                "Nested FCP stage5 original photo CSV file bytes differ from source receipt")
        group_delta={}
        for row in read_csv(b):
            pid=row["photo_id"]
            require(pid not in group_delta and row["incremental_since_depth1"]=="true"
                    and row["U100_source_original_observation_rank"] in {"2","3","4","5"},
                    "Source group FCP rank2to5 was not the nested incremental candidate set")
            group_delta[pid]=row
        group_photos=0
        taxa=set()
        for row in read_csv(a):
            taxon=row["inat_taxon_id"]
            pid=row["photo_id"]
            obs=row["observation_id"]
            require(taxon.isdecimal() and pid.isdecimal() and bool(obs)
                    and (int(taxon)%128)//8==group
                    and pid not in ids and obs not in observations,
                    "Exact global nested stage5 photo/observation/taxon ID duplicated")
            rank=int(row["U100_source_original_observation_rank"])
            require(1<=rank<=DEPTH and row["source_arm"]==ARM
                    and row["prefix_depth"]==str(DEPTH)
                    and row["independent_accepted_species_and_native_stage_verified"]=="false"
                    and row["incremental_since_depth1"]==str(rank>1).lower(),
                    "Original source FCP top-five photo rank/technical contract changed")
            nobs=int(row["full_93m_distinct_original_observations_of_taxon"])
            m=int(row["photos_in_original_observation"])
            require(nobs>=100 and m>=1, "Invalid original FCP selected observation N/M")
            exact1=1/nobs/m
            for key,w in (
                ("depth1_joint_original_photo_inclusion_probability",exact1),
                ("rank2to5_incremental_joint_original_photo_probability",4*exact1),
                ("cumulative_prefix5_joint_original_photo_probability",5*exact1),
                ("parent_U100_joint_original_photo_inclusion_probability",100*exact1)):
                require(math.isclose(float(row[key]),w,rel_tol=1e-12),
                        "Original FCP U100 source probability differs across nested ranks")
            if taxon not in ranks_by_taxon:
                ranks_by_taxon[taxon]={}
            require(rank not in ranks_by_taxon[taxon],
                    "Duplicate observation rank within one original FCP source taxon")
            ranks_by_taxon[taxon][rank]=row
            ids.add(pid);observations.add(obs);group_photos+=1
            taxa.add(taxon)
            technical+=row["photo_source_disposition"]!="supported"
            if rank>1:
                require(pid in group_delta and group_delta.pop(pid)==row,
                        "Nested stage5 rank2to5 source candidates differ from first5")
                delta_ids.add(pid)
            else:
                prior=first.get(taxon)
                require(prior and pid==prior["photo_id"]
                        and obs==prior["observation_id"]
                        and math.isclose(float(row["depth1_joint_original_photo_inclusion_probability"]),
                                         float(prior["depth1_joint_original_photo_inclusion_probability"]),
                                         rel_tol=1e-12),
                        "Nested rank1 is not the immutable original U100 first photo")
            require(row["source_rank1_original_photo_id"]==first[taxon]["photo_id"],
                    "The same FCP species switched initial rank1 original photo")
        require(not group_delta and group_photos==
                receipt["first5_observation_independent_original_photo_candidates"]
                and len(taxa)==receipt["source_group_eligible_nominal_taxa"]
                and group_photos==len(taxa)*DEPTH
                and receipt["new_candidate_ranks2to5_NOT_yet_measured"]==len(taxa)*(DEPTH-1),
                "An original FCP source taxon/photo group lost prefix5 or delta members")
        require(taxa.isdisjoint(old_taxa_groups),
                "Original FCP taxon was allocated to multiple source groups")
        old_taxa_groups.update(taxa)
        group_audits[group]={
            "group":group,"taxa":len(taxa),"first5_photo_candidates":group_photos,
            "incremental_rank2to5_candidates":len(taxa)*4,
            "group_receipt_sha256":sha256(receipt_file),
            "prefix5_csv_sha256":sha256(a),
            "delta_csv_sha256":sha256(b),
        }
    require(set(group_audits)==set(range(expected_groups))
            and len(ranks_by_taxon)==expected_eligible
            and set(ranks_by_taxon)==set(first)
            and len(ids)==len(observations)==expected_prefix5
            and len(delta_ids)==expected_incremental
            and all(set(z)=={1,2,3,4,5} for z in ranks_by_taxon.values()),
            "All worldwide original photo candidate sample groups fail five-observation census")
    output.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(tempfile.mkdtemp(prefix="mois-global-u100-stage5-",dir=output.parent))
    try:
        specs=(
            ("global_FCP_U100_nested_first5_original_photo_candidates.csv.gz",range(1,6)),
            ("global_FCP_U100_incremental_rank2to5_candidates_NEED_DEDUP.csv.gz",range(2,6)),
        )
        shas={}
        for name,ranks in specs:
            dest=tmp/name
            with dest.open("wb") as file:
                with gzip.GzipFile(fileobj=file,mode="wb",filename="",mtime=0) as zipped:
                    with io.TextIOWrapper(zipped,encoding="utf-8",newline="") as stream:
                        writer=csv.DictWriter(stream,fieldnames=FIELDS,lineterminator="\n")
                        writer.writeheader()
                        for taxon in sorted(ranks_by_taxon,key=int):
                            for rank in ranks:
                                writer.writerow(ranks_by_taxon[taxon][rank])
            shas[name]=sha256(dest)
        receipt={
            "schema":SCHEMA,
            "complete_source_original_photo_IDs_METADATA_ONLY":EXPECTED_PHOTOS,
            "complete_nominal_source_taxa":EXPECTED_TAXA,
            "eligible_taxa_with_at_least100_original_observations":expected_eligible,
            "independent_original_photo_observations_per_eligible_taxon":DEPTH,
            "global_first5_nested_source_photo_IDs":len(ids),
            "global_incremental_rank2to5_candidate_photo_IDs_NOT_NEW_INFERENCES":len(delta_ids),
            "global_first1_original_photo_IDs_unchanged":len(first),
            "source_original_photo_technical_unsupported_first5":technical,
            "global_prefix1_original_source_SHA256":expected_first1_sha,
            "global_full93m_U100_parent_source_SHA256":GLOBAL_RECEIPT_SHA,
            "group_original_photo_source_SHA_receipts":
                [group_audits[i] for i in range(expected_groups)],
            "csv_outputs_SHA256":shas,
            "ranks2to5_genuine_never_measured_photo_IDs_not_yet_deduplicated":True,
            "source_probabilities_original_5_over_N_over_M_and_delta_4_over_N_over_M":True,
            "new_image_downloads_and_original_ROI_inferences_performed":0,
            "source_taxa_under100_observations_not_selected_but_retained_in_census":
                EXPECTED_TAXA-expected_eligible,
            "unweighted_global_white_flora_or_wild_native_flowering_inferred":False,
        }
        (tmp/"global_FCP_U100_depth5_source_receipt.json").write_text(
            json.dumps(receipt,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(tmp,output)
    except BaseException:
        shutil.rmtree(tmp,ignore_errors=True)
        raise
    print("MOIS_WORLDWIDE_FCP_U100_NESTED_FIRST5_SOURCE_AUTHENTIC",json.dumps({
        "source_eligible_taxa":len(first),
        "first5_photo_candidates":len(ids),
        "rank2to5_original_photo_candidates_before_prior_ID_reuse":len(delta_ids),
        "original_images_downloaded":0,
        "four_hues_inferred":0,
    },sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--first1",required=True,type=Path)
    p.add_argument("--group",action="append",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    global_source(a.group,a.first1,a.output)


if __name__=="__main__":
    main()
