#!/usr/bin/env python3
"""Audit taxonomic coverage of one million authentic original-image MOIS photos.

The 93,209,351 upstream photo IDs and 42,111 source taxon IDs are source
metadata, not automatically flowering/colour observations. This pass joins
independently SHA-audited 929,138 historic actual ROI photo attempts, all
76,753 completed K1 source photo attempts and the 3,141 NEVER-before-inferred
FCP U100 original photos with the complete 42,111-taxon eligible/rare census.

It prioritizes **missing rare (<100 independent observations) source taxa**
for the NEXT photo-acquisition wave, rather than naively recolouring millions
of already-measured images or equating model technical unknown with white.
No source photo fetch, colour inference or native-flora claim.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path
import os
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, sha, require
from scripts.acquisition.mois_93m_balanced_prior_ROI_reuse_bridge import (
    PRIOR_MAIN,PRIOR_OLD,PRIOR_R1,PRIOR_R2,prior_original_photo_index,
)
from scripts.acquisition.mois_fcp_u100_collect_new_first_photo_ROIs import (
    SCHEMA as FCP_SCHEMA,
)

SCHEMA="mois_42111_taxon_original_photo_true_pixel_depth_coverage_v1"
CENSUS_SHA="6adab400c72932ab1079cdb66632018e90b14f20af2aaf53d6a8521dd8eb8441"
CENSUS_HEADER=(
    "inat_taxon_id","species","full_93m_original_observation_count",
    "full_93m_original_photo_count","FCP_U100_source_eligible",
    "FCP_U100_missing_observations","selected_observed_months",
    "selected_latitude_bands",
)
K1_PARQUET_SHA="0d3de3374a06b1e7b9ad0281e2ab7844632ed139db871488524061191c8f34c1"
K1_RECEIPT_SHA="e86463bc61f9eebd8cbdbe01269f767c140bf243f9e374027008e9f095eb9044"
FCP_PARQUET_SHA="ee79c60671dda65a3d2f4f9319af952dbe81ab805ce269fe32b3bb3448b0d6b8"
FCP_RECEIPT_SHA="5dbd4590ae96729cd9469fc8403a8ee4528d18e1fb5b928d5a108483ea29150b"
OLD_COUNTS={"main":PRIOR_MAIN, "prior_group00":PRIOR_OLD,
            "prior_R1":PRIOR_R1, "prior_R2":PRIOR_R2}

def source_taxon_eligible(nobs,nphotos,flag):
    require(type(nobs) is int and nobs>=1
            and type(nphotos) is int and nphotos>=nobs
            and flag==("true" if nobs>=100 else "false"),
            "Original full93m taxon independent-observation eligibility mutated")
    return nobs>=100


def taxon_coverage_record(taxon,species,nobs,nphotos,flag,attempted,classified):
    source_taxon_eligible(nobs,nphotos,flag)
    require(taxon.isdecimal() and type(attempted) is int
            and type(classified) is int
            and 0<=classified<=attempted,
            "Attempted/classifiable taxon original image denominators invalid")
    return {
        "inat_taxon_id":taxon,"species":species,
        "source_original_distinct_observations":nobs,
        "source_original_photo_IDs":nphotos,
        "U100_eligible_100plus_observations":flag,
        "genuine_original_photo_attempts":attempted,
        "automatic_four_hue_ROIs":classified,
        "acquired_or_evaluated_but_hue_unknown":attempted-classified,
        "has_any_photo_attempt":str(attempted>0).lower(),
        "has_any_four_hue_ROI":str(classified>0).lower(),
        "source_taxon_under100_observations":str(nobs<100).lower(),
    }


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def run(census, historical, k1, fcp, output, *,census_sha=CENSUS_SHA,
        k1_parquet_sha=K1_PARQUET_SHA,k1_receipt_sha=K1_RECEIPT_SHA,
        fcp_parquet_sha=FCP_PARQUET_SHA,fcp_receipt_sha=FCP_RECEIPT_SHA,
        expected_taxa=42111,expected_eligible=20242,expected_photo_source=93209351,
        expected_total=1009032,expected_hues=445805):
    import pandas as pd
    census,historical,k1,fcp,output=map(Path,(census,historical,k1,fcp,output))
    require(not output.exists(),"Refuse overwrite of prior source taxon pixel census")
    require(sha(census)==census_sha,"Full93m original 42111-taxon source population SHA drift")
    taxa={}
    raw_photos=raw_obs=eligible=0
    with census.open(newline="",encoding="utf-8") as fh:
        r=csv.DictReader(fh)
        require(tuple(r.fieldnames or ())==CENSUS_HEADER,
                "Original FCP U100 source census taxonomic grain changed")
        for row in r:
            taxon=row["inat_taxon_id"]
            require(taxon.isdecimal() and taxon not in taxa,
                    "Repeated original source nominal taxon ID")
            nobs=int(row["full_93m_original_observation_count"])
            nphotos=int(row["full_93m_original_photo_count"])
            flag=row["FCP_U100_source_eligible"]
            source_taxon_eligible(nobs,nphotos,flag)
            taxa[taxon]=(row["species"],nobs,nphotos,flag)
            raw_photos+=nphotos
            raw_obs+=nobs
            eligible+=nobs>=100
    require(len(taxa)==expected_taxa and eligible==expected_eligible
            and raw_photos==expected_photo_source,
            "Full globe taxon denominator or complete original photo IDs missing")
    old_sources=[
        ("main",historical/"main"/"combined_actual_photo_checkpoint.parquet",
         historical/"main"/"checkpoint_receipt.json"),
        ("prior_group00",historical/"group00"/"measured_photos.parquet",
         historical/"group00"/"receipt.json"),
        ("prior_R1",historical/"r1"/"measured_photos.parquet",
         historical/"r1"/"receipt.json"),
        ("prior_R2",historical/"r2"/"measured_photos.parquet",
         historical/"r2"/"receipt.json"),
    ]
    old,old_obs,sources,overlap,_,_=prior_original_photo_index(
        old_sources,OLD_COUNTS)
    require(len(old)==929138 and overlap["unique_old_photo_IDs"]==929138
            and not overlap.get("old_cross_arm_same_photo_evidence_conflicts",0),
            "Historical actual photograph identity/source coverage changed")
    seen=set(old)
    per_taxon=defaultdict(lambda:[0,0])
    unknown_taxon_ids=Counter()
    unknown_taxon_images=Counter()
    by_arm={}
    for arm,path,_receipt in old_sources:
        columns=["photo_id","morph"]
        import pyarrow.parquet as pq
        available=set(pq.ParquetFile(path).schema.names)
        if "inat_taxon_id" in available:columns.append("inat_taxon_id")
        df=pd.read_parquet(path,columns=columns)
        require(len(df)==OLD_COUNTS[arm][0] and int(df.morph.isin(HUES).sum())==OLD_COUNTS[arm][1],
                "Historically attempted photo RGB colour count changed")
        if "inat_taxon_id" not in df:
            unknown_taxon_images[arm]=len(df)
        else:
            df["inat_taxon_id"]=df.inat_taxon_id.fillna("").astype(str)
            for taxon,d in df.groupby("inat_taxon_id",dropna=False,sort=False):
                if taxon not in taxa:
                    unknown_taxon_images[arm]+=len(d)
                    unknown_taxon_ids[arm]+=1
                else:
                    per_taxon[taxon][0]+=len(d)
                    per_taxon[taxon][1]+=int(d.morph.isin(HUES).sum())
        by_arm[arm]={"original_photo_attempts":len(df),
                     "automated_four_hue_ROIs":int(df.morph.isin(HUES).sum()),
                     "within_global_42111_source_taxa":len(df)-unknown_taxon_images[arm],
                     "outside_source_taxon_or_id_missing":unknown_taxon_images[arm]}
    def audit_new_arm(root,name,n,hue,pin_parquet,pin_receipt,schema):
        parquet=root/("measured_photos.parquet" if name=="K1" else
                       "new_FCP_depth1_original_photo_ROIs.parquet")
        receipt=root/"receipt.json"
        require(sha(parquet)==pin_parquet and sha(receipt)==pin_receipt,
                "Completed K1/FCP actual original ROI source SHA mutated: "+name)
        report=read_json(receipt)
        if name=="K1":
            require(report.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
                    and report.get("attempted")==n and report.get("classifiable")==hue
                    and report.get("matching_prior_measured_photos_reused")==3280
                    and report.get("measured_photo_parquet_sha256")==pin_parquet,
                    "Full K1 original image photo-ID source receipt altered")
        else:
            require(report.get("schema")==schema
                    and report.get("newly_true_image_measured_original_photo_attempts")==n
                    and report.get("newly_true_image_FCP_four_hue_classified")==hue
                    and report.get("actual_new_image_parquet_SHA256")==pin_parquet,
                    "Original FCP-U100 new photo source receipt altered")
        df=pd.read_parquet(parquet,
            columns=["photo_id","inat_taxon_id","morph"])
        photo=df.photo_id.fillna("").astype(str)
        require(len(photo)==n and photo.is_unique and
                photo.str.fullmatch("[0-9]+").all()
                and not (set(photo)&seen)
                and int(df.morph.isin(HUES).sum())==hue,
                "Original K1/FCP image photo-ID/hue overlaps with historical pixels")
        seen.update(photo)
        group=df.assign(inat_taxon_id=df.inat_taxon_id.fillna("").astype(str))
        for taxon,part in group.groupby("inat_taxon_id",sort=False):
            if taxon not in taxa:
                unknown_taxon_images[name]+=len(part)
                unknown_taxon_ids[name]+=1
            else:
                per_taxon[taxon][0]+=len(part)
                per_taxon[taxon][1]+=int(part.morph.isin(HUES).sum())
        by_arm[name]={"original_photo_attempts":n,"automated_four_hue_ROIs":hue,
                      "within_global_42111_source_taxa":n-unknown_taxon_images[name],
                      "outside_source_taxon_or_id_missing":unknown_taxon_images[name]}
    audit_new_arm(k1,"K1",76753,27689,k1_parquet_sha,k1_receipt_sha,None)
    audit_new_arm(fcp,"FCP_U100_new",3141,1201,fcp_parquet_sha,fcp_receipt_sha,FCP_SCHEMA)
    require(len(seen)==expected_total
            and sum(v["original_photo_attempts"] for v in by_arm.values())==expected_total
            and sum(v["automated_four_hue_ROIs"] for v in by_arm.values())==expected_hues,
            "MOIS real photographic global source-ID/hue expected total changed")
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-42111-pixel-depth-",dir=output.parent))
    try:
        fields=["inat_taxon_id","species","source_original_distinct_observations",
                "source_original_photo_IDs","U100_eligible_100plus_observations",
                "genuine_original_photo_attempts","automatic_four_hue_ROIs",
                "acquired_or_evaluated_but_hue_unknown","has_any_photo_attempt",
                "has_any_four_hue_ROI","source_taxon_under100_observations"]
        full=stage/"all_42111_nominal_source_taxon_actual_original_ROI_coverage.csv"
        missing=stage/"under100_observation_taxa_without_automatic_four_hue_ROI.csv"
        metric=Counter()
        with (full.open("w",encoding="utf-8",newline="") as out,
              missing.open("w",encoding="utf-8",newline="") as out_missing):
            writer=csv.DictWriter(out,fieldnames=fields,lineterminator="\n")
            writer.writeheader()
            gap=csv.DictWriter(out_missing,fieldnames=fields,lineterminator="\n")
            gap.writeheader()
            for taxon,(species,nobs,nphotos,flag) in sorted(taxa.items(),key=lambda z:int(z[0])):
                attempted,classified=per_taxon[taxon]
                row=taxon_coverage_record(
                    taxon,species,nobs,nphotos,flag,attempted,classified)
                writer.writerow(row)
                metric["taxa_in_source"]+=1
                metric["taxa_at_least_one_photo_attempt"]+=attempted>0
                metric["taxa_at_least_one_auto_four_hue_ROI"]+=classified>0
                if nobs<100:
                    metric["under100_observation_source_taxa"]+=1
                    metric["under100_with_any_measured_photo"]+=attempted>0
                    metric["under100_with_any_four_hue_ROI"]+=classified>0
                    if not classified:gap.writerow(row)
                else:
                    metric["U100_eligible_source_taxa"]+=1
                    metric["U100_eligible_with_any_measured_photo"]+=attempted>0
                    metric["U100_eligible_with_any_four_hue_ROI"]+=classified>0
        info={
            "schema":SCHEMA,
            "full_upstream_original_photo_IDs_metadata_only":expected_photo_source,
            "full_nominal_source_taxa":expected_taxa,
            "full_source_distinct_observations":raw_obs,
            "true_genuine_original_photo_attempts":expected_total,
            "true_automatic_four_hue_ROIs":expected_hues,
            "source_taxon_census_counts":dict(metric),
            "actual_source_photo_arm_counts":by_arm,
            "source_taxon_IDs_missing_from_some_original_ROI_arms":dict(unknown_taxon_images),
            "not_assuming_historical_photo_taxa_all_match_current_source":True,
            "original_photo_IDs_deduplicated_within_all_arms":True,
            "source_taxa_with_zero_four_hue_photos_prioritized_without_colour_filter":True,
            "true_open_native_wild_flower_verification_not_performed":True,
            "unweighted_world_white_flora_frequency_NOT_estimated":True,
            "source_census_SHA256":census_sha,
            "actual_historical_photo_source_SHA_receipts":sources,
            "all_taxon_coverage_csv_SHA256":sha(full),
            "under100_zero_classified_taxa_csv_SHA256":sha(missing),
            "photo_bytes_fetched_OR_hue_inferred_by_this_analysis":0,
        }
        (stage/"mois_full42111_original_pixel_taxon_coverage_receipt.json").write_text(
            json.dumps(info,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_FULL93M_SOURCE_TAXON_TRUE_PIXEL_BREADTH_GAPS",json.dumps({
        "source_nominal_taxa":expected_taxa,
        "global_actual_original_photo_attempts":expected_total,
        "global_auto_four_hue_ROIs":expected_hues,
        "taxa_ever_auto_four_hue":metric["taxa_at_least_one_auto_four_hue_ROI"],
        "rare_source_taxa_under100":metric["under100_observation_source_taxa"],
        "rare_source_taxa_without_any_four_hue":
            metric["under100_observation_source_taxa"]-metric["under100_with_any_four_hue_ROI"],
        "source_missing_taxon_join_photos":dict(unknown_taxon_images),
        "original_images_downloaded":0,
    },sort_keys=True),flush=True)
    return info

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("census","historical","k1","fcp","output"):
        p.add_argument("--"+key,required=True,type=Path)
    a=p.parse_args()
    run(a.census,a.historical,a.k1,a.fcp,a.output)

if __name__=="__main__":
    main()
