#!/usr/bin/env python3
"""Target only source taxa with ZERO strict original-photo ROI9, from signed 93m metadata.

One full original-photo metadata group (8 taxon-modulo parts x 14 tiers)
is scanned exactly once, selecting at most K original observations per
target taxon. No image download, colour imputation, or GBIF inference.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import tarfile
from collections import Counter
from pathlib import Path

import pandas as pd
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, EXPECTED_PHOTOS, EXPECTED_TAXA, RESULT_SCHEMA, TIER_ENDS,
    normalized_extension, obs_draw, sha256,
)

SCHEMA="mois_species_first_original93m_zero9_taxon_reverse_lookup_v1"
PHOTO_SHA="3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835"
COVERAGE_SHA="ea602238f3cd152fdd2e36e9ce4bcbc77b795a35dce802e9305f6f6cddd77fb8"
NUM_PARTS=8

def check(ok, reason):
    if not bool(ok):
        raise ValueError("MOIS_BREADTH_FIRST_ZERO9: "+reason)

def image_rank(uuid, pid):
    return hashlib.sha256(("mois-species-first-photo-blind-v1|"+uuid+"|"+pid).encode()).hexdigest()

def source_rank(row):
    return (float(row["observation_inclusion_score"]),row["observation_uuid"])

def keep_best(candidates, row, *, limit):
    """Bounded memory: K lowest observation hash ranks, one source photo/UUID."""
    obs=row["observation_uuid"]
    pid=row["photo_id"]
    specimen=candidates.setdefault(row["inat_taxon_id"], {})
    item=specimen.get(obs)
    if item is not None:
        if image_rank(obs,pid)<image_rank(obs,item["photo_id"]):
            specimen[obs]=row
        return
    if len(specimen)<limit:
        specimen[obs]=row
        return
    worst=max(specimen, key=lambda z:source_rank(specimen[z]))
    if source_rank(row)<source_rank(specimen[worst]):
        del specimen[worst]
        specimen[obs]=row

def select_rows(rows, wanted, seen_photos, seen_obs, *, limit=4):
    """Unit-testable colour-blind source selector, no species fixed-colour claims."""
    result={}
    stats=Counter()
    for row in rows:
        stats["source_rows"]+=1
        taxon=str(row["inat_taxon_id"])
        if taxon not in wanted:continue
        stats["wanted_rows"]+=1
        pid=str(row["photo_id"])
        oid=str(row["observation_uuid"])
        check(pid.isdecimal() and bool(oid),"Malformed original photo/observation ID")
        if pid in seen_photos or oid in seen_obs:
            stats["already_measured_or_observed"]+=1
            continue
        if str(row.get("frozen_loader_supported","")).lower()!="true":
            stats["unsupported_image_extension"]+=1
            continue
        if str(row.get("photo_license","")).strip().lower().startswith("cc") is False:
            stats["missing_nonCC_photo_licence"]+=1
            continue
        ext=normalized_extension(row.get("image_extension",""))
        if not ext:
            stats["unsupported_image_extension"]+=1
            continue
        # An observation in the official signed source has a stable score
        # across all photos of that observation.
        check(abs(float(row["observation_inclusion_score"])-obs_draw(oid))<1e-15,
              "Observation hash rank incompatible with original 93m source")
        payload={k:str(v) for k,v in row.items()}
        payload["inat_taxon_id"]=taxon
        payload["photo_id"]=pid
        payload["original_image_extension"]=ext
        keep_best(result,payload,limit=limit)
        stats["candidate_photo_rows"]+=1
    return result,stats

def group_source(registry,receipt,archives,group,wanted,seen_photo,seen_obs,*,limit=4):
    """SHA-validate every frozen tier/part in a single source group."""
    full=json.loads(Path(registry).read_text(encoding="utf-8"))
    g=json.loads(Path(receipt).read_text(encoding="utf-8"))
    check(full.get("schema")=="mois_complete_full93m_original_photo_target_registry_v1"
          and full.get("exact_original_photo_IDs_in_all_disjoint_tiers")==EXPECTED_PHOTOS
          and full.get("original_named_taxa_from_full_source")==EXPECTED_TAXA
          and full.get("source_only_NO_photo_pixels_or_hue_measured") is True,
          "Original 93m entire source register not verified")
    check(g.get("schema")==RESULT_SCHEMA
          and g.get("source_snapshot_SHA256")==full.get("full_source_manifest_sha256")
          and g.get("all_source_photo_IDs_retained_in_one_and_only_one_tier") is True,
          "Signed taxon group does not match full source")
    parts={int(p["part"]):p for p in g["parts"]}
    expected_parts=set(range(group*NUM_PARTS,(group+1)*NUM_PARTS))
    check(set(parts)==expected_parts, "Group taxon-part identity mismatch")
    candidates={}
    stats=Counter()
    for tier in range(len(TIER_ENDS)):
        name=f"tier-{tier:02d}-group-{group:02d}.tar"
        tar_path=Path(archives)/name
        check(tar_path.is_file() and
              sha256(tar_path)==full["published_tier_archive_SHA256"][name],
              "Signed original source tier TAR SHA mismatch: "+name)
        expected_members={f"tier-{tier:02d}/part-{part:03d}.csv.gz" for part in expected_parts}
        with tarfile.open(tar_path,"r") as source_tar:
            members={m.name:m for m in source_tar.getmembers() if m.isfile()}
            check(set(members)==expected_members and
                  all(m.isfile() or m.isdir() for m in source_tar.getmembers()),
                  "Signed original tier TAR has unexpected/missing source parts")
            for part in sorted(expected_parts):
                member=f"tier-{tier:02d}/part-{part:03d}.csv.gz"
                raw=source_tar.extractfile(members[member]).read()
                check(hashlib.sha256(raw).hexdigest()==
                      parts[part]["ranked_target_csv_sha256"][str(tier)],
                      "Signed original photo CSV part SHA mismatch")
                current=0
                with gzip.GzipFile(fileobj=io.BytesIO(raw)) as z:
                    with io.TextIOWrapper(z,encoding="utf-8",newline="") as f:
                        rd=csv.DictReader(f)
                        check(tuple(rd.fieldnames or ())==COLUMNS,
                              "Signed photo-row schema mutated")
                        for row in rd:
                            check(row["tier"]==str(tier) and
                                  int(row["inat_taxon_id"])%128//8==group,
                                  "Signed taxon/tier metadata contradiction")
                            current+=1
                            stats["source_rows"]+=1
                            taxon=row["inat_taxon_id"]
                            if taxon not in wanted:continue
                            stats["wanted_rows"]+=1
                            if (row["photo_id"] in seen_photo or
                                row["observation_uuid"] in seen_obs):
                                stats["already_measured_or_observed"]+=1
                                continue
                            if (str(row["frozen_loader_supported"]).lower()!="true"
                                    or not normalized_extension(row["image_extension"])):
                                stats["unsupported_image_extension"]+=1
                                continue
                            if not str(row["photo_license"]).strip().lower().startswith("cc"):
                                stats["missing_nonCC_photo_licence"]+=1
                                continue
                            oid=row["observation_uuid"]
                            check(abs(float(row["observation_inclusion_score"])-obs_draw(oid))<1e-15,
                                  "Source observation hash rank incompatible with source")
                            payload={k:str(v) for k,v in row.items()}
                            payload["original_image_extension"]=normalized_extension(row["image_extension"])
                            keep_best(candidates,payload,limit=limit)
                            stats["candidate_photo_rows"]+=1
                check(current==int(parts[part]["counts_by_disjoint_tier"][str(tier)]),
                      "Signed source tier part count mismatch")
    check(stats["source_rows"]==int(g["original_photo_IDs_in_this_group"]),
          "Source group photo row total mismatch")
    return candidates,stats

def output_rows(candidates):
    rows=[]
    for taxon in sorted(candidates,key=lambda t:int(t)):
        chosen=sorted(candidates[taxon].values(),key=source_rank)
        for rank,row in enumerate(chosen,1):
            ext=row["original_image_extension"]
            record={k:row[k] for k in (
                "photo_id","observation_uuid","inat_taxon_id","observer_id",
                "observed_on","latitude","longitude","photo_license",
                "image_extension","observation_inclusion_score","tier")}
            record.update(
                observation_id=row["observation_uuid"],
                source_observation_rank=rank,
                photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{row['photo_id']}/large.{ext}",
                source_arm="full93m_zero_strict9_species_first",
                measured_yet="false")
            rows.append(record)
    return pd.DataFrame(rows,columns=[
        "photo_id","observation_uuid","inat_taxon_id","observer_id",
        "observed_on","latitude","longitude","photo_license","image_extension",
        "observation_inclusion_score","tier","observation_id",
        "source_observation_rank","photo_url_large","source_arm","measured_yet"])

def run(registry,receipt,archives,photos,coverage,output,*,group=0,max_fallback=4,check_sha=True):
    output=Path(output)
    check(not output.exists() and 0<=group<16 and 1<=max_fallback<=10,
          "Unsupported group/selection size or output already exists")
    if check_sha:
        check(sha256(photos)==PHOTO_SHA and sha256(coverage)==COVERAGE_SHA,
              "Frozen actual 1.18m source-photo ROI9 source SHA mismatch")
    frame=pd.read_csv(coverage,compression="infer",dtype={"inat_taxon_id":str})
    check(len(frame)==EXPECTED_TAXA and frame.inat_taxon_id.is_unique,
          "Source 42111 nominal taxon registry incomplete")
    wanted=set(frame.loc[(frame.after_strict9==0)&
            (frame.inat_taxon_id.map(lambda s:int(s)%128//8==group))&
            (frame.n_unique_source_photo_ids>0),"inat_taxon_id"])
    prior=pd.read_parquet(photos,columns=["photo_id","observation_id"])
    prior_photo=set(prior.photo_id.fillna("").astype(str))
    prior_obs=set(prior.observation_id.fillna("").astype(str))-{""}
    check(len(prior)==len(prior_photo)==1181146,
          "Frozen 1,181,146 original photo IDs not unique")
    selections,stats=group_source(registry,receipt,archives,group,wanted,prior_photo,
                                  prior_obs,limit=max_fallback)
    q=output_rows(selections)
    next_wave=q.loc[q.source_observation_rank==1].copy()
    check(next_wave.inat_taxon_id.is_unique and next_wave.photo_id.is_unique
          and next_wave.observation_id.is_unique
          and not set(next_wave.photo_id)&prior_photo,
          "Non-unique new original photo/observation queue")
    output.mkdir(parents=True)
    allpath=output/f"mois_species_zero9_group{group:02d}_source_candidates.csv.gz"
    nextpath=output/f"mois_species_zero9_group{group:02d}_next_first_photo.csv"
    q.to_csv(allpath,index=False,compression="gzip")
    next_wave.to_csv(nextpath,index=False)
    missing=sorted(wanted-set(selections),key=int)
    (output/"no_source_candidate_taxa.json").write_text(json.dumps(missing)+"\n")
    report={"schema":SCHEMA,"group":group,"source_registry_sha":sha256(registry),
            "group_receipt_sha":sha256(receipt),
            "already_measured_photo_ledger_sha":sha256(photos),
            "strict9_taxon_coverage_sha":sha256(coverage),
            "source_zero9_taxa_with_any_photo":len(wanted),
            "source_zero9_taxa_with_eligible_candidate":len(selections),
            "source_zero9_taxa_no_eligible_candidate":len(missing),
            "new_first_wave_distinct_species_photo_opportunities":len(next_wave),
            "total_selected_candidates_including_fallbacks":len(q),
            "observed_original_photo_rows_in_source_group":stats["source_rows"],
            "source_only_technical_measurement":True,
            "new_original_photo_downloads":0,
            "new_measured_9palette_ROIs":0,
            "GBIF_occurrences_downloaded":0,
            "output_sha256":{p.name:sha256(p) for p in (allpath,nextpath)}}
    (output/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_SPECIES_FIRST_ZERO9_FULL_SOURCE",json.dumps(report,sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("registry","group-receipt","archives","photos","coverage","output"):
        p.add_argument("--"+key,required=True)
    p.add_argument("--group",type=int,required=True)
    p.add_argument("--max-fallback",type=int,default=4)
    a=p.parse_args()
    run(a.registry,a.group_receipt,a.archives,a.photos,a.coverage,a.output,
        group=a.group,max_fallback=a.max_fallback)

if __name__=="__main__":
    main()
