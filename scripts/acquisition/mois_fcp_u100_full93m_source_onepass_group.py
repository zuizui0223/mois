#!/usr/bin/env python3
"""FCP-U100 exactly 100 independent original photo observations per taxon from ALL 93m source IDs.

Run ONE independent taxon-modulo source group (eight of 128 original source
partitions). All 14 disjoint pre-colour original photo-ID tiers are SHA checked
and read exactly once. One SQLite table joins photographs to their original
observation UUID and keeps EVERY original observation and number of photos.
For each source-eligible taxon (>=100 original distinct observation UUIDs),
select exactly 100 by the SAME original observation-hash rank that defined the
pre-existing tier partitions, and choose one independently hashed original
photo per selected observation. 100 selected observations really are the
global rank-100 prefix of the entire 93,209,351-photo original source—not
one truncated tier, one calendar season or 100 random photographs possibly
from the same source observation.

Conditional on all original source observation opportunities for taxon i:
  P(original observation in top100) = 100/N_i
  P(original photo | selected observation) = 1/M_observation
  P(original photo in FCP-U100 frame | eligible taxon) = 100/N_i/M_observation
Each taxon's N_i and every selected observation's M are EXACT because all 14
SHA-signed tiers are read. For taxa with fewer than 100 original observations,
no fake 100-photo source is created: they stay in the complete source roster.

FCP already contains 4475 taxa * 100 authentic measured photo outcomes.
Before ANY later full MOIS original ROI-v4 measurement, independently join
FCP/MOIS genuine prior-photo IDs and reuse their existing colour outcomes.
This source step NEVER downloads original images or infers four hues. It is
NOT H9 phenology-photo authentication, a native-flora census or a globally
probability-complete 42,111-taxon sample (noneligible species have 0 in U100).
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import sqlite3
import tarfile
import tempfile

from scripts.acquisition import mois_93m_taxon_month_lat_balanced_photo_sampler as source
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, EXPECTED_PHOTOS, EXPECTED_TAXA, RESULT_SCHEMA, SOURCE_SEED,
    TIER_ENDS, require, sha256, normalized_extension,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import (
    WORKER_COLUMNS, read_taxa,
)

SCHEMA="mois_93m_fcp_u100_full_source_14tier_exact100_unique_original_observation_photos_per_taxon_v1"
FCP_TARGET=100
ALL_TAXA=EXPECTED_TAXA
N_GROUPS=16
PARTS_PER_GROUP=8
TIERS=len(TIER_ENDS)
SOURCE_ARM="mois_full93m_fcp_exact100_original_observation_per_taxon_v1"
PREFSALT="mois-original-observation-rank-photo-first-v1"

def original_observation_priority(uuid):
    # SAME hash key as the source tiers; no resampling a different population.
    return hashlib.sha256((SOURCE_SEED+"|"+str(uuid)).encode()).digest()[:8].hex()

def source_group_taxa(taxon_names, group):
    return {taxon:species for taxon,species in taxon_names.items()
            if int(taxon)%128//PARTS_PER_GROUP==group}

def attach_source_photo(db, row, tier, group):
    uuid=row["observation_uuid"]
    taxon=row["inat_taxon_id"]
    pid=row["photo_id"]
    require(uuid and pid.isdecimal() and taxon.isdecimal()
            and (int(taxon)%128)//PARTS_PER_GROUP==group
            and row["tier"]==str(tier),
            "Original full93m source observation, photo, taxon or tier mismatched")
    rank=source.hash64(source.SALT_PHOTO,uuid+"|"+pid)
    data=(uuid,taxon,source.month_from_recorded_date(row["observed_on"]),
          source.latitude_band(row["latitude"]),
          original_observation_priority(uuid),
          rank,pid,row["photo_license"],row["image_extension"],
          row["observed_on"],row["latitude"],row["longitude"],
          row["observer_id"],row["positional_accuracy"])
    c=db.execute(source.UPSERT,data)
    require(c.rowcount==1,
            "Original source photo attaches conflicting metadata to same observation UUID")

def load_original_group(db, registry, receipt, tiers_dir, group,
                        *, expected_total=EXPECTED_PHOTOS, tiers=TIERS):
    """Audit original 14-tier/8-part SHA coverage and append photos exactly once."""
    reg=json.loads(Path(registry).read_text(encoding="utf-8"))
    g=json.loads(Path(receipt).read_text(encoding="utf-8"))
    require(reg.get("schema")=="mois_complete_full93m_original_photo_target_registry_v1"
            and reg.get("exact_original_photo_IDs_in_all_disjoint_tiers")==expected_total
            and reg.get("original_named_taxa_from_full_source")==EXPECTED_TAXA
            and reg.get("source_only_NO_photo_pixels_or_hue_measured") is True
            and reg.get("native_vs_introduced_not_hard_filtered") is True,
            "Original 93m complete prior-image photo-ID registry changed")
    require(g.get("schema")==RESULT_SCHEMA
            and g.get("source_snapshot_SHA256")==reg["full_source_manifest_sha256"]
            and g.get("source_original_photo_IDs_uniqueness_proved_by_upstream_primary_key") is True
            and g.get("all_source_photo_IDs_retained_in_one_and_only_one_tier") is True
            and g.get("same_observation_uses_same_selection_score") is True
            and g.get("images_downloaded")==0
            and g.get("colour_pixels_measured")==0,
            "Original pre-colour full source taxon-group source photo identity absent")
    expected_parts=set(range(group*PARTS_PER_GROUP,(group+1)*PARTS_PER_GROUP))
    parts={int(x["part"]):x for x in g["parts"]}
    require(set(parts)==expected_parts,
            "Original source taxon-modulo eight source partitions missing or duplicated")
    for n in expected_parts:
        require(int(parts[n]["original_photo_IDs_in_part"])==
                sum(int(v) for v in parts[n]["counts_by_disjoint_tier"].values()),
                "Original signed source photo-ID part count changes across 14 tiers")
    count=0;per_tier=Counter()
    for tier in range(tiers):
        name=f"tier-{tier:02d}-group-{group:02d}.tar"
        archive=Path(tiers_dir)/name
        correct=reg["published_tier_archive_SHA256"].get(name)
        require(correct and sha256(archive)==correct,
                "Original complete93m source tier archive checksum changed: "+name)
        if "archive_sha256_by_tier" in g:
            require(g["archive_sha256_by_tier"].get(name)==correct,
                    "Source group full14tier SHA manifest contradicts immutable master")
        expected_members={f"tier-{tier:02d}/part-{p:03d}.csv.gz" for p in expected_parts}
        with tarfile.open(archive,"r:*") as tf:
            members={}
            for item in tf.getmembers():
                require(item.isfile() or item.isdir(),
                        "Non-original source archive entry cannot become a plant photo")
                if item.isfile():
                    path=item.name.removeprefix("./")
                    require(path in expected_members and path not in members,
                            "Unexpected or duplicated original source photo data member")
                    members[path]=item
            require(set(members)==expected_members,
                    "Original source archive lacks one of eight photo parts")
            for part in sorted(expected_parts):
                filename=f"tier-{tier:02d}/part-{part:03d}.csv.gz"
                raw=tf.extractfile(members[filename]).read()
                require(hashlib.sha256(raw).hexdigest()==
                        parts[part]["ranked_target_csv_sha256"][str(tier)],
                        "Genuine signed original photo source part SHA differs")
                n=0
                with gzip.GzipFile(fileobj=io.BytesIO(raw),mode="rb") as zipped:
                    with io.TextIOWrapper(zipped,encoding="utf-8",newline="") as stream:
                        reader=csv.DictReader(stream)
                        require(tuple(reader.fieldnames or ())==COLUMNS,
                                "Original source photo metadata column order changed")
                        for row in reader:
                            attach_source_photo(db,row,tier,group)
                            n+=1
                require(n==int(parts[part]["counts_by_disjoint_tier"][str(tier)]),
                        "Photo source-ID count mismatch after SHA-locked part decompression")
                count+=n
                per_tier[str(tier)]+=n
            db.commit()
    require(count==int(g["original_photo_IDs_in_this_group"])
            and count==sum(int(x) for x in g["counts_by_disjoint_tier"].values()),
            "Full14 original source photoID group coverage was lost")
    return {"group":group,"original_photos":count,
            "original_tier_counts":dict(per_tier),
            "original_group_receipt_SHA256":sha256(receipt),
            "original_full_93m_source_registry_SHA256":sha256(registry)}

def select_exact100_from_full_source(db, taxon_names):
    """Choose genuine top100 observation-hash ranks for all U100 eligible taxa."""
    inv={}
    for taxon,n,photos in db.execute(
        "SELECT taxon,COUNT(*),SUM(photo_count) FROM obs GROUP BY taxon ORDER BY taxon"):
        require(taxon not in inv and taxon in taxon_names
                and n>0 and photos>=n,
                "Genuine original full 93m per-taxon observation counts invalid")
        inv[taxon]=(n,photos)
    qualified={taxon:n for taxon,(n,_) in inv.items() if n>=FCP_TARGET}
    query="""
    SELECT taxon,uuid,photo_count,best_photo_id,photo_license,
           photo_extension,recorded_on,record_month,lat_band,
           latitude,longitude,observer_id,positional_accuracy,
           obs_rank,
           ROW_NUMBER() OVER(PARTITION BY taxon ORDER BY obs_rank,uuid) AS rn
    FROM obs ORDER BY taxon,obs_rank,uuid
    """
    selected=[];month_support=defaultdict(set);lat_support=defaultdict(set)
    for (taxon,uuid,m,pid,license_code,extension,on,month,band,lat,lon,
         observer_id,accuracy,_rank,rank) in db.execute(query):
        if taxon not in qualified or rank>FCP_TARGET:continue
        n=qualified[taxon]
        require(1<=rank<=FCP_TARGET and m>0 and pid.isdecimal(),
                "A frozen original photo opportunity missing from top100 taxonomy")
        ext=normalized_extension(extension)
        p=FCP_TARGET/n/m
        require(0<p<=1,"Conditional original photo inclusion probability not in (0,1]")
        row={
            "photo_id":pid,"observation_id":uuid,"inat_taxon_id":taxon,
            "species":taxon_names[taxon],"observer_id":observer_id,
            "observed_on":on,"recorded_month":month,"latitude_band":band,
            "latitude":lat,"longitude":lon,"positional_accuracy":accuracy,
            "photo_license":license_code,"original_image_extension":extension,
            "photos_in_original_observation":m,
            "full_93m_distinct_original_observations_of_taxon":n,
            "conditional_species_observation_probability":format(FCP_TARGET/n,".17g"),
            "conditional_one_photo_probability":format(1/m,".17g"),
            "conditional_full_93m_photo_selection_probability":format(p,".17g"),
            "photo_source_disposition":"supported" if ext else "unsupported_original_extension",
            "photo_url_large":(
                f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}"
                if ext else ""),
            "source_arm":SOURCE_ARM,
        }
        selected.append(row)
        month_support[taxon].add(month)
        lat_support[taxon].add(band)
    require(len(selected)==FCP_TARGET*len(qualified)
            and len({x["photo_id"] for x in selected})==len(selected)
            and len({x["observation_id"] for x in selected})==len(selected),
            "FCP source 100 original observation-independent photos per eligible species not proven")
    return inv,qualified,selected,month_support,lat_support

def execute_group(registry,receipt,tiers_dir,taxa_file,output,group,
                  *, expected_total=EXPECTED_PHOTOS,check_taxa_sha=True):
    require(type(group) is int and 0<=group<N_GROUPS,
            "Original source taxon-modulo group must be 0..15")
    output=Path(output)
    require(not output.exists(),"FCP100 original source group cannot be overwritten")
    taxa=read_taxa(taxa_file,check_sha=check_taxa_sha)
    require(len(taxa)==EXPECTED_TAXA,
            "Authentic original all42111 named taxon roster missing")
    owned=source_group_taxa(taxa,group)
    require(len(owned)>0,"Original source group has no matching named taxa")
    output.parent.mkdir(parents=True,exist_ok=True)
    temp=Path(tempfile.mkdtemp(prefix=f"mois-all93m-fcp100-{group:02d}-",dir=output.parent))
    try:
        with sqlite3.connect(temp/"full_photo_observations.sqlite") as db:
            db.execute("PRAGMA journal_mode=OFF")
            db.execute("PRAGMA synchronous=OFF")
            db.execute("PRAGMA cache_size=-65536")
            source.create_database(db)
            group_receipt=load_original_group(db,registry,receipt,tiers_dir,group,
                                              expected_total=expected_total)
            actual_obs,actual_photos=db.execute(
                "SELECT COUNT(*),COALESCE(SUM(photo_count),0) FROM obs"
            ).fetchone()
            require(actual_photos==group_receipt["original_photos"],
                    "Original 93m photo/observation conservation failed for group")
            inventory,eligible,chosen,months,latitudes=select_exact100_from_full_source(db,owned)
        roster=temp/"fcp_original100_complete_nominal_taxa_source_eligibility.csv"
        with roster.open("w",encoding="utf-8",newline="") as f:
            w=csv.writer(f,lineterminator="\n")
            w.writerow(("inat_taxon_id","species","full_93m_original_observation_count",
                        "full_93m_original_photo_count","FCP_U100_source_eligible",
                        "FCP_U100_missing_observations","selected_observed_months",
                        "selected_latitude_bands"))
            for taxon,species in sorted(owned.items(),key=lambda x:int(x[0])):
                n,m=inventory.get(taxon,(0,0))
                w.writerow((taxon,species,n,m,str(n>=FCP_TARGET).lower(),
                           max(0,FCP_TARGET-n),
                           len(months.get(taxon,set())),len(latitudes.get(taxon,set()))))
        csvfile=temp/"exact_100_full93m_original_photo_observation_candidates.csv"
        headers=list(chosen[0]) if chosen else ["photo_id","observation_id","inat_taxon_id","species"]
        with csvfile.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=headers,lineterminator="\n")
            writer.writeheader();writer.writerows(chosen)
        unsupported=temp/"unsupported_photo_original_source_as_technical_unknown.csv"
        unknown=[x for x in chosen if x.get("photo_source_disposition")!="supported"]
        with unsupported.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=headers,lineterminator="\n")
            writer.writeheader();writer.writerows(unknown)
        checked=[x for x in chosen if x.get("photo_source_disposition")=="supported"]
        shards=temp/"shards"
        shards.mkdir()
        # Frozen 256-photo FCP ROI-v4-compatible candidate keys. This is a
        # metadata roster ONLY; prior-photo FCP/MOIS reuse is not yet audited.
        checked.sort(key=lambda r:(source.hash64(
            "mois-FCP-full93m-original-photo-worker-scheduling-v1",
            r["observation_id"]),int(r["photo_id"])))
        hashes={}
        for unit,start in enumerate(range(0,len(checked),256)):
            file=shards/f"shard-{unit:03d}.csv"
            with file.open("w",encoding="utf-8",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=WORKER_COLUMNS,lineterminator="\n")
                writer.writeheader()
                for r in checked[start:start+256]:
                    mid="MOIS-"+hashlib.sha256(
                        ("mois-original-global-RGB-v1|"+r["photo_id"]).encode()
                    ).hexdigest()[:24].upper()
                    writer.writerow({
                        "measurement_id":mid,"photo_id":r["photo_id"],
                        "inat_taxon_id":r["inat_taxon_id"],"species":r["species"],
                        "observation_id":r["observation_id"],
                        "observer_id":r["observer_id"],"observed_on":r["observed_on"],
                        "latitude":r["latitude"],"longitude":r["longitude"],
                        "positional_accuracy":r["positional_accuracy"],
                        "photo_license":r["photo_license"],
                        "photo_url_large":r["photo_url_large"],
                        "image_filename":mid+".jpg","source_arm":SOURCE_ARM,
                        "selection_probability":r["conditional_full_93m_photo_selection_probability"],
                        "source_snapshot_date":"2026-09-27","shard":str(unit),
                    })
            hashes[str(unit)]=sha256(file)
        available=Counter("at_least_100_observations" if n>=100
                          else "below_100_observations" if n>0
                          else "zero_observations"
                          for n,_ in (inventory.get(t,(0,0)) for t in owned))
        report={
            "schema":SCHEMA,
            "full_original_photo_source_IDs":expected_total,
            "full_original_named_plant_taxa":EXPECTED_TAXA,
            "original_taxon_modulo_source_group":group,
            "source_nominal_taxa":len(owned),
            "source_distinct_observations":actual_obs,
            "source_original_photo_IDs":actual_photos,
            "source_tier_archives_read_exactly_once":TIERS,
            "source_group_original_photo_ID_SHA":group_receipt,
            "FCP_U100_source_eligible_nominal_taxa":len(eligible),
            "FCP_U100_source_original_photo_opportunities":len(chosen),
            "FCP_U100_original_photo_ROI_targets_loader_supported":len(checked),
            "FCP_U100_photo_unsupported_format_technical_unknown":len(unknown),
            "source_taxa_eligibility_breakdown":dict(available),
            "original_source_eligibility_and_all_taxa_preserved":True,
            "FCP_original_photo_observation_target_per_taxon":100,
            "original_observation_hash_ranking_same_as_14_source_tiers":True,
            "source_selection_uniform_original_observations_within_taxon":True,
            "one_photo_uniform_per_original_observation_independent_of_observation_hash":True,
            "conditional_FCP_U100_photo_probability":"100/N_taxon/attached_photo_count",
            "eligible_source_taxa_are_not_an_unbiased_census_of_world_flora":True,
            "minimum_classifiable_original_ROIs_for_FCP_intraspecific_statistics":40,
            "full_WCVP_taxonomy_and_native_wild_flowering_photo_status_verified":False,
            "FCP_original_iNaturalist_H9_phenology_annotation_verified_for_all_new_MOIS_photos":False,
            "no_original_image_downloads":True,
            "no_original_ROI_or_hue_inferences":True,
            "all_old_FCP_MOIS_original_photo_IDs_must_be_reused_before_any_new_ROI":True,
            "taxon_source_eligibility_SHA256":sha256(roster),
            "U100_complete_original_photo_candidate_SHA256":sha256(csvfile),
            "unsupported_original_photo_unknown_SHA256":sha256(unsupported),
            "precolour_ROI_worker_shard_SHA256":hashes,
        }
        (temp/"fcp_U100_full93m_group_source_receipt.json").write_text(
            json.dumps(report,indent=2,sort_keys=True)+"\n")
        (temp/"full_photo_observations.sqlite").unlink(missing_ok=True)
        os.replace(temp,output)
        print("MOIS_FULL93M_FCP_U100_EXACT_ORIGINAL_PHOTO_TAXON_GROUP",
              json.dumps({k:report[k] for k in (
                  "original_taxon_modulo_source_group",
                  "source_nominal_taxa","source_original_photo_IDs",
                  "source_distinct_observations",
                  "FCP_U100_source_eligible_nominal_taxa",
                  "FCP_U100_source_original_photo_opportunities",
                  "FCP_U100_photo_unsupported_format_technical_unknown")},sort_keys=True),flush=True)
        return report
    except BaseException:
        shutil.rmtree(temp,ignore_errors=True)
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("registry","receipt","tiers-dir","taxa","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--group",type=int,required=True)
    a=p.parse_args()
    execute_group(a.registry,a.receipt,a.tiers_dir,a.taxa,a.output,a.group)


if __name__=="__main__":
    main()
