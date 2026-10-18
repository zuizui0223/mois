#!/usr/bin/env python3
"""MOIS worldwide flowering-plant breadth source from official iNaturalist Open Data.

This is NOT FCP's 1000-taxon within-species/high-depth photo source.
It scans the *monthly open-licensed metadata* taxonomy -> observation ->
photo snapshot once, caches many photo-ID-first candidate taxa, and seals
the FIRST 4096 original photo IDs for genuine blind RGB ROI-v4 processing.
It never treats source-only photo IDs as a measured flower-colour phenotype.

The 3 publicly released gzipped tables are each streamed once using source
ETag, Content-Length and decompressed gzip CRC; photo image bytes are never
downloaded by THIS script. Snapshot metadata may include non-flower photos.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import re
from collections import Counter,defaultdict
from datetime import date,datetime,timezone
from pathlib import Path

import pandas as pd
import requests

from scripts.acquisition.mois_probe_inat_s3_monthly_metadata_20261017 import (
    BUCKET, probe,
)
from scripts.acquisition.mois_full_historical_source_dates_official_inat_s3 import (
    scan_table,
)

ANGIOSPERMAE="47125"
EXPECTED_PREVIOUS_PHOTOS=313535
EXPECTED_PREVIOUS_COLOUR=151239
MAX_OBS_PER_SPECIES=4
POOL_GOAL=32768
FIRST_WAVE=4096
N_SHARDS=16
PER_SHARD=256
SALT="MOIS_GLOBAL_FLOWERING_PLANTS_BROAD_FIRST_PER_SPECIES_OPEN_SOURCE_20261018_V1"
POOL="mois_global_angiosperms_broad_open_source_photo_ID_pool.csv"
FIRST="mois_broad_first4096_original_open_flower_photo_ID_source.csv"
RECEIPT="mois_broad_angiosperm_new_taxa_32768_pool_first4096_receipt.json"
PLAN_PREFIX="mois_broad_source_only_shard_"
CORE_MANIFEST="primary313535_fullfourth4096_featurelake_manifest.json"
CORE_PHOTO="primary313535_fullfourth4096_original_photo_hues.parquet"
EXCLUDED_OLD="mois_new4096_original_photo_source_after313535_outcome_blind.csv"
EXCLUDED_OLD_RECEIPT="mois_new4096_photo_source_after313535_frozen_receipt.json"
EXCLUDED_SI="mois_targeted4096_photo_original_census_blind_selection.csv"
EXCLUDED_SI_RECEIPT="mois_first4096_original_source_only_freeze.json"
COLUMNS=("photo_id","inat_taxon_id","species","observation_id",
         "observer_id","observed_on","latitude","longitude",
         "positional_accuracy","quality_grade","photo_license","photo_url_large",
         "source_snapshot_date","measurement_id","image_filename",
         "source_arm","source_sampling_priority","genus","grid10",
         "wave_sequence","shard")
PHOTO_LICENSES={"cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa",
                "cc-by-nd","cc-by-nc-nd"}
PHOTO_EXTENSIONS={"jpg","jpeg","png"}
RANK_ALLOWED={"species"}
POOL_PRIORITY=("new_genus_new_species","known_genus_new_species",
               "existing_species_new_region","existing_species_replicate")


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as inp:
        for chunk in iter(lambda:inp.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def priority(*elements):
    b=(SALT+"|"+"|".join(str(t) for t in elements)).encode()
    return int(hashlib.sha256(b).hexdigest()[:16],16)


def norm_genus(name):
    parts=name.strip().split()
    if len(parts)>=2 and re.fullmatch(r"[A-Z][a-zA-Z-]+",parts[0]):
        return parts[0]
    return None


def grid10(lat,lon):
    return f"{int(math.floor((lat+90)/10)):02d}:{int(math.floor((lon+180)/10)):02d}"


def old_photo_source(path):
    p=Path(path)
    source=p/CORE_PHOTO
    receipt=json.loads((p/CORE_MANIFEST).read_text())
    require(receipt["schema"]=="mois_313535_genuine_original_ROI_feature_lake_v2"
            and receipt["actual_original_photo_attempts"]==EXPECTED_PREVIOUS_PHOTOS
            and receipt["four_hue_labels"]==EXPECTED_PREVIOUS_COLOUR
            and receipt["files_sha256"][CORE_PHOTO]==sha(source),
            "Original measured worldwide photograph source receipt changed")
    photos=pd.read_parquet(source,columns=[
        "photo_id","inat_taxon_id","species","latitude","longitude","classified"])
    require(len(photos)==EXPECTED_PREVIOUS_PHOTOS and photos.photo_id.is_unique,
            "Original real-pixel photo IDs or sampling support invalid")
    ids=set(photos.photo_id.astype(str))
    taxa=set(photos.inat_taxon_id.astype(str))
    canon=set(str(v).strip().lower() for v in photos.species.dropna()
              if len(str(v).strip().split())>=2)
    genera=set(x for x in (norm_genus(str(v)) for v in photos.species.dropna())
               if x)
    cells=Counter()
    for lat,lon in zip(photos.latitude,photos.longitude):
        try:
            a,b=float(lat),float(lon)
            if math.isfinite(a) and math.isfinite(b):
                cells[grid10(a,b)]+=1
        except (ValueError,TypeError):
            continue
    return ids,taxa,canon,genera,cells,receipt


def reservation(directory,csvname,receiptname,key,expected):
    t=Path(directory)/csvname
    r=json.loads((Path(directory)/receiptname).read_text())
    require(r["status"]=="SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES"
            and r["new_photo_pixels_actually_measured"] is False
            and r[key]==sha(t),
            "Previously reserved source-photo ID manifest changed")
    ids=set()
    with t.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            pid=row["photo_id"]
            require(pid.isdecimal() and pid not in ids,
                    "Reserved original photo ID duplicated")
            ids.add(pid)
    require(len(ids)==expected,"Existing fifth/SI original photo reservation incomplete")
    return ids,r


def new_eligible_taxon(taxon_id,ancestry,rank,name,active):
    """Only source species under Angiospermae, not all plants or FCP taxa."""
    ancestors=set(re.split(r"[/\\]",ancestry))
    genus=norm_genus(name)
    return (taxon_id.isdecimal() and ANGIOSPERMAE in ancestors
            and rank.lower() in RANK_ALLOWED
            and active.lower() in {"t","true","yes","1"}
            and genus is not None)


def scan_taxa(client,entry):
    taxa={}
    def callback(row,fields):
        tid=row[fields["taxon_id"]].decode("ascii")
        ancestry=row[fields["ancestry"]].decode("ascii")
        rank=row[fields["rank"]].decode("ascii")
        name=row[fields["name"]].decode("utf-8")
        active=row[fields["active"]].decode("ascii")
        if new_eligible_taxon(tid,ancestry,rank,name,active):
            taxa[tid]=(name,norm_genus(name))
    r=scan_table(client,entry,("taxon_id","ancestry","rank","name","active"),
                 callback,"taxa.csv.gz")
    require(len(taxa)>25000,"Official source has unexpectedly few active flowering species")
    return taxa,r


def eligible_observation(taxon_id,grade,observer,uuid,day,lat,lon,accuracy,taxa):
    if (taxon_id not in taxa or grade.lower()!="research"
        or not observer.isdecimal() or not uuid or not day):
        return None
    try:
        when=date.fromisoformat(day)
        x,y=float(lat),float(lon)
        if not (date(1800,1,1)<=when<=date.today()
                and math.isfinite(x) and math.isfinite(y)
                and -90<=x<=90 and -180<=y<=180):
            return None
        if accuracy:
            z=float(accuracy)
            if not math.isfinite(z) or z<0 or z>20000:
                return None
    except (ValueError,TypeError,OverflowError):
        return None
    return {
        "taxon_id":taxon_id,"observation_uuid":uuid,
        "observer_id":observer,"observed_on":day,
        "latitude":lat,"longitude":lon,
        "positional_accuracy":accuracy,
        "grid10":grid10(x,y),"quality_grade":"research",
    }


def select_observations(client,entry,taxa):
    """Top four outcome-blind observation hashes per original named species."""
    heaps=defaultdict(list)
    counters=Counter()
    def one(row,idx):
        counters["all_provider_observation_rows"]+=1
        tid=row[idx["taxon_id"]].decode("ascii")
        if tid not in taxa:
            return
        b={k:row[idx[k]].decode("utf-8") for k in
           ("quality_grade","observer_id","observation_uuid",
            "observed_on","latitude","longitude","positional_accuracy")}
        ob=eligible_observation(tid,b["quality_grade"],b["observer_id"],
            b["observation_uuid"],b["observed_on"],b["latitude"],
            b["longitude"],b["positional_accuracy"],taxa)
        if ob is None:
            counters["species_observation_metadata_excluded"]+=1
            return
        counters["eligible_species_original_observations"]+=1
        # Highest observation_uuid hash is NOT a hue/flower predictor.
        key=priority("observation",ob["observation_uuid"])
        heap=heaps[tid]
        item=(-key,ob["observation_uuid"],ob)
        if len(heap)<MAX_OBS_PER_SPECIES:
            heapq.heappush(heap,item)
        elif item>heap[0]:
            heapq.heapreplace(heap,item)
    receipt=scan_table(client,entry,
        ("taxon_id","quality_grade","observer_id","observation_uuid",
         "observed_on","latitude","longitude","positional_accuracy"),
        one,"observations.csv.gz")
    obs={}
    for taxon,heap in heaps.items():
        by_user=set()
        by_cell=set()
        for _,uuid,val in sorted(heap,key=lambda q:(-q[0],q[1])):
            # Prefer distinct photographers and cells for fallback candidates.
            if val["observer_id"] in by_user:
                continue
            by_user.add(val["observer_id"])
            by_cell.add(val["grid10"])
            obs[uuid]=val
    require(len(obs)>=POOL_GOAL and len(heaps)>=FIRST_WAVE,
            "Source angiosperm observations fail worldwide taxon-breadth threshold")
    return obs,dict(counters),receipt


def choose_photos(client,entry,observations,old_ids,reserved):
    """Prefer earliest positional OPEN PHOTO on each candidate observation."""
    winner={}
    stats=Counter()
    def one(row,idx):
        uuid=row[idx["observation_uuid"]].decode("ascii")
        base=observations.get(uuid)
        if base is None:return
        pid=row[idx["photo_id"]].decode("ascii")
        ext=row[idx["extension"]].decode("ascii").lower().lstrip(".")
        license_name=row[idx["license"]].decode("ascii").lower()
        if pid in old_ids or pid in reserved:
            stats["already_measured_or_reserved_photo_IDs"]+=1
            return
        if not pid.isdecimal() or ext not in PHOTO_EXTENSIONS or (
           license_name not in PHOTO_LICENSES):
            stats["photo_extension_license_or_ID_invalid"]+=1
            return
        pos=row[idx["position"]].decode("ascii")
        try:position=int(pos) if pos else 999
        except ValueError:position=999
        if position<0:position=999
        width=row[idx["width"]].decode("ascii")
        height=row[idx["height"]].decode("ascii")
        if width.isdecimal() and height.isdecimal():
            if min(int(width),int(height))<360:
                stats["photo_below_size_quality_threshold"]+=1
                return
        sort=(position,priority("photo",pid))
        current=winner.get(uuid)
        if current is not None and current[0]<=sort:
            return
        winner[uuid]=(sort,{
            "photo_id":pid,
            "photo_license":license_name,
            "photo_url_large":
               "https://inaturalist-open-data.s3.amazonaws.com/photos/"+
               pid+"/large."+ext,
            **base,
        })
    receipt=scan_table(client,entry,
        ("photo_id","observation_uuid","extension","license","position",
         "width","height"),one,"photos.csv.gz")
    result=[z[1] for z in winner.values()]
    require(len(result)>=FIRST_WAVE,"Unable to identify broad open photo corpus")
    stats["unique_observations_with_valid_photo"]=len(result)
    return result,dict(stats),receipt


def rank_by_breadth(candidates,taxa,old_taxa,old_canon,old_genera,old_cells,
                    *,limit=POOL_GOAL):
    """One taxon per genus per pass BEFORE any within-taxon photo replicate."""
    by_taxon=defaultdict(list)
    for c in candidates:
        by_taxon[c["taxon_id"]].append(c)
    bands={name:defaultdict(list) for name in POOL_PRIORITY}
    for taxon,rows in by_taxon.items():
        if taxon not in taxa:continue
        name,genus=taxa[taxon]
        unseen=(taxon not in old_taxa and name.lower() not in old_canon)
        row_candidates=sorted(rows,key=lambda r:(
            old_cells.get(r["grid10"],0),
            priority("candidate",r["photo_id"])))
        if unseen and genus not in old_genera:
            band="new_genus_new_species"
        elif unseen:
            band="known_genus_new_species"
        else:
            band="existing_species_new_region"
        for k,source in enumerate(row_candidates):
            actual_band=band if k==0 else "existing_species_replicate"
            bands[actual_band][genus].append({
                "taxon_id":taxon,
                "species":name,"genus":genus,
                "source_sampling_priority":actual_band,
                **source,
            })
    chosen=[];used_photo=set();used_obs=set()
    seen_taxa=set();source_counts=Counter()
    for stage in POOL_PRIORITY:
        by_genus=bands[stage]
        keys=sorted(by_genus,key=lambda genus:priority("genus",genus))
        queues={}
        for genus in keys:
            queues[genus]=sorted(by_genus[genus],key=lambda q:(
                old_cells.get(q["grid10"],0),
                priority("rank",q["taxon_id"],q["photo_id"])))
        for turn in range(8):
            for genus in keys:
                if len(chosen)>=limit:break
                pool=queues[genus]
                if turn>=len(pool):continue
                q=pool[turn]
                if q["photo_id"] in used_photo or q["observation_uuid"] in used_obs:
                    continue
                # No more than 1 image per taxon on first pass until exhausted.
                if stage!="existing_species_replicate" and q["taxon_id"] in seen_taxa:
                    continue
                used_photo.add(q["photo_id"]);used_obs.add(q["observation_uuid"])
                seen_taxa.add(q["taxon_id"])
                chosen.append(q);source_counts[stage]+=1
            if len(chosen)>=limit:break
        if len(chosen)>=limit:break
    require(len(chosen)>=FIRST_WAVE,
            "Globally broadened source failed to supply 4096 original photos")
    return chosen,dict(source_counts)


def freeze(source,core,old_plan,season_plan,out,*,client=None):
    require(not out.exists() or not any(out.iterdir()),
            "Never overwrite frozen broad photographic source")
    old_ids,old_taxa,old_canon,old_genera,old_cells,core_receipt=old_photo_source(core)
    fifth,receipt5=reservation(old_plan,EXCLUDED_OLD,EXCLUDED_OLD_RECEIPT,
                             "first_wave_original_4096_csv_sha256",4096)
    season,receipts=reservation(season_plan,EXCLUDED_SI,EXCLUDED_SI_RECEIPT,
                              "first_wave_original_4096_csv_sha256",4096)
    require(fifth.isdisjoint(old_ids) and season.isdisjoint(old_ids)
            and fifth.isdisjoint(season),
            "Existing measured source and separately reserved photo IDs overlap")
    snapshot={}
    with (client or requests.Session()) as session:
        for k in ("taxa","observations","photos"):
            snapshot[k]=probe(BUCKET+"/"+k+".csv.gz",session)
        taxa,taxarec=scan_taxa(session,snapshot["taxa"])
        obs,counts,obsrec=select_observations(
            session,snapshot["observations"],taxa)
        pictures,photo_stats,photorec=choose_photos(
            session,snapshot["photos"],obs,old_ids,fifth|season)
    selected,stage=rank_by_breadth(
        pictures,taxa,old_taxa,old_canon,old_genera,old_cells)
    snapshot_date=datetime.now(timezone.utc).date().isoformat()
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for n,v in enumerate(selected):
        tid=v["taxon_id"]
        photo=v["photo_id"]
        mid="MOIS-"+hashlib.sha256(
            (SALT+"|"+photo).encode()).hexdigest()[:24].upper()
        rows.append({
            "photo_id":photo,
            "inat_taxon_id":tid,"species":v["species"],
            "observation_id":v["observation_uuid"],
            "observer_id":v["observer_id"],
            "observed_on":v["observed_on"],
            "latitude":v["latitude"],"longitude":v["longitude"],
            "positional_accuracy":v["positional_accuracy"],
            "quality_grade":"research","photo_license":v["photo_license"],
            "photo_url_large":v["photo_url_large"],
            "source_snapshot_date":snapshot_date,
            "measurement_id":mid,"image_filename":mid+".jpg",
            "source_arm":"mois_broad_angiosperm_species_first_official_open_S3",
            "source_sampling_priority":v["source_sampling_priority"],
            "genus":v["genus"],"grid10":v["grid10"],
            "wave_sequence":str(n),
            "shard":str(n%N_SHARDS) if n<FIRST_WAVE else "-1",
        })
    def write(path,group):
        with path.open("x",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
            writer.writeheader()
            for item in group:writer.writerow(item)
    write(out/POOL,rows)
    subset=rows[:FIRST_WAVE]
    write(out/FIRST,subset)
    shards={}
    for i in range(N_SHARDS):
        sub=[r for r in subset if int(r["shard"])==i]
        require(len(sub)==PER_SHARD,"Missing planned original broad photo worker")
        filename=f"{PLAN_PREFIX}{i:02d}_256.csv"
        write(out/filename,sub)
        shards[str(i)]=sha(out/filename)
    new_taxa={x["inat_taxon_id"] for x in subset
              if x["inat_taxon_id"] not in old_taxa}
    new_genera={x["genus"] for x in subset if x["genus"] not in old_genera}
    r={
        "schema":"mois_open_global_angiosperm_breadth_first4096_v1",
        "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
        "source_snapshot_date":snapshot_date,
        "angiosperm_ancestor_taxon_id":ANGIOSPERMAE,
        "original_provider_taxon_rank":"species",
        "main_original_real_photo_attempts":EXPECTED_PREVIOUS_PHOTOS,
        "main_original_four_hue_photo_ROIs":EXPECTED_PREVIOUS_COLOUR,
        "first_wave_original_source_photo_ids":FIRST_WAVE,
        "new_photo_pixels_actually_measured":False,
        "no_photo_colour_or_fitness_or_bloom_label_used":True,
        "first_wave_all_source_photo_IDs_are_unprocessed":True,
        "source_photo_IDs_disjoint_from_both_running_and_SI_reservations":True,
        "official_public_open_photo_license_required":True,
        "independently_verified_flowering_petals":False,
        "taxon_population_breadth_is_main_goal_not_within_species_photos":True,
        "source_taxa_eligible":len(taxa),
        "source_observations_with_eligible_species_metadata":len(obs),
        "source_photo_associated_original_candidates":len(pictures),
        "pool_original_photo_IDs":len(rows),
        "first_wave_distinct_photo_ids":len(subset),
        "first_wave_named_taxa":len({x["inat_taxon_id"] for x in subset}),
        "first_wave_named_genera":len({x["genus"] for x in subset}),
        "new_taxa_not_in_original_measurements_by_taxon_id":len(new_taxa),
        "new_genera_not_in_original_measurements_by_name":len(new_genera),
        "first_wave_distinct_observation_UUIDs":len({x["observation_id"] for x in subset}),
        "first_wave_distinct_original_observers":len({x["observer_id"] for x in subset}),
        "first_wave_latlon_10degree_cells":len({x["grid10"] for x in subset}),
        "eligible_species_photo_priority_counts":stage,
        "observation_metadata_status_counts":counts,
        "photo_metadata_filter_status_counts":photo_stats,
        "official_photo_vs_observation_source_SHA256":{
            "taxa":taxarec["compressed_original_SHA256"],
            "observations":obsrec["compressed_original_SHA256"],
            "photos":photorec["compressed_original_SHA256"],
        },
        "source_ETag_first_HEAD_and_gzip_length":{
            key:{"etag":v["etag"],"content_length_bytes":v["content_length_bytes"]}
            for key,v in snapshot.items()},
        "first_wave_original_4096_csv_sha256":sha(out/FIRST),
        "original_broad_all_future_source_pool_SHA256":sha(out/POOL),
        "per_shard_sha256":shards,
        "main_previous_313535_source_sha256":core_receipt["files_sha256"][CORE_PHOTO],
        "fifth4096_source_only_reservation_sha256":
            receipt5["first_wave_original_4096_csv_sha256"],
        "SI4096_source_only_reservation_sha256":
            receipts["first_wave_original_4096_csv_sha256"],
        "quality_grade_research_and_positional_accuracy_at_most_20km":True,
        "single_photo_per_species_has_no_guaranteed_open_bloom":True,
        "not_a_probability_representative_world_flora_sample":True,
        "no_FCP_spatial_polymorphism_W_D_50km_claim":True,
        "licensed_raw_original_RGB_bytes_not_downloaded_or_redistributed":True,
    }
    (out/RECEIPT).write_text(json.dumps(
        r,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("MOIS_BROAD_ANGIOSPERM_SOURCE_CANDIDATES_NOT_IMAGE_PIXEL_HUES",
          json.dumps({k:r[k] for k in (
              "pool_original_photo_IDs","first_wave_named_taxa",
              "first_wave_named_genera",
              "new_taxa_not_in_original_measurements_by_taxon_id",
              "new_genera_not_in_original_measurements_by_name",
              "first_wave_latlon_10degree_cells")},sort_keys=True),flush=True)
    return r


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--previous313535-photo-release",type=Path,required=True)
    ap.add_argument("--fifth4096-reserved-source",type=Path,required=True)
    ap.add_argument("--season4096-reserved-source",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    freeze(None,args.previous313535_photo_release,
           args.fifth4096_reserved_source,args.season4096_reserved_source,
           args.output)
