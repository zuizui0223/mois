#!/usr/bin/env python3
"""Select the NEXT 4096 original angiosperm species PHOTO IDS from cached pool.

No additional iNaturalist 500-million-photo metadata pass is performed.
The first source-only 4096 are excluded, as are their source taxa when
possible. Select new genera/taxa and photo-observer diversity before any RGB,
flower/colour/bloom image signal is read. This code DOES NOT measure pixels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path

import pandas as pd

from scripts.acquisition import mois_global_angiosperm_breadth_source as src

BASE_SOURCE_TAG="mois-angiosperm-broad-species-first-source-cache-v1-20261018"
FIRST_N=4096
SECOND_N=4096
SHARDS=16
PER_SHARD=256
POOL_NAME=src.POOL
FIRST_NAME=src.FIRST
FIRST_RECEIPT=src.RECEIPT
SECOND_NAME="mois_global_angiosperm_second4096_source_only_original_photoids.csv"
SECOND_RECEIPT="mois_global_angiosperm_cached_second4096_source_only_receipt.json"
SHARD_STEM="mois_broad_secondwave_source_shard_"
SELECTION_SALT="MOIS_SECOND_WAVE_AFTER_GLOBAL_SPECIES_FIRST4096_SOURCE_OUTCOME_BLIND_20261018"


def require(condition,message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):
            h.update(data)
    return h.hexdigest()


def rank(*args):
    return hashlib.sha256(
        (SELECTION_SALT+"|"+"|".join(map(str,args))).encode()).hexdigest()


def load_original_pool(path):
    original=json.loads((path/FIRST_RECEIPT).read_text())
    require(
        original["schema"]=="mois_open_global_angiosperm_breadth_first4096_v1"
        and original["status"]=="SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES"
        and original["main_original_real_photo_attempts"]==313535
        and original["first_wave_distinct_photo_ids"]==FIRST_N
        and original["new_photo_pixels_actually_measured"] is False
        and original["source_photo_IDs_disjoint_from_both_running_and_SI_reservations"] is True
        and original["original_broad_all_future_source_pool_SHA256"]==
            sha(path/POOL_NAME)
        and original["first_wave_original_4096_csv_sha256"]==
            sha(path/FIRST_NAME),
        "Cached global angiosperm first source pool is not signed")
    pool=pd.read_csv(path/POOL_NAME,dtype=str,keep_default_na=False)
    first=pd.read_csv(path/FIRST_NAME,dtype=str,keep_default_na=False)
    require(
        len(pool)==original["pool_original_photo_IDs"]==32768
        and len(first)==FIRST_N
        and pool.photo_id.is_unique and first.photo_id.is_unique
        and pool.measurement_id.is_unique
        and set(pool.columns)==set(src.COLUMNS)
        and first.reset_index(drop=True).equals(
            pool.iloc[:FIRST_N].reset_index(drop=True)),
        "Original 32,768 candidate pool or first-wave identity changed")
    return pool,first,original


def choose_next(pool,first,goal=SECOND_N):
    require(
        pool.photo_id.is_unique and first.photo_id.is_unique and
        set(first.photo_id).issubset(set(pool.photo_id)),
        "Second wave requires source-locked first-wave original PHOTO IDs")
    first_ids=set(first.photo_id)
    first_taxa=set(first.inat_taxon_id)
    first_genera=set(first.genus)
    first_regions=Counter(first.grid10)
    # The second source is still not an image or coloured-petal census.
    fresh=pool.loc[~pool.photo_id.isin(first_ids)
                   & ~pool.inat_taxon_id.isin(first_taxa)].copy()
    require(len(fresh)>=goal,
            "Cached pool lacks enough previously unselected angiosperm species")
    groups={"novel_genus":defaultdict(list),"additional_species":defaultdict(list)}
    for r in fresh.to_dict("records"):
        typ="novel_genus" if r["genus"] not in first_genera else "additional_species"
        groups[typ][r["genus"]].append(r)
    selected=[]
    seen_photos=set(first_ids)
    seen_taxa=set(first_taxa)
    seen_observations=set(first.observation_id)
    seen_observers=set()
    counts=Counter()
    # First diversify new genera, then distinct species within other genera.
    # Within each genus use less-photographed 10-degree source cells before
    # source-blind deterministic hashes; never condition on unseen hue.
    for priority,genus_map in groups.items():
        genera=sorted(genus_map,key=lambda g:rank("genus",g))
        for g in genera:
            genus_map[g].sort(key=lambda r:(
                first_regions.get(r["grid10"],0),
                rank("photo",r["photo_id"])))
        cursor={g:0 for g in genera}
        # Two passes: strongly prefer distinct photographers. Do not allow
        # photo or observation duplicates even on the fallback pass.
        for enforce_new_observer in (True,False):
            if not enforce_new_observer:
                # Reconsider candidates rejected only because their
                # photographer was already selected on the first pass.
                # Selected photo/taxon/observation IDs still stay forbidden.
                cursor={g:0 for g in genera}
            while len(selected)<goal:
                round_add=0
                for g in genera:
                    q=genus_map[g]
                    k=cursor[g]
                    while k<len(q):
                        v=q[k]
                        k+=1
                        if (v["photo_id"] in seen_photos or
                            v["inat_taxon_id"] in seen_taxa or
                            v["observation_id"] in seen_observations or
                            (enforce_new_observer and
                             v["observer_id"] in seen_observers)):
                            continue
                        break
                    else:
                        cursor[g]=k
                        continue
                    cursor[g]=k
                    # Keep only published original photo ID metadata columns;
                    # never copy accidental hue or model-result fields.
                    v={key:v.get(key,"") for key in src.COLUMNS}
                    v["source_arm"]="mois_broad_angiosperm_cached_wave2_new_species"
                    v["source_sampling_priority"]=(
                        priority+"_new_observer" if enforce_new_observer
                        else priority+"_observer_reuse")
                    v["wave_sequence"]=str(len(selected))
                    v["shard"]=str(len(selected)%SHARDS)
                    selected.append(v)
                    seen_photos.add(v["photo_id"])
                    seen_taxa.add(v["inat_taxon_id"])
                    seen_observations.add(v["observation_id"])
                    seen_observers.add(v["observer_id"])
                    counts[v["source_sampling_priority"]]+=1
                    round_add+=1
                    if len(selected)==goal:
                        break
                if len(selected)==goal or round_add==0:
                    break
            if len(selected)==goal:
                break
        if len(selected)==goal:
            break
    require(
        len(selected)==goal and
        len({r["photo_id"] for r in selected})==goal and
        len({r["inat_taxon_id"] for r in selected})==goal and
        len({r["observation_id"] for r in selected})==goal and
        not ({r["photo_id"] for r in selected}&first_ids) and
        not ({r["inat_taxon_id"] for r in selected}&first_taxa),
        "Unable to reserve enough globally breadth-increasing ORIGINAL species/photo IDs")
    return selected,dict(counts)


def save(path,rows):
    with path.open("x",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(src.COLUMNS),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(cached,out):
    require(not out.exists() or not any(out.iterdir()),
            "Never overwrite a previously frozen original photo-ID source")
    pool,first,rec=load_original_pool(cached)
    selected,counts=choose_next(pool,first)
    out.mkdir(parents=True,exist_ok=True)
    save(out/SECOND_NAME,selected)
    subtotals={}
    for shard in range(SHARDS):
        rows=[r for r in selected if r["shard"]==str(shard)]
        require(len(rows)==PER_SHARD,
                "Next wave must contain 16 complete 256-original-photo-ID subshards")
        fn=f"{SHARD_STEM}{shard:02d}_256.csv"
        save(out/fn,rows)
        subtotals[str(shard)]=sha(out/fn)
    original_genera=set(first.genus)
    novel_genera={r["genus"] for r in selected}-original_genera
    source={
        "schema":"mois_angiosperm_cached_second4096_new_species_source_only_v1",
        "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
        "parent_original_metadata_snapshot_date":rec["source_snapshot_date"],
        "parent_first_wave_original_4096_csv_sha256":rec[
            "first_wave_original_4096_csv_sha256"],
        "parent_cached_pool_original_32768_photo_SHA256":rec[
            "original_broad_all_future_source_pool_SHA256"],
        "no_repeat_worldwide_iNaturalist_photo_metadata_scan":True,
        "cached_candidate_photo_IDs":len(pool),
        "frozen_first4096_original_photo_IDs_excluded":FIRST_N,
        "second4096_photo_candidates":SECOND_N,
        "second4096_unique_named_taxa":len({r["inat_taxon_id"] for r in selected}),
        "second4096_unique_genera":len({r["genus"] for r in selected}),
        "second4096_additional_genera_not_in_first_wave":len(novel_genera),
        "second4096_10degree_geographic_source_cells":len({
            r["grid10"] for r in selected}),
        "second4096_distinct_original_observation_UUIDs":len({
            r["observation_id"] for r in selected}),
        "second4096_distinct_original_observers":len({
            r["observer_id"] for r in selected}),
        "source_priority_counts":counts,
        "first4096_taxa_and_photo_IDs_excluded":True,
        "first4096_and_new4096_disjoint_photo_IDs":True,
        "original313535_measured_photo_IDs_excluded_via_parent_signed_source":True,
        "other_fifth4096_and_SI4096_reservations_excluded_via_parent":True,
        "only_open_licensed_research_grade_source_metadata":True,
        "new_photo_pixels_actually_measured":False,
        "true_open_petals_photo_specifically_verified":False,
        "no_new_hue_assigned_to_a_source_only_photo":True,
        "first_wave_remaining_RGB_photo_processing_runs_independently":True,
        "second_wave_source_colours_NOT_known_until_actual_RGB_model_runs":True,
        "no_maintained_global_50_50_equilibrium_pigment_benefit_fitness_claim":True,
        "source_only_second4096_csv_sha256":sha(out/SECOND_NAME),
        "source_only_shard_SHA256":subtotals,
    }
    (out/SECOND_RECEIPT).write_text(
        json.dumps(source,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("MOIS_CACHED_SECOND4096_BROAD_TAXA_NOT_RGB_MEASURED",
          json.dumps({k:source[k] for k in (
              "second4096_photo_candidates","second4096_unique_named_taxa",
              "second4096_unique_genera",
              "second4096_additional_genera_not_in_first_wave",
              "second4096_10degree_geographic_source_cells",
              "second4096_distinct_original_observers")},sort_keys=True),
          flush=True)
    return source


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--previously-released-broad-source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.previously_released_broad_source,a.output)
