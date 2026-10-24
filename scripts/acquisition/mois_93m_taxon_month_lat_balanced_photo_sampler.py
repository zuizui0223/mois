#!/usr/bin/env python3
"""One-pass source-metadata-balanced original PHOTO selection from the 93m registry.

All 93,209,351 original photo IDs remain in the upstream immutable target
register. This staged design only reads ONE signed random observation-hash tier
(tier06: probability width 0.0064; 597,070 actual original-photo ID records).
It opens ZERO images. Within each original biological observation, ONE photo is
chosen by an independently salted pseudo-random photo rank, not by colour or
loader format. Observations are then sampled by lowest independently salted
uniform ranks, at most K per (nominal iNat taxon ID x recorded month x latitude
band). Every observation/photo has positive, auditable inclusion probability.

For one photo j from an original observation with m photos, in source stratum s
with N_s original observations in tier06, design probabilities are:
  π_source observation = 0.0128 - 0.0064 = 0.0064;
  π_phase2 observation | tier06 = min(K,N_s)/N_s;
  π_photo | original observation = 1/m;
  π_joint photo = π_source * π_phase2 * π_photo.

These probabilities are true randomization-design probabilities under
independent seeded uniform hash randomizations, NOT botany/photo classifier
recall. Unsupported original image formats remain in a signed technical-unknown
sidecar rather than quietly being declared non-flowering or coloured.

Uses SQLite for observation/photo multiplicity and sorting. Unlike a process
collecting Python sets of 93m IDs, the full photo-source population is scanned
once, with no loaded image and bounded Python-row state. Disk/CPU costs of
the source metadata SQLite index still scale with observations in the tier.

WARNING: This is a new DESIGN arm with different probabilities from existing
p0032, geo-enrichment, and tier06 compute-bucket arms. Never pool their photo
white fractions unweighted. Individual wild/native/anthesis remains unverified.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import date
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import sqlite3
import tarfile
import tempfile

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, EXPECTED_PHOTOS, EXPECTED_TAXA, TIER_ENDS, obs_draw, sha256,
    normalized_extension, require, usable_date,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import (
    WORKER_COLUMNS, read_taxa,
)
from scripts.acquisition.mois_singlepass_signed_tier_bucket_router import ROUTE_SCHEMA
from scripts.acquisition.mois_tier06_all16_second8_broad_warm_photo_queue import PINNED

SCHEMA = "mois_93m_tier06_taxon_month_lat_stratified_onephoto_v1"
SOURCE_ARM = "mois_global93m_tier06_taxon_month_lat_probability_onephoto_K1_v1"
SOURCE_TIER = 6
SOURCE_GROUPS = 16
SOURCE_PHOTOS = 597_070
SOURCE_PROBABILITY = TIER_ENDS[SOURCE_TIER] - TIER_ENDS[SOURCE_TIER - 1]
TIER_REGISTRY_SCHEMA = "mois_complete_full93m_original_photo_target_registry_v1"
SALT_OBSERVATION = "mois-stratified-original-observation-priority-v1"
SALT_PHOTO = "mois-uniform-one-photo-per-original-observation-v1"
MAX_CAP = 120_000

# These are the SAME original zero-pixel source archive signatures that
# underwrite 597,070 tier06 original IDs, not alternative iNaturalist API pulls.


def hash64(salt, value):
    seed = f"{salt}|{value}".encode("utf-8")
    return hashlib.sha256(seed).hexdigest()[:16]


def month_from_recorded_date(value):
    s = str(value or "").strip()
    if usable_date(s):
        return s[5:7]
    return "unknown"


def latitude_band(value):
    try:
        lat = float(value)
    except (TypeError, ValueError):
        return "unknown"
    if not math.isfinite(lat) or lat < -90 or lat > 90:
        return "unknown"
    if -23.44 <= lat <= 23.44:
        return "tropical"
    return "south_extratropical" if lat < -23.44 else "north_extratropical"


UPSERT = """
INSERT INTO obs (
    uuid,taxon,record_month,lat_band,obs_rank,photo_count,
    best_photo_rank,best_photo_id,photo_license,photo_extension,
    recorded_on,latitude,longitude,observer_id,positional_accuracy
) VALUES (?,?,?,?,?,1,?,?,?,?,?,?,?,?,?)
ON CONFLICT(uuid) DO UPDATE SET
    photo_count=obs.photo_count+1,
    best_photo_rank=CASE WHEN excluded.best_photo_rank < obs.best_photo_rank
                         THEN excluded.best_photo_rank ELSE obs.best_photo_rank END,
    best_photo_id=CASE WHEN excluded.best_photo_rank < obs.best_photo_rank
                       THEN excluded.best_photo_id ELSE obs.best_photo_id END,
    photo_license=CASE WHEN excluded.best_photo_rank < obs.best_photo_rank
                       THEN excluded.photo_license ELSE obs.photo_license END,
    photo_extension=CASE WHEN excluded.best_photo_rank < obs.best_photo_rank
                         THEN excluded.photo_extension ELSE obs.photo_extension END
WHERE obs.taxon=excluded.taxon AND obs.record_month=excluded.record_month
  AND obs.lat_band=excluded.lat_band
  AND obs.recorded_on=excluded.recorded_on
  AND obs.latitude=excluded.latitude AND obs.longitude=excluded.longitude
  AND obs.observer_id=excluded.observer_id
  AND obs.positional_accuracy=excluded.positional_accuracy
"""


def import_photo(db, row):
    uuid, pid, taxon = row["observation_uuid"], row["photo_id"], row["inat_taxon_id"]
    require(bool(uuid) and pid.isdecimal() and taxon.isdecimal(),
            "Original photo/observation/taxon ID missing or invalid")
    require(str(row["tier"]) == str(SOURCE_TIER)
            and float(row["observation_inclusion_score"]) == obs_draw(uuid),
            "Original source observation rank or tier was changed")
    original_date = row["observed_on"]
    rank = hash64(SALT_PHOTO, uuid + "|" + pid)
    values = (
        uuid, taxon, month_from_recorded_date(original_date),
        latitude_band(row["latitude"]), hash64(SALT_OBSERVATION, uuid),
        rank, pid, row["photo_license"], row["image_extension"],
        original_date, row["latitude"], row["longitude"], row["observer_id"],
        row["positional_accuracy"],
    )
    cursor = db.execute(UPSERT, values)
    require(cursor.rowcount == 1, "Conflicting source metadata for original observation")


def load_signed_group(db, root, group, *, pinned=PINNED):
    """Verify all 512 original source-bucket bytes before counting any photos."""
    number = f"{group:02d}"
    receipt_path = Path(root) / f"group{number}" / "routing_receipt.json"
    archive = Path(root) / f"group{number}" / f"tier06-group{number}-512-source-buckets.tar"
    pinned_receipt, pinned_tar = pinned[number]
    require(sha256(receipt_path) == pinned_receipt
            and sha256(archive) == pinned_tar,
            "Immutable 93m original source routing group SHA changed: " + number)
    r = json.loads(receipt_path.read_text(encoding="utf-8"))
    require(r.get("schema") == ROUTE_SCHEMA
            and r.get("source_tier") == SOURCE_TIER
            and r.get("source_group") == group
            and r.get("compute_bucket_count") == 512
            and r.get("full_original_source_photo_IDs") == EXPECTED_PHOTOS
            and r.get("native_status_not_hard_filtered") is True
            and r.get("images_downloaded") == 0
            and r.get("flower_pixels_inferred") == 0
            and r.get("source_group_read_exactly_once") is True
            and r.get("no_repeated_source_tier_scan_per_bucket") is True,
            "Signed source group does not cover original photo metadata")
    catalog = {int(x["bucket"]): x for x in r["routed_buckets"]}
    require(len(catalog) == len(r["routed_buckets"])
            and all(0 <= n < 512 for n in catalog),
            "Original observation processing bucket manifest duplicated or malformed")
    count = 0
    with tarfile.open(archive, "r:*") as tf:
        names = {}
        expected = {"routing_receipt.json"} | {
            "bucket_csv/" + x["source_file"] for x in catalog.values()
        }
        for member in tf.getmembers():
            if not member.isfile():
                require(member.isdir(), "Untrusted original source symlink/special tar member")
                continue
            name = member.name.removeprefix("./")
            require(name not in names and name in expected,
                    "Original photo bucket source tar contains an unexpected or duplicate member")
            names[name] = member
        require(set(names) == expected
                and hashlib.sha256(tf.extractfile(names["routing_receipt.json"]).read()).hexdigest()
                    == pinned_receipt,
                "Routed source tar manifest differs from its signed group receipt")
        for bucket, item in sorted(catalog.items()):
            data = tf.extractfile(names["bucket_csv/" + item["source_file"]]).read()
            require(hashlib.sha256(data).hexdigest() == item["source_csv_gzip_sha256"],
                    "One original source observation-bucket gzip SHA changed")
            n = 0
            with gzip.GzipFile(fileobj=io.BytesIO(data), mode="rb") as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as stream:
                    reader = csv.DictReader(stream)
                    require(tuple(reader.fieldnames or ()) == COLUMNS,
                            "Original 93m photo source metadata columns changed")
                    for row in reader:
                        import_photo(db, row)
                        n += 1
            require(n == int(item["original_photo_IDs"]),
                    "Original photo source bucket row count did not match SHA receipt")
            count += n
    require(count == int(r["original_source_photo_IDs_routed"]),
            "Signed source group original photo-ID coverage changed")
    return {
        "source_group": group,
        "original_photos": count,
        "routing_receipt_SHA256": pinned_receipt,
        "routing_archive_SHA256": pinned_tar,
    }


def create_database(db):
    db.execute("PRAGMA cache_size=-32768")
    db.execute("PRAGMA temp_store=FILE")
    db.execute("""
    CREATE TABLE obs (
      uuid TEXT PRIMARY KEY,
      taxon TEXT NOT NULL,
      record_month TEXT NOT NULL,
      lat_band TEXT NOT NULL,
      obs_rank TEXT NOT NULL,
      photo_count INTEGER NOT NULL,
      best_photo_rank TEXT NOT NULL,
      best_photo_id TEXT NOT NULL,
      photo_license TEXT NOT NULL,
      photo_extension TEXT NOT NULL,
      recorded_on TEXT NOT NULL,
      latitude TEXT NOT NULL,
      longitude TEXT NOT NULL,
      observer_id TEXT NOT NULL,
      positional_accuracy TEXT NOT NULL
    )
    """)


def select_candidates(db, taxon_names, k, max_candidates):
    require(1 <= k <= 20, "Balanced tier06 fixed observation count per taxon/month/latitude invalid")
    selected = []
    audit = Counter()
    # All rows visited only from local SQLite. Observations from the same
    # taxon × month × lat band have identical positive conditional probability.
    # Photos inside one observation have known uniform 1/photo_count chance.
    query = """
      SELECT * FROM (
        SELECT uuid,taxon,record_month,lat_band,photo_count,best_photo_id,
               photo_license,photo_extension,recorded_on,latitude,longitude,
               observer_id,positional_accuracy,obs_rank,
               COUNT(*) OVER (
                 PARTITION BY taxon,record_month,lat_band
               ) AS stratum_n,
               ROW_NUMBER() OVER (
                 PARTITION BY taxon,record_month,lat_band
                 ORDER BY obs_rank,uuid
               ) AS rn
        FROM obs
      )
      WHERE rn <= ?
      ORDER BY taxon,record_month,lat_band,obs_rank,uuid
    """
    for row in db.execute(query,(k,)):
        (uuid, taxon, month, band, m, pid, license_name, extension, on,
         lat, lon, observer_id, accuracy, obs_rank, n_stratum, _rank) = row
        require(taxon in taxon_names and m > 0 and n_stratum > 0
                and m <= SOURCE_PHOTOS and n_stratum <= SOURCE_PHOTOS,
                "Missing original source taxon, photo count, or stratum population")
        q_obs = min(k,n_stratum)/n_stratum
        q_photo = 1/m
        joint = SOURCE_PROBABILITY * q_obs * q_photo
        ext = normalized_extension(extension)
        audited = {
            "photo_id": pid, "observation_id": uuid,
            "inat_taxon_id": taxon, "species": taxon_names[taxon],
            "observed_on": on, "record_month": month,
            "latitude_band": band,"latitude":lat,"longitude":lon,
            "observer_id":observer_id,"positional_accuracy":accuracy,
            "photo_license":license_name,"image_extension":extension,
            "photos_in_original_observation":m,
            "observations_in_original_source_stratum":n_stratum,
            "phase1_observation_inclusion_probability":format(SOURCE_PROBABILITY, ".17g"),
            "phase2_observation_inclusion_probability":format(q_obs, ".17g"),
            "conditional_single_photo_probability":format(q_photo, ".17g"),
            "joint_photo_selection_probability":format(joint, ".17g"),
            "source_arm":SOURCE_ARM,
            "photo_source_disposition":"supported" if ext else "unsupported_original_extension",
            "photo_url_large": (
                f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}"
                if ext else ""
            ),
            "rank_for_worker_mixing":hash64("mois-balanced-source-worker-shuffle-v1",uuid),
        }
        selected.append(audited)
        audit["selected_observations"] += 1
        audit["supported_original_images_selected"] += bool(ext)
        audit["unsupported_original_images_selected"] += not bool(ext)
        audit["selected_observations_in_unknown_month"] += month=="unknown"
        audit["selected_observations_in_unknown_latitude"] += band=="unknown"
        if len(selected)>max_candidates:
            raise ValueError("Too many stratified original image candidates: reduce K or increase documented capacity")
    # Exactly one selected photo per original observation, independently ranked.
    require(len({row["observation_id"] for row in selected}) == len(selected),
            "Original observation repeated in balanced source sample")
    selected.sort(key=lambda row:(row["rank_for_worker_mixing"],int(row["photo_id"])))
    return selected,audit


def plan(registry, route, taxa_file, output, *, k=1,
         max_candidates=MAX_CAP, expected_photos=EXPECTED_PHOTOS,
         expected_taxa=EXPECTED_TAXA,
         expected_tier=SOURCE_PHOTOS, pinned=PINNED,
         check_taxa_sha=True, original_batch_size=256):
    import math

    registry,route,taxa_file,output=map(Path,(registry,route,taxa_file,output))
    require(not output.exists(),"Original 93m metadata probability sample may not be overwritten")
    require(original_batch_size==256,
            "Frozen original FCP process uses exactly 256 photo slots per unit")
    full=json.loads(registry.read_text(encoding="utf-8"))
    require(full.get("schema")==TIER_REGISTRY_SCHEMA
            and full.get("exact_original_photo_IDs_in_all_disjoint_tiers")==expected_photos
            and full.get("original_named_taxa_from_full_source")==expected_taxa
            and int(full.get("new_target_counts_by_tier",{}).get(str(SOURCE_TIER),-1))
                ==expected_tier
            and full.get("native_vs_introduced_not_hard_filtered") is True
            and full.get("source_only_NO_photo_pixels_or_hue_measured") is True,
            "Missing original full 93m, non-capped, native-status-neutral photo target frame")
    taxa=read_taxa(taxa_file,check_sha=check_taxa_sha)
    require(len(taxa)==expected_taxa,
            "Missing exact 42,111 original named plant taxa; no implicit taxon filter allowed")
    output.parent.mkdir(parents=True,exist_ok=True)
    temp=Path(tempfile.mkdtemp(prefix="mois-93m-balanced-source-",dir=output.parent))
    try:
        with sqlite3.connect(temp/"source_observation_multiplicity.sqlite") as db:
            create_database(db)
            groups=[]
            for group in range(SOURCE_GROUPS):
                groups.append(load_signed_group(db,route,group,pinned=pinned))
                db.commit()
            photo_total=sum(x["original_photos"] for x in groups)
            require(photo_total==expected_tier,
                    "Original tier06 photo IDs lost across 16 signed taxon groups")
            observations,multiplicity=db.execute(
                "SELECT COUNT(*), SUM(photo_count) FROM obs"
            ).fetchone()
            require(int(multiplicity)==expected_tier and observations>0,
                    "Original photos per observation do not exhaust immutable source")
            stratum_total=db.execute("""
                SELECT COUNT(*) FROM (
                  SELECT 1 FROM obs GROUP BY taxon,record_month,lat_band
                )
            """).fetchone()[0]
            sampled,counts=select_candidates(db,taxa,k,max_candidates)
        require(len(sampled)==counts["selected_observations"]>0,
                "Empty original photo sample from signed complete source")
        header=[x for x in sampled[0] if x!="rank_for_worker_mixing"]
        all_file=temp/"selected_observation_photo_probability_sample.csv"
        unsupported_file=temp/"unsupported_original_photo_source_outcomes.csv"
        with all_file.open("w",encoding="utf-8",newline="") as f,\
             unsupported_file.open("w",encoding="utf-8",newline="") as u:
            writer=csv.DictWriter(f,fieldnames=header,lineterminator="\n")
            unknown=csv.DictWriter(u,fieldnames=header,lineterminator="\n")
            writer.writeheader()
            unknown.writeheader()
            for row in sampled:
                projected={key:row[key] for key in header}
                writer.writerow(projected)
                if row["photo_source_disposition"]!="supported":
                    unknown.writerow(projected)
        supported=[r for r in sampled if r["photo_source_disposition"]=="supported"]
        require(len(supported)>0,"No supported original photo image formats selected")
        shards=temp/"shards"
        shards.mkdir()
        sha_by_unit={}
        for unit,start in enumerate(range(0,len(supported),original_batch_size)):
            file=shards/f"shard-{unit:03d}.csv"
            with file.open("w",encoding="utf-8",newline="") as o:
                writer=csv.DictWriter(o,fieldnames=WORKER_COLUMNS,lineterminator="\n")
                writer.writeheader()
                for row in supported[start:start+original_batch_size]:
                    pid=row["photo_id"]
                    mid="MOIS-"+hashlib.sha256(
                        ("mois-original-global-RGB-v1|"+pid).encode()
                    ).hexdigest()[:24].upper()
                    writer.writerow({
                        "measurement_id":mid,"photo_id":pid,
                        "inat_taxon_id":row["inat_taxon_id"],"species":row["species"],
                        "observation_id":row["observation_id"],
                        "observer_id":row["observer_id"],"observed_on":row["observed_on"],
                        "latitude":row["latitude"],"longitude":row["longitude"],
                        "positional_accuracy":row["positional_accuracy"],
                        "photo_license":row["photo_license"],
                        "photo_url_large":row["photo_url_large"],
                        "image_filename":mid+".jpg",
                        "source_arm":SOURCE_ARM,
                        "selection_probability":row["joint_photo_selection_probability"],
                        "source_snapshot_date":"2026-09-27","shard":str(unit),
                    })
            sha_by_unit[str(unit)]=sha256(file)
        work_units=len(sha_by_unit)
        summary={
            "schema":SCHEMA,
            "all_upstream_original_photo_IDs":expected_photos,
            "original_named_plant_taxa":expected_taxa,
            "source_original_tier":SOURCE_TIER,
            "original_tier06_source_photo_IDs_accounted":photo_total,
            "original_tier06_source_observations":int(observations),
            "original_tier06_taxon_month_latitude_strata":int(stratum_total),
            "source_probability":SOURCE_PROBABILITY,
            "sampled_per_stratum_K":int(k),
            "second_stage_design":"uniform_ranked_observations_without_replacement",
            "within_observation_photo_design":"independent_uniform_minimum_hash_one_original_photo",
            "original_photos_per_observation_counted_across_all_512_source_buckets":True,
            "source_metadata_group_signatures":groups,
            "source_registry_SHA256":sha256(registry),
            "source_taxon_names_SHA256":sha256(taxa_file),
            "sampled_distinct_original_observations":len(sampled),
            "supported_photo_attempt_targets":len(supported),
            "unsupported_original_extension_selected_as_technical_unknown":
                len(sampled)-len(supported),
            "selected_original_observations_missing_recorded_month":
                counts["selected_observations_in_unknown_month"],
            "selected_original_observations_missing_latitude":
                counts["selected_observations_in_unknown_latitude"],
            "expected_total_observation_sample_if_K1":int(stratum_total) if k==1 else None,
            "original_source_sample_full_table_SHA256":sha256(all_file),
            "unsupported_original_source_disposition_SHA256":sha256(unsupported_file),
            "photo_shard_sha256":sha_by_unit,
            "work_units":work_units,
            "group_ids":list(range(math.ceil(work_units/4))),
            "worker_ids":list(range(work_units)),
            "original_model_batch_size":original_batch_size,
            "units_per_group":4,
            "source_arm":SOURCE_ARM,
            "source_sample_decided_without_photo_pixels":True,
            "all_source_observations_have_strictly_positive_phase2_probability":True,
            "all_photo_ids_have_strictly_positive_conditional_photo_selection_probability":True,
            "probability_multiplication_valid_for_nested_observation_photo_design":True,
            "all_original_photo_opportunities_retained_upstream_without_native_filter":True,
            "no_original_image_downloads":True,
            "no_ROI_colour_inference":True,
            "no_botanical_anthesis_or_individual_wild_native_validation":True,
            "not_a_complete_world_flora_photo_hue_dataset":True,
            "native_and_introduced_species_strata_must_be_added_non_destructively":True,
            "missing_image_ROI_nonresponse_must_be_adjusted_separately":True,
            "same_source_arm_photos_cannot_be_naively_pooled_with_other_designs":True,
        }
        (temp/"sampling_receipt.json").write_text(
            json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        # The SQLite source photo/observation index is a temporary calculation
        # aid: outputs already conserve every source photo count and strata.
        (temp/"source_observation_multiplicity.sqlite").unlink(missing_ok=True)
        os.replace(temp,output)
        print("MOIS_ORIGINAL93M_BALANCED_PROBABILITY_PHOTO_SAMPLE",
              json.dumps({k:summary[k] for k in (
                  "original_tier06_source_photo_IDs_accounted",
                  "original_tier06_source_observations",
                  "original_tier06_taxon_month_latitude_strata",
                  "sampled_distinct_original_observations",
                  "supported_photo_attempt_targets",
                  "unsupported_original_extension_selected_as_technical_unknown",
                  "work_units")},sort_keys=True),flush=True)
        return summary
    except BaseException:
        shutil.rmtree(temp,ignore_errors=True)
        raise


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ("registry","route","taxa-file","output"):
        ap.add_argument("--"+key,type=Path,required=True)
    ap.add_argument("--k",type=int,default=1)
    ap.add_argument("--max-candidates",type=int,default=MAX_CAP)
    opts=ap.parse_args()
    plan(opts.registry,opts.route,opts.taxa_file,opts.output,
         k=opts.k,max_candidates=opts.max_candidates)


if __name__=="__main__":
    main()
