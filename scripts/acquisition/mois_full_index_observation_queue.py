#!/usr/bin/env python3
"""Original-source-wide observation probability draw; no taxa/photo inclusion cap.

Use the independently signed 93,209,351 source PHOTO-ID and observation
catalogue. A 0<p<=1 nested observation-level hash threshold determines the
current batch; all unselected source photos remain eligible for higher p.
Every original photo from each selected observation enters the pre-colour
source (unless the photo ID was already genuinely measured). The output is
a reproducible image-acquisition queue, NOT any flower-colour data.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math
from pathlib import Path

SCHEMA="mois_full_42111_93209351_uncapped_original_photo_observation_index_v1"
OUT_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
COLUMNS=["measurement_id","photo_id","inat_taxon_id","species","observation_id",
         "observer_id","observed_on","latitude","longitude","positional_accuracy",
         "photo_license","photo_url_large","image_filename","source_arm",
         "selection_probability","source_snapshot_date","shard"]
SEED="mois-global-observation-v1"
N_PHOTOS=93209351
N_TAXA=42111
N_SHARDS=128
def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()
SUPPORTED_FROZEN_IMAGE_EXTENSIONS=frozenset({"jpg","jpeg","png"})
def frozen_original_extension(raw):
    """Accept only formats supported by the SHA-pinned FCP original image loader.

    Case/dot/whitespace differences may be normalized only to the original
    allowed extensions. Unsupported source rows are sampled and retained in
    a SOURCE-ONLY sidecar, not silently dropped or given invented RGB hues.
    """
    value=str(raw or "").strip().lower().lstrip(".")
    return value if value in SUPPORTED_FROZEN_IMAGE_EXTENSIONS else None

def obs_draw(uuid):
    value=hashlib.sha256((SEED+"|"+str(uuid)).encode()).digest()[:8]
    return int.from_bytes(value,"big")/2**64
def csvrows(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        yield from csv.DictReader(f)
def gzrows(path):
    with gzip.open(path,"rt",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f)
def known_ids(paths):
    ids=set()
    for path in paths:
        for r in csvrows(path):
            pid=str(r["photo_id"])
            require(pid.isdecimal(),"Previously processed photo ID malformed")
            ids.add(pid)
    return ids
def prepare(index,counts,prior_files,out,*,p=0.0002,batch_size=256):
    index=Path(index);out=Path(out)
    require(0<p<=1,"Observation inclusion probability must be in (0,1]")
    require(1<=batch_size<=256,"Original ROI size must be <=256")
    require(not out.exists() or not any(out.iterdir()),"Source output must not overwrite frozen measurements")
    source=json.loads((index/"manifest.json").read_text())
    require(source["schema"]==SCHEMA and source["photo_caps"] is None
            and source["photos_emitted"]==N_PHOTOS
            and source["taxa_in_frame"]==N_TAXA
            and not source["photo_colours_measured"]
            and len(source["shards"])==N_SHARDS,
            "Must use the complete authenticated original PHOTO IDs, not a count-only or selected source")
    names={int(x["inat_taxon_id"]):x["species"] for x in csvrows(counts)}
    require(len(names)==N_TAXA,"All original source taxon names required")
    previous=known_ids(prior_files)
    seen=set(); selected=[]; unsupported=[]; reused=0; total=0; candidate_n=0; nobs=0
    for item in source["shards"]:
        i=int(item["part"])
        photos=index/"photos"/f"part-{i:03d}.csv.gz"
        observations=index/"observations"/f"part-{i:03d}.csv.gz"
        require(sha(photos)==item["photos_sha256"] and sha(observations)==item["observations_sha256"],
                "Original photo/observation archive shard SHA mismatch")
        picks=[];wanted=set(); scanned=0
        for row in gzrows(photos):
            scanned+=1
            # Deterministic per-OBSERVATION selection (not one image per species).
            if obs_draw(row["observation_uuid"])>=p:continue
            candidate_n+=1
            pid=row["photo_id"];tax=int(row["inat_taxon_id"])
            require(pid.isdecimal() and tax in names and pid not in seen,
                    "Duplicate original photo ID or unknown taxon in selected source")
            seen.add(pid)
            if pid in previous:
                reused+=1
                continue
            # The source population is still the ENTIRE signed 93m catalogue.
            # The prior frozen FCP image loader can decode only jpg/jpeg/png.
            # Never pretend that an unsupported (or absent) file extension
            # is a photographed flower colour or abort all compatible photos.
            normalised=frozen_original_extension(row["image_extension"])
            if normalised is None:
                unsupported.append({
                    "photo_id":pid,
                    "observation_id":row["observation_uuid"],
                    "inat_taxon_id":str(tax),
                    "source_image_extension":row["image_extension"],
                    "reason":"unsupported_by_frozen_original_FCP_image_loader",
                })
                continue
            wanted.add(row["observation_uuid"])
            picks.append({**row,"image_extension":normalised})
        require(scanned==item["photo_count"],"Original photo-count shard mismatch")
        total+=scanned
        if not picks:continue
        obs={};obs_scanned=0
        for row in gzrows(observations):
            obs_scanned+=1
            if row["observation_uuid"] in wanted:
                require(row["observation_uuid"] not in obs,"Duplicate observation in original source")
                obs[row["observation_uuid"]]=row
        require(obs_scanned==item["observation_count"] and set(obs)==wanted,
                "Original selected photos lack source observation date/site")
        nobs+=len(obs)
        for row in picks:
            o=obs[row["observation_uuid"]]
            require(row["inat_taxon_id"]==o["inat_taxon_id"],
                    "Original photo taxon and observation taxon mismatch")
            pid=row["photo_id"];tax=row["inat_taxon_id"]
            mid="MOIS-"+hashlib.sha256(("mois-original-global-RGB-v1|"+pid).encode()).hexdigest()[:24].upper()
            selected.append({
                "measurement_id":mid,"photo_id":pid,
                "inat_taxon_id":tax,"species":names[int(tax)],
                "observation_id":row["observation_uuid"],
                "observer_id":o["observer_id"],"observed_on":o["observed_on"],
                "latitude":o["latitude"],"longitude":o["longitude"],
                "positional_accuracy":o["positional_accuracy"],
                "photo_license":row["photo_license"],
                "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{row['image_extension']}",
                "image_filename":mid+".jpg",
                "source_arm":"global93m_probability_sample_observations",
                "selection_probability":str(float(p)),
                "source_snapshot_date":"2026-09-27",
                "shard":"",
            })
    require(total==source["photos_emitted"]
            and len(selected)+reused+len(unsupported)==candidate_n,
            "Original indexed photo source sampled/reused/unsupported accounting failed")
    require(len({r["photo_id"] for r in unsupported})==len(unsupported)
            and not ({r["photo_id"] for r in unsupported} &
                     {r["photo_id"] for r in selected}),
            "Source photo ID repeated between supported and unsupported image formats")
    # Randomise processing order independent of input filename/species grouping.
    selected.sort(key=lambda x:hashlib.sha256(("mois-prioritise-observation|"+x["photo_id"]).encode()).digest())
    out.mkdir(parents=True,exist_ok=True)
    (out/"shards").mkdir()
    # Keep the complete original observation-level draw auditable even where
    # FCP's original frozen network/model cannot process that source format.
    # JSONL is source-only (no image bytes or inferred flower colour).
    unsupported_file=out/"unsupported_selected_original_photo_formats.jsonl"
    with unsupported_file.open("w",encoding="utf-8",newline="\n") as f:
        for row in sorted(unsupported,key=lambda r:int(r["photo_id"])):
            f.write(json.dumps(row,sort_keys=True,ensure_ascii=False)+"\n")
    units=math.ceil(len(selected)/batch_size)
    require(units>0,"Initial precolour observation source contains no new photographs")
    for n,r in enumerate(selected):r["shard"]=str(n//batch_size)
    checksums={}
    for unit in range(units):
        dest=out/"shards"/f"shard-{unit:03d}.csv"
        with dest.open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
            writer.writeheader()
            writer.writerows(selected[unit*batch_size:(unit+1)*batch_size])
        checksums[str(unit)]=sha(dest)
    receipt={
        "schema":OUT_SCHEMA,
        "full_source_manifest_sha256":sha(index/"manifest.json"),
        "full_frame_photo_ids":total,
        "sampled_photo_ids_before_prior_dedup":candidate_n,
        "previously_measured_photo_ids_skipped":reused,
        "selected_new_photos_unsupported_by_frozen_loader":len(unsupported),
        "selected_new_photo_formats_unsupported_sha256":sha(unsupported_file),
        "unsupported_original_source_photo_ext_counts":{
            ext:n for ext,n in sorted(__import__("collections").Counter(
                r["source_image_extension"] for r in unsupported).items())
        },
        "supported_frozen_source_image_extensions":sorted(SUPPORTED_FROZEN_IMAGE_EXTENSIONS),
        "new_photo_attempts_planned":len(selected),
        "selected_observations_with_new_photos":nobs,
        "selection_probability":p,
        "selection_unit":"observation_uuid",
        "sampled_before_any_hue_measurement":True,
        "work_units":units,"worker_ids":list(range(units)),
        "group_ids":list(range(math.ceil(units/4))),
        "units_per_group":4,"original_model_batch_size":batch_size,
        "photo_shard_sha256":checksums,
        "all_unselected_source_photos_eligible_at_higher_probability":True,
        "no_per_taxon_or_observation_photo_count_cap":True,
        "source_original_images_measured":False,
        "colour_measured":False,"ecological_analysis_run":False,
    }
    (out/"manifest.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_GLOBAL_FULL_INDEX_PROBABILITY_SAMPLE_SOURCE_ONLY",json.dumps({
        k:receipt[k] for k in ("full_frame_photo_ids","sampled_photo_ids_before_prior_dedup",
                               "previously_measured_photo_ids_skipped",
                               "selected_new_photos_unsupported_by_frozen_loader",
                               "unsupported_original_source_photo_ext_counts",
                               "new_photo_attempts_planned","work_units","selection_probability")
    },sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--index",type=Path,required=True)
    p.add_argument("--counts",type=Path,required=True)
    p.add_argument("--previous",type=Path,action="append",default=[])
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--probability",type=float,default=0.0002)
    a=p.parse_args()
    prepare(a.index,a.counts,a.previous,a.output,p=a.probability)
