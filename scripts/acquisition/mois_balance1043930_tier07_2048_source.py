#!/usr/bin/env python3
"""MOIS acquisition designed to REDUCE existing global flower-colour sample bias.

Use the COMPLETE source-signed 93,209,351 original-photo index, but read only
its one-time routed original tier07 (~1.19m metadata candidates). NEVER again
re-download 93m source archives or over-sample taxa merely because they lack a
colour result. Pre-colour stage uses the 400,710 source-authentic ROI9 photos
to allocate NEW image measurements towards poorly measured hemisphere ×
latitude/biogeographic band × local recorded season groups.

Two SEPARATE source cohorts:
  - A 256-photo independent original-observation-uniform REFERENCE cohort
    drawn from distinct, decodable, previously unattempted observations,
    to evaluate photo-source and detector nonresponse without choosing colour.
  - A 1,792-photo geographic/season gap-targeted AUGMENTATION with one
    independently hashed original observation per taxon×zone×local season.
    Allocate strata by source support, prior strict ROI9 coverage and original
    ROI success probability; then select species uniformly within each stratum.
    These are NOT a probability sample of native world plant abundance and
    MUST NOT be naively pooled with original 1043930 source colour ratios.

All original photo IDs/observation UUIDs are disjoint from the 1,043,930 photo
attempts and mutually unique in the new source. Keep all original 93m photo
IDs eligible for later processing. Predeclare hashes, sampling-stage sizes,
cell counts/quotas and conditional probability labels BEFORE image fetch.
The source archive is pre-colour metadata only. Native/island/flower organ
status, CIELAB, 9-palette values and future ROI success are not inferred.

One 93m-source tier scan, streaming old-photo membership and temporary
SQLite observation aggregation; no source photo pixels fetched in this step.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,math,sqlite3,tarfile
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS,RESULT_SCHEMA,TIER_ENDS,obs_draw
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS

SCHEMA="mois_tier07_source_actual1043930_regional_seasonal_bias_repair_twoarm_v1"
OLD_PHOTO_SHA="df0f562b58bdc7ba06e901530eb809f513a763e2b0524d3b195e99ee56e930ef"
OLD_OBS_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
PARENT_N=93209351
N_OLD=1043930
N_OLD_STRICT=400710
TIER=7
SOURCE_P=TIER_ENDS[TIER]-TIER_ENDS[TIER-1]
SEED="mois-20261026-season-and-hemisphere-gap-repair-tier07-independent-v1"
REFERENCE=256
AUGMENTATION=1792
TOTAL=REFERENCE+AUGMENTATION
SOURCE_ARMS={
 "reference":"mois_tier07_geo_season_photo_reference_original_observation_v1",
 "targeted":"mois_tier07_1043930_photo_bias_gap_adaptive_species_season_v1",
}
ZONE=["N_tropical","N_temperate","N_subpolar","N_arctic",
      "S_tropical","S_temperate","S_subpolar","S_antarctic"]
SEASONS=["spring","summer","autumn","winter"]
FRAC=("white","yellow","orange","red","pink","magenta",
      "purple","blue","bronze")

def require(condition,msg):
    if not condition:raise ValueError("MOIS_BIAS_REPAIR_SOURCE: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for c in iter(lambda:f.read(4<<20),b""):h.update(c)
    return h.hexdigest()

def rank(kind,*keys):
    return hashlib.sha256(("|".join((SEED,kind,*map(str,keys)))).encode()).hexdigest()

def zone(lat):
    try: x=float(lat)
    except (ValueError,TypeError):return "unknown"
    if not math.isfinite(x) or abs(x)>90:return "unknown"
    if abs(x)>=66.5622:suffix="arctic" if x>=0 else "antarctic"
    elif abs(x)>=55:suffix="subpolar"
    elif abs(x)>=23.5:suffix="temperate"
    else:suffix="tropical"
    return ("N_" if x>=0 else "S_")+suffix

def season(lat,date):
    hemisphere=zone(lat)
    if hemisphere=="unknown":return "unknown"
    value=str(date or "")
    try:
        import datetime
        d=datetime.date.fromisoformat(value)
        if len(value)!=10:return "unknown"
    except (TypeError,ValueError):return "unknown"
    month=d.month
    north=hemisphere.startswith("N_")
    if month in (12,1,2):return "winter" if north else "summer"
    if month in (3,4,5):return "spring" if north else "autumn"
    if month in (6,7,8):return "summer" if north else "winter"
    return "autumn" if north else "spring"

def numbers(frame):
    old=frame.copy()
    old["region"]=[zone(x) for x in old.latitude.to_numpy()]
    old["local_season"]=[season(lat,day) for lat,day in zip(
        old.latitude.to_numpy(),old.date.fillna("").to_numpy())]
    groups=old.groupby(["region","local_season"],dropna=False,observed=True).agg(
        prior_photo_attempts=("photo_id","size"),
        prior_strict9_ROIs=("continuous9_strict","sum"))
    return {
        (str(z),str(t)):{"attempts":int(row.prior_photo_attempts),
                        "strict":int(row.prior_strict9_ROIs)}
        for (z,t),row in groups.iterrows()
    }

def quotas(source_species_counts,old_counts,*,total=AUGMENTATION,cap=160,
           floor=8):
    """Bounded coverage allocation; NO outcomes of the NEW photos are used."""
    require(total>0 and cap>=floor>0,"Invalid acquisition budget")
    available={k:int(v) for k,v in source_species_counts.items() if
               k[0] in ZONE and k[1] in SEASONS and int(v)>0}
    require(sum(min(cap,n) for n in available.values())>=total,
            "The original source does not support requested balancing budget")
    allot={k:min(n,floor) for k,n in available.items()}
    remaining=total-sum(allot.values())
    require(remaining>=0,"Overlarge per-stratum guaranteed lower-bound quota")
    # Number of eligible distinct source taxa is used instead of source images
    # because one targeted source photograph/taxon×zone×season is admitted.
    # Existing strict9 support and source-photo nonresponse affect the gap
    # score, with bounded pseudocounts to avoid 0/0 and runaway rare strata.
    weights={}
    for k,n in available.items():
        observed=old_counts.get(k,{})
        attempts=float(observed.get("attempts",0))
        strict=float(observed.get("strict",0))
        estimated_p=max(.075,min(.9,(strict+8.)/(attempts+20.)))
        weights[k]=(n**.55)/(max(strict,20.)**.35*math.sqrt(estimated_p))
    # Deterministic bounded proportional water-filling; capped strata stop
    # collecting rather than redirecting into already well measured north.
    while remaining:
        choices=[k for k in allot if allot[k]<min(cap,available[k])]
        require(bool(choices),"No supported source stratum for remaining quota")
        k=max(choices,key=lambda a:(weights[a]/(1+allot[a]),-ZONE.index(a[0]),
                                      -SEASONS.index(a[1])))
        allot[k]+=1
        remaining-=1
    require(sum(allot.values())==total,
            "Source metadata allocation did not exhaust the fixed actual budget")
    return allot,weights

UPSERT="""
INSERT INTO obs (
  uuid,taxon,observed_on,latitude,longitude,accuracy,observer,license,
  photo_id,extension,photo_rank,photo_n,zone,season,reference_rank,observation_rank,
  part
) VALUES (?,?,?,?,?,?,?,?,?,?,?,1,?,?,?,?,?)
ON CONFLICT(uuid) DO UPDATE SET
  photo_n=obs.photo_n+1,
  photo_rank=CASE WHEN excluded.photo_rank<obs.photo_rank
                  THEN excluded.photo_rank ELSE obs.photo_rank END,
  photo_id=CASE WHEN excluded.photo_rank<obs.photo_rank
                THEN excluded.photo_id ELSE obs.photo_id END,
  extension=CASE WHEN excluded.photo_rank<obs.photo_rank
                 THEN excluded.extension ELSE obs.extension END,
  license=CASE WHEN excluded.photo_rank<obs.photo_rank
               THEN excluded.license ELSE obs.license END
WHERE obs.taxon=excluded.taxon
  AND obs.observed_on=excluded.observed_on
  AND obs.latitude=excluded.latitude AND obs.longitude=excluded.longitude
  AND obs.accuracy=excluded.accuracy AND obs.observer=excluded.observer
  AND obs.zone=excluded.zone AND obs.season=excluded.season
"""
def create_db(db):
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("PRAGMA temp_store=FILE")
    db.execute("PRAGMA cache_size=-262144")
    db.execute("""CREATE TABLE obs(
      uuid TEXT PRIMARY KEY,taxon TEXT NOT NULL,observed_on TEXT,
      latitude TEXT,longitude TEXT,accuracy TEXT,observer TEXT,
      license TEXT,photo_id TEXT,extension TEXT,photo_rank TEXT,
      photo_n INTEGER,zone TEXT,season TEXT,reference_rank TEXT,
      observation_rank TEXT,part INTEGER)""")

def read_source_group(db,root,index,master,used_photos,used_obs,counters):
    tar_path=root/f"tier-07-group-{index:02d}.tar"
    receipt_path=root/f"group-{index:02d}-receipt.json"
    doc=json.loads(receipt_path.read_text())
    require(doc.get("schema")==RESULT_SCHEMA
            and doc.get("source_snapshot_SHA256")==master["full_source_manifest_sha256"]
            and sha(tar_path)==master["published_tier_archive_SHA256"][tar_path.name]
            ==doc["archive_sha256_by_tier"][tar_path.name],
            "Signed 93m tier07 source archive or original group receipt mutated")
    parts={int(a["part"]):a for a in doc["parts"]}
    require(len(parts)==8 and all(i//8==index for i in parts),
            "The source group must cover exactly eight disjoint original parts")
    c=0
    with tarfile.open(tar_path) as archive:
        members={m.name.removeprefix("./"):m for m in archive.getmembers() if m.isfile()}
        expected={f"tier-07/part-{k:03d}.csv.gz" for k in parts}
        require(set(members)==expected,"Tier07 archive contains unexpected source members")
        for part,proof in sorted(parts.items()):
            data=archive.extractfile(members[f"tier-07/part-{part:03d}.csv.gz"]).read()
            require(hashlib.sha256(data).hexdigest()==
                    proof["ranked_target_csv_sha256"][str(TIER)],
                    "Original source photo/observation byte SHA mismatch")
            lines=0
            with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
                with io.TextIOWrapper(gz,encoding="utf-8",newline="") as stream:
                    reader=csv.DictReader(stream)
                    require(tuple(reader.fieldnames or [])==COLUMNS,
                            "Original source photo/observation schema changed")
                    for row in reader:
                        lines+=1
                        pid,uuid,tax=row["photo_id"],row["observation_uuid"],row["inat_taxon_id"]
                        require(pid.isdecimal() and tax.isdecimal() and uuid
                                and row["tier"]==str(TIER)
                                and float(row["observation_inclusion_score"])==obs_draw(uuid),
                                "Original 93m source ID or observation randomization corrupted")
                        counters["source_tier07_original_photos"]+=1
                        if pid in used_photos:
                            counters["previous_photo_IDs_excluded"]+=1
                            continue
                        if uuid in used_obs:
                            counters["previous_observation_IDs_excluded"]+=1
                            continue
                        ext=row["image_extension"].lower().strip().lstrip(".")
                        if row["frozen_loader_supported"]!="true" or ext not in ("jpeg","jpg","png"):
                            counters["unsupported_image_extension"]+=1
                            continue
                        z=zone(row["latitude"])
                        t=season(row["latitude"],row["observed_on"])
                        ph=rank("photo",uuid,pid)
                        vals=(uuid,tax,row["observed_on"],row["latitude"],row["longitude"],
                              row["positional_accuracy"],row["observer_id"],row["photo_license"],
                              pid,ext,ph,z,t,rank("reference",uuid),
                              rank("independent_observation",uuid),part)
                        db.execute(UPSERT,vals)
                        counters["new_supported_original_photo_IDs"]+=1
            require(lines==proof["counts_by_disjoint_tier"][str(TIER)],
                    "Original tier07 part metadata photo count mismatch")
            c+=lines
    require(c==doc["counts_by_disjoint_tier"][str(TIER)],
            "Original tier07 group signed photo count mismatch")
    db.commit()
    return {"group":index,"source_photos":c,"source_tar_SHA256":sha(tar_path)}

FIELDS=("uuid","taxon","observed_on","latitude","longitude","accuracy",
        "observer","license","photo_id","extension","photo_rank",
        "photo_n","zone","season","reference_rank","observation_rank","part")
def obs_record(row):
    return dict(zip(FIELDS,row))

def sample(db,prior,*,reference=REFERENCE,targeted=AUGMENTATION):
    require(reference>=32 and targeted>=64,
            "A balanced study needs a reference group and meaningful supplementation")
    db.execute("CREATE INDEX obs_taxon_zone_season ON obs(taxon,zone,season)")
    n_obs=int(db.execute("SELECT COUNT(*) FROM obs").fetchone()[0])
    require(n_obs>=reference+targeted,"Original 93m source has insufficient new observations")
    # Stage A: uniform reference original OBSERVATIONS from all remaining supported
    # photo-bearing observations, regardless of their colour/region/season.
    base=list(db.execute(
        "SELECT "+",".join(FIELDS)+" FROM obs ORDER BY reference_rank,uuid LIMIT ?",
        (reference,)))
    chosen_ref=[obs_record(x) for x in base]
    n_ref=len(chosen_ref)
    require(n_ref==reference,"Reference photo cohort truncated")
    db.execute("CREATE TEMP TABLE ref (uuid TEXT PRIMARY KEY)")
    db.executemany("INSERT INTO ref VALUES (?)",[(v["uuid"],) for v in chosen_ref])
    # Stage B: one independent original observation per taxon×latitude region×
    # local season, then one hashed taxon per cell with no outcome-based choice.
    query="""SELECT * FROM (
        SELECT """+",".join("o."+x for x in FIELDS)+""",
        COUNT(*) OVER(PARTITION BY o.taxon,o.zone,o.season) AS taxon_cell_observations,
        ROW_NUMBER() OVER(PARTITION BY o.taxon,o.zone,o.season
                          ORDER BY o.observation_rank,o.uuid) AS ord
        FROM obs o WHERE NOT EXISTS(SELECT 1 FROM ref r WHERE r.uuid=o.uuid)
          AND o.zone!='unknown' AND o.season!='unknown'
    ) WHERE ord=1"""
    pool=defaultdict(list)
    for v in db.execute(query):
        item=obs_record(v[:len(FIELDS)])
        item["taxon_cell_observations"]=int(v[len(FIELDS)])
        pool[(item["zone"],item["season"])].append(item)
    n_source={key:len(values) for key,values in pool.items()}
    quota,weights=quotas(n_source,prior,total=targeted)
    selected_target=[]
    cell_counts=[]
    for z in ZONE:
        for ss in SEASONS:
            key=(z,ss);items=pool.get(key,[])
            quota_n=quota.get(key,0)
            items.sort(key=lambda a:rank("taxon_in_stratum",z,ss,a["taxon"]))
            chosen=items[:quota_n]
            for item in chosen:
                item["cohort"]="targeted"
                item["stratum_species_source_size"]=len(items)
                item["stratum_photo_quota"]=quota_n
                selected_target.append(item)
            pre=prior.get(key,{"attempts":0,"strict":0})
            cell_counts.append({
                "region":z,"recorded_local_season":ss,
                "source_remaining_photo_bearing_taxa":len(items),
                "prior_original_photo_attempts":pre["attempts"],
                "prior_strict_roi9_measured":pre["strict"],
                "prior_original_photo_roi9_success_fraction":(
                     pre["strict"]/pre["attempts"] if pre["attempts"] else None),
                "targeted_sample_photo_quota":quota_n,
                "source_balance_coverage_weight":weights.get(key),
                "conditional_taxon_draw_prob_if_stage2_source_fixed":(
                    quota_n/len(items) if items else None),
                "selected_original_observation_photo_count":len(chosen),
            })
    for item in chosen_ref:
        item["cohort"]="reference"
        item["stratum_species_source_size"]=None
        item["stratum_photo_quota"]=None
    selected=chosen_ref+selected_target
    require(len(selected)==reference+targeted
            and len({x["uuid"] for x in selected})==len(selected)
            and len({x["photo_id"] for x in selected})==len(selected),
            "Reference + targeted original photos are not disjoint")
    source_log={
        "eligible_disjoint_original_photo_observations":n_obs,
        "reference_uniform_original_observation_count":reference,
        "reference_conditional_probability_per_eligible_observation":reference/n_obs,
        "targeted_original_photo_count":targeted,
        "targeted_minimum_8_and_maximum_160_per_supported_region_season":True,
        "targeted_strata_actual_source_and_old9_roi9_support":cell_counts,
        "targeted_one_observation_per_taxon_zone_recorded_season":True,
        "source_original_strata_without_candidates_are_not_assigned_fake_photos":True,
    }
    return selected,source_log

def plan(root,registry,atlas,atlas_receipt,union,union_receipt,taxa,out,
         reference=REFERENCE,targeted=AUGMENTATION):
    root,registry,atlas,atlas_receipt,union,union_receipt,taxa,out=map(
        Path,(root,registry,atlas,atlas_receipt,union,union_receipt,taxa,out))
    require(not out.exists(),"Refuse overwrite of prior bias-correcting source design")
    meta=json.loads(registry.read_text(encoding="utf-8"))
    newer=json.loads(atlas_receipt.read_text(encoding="utf-8"))
    prior=json.loads(union_receipt.read_text(encoding="utf-8"))
    require(meta.get("schema")==
            "mois_complete_full93m_original_photo_target_registry_v1"
            and meta.get("exact_original_photo_IDs_in_all_disjoint_tiers")==PARENT_N
            and meta.get("original_named_taxa_from_full_source")==42111
            and len(meta.get("published_tier_archive_SHA256",{}))>=16,
            "No signed complete 93m original-photo source")
    require(newer.get("new_source_photo9_lake_SHA256")==OLD_PHOTO_SHA==sha(atlas)
            and newer.get("new_original_strict9_photo_count")==N_OLD_STRICT
            and newer.get("full_original_photo_attempts")==N_OLD
            and prior.get("source_pixel_attempts")==N_OLD
            and prior.get("photo_parquet_sha256")==OLD_OBS_SHA==sha(union),
            "Unrecognised original measured photo and original observation sources")
    history=pd.read_parquet(atlas,columns=[
        "photo_id","inat_taxon_id","latitude","date","continuous9_strict"])
    known=pd.read_parquet(union,columns=["photo_id","observation_id"])
    history.photo_id=history.photo_id.astype(str)
    known.photo_id=known.photo_id.astype(str)
    require(len(history)==len(known)==N_OLD
            and history.photo_id.is_unique and known.photo_id.is_unique
            and history.photo_id.tolist()==known.photo_id.tolist()
            and int(history.continuous9_strict.sum())==N_OLD_STRICT,
            "Original photo IDs and strict real ROI9 measures disagree")
    original_photos=set(known.photo_id)
    original_observations=set(known.observation_id.fillna("").astype(str))
    original_observations.discard("")
    frame=pd.read_csv(taxa,compression="gzip",dtype={"inat_taxon_id":str})
    require(len(frame)==42111 and frame.inat_taxon_id.is_unique,
            "Use full unfiltered 42,111 source plant taxa, NOT only colour-failed species")
    species=frame.set_index("inat_taxon_id").species.to_dict()
    old= numbers(history)
    # SQLite temp is in work directory and never uploaded in release or photo
    # pixels; streaming tar entry validation avoids 93m repeated rescans.
    scratch=out.parent/(out.name+"_temporary.sqlite")
    require(not scratch.exists(),"Cannot reuse unfinished old photo source sqlite")
    receipts=[];audit=Counter()
    try:
        with sqlite3.connect(scratch) as db:
            create_db(db)
            for group in range(16):
                receipts.append(read_source_group(
                    db,root,group,meta,original_photos,original_observations,audit))
            require(audit["source_tier07_original_photos"]==
                    meta["new_target_counts_by_tier"][str(TIER)],
                    "Not all original tier07 source photo opportunities visited")
            selected,metrics=sample(db,old,reference=reference,targeted=targeted)
        require(len(selected)==reference+targeted,"No full source balance sample")
    finally:
        scratch.unlink(missing_ok=True)
    # Colour-blind IDs/stage2 selection; record conditional probabilities and
    # source roles. 'selection_probability' inside legacy ROI collector is
    # conditional on the specific previous-photo exclusions, NOT a known
    # unconditional 93m photo inclusion rate.
    output=[];detail=[]
    for p in selected:
        role=p["cohort"]
        taxon=p["taxon"];pid=p["photo_id"]
        sregion=p["zone"];sseason=p["season"]
        if role=="reference":
            conditional=metrics["reference_conditional_probability_per_eligible_observation"]
        else:
            conditional=(p["stratum_photo_quota"]/
                        p["stratum_species_source_size"]/
                        p["taxon_cell_observations"])
        # One photo drawn uniformly from all supported attached photos within
        # the already sampled source-tier original observation.
        cp=conditional/p["photo_n"]
        require(0<cp<=1 and taxon in species
                and pid not in original_photos and p["uuid"] not in original_observations,
                "Outcome blind original photo/observation candidate invalid")
        mhash=rank("measurement_ID",pid).upper()[:24]
        photo_url=(
            f"https://inaturalist-open-data.s3.amazonaws.com/photos/"
            f"{pid}/large.{p['extension']}")
        result=dict(
            measurement_id="MOIS-"+mhash,photo_id=pid,
            inat_taxon_id=taxon,species=species[taxon],
            observation_id=p["uuid"],observer_id=p["observer"],
            observed_on=p["observed_on"],latitude=p["latitude"],
            longitude=p["longitude"],positional_accuracy=p["accuracy"],
            photo_license=p["license"],photo_url_large=photo_url,
            image_filename="MOIS-"+mhash+".jpg",
            source_arm=SOURCE_ARMS[role],
            selection_probability=format(cp,".17g"),
            source_snapshot_date="2026-09-27",shard="",
        )
        output.append(result)
        detail.append({
            "photo_id":pid,"original_observation_id":p["uuid"],
            "inat_taxon_id":taxon,
            "design_cohort":role,
            "latitude_region":sregion,"recorded_local_season":sseason,
            "source_5degree_grid_cell":(
                 f"{int(math.floor((float(p['latitude'])+90)/5))}:"+
                 f"{int(math.floor((float(p['longitude'])+180)/5))}"
                 if sregion!="unknown" and p["longitude"] not in ("","nan") else "unknown"),
            "source_photos_per_original_observation":p["photo_n"],
            "source_original_tier07_probability_width":SOURCE_P,
            "conditional_stage2_observation_and_photo_probability":cp,
            "conditional_source_design_prob_NOT_unconditional_population_IPW":True,
            "source_previous_1043930_photo_and_observation_overlap":False,
            "targeted_source_taxa_in_region_season":p["stratum_species_source_size"],
            "targeted_photo_quota":p["stratum_photo_quota"],
            "source_taxon_region_season_independent_observations":(
                p.get("taxon_cell_observations")),
        })
    # Reference is a separate outcome-independent diagnostic source and must
    # not be pooled into the targeted ecological white-frequency estimator.
    output_and_details=list(zip(output,detail))
    output_and_details.sort(key=lambda q:rank("warm_worker_shuffle",q[0]["photo_id"]))
    output=[x for x,_ in output_and_details]
    detail=[x for _,x in output_and_details]
    require(len({x["measurement_id"] for x in output})==len(output)
            and len({x["photo_id"] for x in output})==len(output)
            and len({x["observation_id"] for x in output})==len(output),
            "Two distinct source images share a unique photo/observation identity")
    out.mkdir(parents=True)
    (out/"shards").mkdir()
    n_units=math.ceil(len(output)/256)
    SHAs={}
    for n,start in enumerate(range(0,len(output),256)):
        dest=out/"shards"/f"shard-{n:03d}.csv"
        with dest.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=WORKER_COLUMNS,lineterminator="\n")
            writer.writeheader()
            for v in output[start:start+256]:
                writer.writerow({**v,"shard":str(n)})
        SHAs[str(n)]=sha(dest)
    with (out/"selected_original_photo_design_provenance.csv").open(
           "w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(detail[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(detail)
    counts=Counter(x["source_arm"] for x in output)
    req={
        "schema":"mois_full_source_observation_sampling_cohort_v1",
        "sampling_schema":SCHEMA,
        "colour_measured":False,
        "new_photo_attempts_planned":len(output),
        "original_photo_source_total_candidate_photos":PARENT_N,
        "full_frame_photo_ids":PARENT_N,
        "source_tier":TIER,
        "source_probability_width":SOURCE_P,
        "all_original_source_ids_remain_eligible_for_future_photo_runs":True,
        "source_region_season_balanced_but_not_random_global_plants":True,
        "original_1043930_source_colour_composition_used_FOR_STRATUM_ALLOCATION_NOT_LABELS":True,
        "source_previous_photo_ID_and_observation_ID_dedup":N_OLD,
        "prior_real_photo_strict9_source":N_OLD_STRICT,
        "original_source_tier07_authentic_photo_opportunities":audit[
            "source_tier07_original_photos"],
        "previous_photo_IDs_skipped_in_tier07":audit["previous_photo_IDs_excluded"],
        "previous_observation_IDs_skipped_in_tier07":audit["previous_observation_IDs_excluded"],
        "unsupported_original_source_photo_extension_count":audit[
            "unsupported_image_extension"],
        "remaining_new_supported_original_photo_candidates":audit[
            "new_supported_original_photo_IDs"],
        "photo_shard_sha256":SHAs,
        "original_model_batch_size":256,
        "worker_ids":list(range(n_units)),
        "work_units":n_units,
        "units_per_group":4,
        "group_ids":list(range(math.ceil(n_units/4))),
        "source_arm_counts":dict(counts),
        "reference_uniform_source_and_targeted_stratum_independently_separated":True,
        "source_new_image_IDs_and_observations_unique_and_prior_disjoint":True,
        "same_recorded_local_season_defined_with_hemisphere":True,
        "sampling_stage_source_details":metrics,
        "original_tier07_tar_verification":receipts,
        "photo_shard_SHA256_used_by_worker_before_any_photo_fetch":True,
        "sampling_conditional_probabilities_in_photo_shards_are_NOT_93m_population_inclusion":True,
        "provider_original_image_bytes_downloaded":0,
        "no_colour_inferred_from_source_photo_metadata":True,
        "new_images_sex_or_photo_native_status_not_guessed":True,
        "all_previous_photo_attempts_sha256":OLD_OBS_SHA,
        "existing_strict9_actual_original_roi_SHA256":OLD_PHOTO_SHA,
        "source_register_SHA256":sha(registry),
        "original_sampling_provenance_SHA256":sha(
            out/"selected_original_photo_design_provenance.csv"),
    }
    (out/"manifest.json").write_text(json.dumps(req,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL93M_BIAS_REPAIR_ORIGINAL_2048_PHOTO_SELECTION",json.dumps({
        "new_photo_ids":len(output),"source_arm_counts":dict(counts),
        "fresh_source_observations":metrics["eligible_disjoint_original_photo_observations"],
        "new_supported_source_photos":audit["new_supported_original_photo_IDs"],
        "targeted_source_strata_nonempty":sum(
            bool(x["source_remaining_photo_bearing_taxa"])
            for x in metrics["targeted_strata_actual_source_and_old9_roi9_support"]),
        "highest_new_photo_allocations":sorted(
            [{"zone":a["region"],"season":a["recorded_local_season"],
              "quota":a["targeted_sample_photo_quota"],
              "prior_strict":a["prior_strict_roi9_measured"],
              "source_taxa":a["source_remaining_photo_bearing_taxa"]}
             for a in metrics["targeted_strata_actual_source_and_old9_roi9_support"]],
             key=lambda d:d["quota"],reverse=True)[:12],
        "no_new_images_downloaded":True
    },sort_keys=True),flush=True)
    return req

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("groups","registry","atlas","atlas-receipt",
                 "union","union-receipt","taxa","output"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--reference",type=int,default=REFERENCE)
    p.add_argument("--targeted",type=int,default=AUGMENTATION)
    x=p.parse_args()
    plan(x.groups,x.registry,x.atlas,x.atlas_receipt,
         x.union,x.union_receipt,x.taxa,x.output,
         reference=x.reference,targeted=x.targeted)

if __name__=="__main__":
    main()
