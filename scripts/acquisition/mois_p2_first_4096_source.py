#!/usr/bin/env python3
"""MOIS first practical taxon-gap P2-only true source-photo refill (up to 4096).

READ ONLY the 16 SHA-verified tier07 archives (~1.19m original photo IDs)
of the already completed 93,209,351-original-photo registry. Never enumerate
93m photo IDs again. Before colour download, restrict the SQLite candidates
to P2 taxa with >=10 source photo opportunities and *zero* valid ROI9 in the
updated 1,045,978-photo source. Exactly one independent observation and one
source photo per selected taxon. Balance across observed local seasons and
latitude bands *within genuinely available P2 candidates*.

Data-only source selection: no photo pixels, no botanical petal/native truth,
no worldwide flora inclusion probability, no re-inference of earlier images.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv,gzip,hashlib,io,json,math,sqlite3,tarfile
from pathlib import Path

import pandas as pd

from scripts.acquisition.mois_balance1043930_tier07_2048_source import (
    PARENT_N, TIER, SOURCE_P, ZONE, SEASONS,
    OLD_OBS_SHA, UPSERT, FIELDS, create_db, zone, season, sha, require,
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS, RESULT_SCHEMA, obs_draw
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS

SCHEMA="mois_tier07_p2_zeroROI9_taxon_once_geoseason_original_photo_v1"
ARM="mois_tier07_P2_zero_strict9_one_taxon_one_original_photo_v1"
SEED="mois-p2-first-20261026-t07-4096-strict9-gap-v1"
OLD_PRIORITY_SHA="2dcc365c3b771a83477f1885703a2f573568bb3a7c996dcc11a194182068e8cc"
NEW_ROI2048_SHA="dea4699566430aa8fe33802e20d997646357ef46fd3429454959ba54fbd535a9"
PILOT64_SHA="c8b9be91bb27edbe66b52725e34fba0eab8dd5ec5617bca1b99bdaf893ab36ed"
P2_NAME="P2_source_has_photos_no_colour_10plus"
TARGET=4096
BATCH_SIZE=128
HUES=("white","yellow_orange","red_pink","blue_purple")
REGION_SUCCESS={
    "N_arctic":.361,"N_subpolar":.410,"N_temperate":.418,
    "N_tropical":.252,"S_temperate":.321,"S_tropical":.159,
    "S_subpolar":.138,"S_antarctic":.10,
}

def stable_rank(kind,*keys):
    return hashlib.sha256(
        ("|".join((SEED,kind,*map(str,keys)))).encode()
    ).hexdigest()

def effort_bin(n):
    n=int(n)
    return "E0_1" if n<=1 else ("E2_4" if n<=4 else "E5plus")

def frozen_prior(priority_path, new_path, pilot_path, union_path):
    priority_path,new_path,pilot_path,union_path=map(
        Path,(priority_path,new_path,pilot_path,union_path))
    require(sha(priority_path)==OLD_PRIORITY_SHA
            and sha(new_path)==NEW_ROI2048_SHA
            and sha(pilot_path)==PILOT64_SHA
            and sha(union_path)==OLD_OBS_SHA,
            "Original source taxon, old photo union, 2048 or 64-photo proof SHA changed")
    f=pd.read_csv(priority_path,compression="gzip",
                  dtype={"inat_taxon_id":str},low_memory=False)
    d=pd.read_parquet(new_path,columns=[
        "photo_id","observation_id","inat_taxon_id","continuous9_strict"])
    old=pd.read_parquet(union_path,columns=["photo_id","observation_id"])
    p=pd.read_csv(pilot_path,dtype={
        "photo_id":str,"observation_id":str,"taxon_id":str},low_memory=False)
    require(len(f)==42111 and f.inat_taxon_id.is_unique
            and int(f.n_strict_roi9.sum())>=390541
            and len(old)==1043930 and old.photo_id.astype(str).is_unique
            and len(d)==2048 and d.photo_id.astype(str).is_unique
            and len(p)==64 and p.photo_id.is_unique
            and d.observation_id.astype(str).is_unique
            and p.observation_id.is_unique
            and int(d.continuous9_strict.sum())==682,
            "Actual source input census or ROI9 results do not match frozen receipt")
    old_photo=set(old.photo_id.astype(str))
    old_obs=set(old.observation_id.fillna("").astype(str))-{""}
    new_photo=set(d.photo_id.astype(str)); new_obs=set(d.observation_id.astype(str))
    pilot_photo=set(p.photo_id.astype(str)); pilot_obs=set(p.observation_id.astype(str))
    require(not (old_photo & new_photo)
            and not (old_photo & pilot_photo)
            and not (new_photo & pilot_photo)
            and not (new_obs & pilot_obs),
            "An already attempted original photo/observation is duplicated across source arms")
    prior_zero=f.loc[f.priority.eq(P2_NAME)].copy()
    newly_coloured=set(d.loc[d.continuous9_strict,"inat_taxon_id"].astype(str))
    eligible=prior_zero.loc[~prior_zero.inat_taxon_id.isin(newly_coloured)]
    require(len(prior_zero)==17812 and len(eligible)==17727
            and int(f.n_strict_roi9.sum())<=400710
            and int(eligible.n_strict_roi9.sum())==0,
            "Real prior 2048 acquired nine-palette colours did not reduce P2 by 85 taxa")
    previous_attempts=Counter(d.inat_taxon_id.astype(str))
    previous_attempts.update(p.taxon_id.astype(str))
    info={}
    for row in eligible.itertuples(index=False):
        n=str(row.inat_taxon_id)
        info[n]={
            "species":str(row.species),
            "prior_image_attempts":int(row.n_image_photo_attempts)+previous_attempts[n],
            "frame_photo_opportunities":int(row.n_unique_source_photo_ids),
        }
    require(all(v["frame_photo_opportunities"]>=10 for v in info.values()),
            "Non-P2 or source-empty taxon passed the P2-only preselection")
    return info,old_photo|new_photo|pilot_photo,old_obs|new_obs|pilot_obs,{
        "prior_frame_taxa":len(f),
        "p2_taxa_before_new2048":len(prior_zero),
        "p2_taxa_remaining":len(info),
        "prior_photo_attempts":len(old_photo),
        "measured_new2048_ids":len(new_photo),
        "separate_64_lab_ids":len(pilot_photo),
        "prior_source_success_taxa_excluded":len(prior_zero)-len(info),
        "all_unique_prior_photo_IDs":len(old_photo|new_photo|pilot_photo),
        "old_observation_IDs_present":len(old_obs),
        "all_unique_prior_observations":len(old_obs|new_obs|pilot_obs),
    }

def read_filtered_tier07(db,root,master,allowed,used_photo,used_obs):
    """Visit signed tier07 source once; store P2 observations ONLY."""
    base=Path(root)
    audit=Counter()
    proofs=[]
    for group in range(16):
        tar_path=base/f"tier-07-group-{group:02d}.tar"
        receipt_path=base/f"group-{group:02d}-receipt.json"
        doc=json.loads(receipt_path.read_text())
        archive_sha=sha(tar_path)
        require(doc.get("schema")==RESULT_SCHEMA
                and doc.get("source_snapshot_SHA256")==master["full_source_manifest_sha256"]
                and archive_sha==master["published_tier_archive_SHA256"][tar_path.name]
                ==doc["archive_sha256_by_tier"][tar_path.name],
                "Tier07 signature changed: "+str(group))
        parts={int(a["part"]):a for a in doc["parts"]}
        require(len(parts)==8 and all(i//8==group for i in parts),
                "Original signed group does not cover exactly 8 source parts")
        group_rows=0
        with tarfile.open(tar_path) as arc:
            members={m.name.removeprefix("./"):m for m in arc.getmembers() if m.isfile()}
            expected={f"tier-07/part-{part:03d}.csv.gz" for part in parts}
            require(set(members)==expected,"Unexpected source tier07 archive contents")
            for part,proof in sorted(parts.items()):
                raw=arc.extractfile(members[f"tier-07/part-{part:03d}.csv.gz"]).read()
                require(hashlib.sha256(raw).hexdigest()==
                        proof["ranked_target_csv_sha256"][str(TIER)],
                        "Signed source gz data SHA changed")
                seen=0
                with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
                    with io.TextIOWrapper(gz,encoding="utf-8",newline="") as fh:
                        reader=csv.DictReader(fh)
                        require(tuple(reader.fieldnames or [])==COLUMNS,
                                "Original full93m photo metadata schema changed")
                        for row in reader:
                            seen+=1
                            pid,uuid,taxon=row["photo_id"],row["observation_uuid"],row["inat_taxon_id"]
                            require(pid.isdecimal() and taxon.isdecimal() and uuid
                                    and row["tier"]==str(TIER)
                                    and float(row["observation_inclusion_score"])==obs_draw(uuid),
                                    "Signed observation/photo identity or tier hash changed")
                            audit["photo_IDs_verified"]+=1
                            if taxon not in allowed:
                                audit["not_target_P2_taxon"]+=1;continue
                            audit["P2_source_photo_opportunities_in_tier07"]+=1
                            if pid in used_photo or uuid in used_obs:
                                audit["P2_prior_attempt_or_observation_excluded"]+=1
                                continue
                            ext=row["image_extension"].lower().strip().lstrip(".")
                            if row["frozen_loader_supported"]!="true" or ext not in ("jpeg","jpg","png"):
                                audit["P2_unsupported_source_format"]+=1
                                continue
                            lat,lon=row["latitude"],row["longitude"]
                            z,t=zone(lat),season(lat,row["observed_on"])
                            if z=="unknown" or t=="unknown":
                                audit["P2_missing_coordinate_or_observed_day"]+=1
                                continue
                            values=(
                                uuid,taxon,row["observed_on"],lat,lon,
                                row["positional_accuracy"],row["observer_id"],
                                row["photo_license"],pid,ext,
                                stable_rank("photo",uuid,pid),z,t,
                                stable_rank("reference",uuid),
                                stable_rank("anchor_observation",uuid),part,
                            )
                            db.execute(UPSERT,values)
                            audit["P2_target_original_photo_rows_indexed"]+=1
                require(seen==proof["counts_by_disjoint_tier"][str(TIER)],
                        "Missing original photo IDs in source tier07 part")
                group_rows+=seen
        require(group_rows==doc["counts_by_disjoint_tier"][str(TIER)],
                "Original full signed source group photo count changed")
        db.commit()
        proofs.append({"group":group,"original_tier07_photo_IDs":group_rows,
                       "source_tar_SHA256":archive_sha})
    require(audit["photo_IDs_verified"]==
            master["new_target_counts_by_tier"][str(TIER)],
            "Full original tier07 16-group source not completely scanned")
    return audit,proofs

def gap_quota(counts,target):
    """Weighted balanced quotas on SUPPORTED source cells, not geographic 50:50."""
    counts={k:int(v) for k,v in counts.items() if int(v)>0}
    require(counts and target>0 and sum(counts.values())>=target,
            "Not enough supported independent P2 source taxa for stated quota")
    # Soft maximum of 18% per cell; relax only when source availability requires.
    cap=max(1,math.ceil(.18*target))
    while sum(min(cap,v) for v in counts.values())<target:
        cap=min(target,cap*2)
    quota={key:min(2,n) for key,n in counts.items()}
    if sum(quota.values())>target:
        quota={key:0 for key in counts}
    weight={key:(math.sqrt(n)/math.sqrt(REGION_SUCCESS.get(key[0],.3)))
            for key,n in counts.items()}
    left=target-sum(quota.values())
    for _ in range(left):
        choices=[key for key,n in counts.items() if quota[key]<min(cap,n)]
        require(choices,"Available P2 source cells could not fill targeted quota")
        key=max(choices,key=lambda x:(weight[x]/(1+quota[x]),
                -ZONE.index(x[0]),-SEASONS.index(x[1])))
        quota[key]+=1
    require(sum(quota.values())==target and
            all(0<=quota[k]<=counts[k] for k in quota),
            "P2 allocation did not preserve finite cell availability")
    return quota

def select_one_per_taxon(db,allowed,target=TARGET):
    """Random hashed observation per taxon, then cell/effort stratified ranks."""
    n_db=int(db.execute("SELECT COUNT(*) FROM obs").fetchone()[0])
    require(n_db>0,"No P2 original observations in signed source")
    db.execute("CREATE INDEX obs_by_taxon_rank ON obs(taxon,observation_rank)")
    query="""SELECT * FROM (
        SELECT """+",".join(FIELDS)+""",
        COUNT(*) OVER (PARTITION BY taxon) AS source_taxon_observations,
        ROW_NUMBER() OVER (PARTITION BY taxon ORDER BY observation_rank,uuid) AS rank_in_taxon
        FROM obs
    ) WHERE rank_in_taxon=1"""
    pool=defaultdict(list)
    for row in db.execute(query):
        obs=dict(zip(FIELDS,row[:len(FIELDS)]))
        obs["source_taxon_observations"]=int(row[len(FIELDS)])
        tax=obs["taxon"]
        require(tax in allowed and obs["photo_n"]>0,
                "P2 anchor observation status or multiplicity invalid")
        obs.update(allowed[tax])
        obs["effort_bin"]=effort_bin(obs["prior_image_attempts"])
        pool[(obs["zone"],obs["season"])].append(obs)
    total_available=sum(map(len,pool.values()))
    require(total_available>=min(target,512),
            "Fewer than 512 independent never-measured P2 taxa are available")
    selected_n=min(target,total_available)
    counts={k:len(items) for k,items in pool.items()}
    quota=gap_quota(counts,selected_n)
    selected=[];summary=[];selcounts=Counter()
    for (z,s),arr in sorted(pool.items()):
        k=quota[(z,s)]
        # Each effort tier has its own seeded uniform taxon rank; P2 taxa with
        # 0/1 previous failures are exhausted before repeatedly failing taxa.
        arr.sort(key=lambda o:(
            0 if o["effort_bin"]=="E0_1" else
            (1 if o["effort_bin"]=="E2_4" else 2),
            stable_rank("taxon_within_effort_cell",z,s,o["taxon"])))
        chosen=arr[:k]
        eligible_effort=Counter(o["effort_bin"] for o in arr)
        selected_effort=Counter(o["effort_bin"] for o in chosen)
        for item in chosen:
            m=int(item["photo_n"]);N=int(item["source_taxon_observations"])
            effort=item["effort_bin"]
            q=selected_effort[effort]/eligible_effort[effort]
            # These are probabilities *conditional on prior adaptive source
            # strata*, never global iNaturalist or botanical flora IPW.
            conditional=SOURCE_P*q/N/m
            require(0<conditional<=1,"P2 nested photo selection chance invalid")
            item["conditional_design_probability"]=conditional
            item["eligible_taxa_in_effort_cell"]=eligible_effort[effort]
            item["selected_taxa_in_effort_cell"]=selected_effort[effort]
            selected.append(item)
            selcounts[effort]+=1
        summary.append({
            "latitude_region":z,"recorded_local_season":s,
            "supported_eligible_P2_taxa":len(arr),"selected_photo_quota":k,
            "selected_P2_effort_0to1":selected_effort["E0_1"],
            "selected_P2_effort_2to4":selected_effort["E2_4"],
            "selected_P2_effort_5plus":selected_effort["E5plus"],
        })
    require(len(selected)==selected_n and
            len(set(v["taxon"] for v in selected))==selected_n and
            len(set(v["uuid"] for v in selected))==selected_n and
            len(set(v["photo_id"] for v in selected))==selected_n,
            "At most one original photo, observation and accepted source taxon per P2 slot")
    require(not any(v["prior_image_attempts"]>=5 for v in selected)
            or selcounts["E0_1"]+selcounts["E2_4"]<selected_n,
            "Repeated ROI-failure taxa received photos before scarce easy P2 strata")
    selected.sort(key=lambda o:stable_rank("worker_shuffle",o["uuid"]))
    return selected,summary,dict(selcounts),total_available,n_db

def plan(source,registry,priority,new2048,pilot64,old_union,out,
         *,target=TARGET,batch=BATCH_SIZE):
    source,registry,priority,new2048,pilot64,old_union,out=map(Path,(
        source,registry,priority,new2048,pilot64,old_union,out))
    require(not out.exists() and target>=512 and batch>0 and batch<=256,
            "Refuse overwrite or invalid original ROI capacity")
    master=json.loads(registry.read_text(encoding="utf-8"))
    require(master.get("schema")=="mois_complete_full93m_original_photo_target_registry_v1"
            and master.get("exact_original_photo_IDs_in_all_disjoint_tiers")==PARENT_N
            and master.get("original_named_taxa_from_full_source")==42111
            and len(master.get("published_tier_archive_SHA256",{}))>=16,
            "Original full 93m photo source metadata registry absent")
    allowed,used_photos,used_obs,prior_audit=frozen_prior(
        priority,new2048,pilot64,old_union)
    scratch=out.parent/(out.name+"_p2_source_tmp.sqlite")
    require(not scratch.exists(),"Refuse unsafe reuse of mutable SQLite input")
    try:
        with sqlite3.connect(scratch) as db:
            create_db(db)
            audit,proof=read_filtered_tier07(
                db,source,master,allowed,used_photos,used_obs)
            selected,cells,effort,possible,source_obs=select_one_per_taxon(
                db,allowed,target)
    finally:
        scratch.unlink(missing_ok=True)
    require(all(v["photo_id"] not in used_photos and
                v["uuid"] not in used_obs for v in selected),
            "P2 future queue photo/observation overlaps a measured photo")
    out.mkdir(parents=True)
    (out/"shards").mkdir()
    workers=[];details=[]
    for obs in selected:
        pid=obs["photo_id"]
        mid="MOIS-"+stable_rank("measurement_ID",pid).upper()[:24]
        lat,lon=obs["latitude"],obs["longitude"]
        grid=f"{int(math.floor((float(lat)+90)/5))}:{int(math.floor((float(lon)+180)/5))}"
        workers.append({
            "measurement_id":mid,"photo_id":pid,
            "inat_taxon_id":obs["taxon"],"species":obs["species"],
            "observation_id":obs["uuid"],"observer_id":obs["observer"],
            "observed_on":obs["observed_on"],"latitude":lat,"longitude":lon,
            "positional_accuracy":obs["accuracy"],"photo_license":obs["license"],
            "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{obs['extension']}",
            "image_filename":mid+".jpg","source_arm":ARM,
            "selection_probability":format(obs["conditional_design_probability"],".17g"),
            "source_snapshot_date":"2026-09-27","shard":"",
        })
        details.append({
            "photo_id":pid,"original_observation_id":obs["uuid"],
            "inat_taxon_id":obs["taxon"],"source_gap_priority":"P2",
            "prior_strict9":0,"prior_image_attempts":obs["prior_image_attempts"],
            "original_source_frame_photo_opportunities":obs["frame_photo_opportunities"],
            "latitude_region":obs["zone"],"recorded_local_season":obs["season"],
            "source_5degree_grid_cell":grid,"prior_attempt_effort_bin":obs["effort_bin"],
            "source_taxon_eligible_observations":obs["source_taxon_observations"],
            "source_eligible_supported_photos_in_observation":obs["photo_n"],
            "source_taxa_in_effort_zone_season":obs["eligible_taxa_in_effort_cell"],
            "selected_taxa_in_effort_zone_season":obs["selected_taxa_in_effort_cell"],
            "conditional_design_probability_NOT_global_botanical_inclusion":obs["conditional_design_probability"],
            "source_previous_all_measured_photo_and_observation_ID_overlap":False,
        })
    n_units=math.ceil(len(workers)/batch)
    shard_sha={}
    for unit,start in enumerate(range(0,len(workers),batch)):
        path=out/"shards"/f"shard-{unit:03d}.csv"
        with path.open("w",newline="",encoding="utf-8") as fh:
            w=csv.DictWriter(fh,fieldnames=WORKER_COLUMNS,lineterminator="\n")
            w.writeheader()
            for row in workers[start:start+batch]:
                w.writerow({**row,"shard":str(unit)})
        shard_sha[str(unit)]=sha(path)
    prov=out/"selected_P2_photo_source_provenance.csv"
    with prov.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(details[0]),lineterminator="\n")
        w.writeheader();w.writerows(details)
    original={
        "schema":"mois_full_source_observation_sampling_cohort_v1",
        "sampling_schema":SCHEMA,
        "colour_measured":False,
        "new_photo_attempts_planned":len(workers),
        "full_frame_photo_ids":PARENT_N,
        "original_upstream_photo_ID_count":PARENT_N,
        "source_tier":TIER,
        "source_tier_width":SOURCE_P,
        "source_arm_counts":{ARM:len(workers)},
        "photo_shard_sha256":shard_sha,
        "original_model_batch_size":batch,
        "worker_ids":list(range(n_units)),"work_units":n_units,
        "units_per_group":4,"group_ids":list(range(math.ceil(n_units/4))),
        "P2_only_no_existing_strict9_source_taxa":True,
        "distinct_original_plant_taxa":len(workers),
        "distinct_original_observations":len(workers),
        "earlier_reference_cohort_256_kept_separately":True,
        "source_P5_high_colour_taxa_targeted":0,
        "remaining_original_P2_source_taxa_before_selection":len(allowed),
        "candidate_distinct_P2_source_taxa":possible,
        "candidate_P2_source_observations_indexed":source_obs,
        "source_photo_registry_shas":proof,
        "full_tier07_source_audit":dict(audit),
        "prior_image_and_taxon_audit":prior_audit,
        "frozen_source_region_local_season_quota":cells,
        "selection_by_previous_photo_effort":effort,
        "photo_source_provenance_SHA256":sha(prov),
        "source_registry_SHA256":sha(registry),
        "old_1043930_attempted_photo_union_SHA256":sha(old_union),
        "old_42111_priority_after_400710_SHA256":sha(priority),
        "new_2048_actual_original_photo_ROI_SHA256":sha(new2048),
        "prior_64_separate_boxLab_pilot_SHA256":sha(pilot64),
        "all_source_photo_ID_and_observation_IDs_new_and_unique":True,
        "one_photo_per_source_taxon_for_priority_P2":True,
        "ROI9_from_any_P2_photo_before_new_attempt":0,
        "selection_probabilities_are_conditional_NOT_93m_population_IPW":True,
        "native_wild_and_petal_anthesis_unverified":True,
        "no_original_photo_image_download_or_inference":True,
        "all_original_93m_photo_opportunities_retained_unmodified":True,
    }
    (out/"manifest.json").write_text(
        json.dumps(original,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_REAL93M_P2_ONLY_TAXON_BREADTH_EFFICIENT_SOURCE_QUEUE",
          json.dumps({
            "selected_photo_and_unique_taxon_IDs":len(workers),
            "unique_P2_taxa_in_tier07":possible,
            "prior_P2_taxa_without_strict9":len(allowed),
            "upstream_source_photo_IDs_verified":audit["photo_IDs_verified"],
            "only_P2_photo_IDs_indexed":audit["P2_target_original_photo_rows_indexed"],
            "prior_photo_ids_deduped":len(used_photos),
            "effort":effort,"work_units":n_units,
            "group_ids":original["group_ids"],
            "new_images_processed":0
          },sort_keys=True),flush=True)
    return original

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("source","registry","priority","new2048","pilot64",
                 "old-union","output"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--target",type=int,default=TARGET)
    p.add_argument("--batch",type=int,default=BATCH_SIZE)
    x=p.parse_args()
    plan(x.source,x.registry,x.priority,x.new2048,x.pilot64,
         x.old_union,x.output,target=x.target,batch=x.batch)

if __name__=="__main__":
    main()
