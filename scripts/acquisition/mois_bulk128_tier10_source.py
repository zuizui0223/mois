#!/usr/bin/env python3
"""MOIS one-shot 128-shard / up-to-131072 original-photo bulk source.

Exactly one SHA-verified 93m observation-hash tier (tier10, ~9.5m original
source photo records, disjoint from earlier tier07 queues) is streamed ONCE.
Original source has 128 immutable taxon-modulo parts packed in 16 archives.
Each part yields up to 1024 independent new original photo/observation IDs:
256 source-frame random reference + up to 768 taxon/region/local-season gap
targets (P1/P2 first; P3/P4 next; P5 *never* targeted).
One model initialization per 1024-photo shard; four <=256-photo FCP packets.
No images are fetched by this selection script. Do NOT treat this deliberately
outcome-dependent design as an unbiased world-flora flower-colour sample.
"""
from __future__ import annotations

import argparse,bisect,collections,csv,gzip,hashlib,io,json,math,tarfile
from pathlib import Path
import pandas as pd

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, RESULT_SCHEMA, TIER_ENDS, obs_draw, require as verify,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS
from scripts.acquisition.mois_balance1043930_tier07_2048_source import zone,season

SCHEMA="mois_original93m_tier10_bulk128_probability_reference_plus_gap_v1"
TIER=10
SHARDS=128
SOURCE_GROUPS=16
SHARDS_PER_GROUP=8
TARGET_PER_SHARD=1024
REFERENCE_PER_SHARD=256
WORK_UNITS_PER_SHARD=4
SEED="mois-bulk128-131072-geoseason-taxongap-tier10-v1"
PARENT_PHOTOS=93_209_351
PRIOR_PHOTOS=1_043_930
NO_COLOUR_PRIORITY=("P1_source_rare_no_colour_1to9",
                    "P2_source_has_photos_no_colour_10plus")
THIN_PRIORITY=("P3_thin_colour_1to4","P4_thin_colour_5to9")
P5="P5_existing_10plus_strict9"
REFERENCE_ARM="mois_bulk128_tier10_photo_reference_v1"
GAP_ARM="mois_bulk128_tier10_taxon_region_season_gap_v1"
SOURCE_WIDTH=TIER_ENDS[TIER]-TIER_ENDS[TIER-1]
HUES=("white","yellow_orange","red_pink","blue_purple")
PHOTO_PRIORITY_SHA="2dcc365c3b771a83477f1885703a2f573568bb3a7c996dcc11a194182068e8cc"
OLD_PHOTO_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
NEW2048_SHA="dea4699566430aa8fe33802e20d997646357ef46fd3429454959ba54fbd535a9"
PILOT64_SHA="c8b9be91bb27edbe66b52725e34fba0eab8dd5ec5617bca1b99bdaf893ab36ed"
P2_RESERVE_SHA="1f7b146a98d31b3f465f7b8deec6ac42cc644f877372885833461aedf8f1b806"
NEWCOLOUR_SHA="c47b5ecfd66251ecf3051833d03a4137bf9c8087857a1722970db87f4c943195"
ZONE_WEIGHT={"N_arctic":1.0,"N_subpolar":1.0,"N_temperate":0.7,
             "N_tropical":2.0,"S_subpolar":1.3,
             "S_temperate":1.6,"S_tropical":2.7}
SEASON_WEIGHT={"winter":1.35,"autumn":1.2,"spring":1.0,"summer":.9}
CAP={NO_COLOUR_PRIORITY[0]:4,NO_COLOUR_PRIORITY[1]:8,
     THIN_PRIORITY[0]:8,THIN_PRIORITY[1]:5}

def require(cond,msg):
    if not cond:
        raise ValueError("MOIS_BULK128_SOURCE: "+msg)

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):
            digest.update(chunk)
    return digest.hexdigest()

def rank(phase,*components):
    raw="|".join((SEED,phase,*map(str,components)))
    return hashlib.sha256(raw.encode()).hexdigest()

def load_prior(priority,old_union,new2048,pilot64,p2_reserved,recovered):
    """Exact photo+observation ID reuse gates; no guessed native/flower state."""
    paths=[Path(x) for x in (priority,old_union,new2048,pilot64,p2_reserved,recovered)]
    checks=[PHOTO_PRIORITY_SHA,OLD_PHOTO_SHA,NEW2048_SHA,
            PILOT64_SHA,P2_RESERVE_SHA,NEWCOLOUR_SHA]
    for path,digest in zip(paths,checks):
        require(sha(path)==digest,"Frozen baseline changed: "+path.name)
    base=pd.read_csv(paths[0],compression="gzip",dtype={"inat_taxon_id":str},low_memory=False)
    old=pd.read_parquet(paths[1],columns=["photo_id","observation_id"])
    new=pd.read_parquet(paths[2],columns=["photo_id","observation_id","inat_taxon_id","continuous9_strict"])
    sixtyfour=pd.read_csv(paths[3],dtype={"photo_id":str,"observation_id":str,"taxon_id":str})
    reserved=pd.read_csv(paths[4],dtype={"photo_id":str,"original_observation_id":str,"inat_taxon_id":str})
    covered=pd.read_csv(paths[5],dtype={"inat_taxon_id":str})
    require(len(base)==42111 and base.inat_taxon_id.is_unique
            and len(old)==PRIOR_PHOTOS and old.photo_id.astype(str).is_unique
            and len(new)==2048 and new.photo_id.astype(str).is_unique
            and int(new.continuous9_strict.sum())==682
            and len(sixtyfour)==64 and sixtyfour.photo_id.is_unique
            and len(reserved)==4096 and reserved.photo_id.is_unique
            and len(covered)==85 and covered.inat_taxon_id.is_unique,
            "Before-colour original source frame and 2048+4096 reservations changed")
    ids=[
        set(old.photo_id.astype(str)),
        set(new.photo_id.astype(str)),
        set(sixtyfour.photo_id.astype(str)),
        set(reserved.photo_id.astype(str)),
    ]
    require(len(set.union(*ids))==sum(map(len,ids)),
            "Overlap in independent frozen source-photo ID reservations")
    obs=[
        set(old.observation_id.fillna("").astype(str))-{""},
        set(new.observation_id.astype(str)),
        set(sixtyfour.observation_id.astype(str)),
        set(reserved.original_observation_id.astype(str)),
    ]
    # Legacy original observed IDs can be missing/aliased; photo-ID exact
    # provenance is stronger and is enforced independently above.
    require(not set.intersection(obs[1],obs[2])
            and not set.intersection(obs[1],obs[3])
            and not set.intersection(obs[2],obs[3]),
            "New source observation UUID reservations conflict")
    known_taxa=set(covered.inat_taxon_id.astype(str))
    require(set(new.loc[new.continuous9_strict & new.inat_taxon_id.astype(str).isin(base.inat_taxon_id)
                        & new.inat_taxon_id.astype(str).isin(
                            base.loc[base.priority.isin(NO_COLOUR_PRIORITY),"inat_taxon_id"]),
                        "inat_taxon_id"].astype(str)).issubset(known_taxa),
            "True newly coloured original taxa absent from delta receipt")
    priorities=base.set_index("inat_taxon_id")["priority"].to_dict()
    for taxon in known_taxa:
        if priorities.get(taxon) in NO_COLOUR_PRIORITY:
            priorities[taxon]=THIN_PRIORITY[0]
    reserved_no_colour=set(reserved.inat_taxon_id)
    for taxon in reserved_no_colour:
        if priorities.get(taxon) in NO_COLOUR_PRIORITY:
            priorities[taxon]="RESERVED_PREVIOUS_P2_4096_NOT_FOR_TARGETED"
    counts={
        "original_prior_photo_IDs":len(ids[0]),
        "original_2048_new_photo_IDs":len(ids[1]),
        "original_64_separate_Lab_photo_IDs":len(ids[2]),
        "original_4096_P2_reserved_photo_IDs":len(ids[3]),
        "exact_photo_IDs_reserved_or_measured":len(set.union(*ids)),
        "original_nominal_source_taxa":len(base),
        "prior_P2_photo_IDs_reserved_by_different_running_job":len(reserved_no_colour),
        "originally_zero9_taxa_now_known":len(known_taxa),
    }
    return priorities,set.union(*ids),set.union(*obs),counts

def add_original_photo(table,row,part):
    pid,uuid,taxon=row["photo_id"],row["observation_uuid"],row["inat_taxon_id"]
    require(pid.isdecimal() and uuid and taxon.isdecimal() and int(taxon)%SHARDS==part,
            "Unexpected original taxon modulo128 partition or missing photo/observation ID")
    require(row["tier"]==str(TIER)
            and float(row["observation_inclusion_score"])==obs_draw(uuid),
            "Signed prior observation hash or disjoint tier09+ source changed")
    z=zone(row["latitude"])
    ss=season(row["latitude"],row["observed_on"])
    photorank=rank("random_photo_per_observation",uuid,pid)
    if uuid not in table:
        table[uuid]={
            "uuid":uuid,"taxon":taxon,"date":row["observed_on"],
            "latitude":row["latitude"],"longitude":row["longitude"],
            "observer":row["observer_id"],"accuracy":row["positional_accuracy"],
            "zone":z,"season":ss,"photo_id":pid,
            "extension":row["image_extension"],"license":row["photo_license"],
            "photorank":photorank,"nphotos":1,
            "ref_rank":rank("uniform_reference_observation",uuid),
            "target_rank":rank("independent_target_observation",uuid),
            "part":part
        }
    else:
        o=table[uuid]
        require((taxon,row["observed_on"],row["latitude"],
                 row["longitude"],row["observer_id"],row["positional_accuracy"])
                ==(o["taxon"],o["date"],o["latitude"],o["longitude"],
                   o["observer"],o["accuracy"]),
                "Same original observation UUID carries conflicting botany/geography/date")
        o["nphotos"]+=1
        if photorank<o["photorank"]:
            o.update(photo_id=pid,extension=row["image_extension"],
                     license=row["photo_license"],photorank=photorank)

def take_priority(candidates,already,limit,served,weight,per_taxon_cap):
    """Prefer 1st observation of a taxon before 2nd, across supported cells."""
    if limit<=0:return []
    by_taxon=collections.defaultdict(list)
    for o in candidates:
        if o["uuid"] not in already:
            by_taxon[o["taxon"]].append(o)
    bins=collections.defaultdict(list)
    for taxon,items in by_taxon.items():
        items.sort(key=lambda x:(x["target_rank"],x["uuid"]))
        for depth,o in enumerate(items[:per_taxon_cap]):
            cell=(o["zone"],o["season"])
            bins[(depth,cell)].append(o)
    for (depth,cell),items in bins.items():
        items.sort(key=lambda o:rank("taxon_cell",o["taxon"],o["uuid"]))
    selected=[]
    for depth in range(per_taxon_cap):
        active={cell:iter(items) for (d,cell),items in bins.items() if d==depth}
        nexts={}
        for cell,it in active.items():
            try:nexts[cell]=next(it)
            except StopIteration:pass
        while nexts and len(selected)<limit:
            # Cells with genuinely scarce measured colour and valid candidates
            # receive more turns; empty cells never receive fabricated photos.
            cell=min(nexts,key=lambda c:(
                (served[c]+1)/(weight.get(c,1.0)),
                rank("stable_cell_tiebreak",*c)))
            o=nexts.pop(cell)
            if o["uuid"] not in already:
                selected.append(o);already.add(o["uuid"]);served[cell]+=1
            try:nexts[cell]=next(active[cell])
            except StopIteration:pass
        if len(selected)>=limit:break
    return selected

def choose_part(observations,priority,part,*,quota=TARGET_PER_SHARD,
                reference=REFERENCE_PER_SHARD):
    require(0<=part<SHARDS and 1<=reference<quota,
            "Invalid immutable 128-shard photo budget")
    nobs=len(observations)
    require(nobs>=reference,"Source taxon partition lacks reference observations")
    allobs=list(observations.values())
    references=sorted(allobs,key=lambda o:(o["ref_rank"],o["uuid"]))[:reference]
    chosen={o["uuid"] for o in references}
    served=collections.Counter()
    weight={(z,s):ZONE_WEIGHT.get(z,.9)*SEASON_WEIGHT.get(s,1)
            for z in ZONE_WEIGHT for s in SEASON_WEIGHT}
    targeted_budget=quota-reference
    bypriority=collections.defaultdict(list)
    for o in allobs:
        if o["uuid"] in chosen or o["zone"]=="unknown" or o["season"]=="unknown":
            continue
        status=priority.get(o["taxon"],"unresolved")
        if status in CAP:
            bypriority[status].append(o)
    targets=[]
    # Force source-gap breadth as the first stage, but never pretend absent
    # P1/P2 photo/observations are present. P5 never enters this pool.
    urgent_cap=math.ceil(targeted_budget*.65)
    for status in NO_COLOUR_PRIORITY:
        part_candidates=bypriority[status]
        new=take_priority(part_candidates,chosen,urgent_cap-len(targets),
                          served,weight,CAP[status])
        targets.extend((status,o) for o in new)
    for status in THIN_PRIORITY:
        space=targeted_budget-len(targets)
        new=take_priority(bypriority[status],chosen,space,
                          served,weight,CAP[status])
        targets.extend((status,o) for o in new)
    if len(targets)<targeted_budget:
        # Use supported P1/P2 observations beyond their initial 65% quota
        # before expanding to already well-measured taxa.
        for status in (*NO_COLOUR_PRIORITY,*THIN_PRIORITY):
            space=targeted_budget-len(targets)
            new=take_priority(bypriority[status],chosen,space,
                              served,weight,CAP[status])
            targets.extend((status,o) for o in new)
            if len(targets)==targeted_budget:break
    result=[("reference",priority.get(o["taxon"],"unknown"),o) for o in references]
    result +=[("targeted",kind,o) for kind,o in targets]
    require(len(result)<=quota
            and len({o["uuid"] for _,_,o in result})==len(result)
            and len({o["photo_id"] for _,_,o in result})==len(result)
            and all(label!=P5 for role,label,_ in result if role=="targeted"),
            "Target photo IDs or source priorities duplicated/polluted")
    return result,{
        "original_source_eligible_observations":nobs,
        "selected_reference_photo_attempts":len(references),
        "selected_targeted_photo_attempts":len(targets),
        "targeted_unfilled_from_unrepresented_source_support":
            targeted_budget-len(targets),
        "source_targeted_priority_counts":
            dict(collections.Counter(status for status,_ in targets)),
        "source_targeted_region_season_counts":{
            "|".join(k):v for k,v in sorted(served.items())},
    }

def photo_worker(o,role,part,unit,source_count):
    pid=o["photo_id"]
    image_id="MOIS-"+rank("measurement",pid).upper()[:24]
    # This number is for technical worker compatibility. For TARGETED rows
    # it reflects only phase-1 observation tier + within-observation photo
    # randomization; the adaptive phase-2 inclusion probability is UNKNOWN.
    # No Horvitz-Thompson global-flora weighting may use this field.
    sampled_prob=SOURCE_WIDTH/o["nphotos"]
    if role=="reference":
        sampled_prob*=min(REFERENCE_PER_SHARD,source_count)/source_count
    ext=o["extension"].lower().lstrip(".")
    return {
       "measurement_id":image_id,"photo_id":pid,
       "inat_taxon_id":o["taxon"],"species":o["species"],
       "observation_id":o["uuid"],"observer_id":o["observer"],
       "observed_on":o["date"],"latitude":o["latitude"],
       "longitude":o["longitude"],"positional_accuracy":o["accuracy"],
       "photo_license":o["license"],
       "photo_url_large":(
          f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}"),
       "image_filename":image_id+".jpg",
       "source_arm":REFERENCE_ARM if role=="reference" else GAP_ARM,
       "selection_probability":format(sampled_prob,".17g"),
       "source_snapshot_date":"2026-09-27","shard":str(unit),
    }

def select_and_write_part(part,source_gz,proof,priorities,used_photo,used_obs,out,
                          *,quota=TARGET_PER_SHARD,reference=REFERENCE_PER_SHARD):
    require(int(proof["part"])==part
            and sha(source_gz)==proof["ranked_target_csv_sha256"][str(TIER)],
            "Original source partition's immutable photo bytes changed")
    frame=collections.Counter()
    obs={}
    with gzip.open(source_gz,"rt",encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        require(tuple(reader.fieldnames or ())==COLUMNS,
                "Original source photo-column contract altered")
        for row in reader:
            frame["tier_source_photo_IDs"]+=1
            pid,uuid=row["photo_id"],row["observation_uuid"]
            # Even excluded photos are SHA-validated at gz/source-group level.
            require(pid.isdecimal() and uuid and row["tier"]==str(TIER),
                    "Missing source original-photo/observation identity")
            if pid in used_photo or uuid in used_obs:
                frame["prior_photo_or_observation_reserved"]+=1
                continue
            ext=row["image_extension"].lower().lstrip(".")
            if row["frozen_loader_supported"]!="true" or ext not in ("jpg","jpeg","png"):
                frame["unsupported_original_image"]+=1;continue
            add_original_photo(obs,row,part)
            frame["valid_new_source_photo_rows"]+=1
    require(frame["tier_source_photo_IDs"]==proof["counts_by_disjoint_tier"][str(TIER)],
            "Source part original photo-ID count lost or inflated")
    for o in obs.values():
        o["species"]=priorities.get(o["taxon"],{}).get("species","")
    selection,status=choose_part(obs,{k:v["priority"] for k,v in priorities.items()},
                                part,quota=quota,reference=reference)
    root=Path(out)
    work=root/"shards";work.mkdir(parents=True,exist_ok=True)
    provenance=root/"provenance";provenance.mkdir(parents=True,exist_ok=True)
    total=len(selection);unit_count=WORK_UNITS_PER_SHARD
    require(total>=unit_count,"No distinct original-photo source items in partition")
    order=sorted(selection,key=lambda x:rank("warm_shuffle",x[2]["uuid"]))
    chunks=[];start=0
    for i in range(unit_count):
        take=min(256,math.ceil((total-start)/(unit_count-i)))
        require(take>0,"Invalid original source shard worker packet")
        chunks.append(order[start:start+take]);start+=take
    require(start==total and len(chunks)==4
            and all(0<len(c)<=256 for c in chunks),
            "128-shard work units lost original photo IDs")
    work_checksums={}
    details=[]
    for i,chunk in enumerate(chunks):
        unit=part*4+i
        p=work/f"shard-{unit:03d}.csv"
        with p.open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=WORKER_COLUMNS,lineterminator="\n")
            writer.writeheader()
            for role,pri,o in chunk:
                writer.writerow(photo_worker(o,role,part,unit,len(obs)))
                details.append({
                    "photo_id":o["photo_id"],"observation_id":o["uuid"],
                    "inat_taxon_id":o["taxon"],"species":o["species"],
                    "acquisition_design_cohort":role,
                    "original_source_priority":pri,
                    "latitude_region":o["zone"],"recorded_local_season":o["season"],
                    "source_5degree_cell":(
                       f"{int(math.floor((float(o['latitude'])+90)/5))}:"
                       f"{int(math.floor((float(o['longitude'])+180)/5))}"
                       if o["zone"]!="unknown" and o["longitude"] else "unknown"),
                    "source_photo_group":part,"source_worker_unit":unit,
                    "photos_per_observation_in_new_source_part":o["nphotos"],
                    "frame_photo_probability_component":SOURCE_WIDTH,
                    "technical_worker_probability_not_targeted_full_joint":
                        photo_worker(o,role,part,unit,len(obs))["selection_probability"],
                    "adaptive_phase2_inclusion_probability_IDENTIFIED":
                        role=="reference",
                    "previous_photo_ID_and_observation_overlap":False,
                    "is_botanically_verified_open_native_flower":False,
                })
        work_checksums[str(unit)]=sha(p)
    detail_path=provenance/f"part-{part:03d}.csv"
    with detail_path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(details[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(details)
    receipt={
       "schema":SCHEMA,"part":part,"source_tier":TIER,
       "original_source_part_photo_count":frame["tier_source_photo_IDs"],
       "original_source_part_counts":dict(frame),
       "selected_new_original_photo_IDs":total,
       "selected_distinct_original_observations":total,
       "original_source_partition_zero_duplication":True,
       "source_photo_worker_unit_SHA256":work_checksums,
       "source_photo_provenance_SHA256":sha(detail_path),
       "selected_photo_attempts_by_design":
           dict(collections.Counter(role for role,_,_ in selection)),
       "selected_photo_attempts_by_original_priority":
           dict(collections.Counter(label for _,label,_ in selection)),
       "selected_photo_attempts_by_zone":
           dict(collections.Counter(o["zone"] for _,_,o in selection)),
       "selected_photo_attempts_by_recorded_local_season":
           dict(collections.Counter(o["season"] for _,_,o in selection)),
       "P5_targeted_original_photo_count":0,
       "target_photo_count_if_source_fully_supports":quota,
       "source_selection_details":status,
       "no_real_image_download_or_colour_classifier_in_this_process":True,
    }
    rp=root/f"part-{part:03d}.receipt.json"
    rp.write_text(json.dumps(receipt,sort_keys=True,indent=2,allow_nan=False)+"\n")
    return receipt

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("source","registry","priority","old-union","new2048",
                 "pilot64","p2-reserved","newcovered","output"):
        parser.add_argument("--"+name,required=True,type=Path)
    parser.add_argument("--target-per-shard",type=int,default=TARGET_PER_SHARD)
    parser.add_argument("--reference-per-shard",type=int,default=REFERENCE_PER_SHARD)
    args=parser.parse_args()
    require(not args.output.exists(),"Source plan may not silently overwrite previous selection")
    master=json.loads(args.registry.read_text())
    require(master.get("schema")=="mois_complete_full93m_original_photo_target_registry_v1"
            and master.get("exact_original_photo_IDs_in_all_disjoint_tiers")==PARENT_PHOTOS
            and master.get("original_named_taxa_from_full_source")==42111
            and sha(args.registry)=="4ed5139e510ba8442d4cdb03f2b769d94fac942832b1d544d7a7925bc431255b",
            "Not the original immutable 93m metadata source")
    pri,used_photo,used_obs,prior_counts=load_prior(
        args.priority,args.old_union,args.new2048,args.pilot64,
        args.p2_reserved,args.newcovered)
    catalogue=pd.read_csv(args.priority,compression="gzip",
                          dtype={"inat_taxon_id":str},low_memory=False)
    species=catalogue.set_index("inat_taxon_id")["species"].fillna("").astype(str).to_dict()
    priority={taxon:{"priority":label,"species":species.get(taxon,"")}
              for taxon,label in pri.items()}
    args.output.mkdir(parents=True)
    full=[];selected_all=0;total_source=0
    for group in range(SOURCE_GROUPS):
        archive=args.source/f"tier-{TIER:02d}-group-{group:02d}.tar"
        receipt=args.source/f"group-{group:02d}-receipt.json"
        group_doc=json.loads(receipt.read_text())
        require(group_doc.get("schema")==RESULT_SCHEMA
                and group_doc["source_snapshot_SHA256"]==master["full_source_manifest_sha256"]
                and sha(archive)==master["published_tier_archive_SHA256"][archive.name]
                ==group_doc["archive_sha256_by_tier"][archive.name],
                "Frozen source tar group failed authentication: "+str(group))
        proofs={int(p["part"]):p for p in group_doc["parts"]}
        require(set(proofs)==set(range(group*8,(group+1)*8)),
                "Missing one of eight source photo partitions")
        with tarfile.open(archive) as tar:
            entries={m.name.removeprefix("./"):m for m in tar.getmembers() if m.isfile()}
            expected={f"tier-{TIER:02d}/part-{p:03d}.csv.gz" for p in proofs}
            require(set(entries)==expected,
                    "Original tier part archive has missing or unexpected members")
            for part in sorted(proofs):
                packed=tar.extractfile(entries[f"tier-{TIER:02d}/part-{part:03d}.csv.gz"]).read()
                # Only one original shard's compressed metadata lives in memory.
                file=args.output/f"tmp-source-part-{part:03d}.gz"
                file.write_bytes(packed)
                result=select_and_write_part(part,file,proofs[part],
                         priority,used_photo,used_obs,args.output,
                         quota=args.target_per_shard,
                         reference=args.reference_per_shard)
                file.unlink(missing_ok=True)
                selected_all+=result["selected_new_original_photo_IDs"]
                total_source+=result["original_source_part_photo_count"]
                full.append(result)
        require(sum(x["original_source_part_photo_count"] for x in full
                    if x["part"]//8==group)
                ==group_doc["counts_by_disjoint_tier"][str(TIER)],
                "One of original 16 source group photo counts does not reconcile")
        print("MOIS_BULK128_ORIGINAL_SOURCE_GROUP_DONE",
              json.dumps({"group":group,"verified_parts":len(proofs),
                          "bulk_source_new_photo_targets":selected_all},
                         sort_keys=True),flush=True)
    require(len(full)==SHARDS and total_source==
            master["new_target_counts_by_tier"][str(TIER)],
            "Signed 128 original source partitions do not conserve source tier")
    all_ids=set();all_obs=set();all_mids=set()
    for p in range(SHARDS):
        for u in range(p*4,p*4+4):
            f=args.output/"shards"/f"shard-{u:03d}.csv"
            with f.open(newline="",encoding="utf-8") as h:
                for row in csv.DictReader(h):
                    pid,obs,mid=row["photo_id"],row["observation_id"],row["measurement_id"]
                    require(pid not in all_ids and obs not in all_obs and mid not in all_mids,
                            "Original image/observation ID reused between 128 shards")
                    require(pid not in used_photo and obs not in used_obs,
                            "New photo acquisition ID was already attempted/reserved")
                    all_ids.add(pid);all_obs.add(obs);all_mids.add(mid)
    require(len(all_ids)==selected_all==len(all_obs)==len(all_mids)
            and selected_all<=SHARDS*args.target_per_shard,
            "Full frozen source photo/observation target roster mismatch")
    manifest={
       "schema":"mois_full_source_observation_sampling_cohort_v1",
       "sampling_schema":SCHEMA,
       "source_tier":TIER,
       "source_inclusion_probability_width":SOURCE_WIDTH,
       "full_frame_photo_ids":PARENT_PHOTOS,
       "all_upstream_128_source_taxon_shards_verified":True,
       "original_photo_source_tier_IDs_scanned":total_source,
       "source_total_original_taxa":42111,
       "source_128_compute_shards":SHARDS,
       "max_photo_attempts_per_compute_shard":args.target_per_shard,
       "new_photo_attempts_planned":selected_all,
       "planned_maximum_if_all_source_strata_supported":
           SHARDS*args.target_per_shard,
       "baseline_photo_and_observation_IDs_checked":prior_counts,
       "colour_measured":False,
       "worker_ids":list(range(SHARDS*4)),
       "work_units":SHARDS*4,
       "group_ids":list(range(SHARDS)),
       "units_per_group":4,
       "original_model_batch_size":256,
       "photo_shard_sha256":{
           k:v for r in full for k,v in r["source_photo_worker_unit_SHA256"].items()},
       "selected_reference_photo_attempts":sum(
           x["selected_photo_attempts_by_design"].get("reference",0) for x in full),
       "selected_targeted_photo_attempts":sum(
           x["selected_photo_attempts_by_design"].get("targeted",0) for x in full),
       "P5_targeted_photo_attempts":0,
       "per_128_shard_receipts_sha256":{
           str(i):sha(args.output/f"part-{i:03d}.receipt.json") for i in range(SHARDS)},
       "original_source_photo_IDs_and_observation_UUIDs_all_unique":True,
       "each_selected_photo_ID_and_observation_distinct_from_old_and_reserved":True,
       "source_selection_exclusively_before_original_image_fetch":True,
       "one_random_reference_and_gap_targeted_arm_not_naively_poolable":True,
       "targeted_adaptive_phase2_photo_inclusion_probability_not_identified":True,
       "original_work_shard_selection_probability_for_targeted_IS_NOT_joint_IPW":True,
       "unknown_images_cannot_be_promoted_to_colours_from_metadata":True,
       "native_wild_anthesis_petal_status_not_independently_verified":True,
       "no_photo_pixel_download_or_colour_inference":True,
    }
    (args.output/"manifest.json").write_text(
        json.dumps(manifest,sort_keys=True,indent=2,allow_nan=False)+"\n")
    (args.output/"bulk128_source_summary.csv").write_text(
        ",".join(("part","photo_targets","reference","targeted",
                 "P1","P2","P3","P4","N_tropical","S_tropical"))+"\n"+
        "".join(",".join(str(v) for v in (
             r["part"],r["selected_new_original_photo_IDs"],
             r["selected_photo_attempts_by_design"].get("reference",0),
             r["selected_photo_attempts_by_design"].get("targeted",0),
             *(r["selected_photo_attempts_by_original_priority"].get(p,0)
               for p in (*NO_COLOUR_PRIORITY,*THIN_PRIORITY)),
             r["selected_photo_attempts_by_zone"].get("N_tropical",0),
             r["selected_photo_attempts_by_zone"].get("S_tropical",0)))+"\n"
            for r in full),encoding="utf-8")
    print("MOIS_REAL_GLOBAL_BULK128_ORIGINAL_SOURCE_ONLY_NEW_ROSTER",
          json.dumps({
           "128_compute_shards":SHARDS,
           "source_original_IDs_scanned":total_source,
           "genuinely_new_photo_and_observation_IDs":selected_all,
           "random_reference_photo_targets":manifest["selected_reference_photo_attempts"],
           "taxon_geo_season_priority_photo_targets":manifest["selected_targeted_photo_attempts"],
           "P5_targeted":0,
           "no_raw_original_photo_downloads":True
          },sort_keys=True),flush=True)
    return manifest

if __name__=="__main__":
    main()
