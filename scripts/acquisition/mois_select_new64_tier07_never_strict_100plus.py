#!/usr/bin/env python3
"""MOIS minimal 64 ORIGINAL photo bias-controlled new-image sample from frozen 93m.

All 42,111 taxa with source photos have already had at least 1 real original
image attempted. Consequently "taxa with >=100 original photos and 0 previous
image attempts" is an EMPTY pool. This supplementary pilot samples taxa with
>=100 source photo opportunities, 0 strict original ROI9 colours, and 1-4
actual prior photo attempts. It is an OUTCOME-CONDITIONED supplemental design,
not a global probability sample or unbiased flora color frequency estimator.

Read all 16 original tier07 PHOTO archives (observation hash tier probability
0.0128) authenticated against all-source photo-shard SHA manifest; retain
prior photo-id and observation UUID exclusion; select uniformly by independent
hash from species / observations / attached photos; stratify by source
latitude band × number of prior photo attempts.

Cap: eight species/photos per 4 geographic bands × 2 prior effort bins =
64 distinct new species, source PHOTO IDs and original observation UUIDs.
Log candidate population, conditional selection probabilities and failed
photo bytes individually; never fill quotas from photographed-success
outcomes. The original full 93m photo opportunity frame remains eligible.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,math,tarfile
from collections import Counter,defaultdict
from pathlib import Path
import pandas as pd
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, RESULT_SCHEMA, TIER_ENDS, obs_draw,
)
SCHEMA="mois_minimal_64_unmeasured_100plus_taxa_tier07_geo_effort_sampling_v1"
TAG="mois-complete-93209351-original-photo-target-tiers-v1"
SEED="mois-20261026-minimal64-tier07-source-taxon-geo-phase2-v1"
TIER=7
PARENT_N=93_209_351
TAXA_N=42_111
OLD_N=1_043_930
SOURCE_PHOTO_SHA="df0f562b58bdc7ba06e901530eb809f513a763e2b0524d3b195e99ee56e930ef"
SOURCE_UNION_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
OLD_PRIOR_SOURCE_SHA="b20bbefe15020c21c3b5b0f7edf3f8e08148a4ee7dda8bf772a89bf6b126f447"
ORIGINAL_TAXON_PRIORITY_SHA="ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
STRATA=(
    "north_extratropical",
    "south_extratropical",
    "north_tropical",
    "south_tropical"
)
EFFORT=("one_attempt","two_to_four_attempts")
TARGET_PER_CELL=8

def require(condition,reason):
    if not condition:
        raise ValueError("MOIS_NEW64_SOURCE: "+reason)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for c in iter(lambda:f.read(4<<20),b""):h.update(c)
    return h.hexdigest()

def rank(*parts):
    return hashlib.sha256(("|".join((SEED,*map(str,parts)))).encode()).hexdigest()

def band4(latitude):
    try: v=float(latitude)
    except (ValueError,TypeError): return "unknown"
    if not math.isfinite(v) or abs(v)>90:return "unknown"
    if v>=23.5:return "north_extratropical"
    if v<=-23.5:return "south_extratropical"
    return "north_tropical" if v>=0 else "south_tropical"

def effort_class(attempts):
    if int(attempts)==1:return EFFORT[0]
    if 2<=int(attempts)<=4:return EFFORT[1]
    return None

def read_taxa(register,old_new):
    # per-species source counts are frozen in earlier ID-level original audit.
    with gzip.open(register,"rt",encoding="utf-8") as f:
        taxa=pd.read_csv(f,dtype={"inat_taxon_id":str})
    require(len(taxa)==TAXA_N and taxa.inat_taxon_id.is_unique,
            "Original source taxon opportunity register no longer 42,111 rows")
    old_new.photo_id=old_new.photo_id.astype(str)
    old_new.inat_taxon_id=old_new.inat_taxon_id.fillna("").astype(str)
    real=old_new.groupby("inat_taxon_id",observed=True).continuous9_strict.sum()
    old_attempts=taxa.set_index("inat_taxon_id").n_image_photo_attempts.astype(int)
    counts=taxa.set_index("inat_taxon_id").n_unique_source_photo_ids.astype(int)
    require(len(old_new)==OLD_N and old_new.photo_id.is_unique
            and int(old_new.continuous9_strict.sum())==400710,
            "Current original 400710 strict colour source changed")
    # Actual source priority register originally counted 390541; rebase zero
    # source-colour taxa to the recovered 400710, not the previous count.
    now=real.reindex(counts.index,fill_value=0).astype(int)
    eligible=pd.DataFrame({
        "inat_taxon_id":counts.index.to_numpy(),
        "species":taxa.set_index("inat_taxon_id").species.to_numpy(),
        "n_source_photos":counts.to_numpy(),
        "n_previous_attempts":old_attempts.to_numpy(),
        "n_current_strict9":now.to_numpy(),
    }).set_index("inat_taxon_id")
    none=int((counts.ge(100)&old_attempts.eq(0)).sum())
    require(none==0,"Zero-image species pool unexpectedly newly appeared")
    chosen=eligible.loc[(eligible.n_source_photos>=100)
              &eligible.n_previous_attempts.between(1,4)
              &eligible.n_current_strict9.eq(0)].copy()
    require(len(chosen)>=1000,"Original source 100plus/zero9/1to4 previous attempt pool unavailable")
    return eligible,chosen

def load_one_group(tarfile_path,receipt,master,seen,old_photos,old_obs,eligible,photo_by_obs,counters):
    i=int(tarfile_path.name[len("tier-07-group-"):-4])
    document=json.loads(receipt.read_text(encoding="utf-8"))
    expected=master["published_tier_archive_SHA256"].get(tarfile_path.name)
    require(expected and sha(tarfile_path)==expected==
            document["archive_sha256_by_tier"][tarfile_path.name]
            and document["schema"]==RESULT_SCHEMA,
            "Not a complete source-signed original tier07 archive: "+str(i))
    members={int(x["part"]):x for x in document["parts"]}
    require(len(members)==8 and all(n//8==i for n in members),
            "Source group not eight disjoint authoritative taxon-modulo shards")
    total=0
    with tarfile.open(tarfile_path,"r:") as tf:
        memberset={m.name.removeprefix("./"):m for m in tf.getmembers() if m.isfile()}
        expected_names={f"tier-07/part-{j:03d}.csv.gz" for j in members}
        require(set(memberset)==expected_names,
                "Original source release tier07 tar members incorrect")
        for j,meta in sorted(members.items()):
            member=memberset[f"tier-07/part-{j:03d}.csv.gz"]
            raw=tf.extractfile(member).read()
            require(hashlib.sha256(raw).hexdigest()==
                    meta["ranked_target_csv_sha256"][str(TIER)],
                    "Actual 93m photo-source shard file SHA changed")
            count=0
            with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
                with io.TextIOWrapper(gz,encoding="utf-8",newline="") as stream:
                    reader=csv.DictReader(stream)
                    require(tuple(reader.fieldnames or [])==COLUMNS,
                            "Original tier07 source columns changed")
                    for row in reader:
                        count+=1
                        pid,uuid,taxon=row["photo_id"],row["observation_uuid"],row["inat_taxon_id"]
                        require(pid.isdecimal() and uuid and taxon.isdecimal()
                                and row["tier"]==str(TIER)
                                and float(row["observation_inclusion_score"])==obs_draw(uuid),
                                "Source observation/photo identity or tier hash corrupted")
                        counters["tier07_original_photo_opp"]+=1
                        if taxon not in eligible.index:
                            continue
                        counters["source_species_eligibility_photo_opp"]+=1
                        if pid in old_photos:
                            counters["previously_attempted_same_photo"]+=1
                            continue
                        if uuid in old_obs:
                            counters["previously_attempted_same_observation"]+=1
                            continue
                        b=band4(row["latitude"])
                        if b not in STRATA:
                            counters["missing_or_invalid_geography_source"]+=1
                            continue
                        if row["frozen_loader_supported"]!="true":
                            counters["unsupported_image_extension"]+=1
                            continue
                        ext=row["image_extension"].lower().strip()
                        if ext not in ("jpg","jpeg","png"):
                            counters["unsupported_image_extension"]+=1
                            continue
                        key=uuid
                        order=rank("photo",uuid,pid)
                        prev=photo_by_obs.get(key)
                        if prev is None or order<prev[0]:
                            photo_by_obs[key]=(order,{**row,"source_part":j,
                                                       "source_group":i,"geographic_band":b})
            require(count==meta["counts_by_disjoint_tier"][str(TIER)],
                    "Disjoint original source photo shard count drift")
            total+=count
    require(total==document["counts_by_disjoint_tier"][str(TIER)],
            "Complete original tier07 group count differs from full original source index")
    return total

def design(groups,registry,register,atlas,atlas_receipt,union,union_receipt,out,
           per_cell=TARGET_PER_CELL):
    groups,registry,register,atlas,atlas_receipt,union,union_receipt,out=map(
        Path,(groups,registry,register,atlas,atlas_receipt,union,union_receipt,out))
    require(not out.exists(),"Refuse modifying any prior original source pilot sample")
    require(1<=per_cell<=32,"Only small, predeclared sampling quotas accepted")
    reg=json.loads(registry.read_text())
    require(reg["schema"]=="mois_complete_full93m_original_photo_target_registry_v1"
            and reg["exact_original_photo_IDs_in_all_disjoint_tiers"]==PARENT_N
            and reg["original_named_taxa_from_full_source"]==TAXA_N,
            "Parent 93m source target index is not authentic")
    receipt=json.loads(atlas_receipt.read_text())
    previous=json.loads(union_receipt.read_text())
    require(receipt["new_source_photo9_lake_SHA256"]==SOURCE_PHOTO_SHA==sha(atlas)
            and receipt["new_original_strict9_photo_count"]==400710
            and receipt["full_original_photo_attempts"]==OLD_N
            and previous["source_pixel_attempts"]==OLD_N
            and previous["photo_parquet_sha256"]==SOURCE_UNION_SHA==sha(union),
            "Previous 1043930 genuinely obtained photos or recovered nine-ROI data changed")
    colour=pd.read_parquet(atlas,columns=["photo_id","inat_taxon_id","continuous9_strict"])
    previous_ids=pd.read_parquet(union,columns=["photo_id","observation_id"])
    colour.photo_id=colour.photo_id.astype(str)
    previous_ids.photo_id=previous_ids.photo_id.astype(str)
    require(len(previous_ids)==OLD_N
            and colour.photo_id.tolist()==previous_ids.photo_id.tolist()
            and previous_ids.photo_id.is_unique,"Original ID and observation identity drift")
    original_photo=set(previous_ids.photo_id)
    original_observations=set(previous_ids.observation_id.dropna().astype(str))
    original_observations.discard("")
    taxa,eligible=read_taxa(register,colour)
    candidates={}
    stats=Counter()
    archive_sha={}
    total=0
    for group in range(16):
        path=groups/f"tier-07-group-{group:02d}.tar"
        receiptfile=groups/f"group-{group:02d}-receipt.json"
        require(path.is_file() and receiptfile.is_file(),
                "One signed full tier07 source group absent")
        n=load_one_group(path,receiptfile,reg,None,original_photo,
                         original_observations,eligible,candidates,stats)
        total+=n
        archive_sha[path.name]=sha(path)
    require(total==reg["new_target_counts_by_tier"][str(TIER)],
            "All 16 signed original source tier07 archive IDs not fully visited")
    # One independent original observation and PHOTO per target taxon. Source
    # taxon/photo ID selection uses deterministic salt unrelated to colour.
    eligible_distinct_new_observations=len(candidates)
    by_species=defaultdict(list)
    for _,row in candidates.values():
        by_species[row["inat_taxon_id"]].append(row)
    one_by_species=[]
    for taxon,photos in by_species.items():
        choice=min(photos,key=lambda r:rank("observation",taxon,r["observation_uuid"]))
        choice["source_tier07_eligible_observations_for_this_taxon"]=len(photos)
        one_by_species.append(choice)
    design_cells=defaultdict(list)
    for photo in one_by_species:
        taxon=photo["inat_taxon_id"]
        effort=effort_class(eligible.loc[taxon,"n_previous_attempts"])
        require(effort is not None,"Unexpected photo attempt class outside source selection")
        design_cells[(photo["geographic_band"],effort)].append(photo)
    selected=[];cell_stats=[]
    for band in STRATA:
        for eff in EFFORT:
            candidates=design_cells[(band,eff)]
            order=sorted(candidates,key=lambda row:rank("taxon",band,eff,row["inat_taxon_id"]))
            take=order[:per_cell]
            cell_stats.append({"region":band,"prior_attempt_class":eff,
                "source_tier07_after_existing_photo_observation_dedup_taxa":len(candidates),
                "selected_new_distinct_taxa_and_observations":len(take),
                "conditional_species_sampling_probability_given_eligible_tier07":
                    len(take)/len(candidates) if candidates else None,
                "required_target_photos_not_filled_for_lack_of_candidates":
                    per_cell-len(take)})
            for row in take:
                row["design_effort_stratum"]=eff
                row["design_region_stratum"]=band
                row["source_tier07_eligible_taxa_in_stratum"]=len(candidates)
                row["conditional_tier07_species_probability"]=(len(take)/len(candidates))
                selected.append(row)
    require(len(selected)>=24,
            "No bias-aware small original-source sample available; report source instead")
    require(len({r["photo_id"] for r in selected})==len(selected) and
            len({r["observation_uuid"] for r in selected})==len(selected) and
            len({r["inat_taxon_id"] for r in selected})==len(selected) and
            all(r["photo_id"] not in original_photo and
                r["observation_uuid"] not in original_observations
                for r in selected),
            "New photo is duplicate original image, observation or taxon")
    source_rows=[];ledger=[]
    selected.sort(key=lambda r:rank("measurement_worker",r["photo_id"]))
    for r in selected:
        pid=r["photo_id"];taxon=r["inat_taxon_id"]
        attrs=eligible.loc[taxon]
        ext=r["image_extension"].lower()
        source_rows.append({
            "photo_id":pid,"observation_uuid":r["observation_uuid"],
            "inat_taxon_id":taxon,"observer_id":r["observer_id"],
            "photographer_id":"",
            "observed_on":r["observed_on"],
            "latitude":r["latitude"],"longitude":r["longitude"],
            "positional_accuracy":r["positional_accuracy"],
            "photo_license":r["photo_license"],"image_extension":ext,
            "frozen_loader_supported":"true",
            "valid_observed_date":r["valid_observed_date"],
            "valid_geopoint":r["valid_geopoint"],
            "observation_inclusion_score":r["observation_inclusion_score"],
            "tier":str(TIER)
        })
        ledger.append({
            "photo_id":pid,"observation_id":r["observation_uuid"],
            "taxon_id":taxon,"species":str(attrs.species),
            "original_fullframe_taxon_photo_opportunities":int(attrs.n_source_photos),
            "prior_photo_attempts":int(attrs.n_previous_attempts),
            "prior_valid9_colour_ROIs":int(attrs.n_current_strict9),
            "source_hemisphere_band":r["design_region_stratum"],
            "sampling_prior_attempt_class":r["design_effort_stratum"],
            "original_date":r["observed_on"],
            "source_93m_tier":TIER,"source_tier_presentation_probability":(
                 TIER_ENDS[TIER]-TIER_ENDS[TIER-1]),
            "source_eligible_tier07_observations_for_taxon":
                 r["source_tier07_eligible_observations_for_this_taxon"],
            "source_eligible_tier07_taxa_in_pilot_stratum":
                 r["source_tier07_eligible_taxa_in_stratum"],
            "conditional_species_probability_given_eligible_source_tier":
                 r["conditional_tier07_species_probability"],
            "conditional_one_observation_probability_given_eligible_taxon":
                 1/r["source_tier07_eligible_observations_for_this_taxon"],
            "actual_original_photoID_already_attempted":False,
            "original_observation_already_attempted":False,
            "original_global_93m_photo_probability_UNIDENTIFIED_from_selected_legacy_overlap":True
        })
    out.mkdir(parents=True)
    sourcefile=out/"mois_bounded_bias_aware_tier07_new64_original_source.csv.gz"
    with gzip.open(sourcefile,"wt",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
        writer.writeheader();writer.writerows(source_rows)
    details=out/"mois_bounded_new_taxon_original_photo_sampling_provenance.csv"
    with details.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(ledger[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(ledger)
    report={
        "schema":SCHEMA,"full_93m_parent_photo_ids":PARENT_N,
        "full_source_original_named_taxa":TAXA_N,
        "all_1043930_previously_attempted_original_photo_ids_excluded":OLD_N,
        "previous_9colour_success_source":400710,
        "source_tier":TIER,
        "source_tier_photos_authenticated_and_visited":total,
        "source_tier_probability_width":TIER_ENDS[TIER]-TIER_ENDS[TIER-1],
        "eligible_taxa_100plus_source_photo_and_nine_Roi0_prior1to4":len(eligible),
        "eligible_source_distinct_observations_post_existing_dedup":eligible_distinct_new_observations,
        "eligible_target_taxa_with_one_original_source_observation":len(one_by_species),
        "actual_new_original_photo_IDs_SELECTED_not_yet_image_downloaded":len(source_rows),
        "one_photo_one_source_observation_one_taxon_in_minimal_pilot":True,
        "geographic_and_previous_effort_strata":cell_stats,
        "all_16_original_tier07_source_tar_archives_SHA_verified":True,
        "original_tier07_group_SHA256":archive_sha,
        "genuine_selection_probability_not_global_without_legacy_design":True,
        "all_photo_rows_eligible_regardless_of_flower_bloom_status":True,
        "sampling_conditions_depend_on_PREVIOUS_ROI_success_and_are_not_population_representative":True,
        "image_bytes_newly_downloaded":0,
        "original_ROI_or_photo_Lab_inferred":0,
        "photo_source_SHA256":sha(sourcefile),
        "provenance_file_SHA256":sha(details),
        "source_latest1043930_photo_sha256":SOURCE_PHOTO_SHA,
        "source_original_obs_ID_sha256":SOURCE_UNION_SHA,
    }
    (out/"mois_bounded_bias_aware_newphoto_precolour_source_receipt.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_SOURCE_NEW_TAXA_100PLUS_ZERO9_BIAS_AWARE_MINIMAL_PRECOLOUR",
         json.dumps({k:v for k,v in report.items() if k not in (
                "original_tier07_group_SHA256",)},sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("groups","registry","register","atlas","atlas-receipt",
                 "union","union-receipt","output"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--per-cell",type=int,default=TARGET_PER_CELL)
    a=p.parse_args()
    design(a.groups,a.registry,a.register,a.atlas,a.atlas_receipt,
           a.union,a.union_receipt,a.output,per_cell=a.per_cell)
if __name__=="__main__":
    main()
