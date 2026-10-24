#!/usr/bin/env python3
"""FCP-like U100 original-photo source design for MOIS; never open photo pixels.

Adapt the VERIFIED FCP source protocol (100 distinct original observations per
taxon, 256-slot blind ROI-v4 workers, same four terminal morph states, and
taxon-wise >=40-classified-photo eligibility) rather than measuring all 93m
photographs. For this source-only stage, re-use the already SHA-verified
597070-photo tier06 portion of the full 93209351-photo signed registry.

An exact U100-eligible taxon has >=100 DIFFERENT original observation UUIDs in
this tier. Select 100 observations by SHA rank and within each observation a
single original photo by independent SHA rank (not by hue). Record every
nominal taxon with source availability, including 0/<100 taxa, before opening
ANY images. No replacement following a missing or unclassifiable ROI.

Probability within signed source tier06:
  conditional P(observation selected | tier06 and eligible taxon) = 100/N_taxon;
  conditional P(photo selected | observation selected) = 1/M_observation.
The full 93m phase-1 tier hash width is 0.0064. Multiplying gives a *nominal*
two-phase source probability, but the actual sample is CONDITIONAL on an
eligible taxon and realized tier06 stratum N. Do NOT pretend U100-eligible taxa
are a full-probability-complete worldwide species sample. For species <100
eligible observations this source records missing support; future disjoint
tiers must be used for prospective source completion.

This is an FCP-STYLE source+measurement architecture, NOT an assertion of
FCP H9 flowering annotation for each original MOIS photo. Original source
observed_on is not a botanist-confirmed open flower. No native-only filtering,
image downloads, LO-res image exclusion, ROI segmentation or new colours.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile

from scripts.acquisition import mois_93m_taxon_month_lat_balanced_photo_sampler as core
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    EXPECTED_PHOTOS, EXPECTED_TAXA, sha256, require,
    normalized_extension,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import (
    WORKER_COLUMNS, read_taxa,
)

SCHEMA = "mois_fcp_style_exact100_original_observations_per_taxon_tier06_source_v1"
WORKER_SCHEMA = "mois_full_source_observation_sampling_cohort_v1"
DESIGN = "MOIS_FCP_U100_100_DISTINCT_OBSERVATIONS_PER_NOMINAL_TAXON"
FCP_SOURCE_TARGET = 100
FCP_MIN_CLASSIFIABLE = 40
SOURCE_ARM = "mois_93m_tier06_fcp_u100_100obs_singlephoto_source_v1"
FULL_TIER_PHOTOS = core.SOURCE_PHOTOS


def original_taxon_source_counts(db):
    found={}
    for taxon,n_obs,n_photos in db.execute(
            "SELECT taxon,COUNT(*) AS n_obs,SUM(photo_count) AS n_photos "
            "FROM obs GROUP BY taxon ORDER BY taxon"):
        require(taxon not in found and n_obs>0 and n_photos>=n_obs,
                "Incomplete original taxon observation/photo source inventory")
        found[taxon]=(n_obs,n_photos)
    return found


def sample_fcp_u100(db, names, *, target=FCP_SOURCE_TARGET):
    """One uniform-hash-observation draw of EXACTLY 100 for each eligible taxon.

    An observation's random photo choice is independent of its observation
    selection rank; the final table contains source selection probabilities.
    Every lower-support taxon stays in the availability table with no invented
    100-photo guarantee. No within-source pseudo-replication of 2+ images
    attached to one biological observation is permitted.
    """
    require(target == FCP_SOURCE_TARGET,
            "The FCP U100 protocol is fixed at exactly 100 original photos/taxon")
    inventory=original_taxon_source_counts(db)
    require(len(names)>0 and set(inventory)<=set(names),
            "Taxa in authentic photo metadata absent from 42,111-taxon source")
    eligible={taxon:n for taxon,(n,_) in inventory.items() if n>=target}
    # Sorting uses a photo-independent observation rank. The 1/M within-
    # observation draw was precomputed over ALL photos within that observation
    # during first-pass source ingestion; its photographic colour was unknown.
    sql = """
      SELECT taxon,uuid,photo_count,best_photo_id,photo_license,
             photo_extension,recorded_on,record_month,lat_band,
             latitude,longitude,observer_id,positional_accuracy,
             obs_rank,
             ROW_NUMBER() OVER (
               PARTITION BY taxon ORDER BY obs_rank,uuid
             ) AS within_taxon_rank
      FROM obs ORDER BY taxon,obs_rank,uuid
    """
    rows=[]
    regions=defaultdict(Counter)
    months=defaultdict(Counter)
    seen_obs=set()
    for (taxon,uuid,m,pid,license_code,extension,on,month,band,
         lat,lon,observer_id,accuracy,obsrank,rank) in db.execute(sql):
        if taxon not in eligible or rank>target:
            continue
        require(uuid not in seen_obs and m>=1 and pid.isdecimal()
                and taxon in names and 1<=rank<=target,
                "FCP U100 selected original photo/observation identity invalid")
        seen_obs.add(uuid)
        n=eligible[taxon]
        conditional = target/n/m
        nominal = core.SOURCE_PROBABILITY * conditional
        require(0 < nominal <=core.SOURCE_PROBABILITY,
                "Original observed source tier/100-photo inclusion probability invalid")
        ext=normalized_extension(extension)
        row={
            "photo_id":pid,"observation_id":uuid,"inat_taxon_id":taxon,
            "species":names[taxon],"observer_id":observer_id,
            "observed_on":on,"recorded_month":month,"latitude_band":band,
            "latitude":lat,"longitude":lon,"positional_accuracy":accuracy,
            "photo_license":license_code,"original_image_extension":extension,
            "photos_attached_to_original_observation":m,
            "original_observations_in_taxon_tier06":n,
            "conditional_observation_inclusion_probability":format(target/n,".17g"),
            "conditional_one_photo_per_observation_probability":format(1/m,".17g"),
            "conditional_photo_inclusion_probability_in_tier06":format(conditional,".17g"),
            "nominal_source_tier_plus_photo_probability":format(nominal,".17g"),
            "photo_source_disposition":"supported" if ext else "unsupported_original_extension",
            "photo_url_large":(
                f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}"
                if ext else ""
            ),
            "source_arm":SOURCE_ARM,
        }
        rows.append(row)
        regions[taxon][band]+=1
        months[taxon][month]+=1

    require(len(rows)==len(eligible)*target
            and len(seen_obs)==len(rows)
            and len({r["photo_id"] for r in rows})==len(rows),
            "Original source photo/observation FCP U100 support not exact/unique")
    return rows,inventory,eligible,months,regions


def run(registry, routed, counts, output, *,
        expected_photos=EXPECTED_PHOTOS, expected_taxa=EXPECTED_TAXA,
        expected_tier_photos=FULL_TIER_PHOTOS,
        pinned=core.PINNED, check_taxa_sha=True,
        target=FCP_SOURCE_TARGET):
    registry,routed,counts,output=map(Path,(registry,routed,counts,output))
    require(not output.exists(),"Refuse to overwrite frozen FCP U100 source assignments")
    require(target==100,"Only exact FCP 100-photos-per-taxon sampling supported")
    upstream=json.loads(registry.read_text(encoding="utf-8"))
    require(upstream.get("schema")==core.TIER_REGISTRY_SCHEMA
            and upstream.get("exact_original_photo_IDs_in_all_disjoint_tiers")==expected_photos
            and upstream.get("original_named_taxa_from_full_source")==expected_taxa
            and upstream.get("new_target_counts_by_tier",{}).get("6")==expected_tier_photos
            and upstream.get("source_only_NO_photo_pixels_or_hue_measured") is True
            and upstream.get("native_vs_introduced_not_hard_filtered") is True,
            "Full original 93m-source photo/observation registry not proven")
    names=read_taxa(counts,check_sha=check_taxa_sha)
    require(len(names)==expected_taxa,
            "Source-locked full nominal plant taxon roster changed")
    output.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(tempfile.mkdtemp(prefix="mois-u100-source-",dir=output.parent))
    try:
        with sqlite3.connect(tmp/"original_observations.sqlite") as db:
            core.create_database(db)
            group_records=[]
            for n in range(core.SOURCE_GROUPS):
                group_records.append(core.load_signed_group(
                    db,routed,n,pinned=pinned,expected_full=expected_photos))
                db.commit()
            source_count=sum(int(g["original_photos"]) for g in group_records)
            observed,all_photo_counts=db.execute(
                "SELECT COUNT(*),COALESCE(SUM(photo_count),0) FROM obs"
            ).fetchone()
            require(source_count==expected_tier_photos
                    and all_photo_counts==source_count,
                    "FCP 100-observation source frame lost genuine photo opportunities")
            selected,inventory,eligible,months,regions=sample_fcp_u100(
                db,names,target=target)
        del all_photo_counts
        root=tmp/"original_source"
        root.mkdir()
        # The complete 42,111-taxonomical source roster includes under-100
        # and zero-observation taxa; they are never silently discarded.
        availability=root/"all_nominal_taxa_FCP_U100_source_eligibility.csv"
        with availability.open("w",newline="",encoding="utf-8") as f:
            w=csv.writer(f,lineterminator="\n")
            w.writerow(("inat_taxon_id","species","original_observations_in_tier06",
                        "original_photo_ids_in_tier06","fcp_U100_source_eligible",
                        "source_photo_deficit_to_U100","selected_observation_months",
                        "selected_latitude_bands"))
            for taxon,species in sorted(names.items(),key=lambda x:int(x[0])):
                n,m=inventory.get(taxon,(0,0))
                w.writerow((taxon,species,n,m,str(n>=target).lower(),
                           max(0,target-n),
                           len(months.get(taxon,{})),
                           len(regions.get(taxon,{}))))
        selected_file=root/"exact_100_original_photo_per_eligible_taxon.csv"
        fields=list(selected[0]) if selected else [
            "photo_id","observation_id","inat_taxon_id","species"
        ]
        with selected_file.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
            w.writeheader();w.writerows(selected)
        unsupported=[r for r in selected if r["photo_source_disposition"]!="supported"]
        unsupported_file=root/"selected_photo_unsupported_loader_as_unknown.csv"
        with unsupported_file.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
            w.writeheader();w.writerows(unsupported)

        supported=[r for r in selected if r["photo_source_disposition"]=="supported"]
        # Preserves FCP's source-photo upfront blind assignments; colours,
        # geographies and species are NOT sent to its 3-column mask worker.
        # This is an immutable source-only future ROI candidate, not a job launch.
        sha_per={}
        shards=root/"shards"
        shards.mkdir()
        supported.sort(key=lambda x:(core.hash64(
            "mois-fcp-u100-acquisition-worker-order-v1",
            x["observation_id"]),int(x["photo_id"])))
        for unit,start in enumerate(range(0,len(supported),256)):
            file=shards/f"shard-{unit:03d}.csv"
            with file.open("w",newline="",encoding="utf-8") as f:
                w=csv.DictWriter(f,fieldnames=WORKER_COLUMNS,lineterminator="\n")
                w.writeheader()
                for r in supported[start:start+256]:
                    mid="MOIS-"+hashlib.sha256(
                        ("mois-original-global-RGB-v1|"+r["photo_id"]).encode()
                    ).hexdigest()[:24].upper()
                    w.writerow({
                        "measurement_id":mid,"photo_id":r["photo_id"],
                        "inat_taxon_id":r["inat_taxon_id"],"species":r["species"],
                        "observation_id":r["observation_id"],
                        "observer_id":r["observer_id"],
                        "observed_on":r["observed_on"],
                        "latitude":r["latitude"],"longitude":r["longitude"],
                        "positional_accuracy":r["positional_accuracy"],
                        "photo_license":r["photo_license"],
                        "photo_url_large":r["photo_url_large"],
                        "image_filename":mid+".jpg",
                        "source_arm":SOURCE_ARM,
                        "selection_probability":r["nominal_source_tier_plus_photo_probability"],
                        "source_snapshot_date":"2026-09-27","shard":str(unit),
                    })
            sha_per[str(unit)]=sha256(file)
        n_units=len(sha_per)
        measured_future_count=len(supported)
        occupancy=Counter("FCP_U100_exact100_source" if n>=100 else
                          "under100_at_current_source_tier"
                          if n>0 else "absent_at_this_tier"
                          for n,_ in (inventory.get(t,(0,0)) for t in names))
        summary={
            "schema":SCHEMA,
            "fcp_original_source_protocol":"100 photo opportunities from 100 distinct observations per eligible nominal taxon",
            "fcp_original_four_hue_roi_protocol":"FCP frozen original ROI-v4 / EfficientSAM; no altered thresholds",
            "full_93m_original_photo_source_ids":expected_photos,
            "full_original_nominal_taxa":expected_taxa,
            "stage1_source_tier":6,
            "source_tier_original_photo_ids":source_count,
            "source_tier_distinct_original_observations":int(observed),
            "distinct_original_taxa_with_at_least_one_source_observation":len(inventory),
            "FCP_U100_original_source_eligible_taxa":len(eligible),
            "FCP_U100_exact100_source_photo_opportunities":len(selected),
            "fcp_minimum_classified_ROIs_per_taxon_for_species_conditional_inference":FCP_MIN_CLASSIFIABLE,
            "exactly_100_distinct_original_observations_per_eligible_taxon":True,
            "source_eligibility_counts":dict(occupancy),
            "FCP_U100_selected_supported_photo_ROI_targets":measured_future_count,
            "FCP_U100_selected_unsupported_source_photo_unknowns":len(unsupported),
            "under100_taxa_remain_in_complete_roster_for_later_tier_completion":True,
            "source_eligibility_is_outcome_blind_but_source_abundance_biased":True,
            "conditional_observation_probability_formula":"100/original_observation_count_of_eligible_taxon_in_tier06",
            "conditional_one_photo_probability_formula":"1/photo_count_of_original_observation",
            "source_tier_observation_random_hash_width":core.SOURCE_PROBABILITY,
            "nominal_joint_photo_probability_formula":"0.0064*100/N_taxon/M_observation conditional on source eligibility and realized source tier",
            "not_unconditionally_design_unbiased_across_all_42111_taxa":True,
            "fcp_H9_iNaturalist_flowering_annotation_verified_for_this_MOIS_photo_frame":False,
            "per_photo_individual_native_wild_flowering_status_verified":False,
            "source_group_SHA_receipts":group_records,
            "signed_full_source_registry_SHA256":sha256(registry),
            "signed_taxon_source_SHA256":sha256(counts),
            "original_photo_U100_assignments_SHA256":sha256(selected_file),
            "nominal_taxon_eligibility_complete_roster_SHA256":sha256(availability),
            "unsupported_source_photos_SHA256":sha256(unsupported_file),
            "source_original_photo_shard_SHA256":sha_per,
            "worker_original_photo_units":n_units,
            "original_ROI_photo_pixels_downloaded":0,
            "original_four_hue_image_ROIs_inferred":0,
            "source_arm":SOURCE_ARM,
            "no_prior_actual_photo_reuse_audit_has_been_run_for_this_new_U100_frame":True,
            "full_frame_photo_ids":expected_photos,
        }
        (root/"fcp_U100_source_only_receipt.json").write_text(
            json.dumps(summary,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        (tmp/"original_observations.sqlite").unlink(missing_ok=True)
        os.replace(root,output)
        shutil.rmtree(tmp,ignore_errors=True)
        print("MOIS_FCP_U100_SOURCE_OPPORTUNITY_REAL93M_TIER06",
              json.dumps({k:summary[k] for k in (
                  "source_tier_original_photo_ids",
                  "source_tier_distinct_original_observations",
                  "FCP_U100_original_source_eligible_taxa",
                  "FCP_U100_exact100_source_photo_opportunities",
                  "FCP_U100_selected_supported_photo_ROI_targets",
                  "source_eligibility_counts")},sort_keys=True),flush=True)
        return summary
    except BaseException:
        shutil.rmtree(tmp,ignore_errors=True)
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("registry","routed","taxa","output"):
        p.add_argument("--"+key,type=Path,required=True)
    args=p.parse_args()
    run(args.registry,args.routed,args.taxa,args.output)


if __name__ == "__main__":
    main()
