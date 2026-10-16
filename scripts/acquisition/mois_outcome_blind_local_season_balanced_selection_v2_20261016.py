#!/usr/bin/env python3
"""Prospective alternate Mois source-photo design: local-season-balanced v2.

The FCP 42,111-taxon / 93,209,351 original photo-ID source and the Mois
maximum of forty DISTINCT source observations per taxon do not change.

V1 chose <=40 observations using original region × annual BIO12 strata.
V2 preserves those exact stratum quotas, but distributes each stratum's
allocation over observed local seasons when its quota covers all occupied
seasons. If a stratum has MORE occupied seasons than available quota, V2
falls back to the original unstratified stratum sample, so every original
observation has positive (design) inclusion probability. No photo pixel,
historical colour, ROI result or Flowers annotation influences selection.

This is an explicitly DIFFERENT prospectively frozen photo cohort, not a
replacement or an independent replication of V1. No primary estimand is
chosen retrospectively from hue data, and duplicates between designs count
once if the resulting cohorts are subsequently compared.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.acquisition.mois_stratified_photo_selection_20261016 import (
    CAP, PRIVATE_FIELDS, classify_site, priority, quotas,
)
from scripts.analysis.mois_seasonal_flower_colour_opportunity_20261016 import (
    local_season,
)
from scripts.acquisition.mois_sample_direct_from_fcp_source_20261016 import (
    PinnedBio12, calibrate, run as run_direct,
)

SCHEMA="mois_direct_fcp_source_stratified_photo_selection_v2"
STRATEGY="OUTCOME_BLIND_NESTED_LOCAL_SEASON_WITHIN_ORIGINAL_REGION_BIO12_V2"
# Original historical source: 45,953/100,543 classifiable photos.
# At p=0.457 independently per observation, 10 per each of two seasons
# yields ~0.825 probability that BOTH have >=3 classifiable images.
# This is a design scenario, not a guarantee about new season/region photos.
MIN_TARGET_PER_OCCUPIED_SEASON=10


def nested_season_quotas(opportunities:dict[str,int],quota:int)->dict[str,int]|None:
    """Positive first-order pi for every source observation.

    Return None when the region×BIO12 quota is too small even to cover
    all occupied season strata; use V1 whole-stratum SRS in that case.
    """
    if not opportunities or not 1<=quota<=sum(opportunities.values()):
        raise ValueError("Invalid source stratum annual precipitation quota")
    if any(n<=0 for n in opportunities.values()):
        raise ValueError("Season with no source observations must not be occupied")
    if quota<len(opportunities):
        return None
    chosen={phase:1 for phase in sorted(opportunities)}
    while sum(chosen.values())<quota:
        # First lift occupied seasons toward a 10-observation planning target
        # (or actual available source count), then allocate by sqrt(N)
        # allocation principle (bounded by actual seasonal opportunities).
        support=[
            phase for phase in chosen
            if chosen[phase]<min(MIN_TARGET_PER_OCCUPIED_SEASON,opportunities[phase])
        ]
        if support:
            phase=min(support,key=lambda p:(chosen[p],p))
        else:
            support=[phase for phase in chosen
                     if chosen[phase]<opportunities[phase]]
            if not support:
                raise ValueError("Seasonal source quota capacity contradicts opportunity counts")
            phase=min(support,key=lambda p:(
                chosen[p]/math.sqrt(opportunities[p]),p
            ))
        chosen[phase]+=1
    if (sum(chosen.values())!=quota or
        any(not 1<=k<=opportunities[p] for p,k in chosen.items())):
        raise AssertionError("Nested season quota normalization failed")
    return chosen


def sample_species_season_v2(taxon:int,genus:str,rows,writer):
    observations={}
    photo_seen=set()
    n_photos=0
    for row in rows:
        n_photos+=1
        oid=row["observation_id"].strip()
        photo=row["photo_id"].strip()
        if (not oid or not photo or photo in photo_seen
            or not row["photo_license"].strip() or row["genus"].strip()!=genus):
            raise ValueError("Missing/duplicate source photograph or genus/license drift")
        photo_seen.add(photo)
        group=classify_site(row)
        _,season,_,_=local_season(row["latitude"].strip(),row["observed_on"].strip())
        obs=observations.get(oid)
        if obs is not None:
            if (obs["group"]!=group or obs["season"]!=season
                or obs["observed_on"]!=row["observed_on"].strip()
                or obs["observer"]!=row["observer_id"].strip()):
                raise ValueError("One original observation has inconsistent site/season/observer")
            obs["m"]+=1
        score=priority("representative-photo",photo)
        if obs is None or score<obs["photo_score"]:
            observations[oid]={
                "group":group,"season":season,
                "observed_on":row["observed_on"].strip(),
                "observer":row["observer_id"].strip(),
                "m":1 if obs is None else obs["m"],
                "photo_score":score,"row":row,
            }
    if not observations:
        raise ValueError("A nonzero source taxon has no original observations")
    by_geo=defaultdict(list)
    for oid,obs in observations.items():
        by_geo[obs["group"]].append((
            priority("observation",f"{taxon}/{oid}"),oid,obs
        ))
    quota_geo=quotas({group:len(items) for group,items in by_geo.items()})
    bins=defaultdict(int)
    pilot=int(priority("taxon-pilot",str(taxon))%10==0)
    selected=0
    for group in sorted(by_geo):
        entries=by_geo[group]
        q=quota_geo[group]
        by_season=defaultdict(list)
        for entry in entries:
            by_season[entry[2]["season"]].append(entry)
        seasonal=nested_season_quotas(
            {phase:len(items) for phase,items in by_season.items()},q
        )
        if seasonal is None:
            # Crucial nonzero sampling inclusion for all occupied seasons
            # when there are more seasons than available places.
            picked=[(entry,len(entries),q,"FALLBACK_ALL_SEASONS")
                    for entry in sorted(entries)[:q]]
            bins["__fallback_observations__"]+=len(picked)
        else:
            picked=[]
            for phase in sorted(by_season):
                source=by_season[phase]
                k=seasonal[phase]
                picked.extend(
                    (entry,len(source),k,"NESTED_SEASON")
                    for entry in sorted(source)[:k]
                )
            bins["__nested_observations__"]+=len(picked)
        for (_,oid,obs),N,k,mode in picked:
            row=obs["row"]
            obs_pi=k/N
            photo_pi=obs_pi/obs["m"]
            if not (0<obs_pi<=1 and 0<photo_pi<=1):
                raise ValueError("Invalid positive photo/observation probability")
            writer.writerow({
                **{key:row[key] for key in (
                    "inat_taxon_id","genus","photo_id","observation_id",
                    "photo_license","latitude","longitude","bio12_mm"
                )},
                **{key:str(row.get(key) or "") for key in (
                    "photo_url_large","attribution","observer_id","observed_on",
                    "positional_accuracy_m"
                )},
                "region_stratum":group[0],"precip_stratum":group[1],
                "source_photos_in_observation":obs["m"],
                "available_observations_in_stratum":N,
                "sampled_observations_in_stratum":k,
                "observation_inclusion_probability":f"{obs_pi:.12g}",
                "photo_inclusion_probability":f"{photo_pi:.12g}",
                "photo_design_weight":f"{1/photo_pi:.12g}",
                "technical_pilot_flag":pilot,
            })
            selected+=1
            bins["__unknown_site__"]+=group[0]=="G_UNKNOWN"
            bins["__unknown_climate__"]+=group[1]=="P_UNKNOWN"
            bins[f"{group[0]}:{group[1]}"]+=1
    if selected!=min(CAP,len(observations)):
        raise AssertionError("Nested season design changed taxon-specific sample cap")
    return n_photos,len(observations),selected,pilot*selected,dict(bins)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fcp-source-root",type=Path,required=True)
    ap.add_argument("--complete-counts",type=Path,required=True)
    ap.add_argument("--source-admission",type=Path,required=True)
    ap.add_argument("--bio12-geotiff",type=Path,required=True)
    ap.add_argument("--bio12-sha256",required=True)
    ap.add_argument("--archived-site-calibration",type=Path,required=True)
    ap.add_argument("--archived-site-calibration-sha256",required=True)
    ap.add_argument("--out-private-photo-csv",type=Path,required=True)
    ap.add_argument("--out-aggregate-receipt",type=Path,required=True)
    a=ap.parse_args()
    climate=PinnedBio12(a.bio12_geotiff,a.bio12_sha256)
    calibration=calibrate(climate,a.archived_site_calibration,
                          a.archived_site_calibration_sha256)
    outcome=run_direct(
        a.fcp_source_root,a.complete_counts,a.source_admission,climate,
        a.out_private_photo_csv,a.out_aggregate_receipt,
        calibration=calibration,species_sampler=sample_species_season_v2,
        selection_protocol=STRATEGY
    )
    print(json.dumps({k:outcome[k] for k in (
        "status","source_species","source_photo_ids_authenticated",
        "selected_original_photo_ids","sampling_protocol"
    )},sort_keys=True))


if __name__=="__main__":
    main()
