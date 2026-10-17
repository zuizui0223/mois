#!/usr/bin/env python3
"""MOIS genuine dated photograph × original four-state hue temporal SUPPORT audit.

This does NOT calculate the prevalence of verified open white flowers.
Older 100,543 Mois photos have no frozen observed_on in the 247,032-ID
measured union and cannot be treated as dated. Only FCP's additional
146,489 photo IDs have unaltered source dates. Original photo ROI labels
are not botanical confirmation of anthesis, pigment chemistry, fitness or
weather effects.

Independent observation IDs (when present) are used for temporal support,
not pixel or photo count. Source photo coordinate accuracy AND photographer
IDs are absent from the measured union; fine 0.25-degree cells are merely
reported sampled coordinates, not verified true local populations.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import zipfile
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path

UNION_STATE="REAL_FROZEN_FOUR_HUE_PHOTO_UNION_NO_NEW_IMAGE_PIXELS"
UNION_IDS=247032
PHOTO_SHA="c344e34be5508964e0de092dfcc7df63997332e2288bda257c3f561e937e10aa"
SOURCE_SHA="b961b50ee6efe9ec08a5575a8eb1e34b1edf7b7fcb426226fb11ee901cfd3329"
COLOURS=("white","yellow_orange","red_pink","blue_purple")
SOURCE_YEAR=date(2026,9,27)
METHOD_DATE="source original observed_on exact day, not upload date"
MONTH_COLS=("year_month","source_photo_opportunities","source_automatic_four_state_ROIs",
            "source_unclassifiable_or_mixed_ROIs","source_automatic_white_ROIs",
            "source_automatic_chromatic_ROIs",
            "photographic_white_fraction_among_AUTOMATIC_ROIs_only",
            "botanical_open_petaloid_flower_confirmed","actual_month_climate_join_completed")
SUMMARY_COLS=("grid_width_degrees","years_in_scope","minimum_independent_source_observations_in_each_of_two_seasons",
              "matched_taxon_grid_years","matched_named_taxa","all_matched_season_photographic_colour_observations",
              "source_coordinate_accuracy_verified","source_photographer_independence_verified",
              "flower_open_organ_independently_validated")
MATCH_COLS=("inat_taxon_id","grid_lat_idx","grid_lon_idx","local_meteorological_season_year",
            "season_a","season_b","independent_observations_a","independent_observations_b",
            "automated_photo_white_labels_a","automated_photo_white_labels_b",
            "automatic_4state_observations_other_colours_a",
            "automatic_4state_observations_other_colours_b",
            "source_grid_degrees","source_geoprecision_verified","photo_open_flower_verified",
            "same_calendar_year_climate_available_both_seasons_through_2025")


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()


def season(year:int,month:int,lat:float)->tuple[str,int]:
    if abs(lat)<23.5:
        return "tropical_month_only",year
    north=("winter" if month in (12,1,2) else
           "spring" if month in (3,4,5) else
           "summer" if month in (6,7,8) else "autumn")
    if lat<0:
        north={"winter":"summer","spring":"autumn",
               "summer":"winter","autumn":"spring"}[north]
    return north,year+int(month==12)


def grid(lat:float,lon:float,step:float)->tuple[int,int]:
    if not (-90<=lat<=90 and -180<=lon<=180 and math.isfinite(lat) and math.isfinite(lon)):
        raise ValueError("Original measured source latitude/longitude outside Earth")
    # 180°E and 180°W refer to the same wrapped 360° longitude location.
    return min(int(180/step)-1,int(math.floor((lat+90)/step))),int(math.floor(((lon+180)%360)/step))


def check_inputs(archive:Path,400_source:Path|None=None,*,test_only=False)->tuple[dict,set[str]]:
    with zipfile.ZipFile(archive) as z:
        receipt=json.loads(z.read("mois_fcp_measured_photo_union_receipt.json"))
    if (receipt.get("status")!=UNION_STATE or
        receipt.get("combined_unique_photo_ids")!=(8 if test_only else UNION_IDS) or
        receipt.get("historical_global_genus_moisture_slopes_reestimated") is not False or
        receipt.get("cross_source_measurement_independence") is not False or
        receipt.get("new_fcp_photo_site_month_climate_join_complete") is not False):
        raise ValueError("Input photo metadata is not the exact frozen original MOIS+FCP colour union")
    if not test_only and receipt.get("photo_hue_union_csv_sha256")!=PHOTO_SHA:
        raise ValueError("Original 247032 combined photo hue source SHA changed")
    selected=set()
    if 400_source is not None:
        if not test_only and sha(400_source)!=SOURCE_SHA:
            raise ValueError("Source-only 400 original photo selection SHA invalid")
        with 400_source.open(newline="",encoding="utf-8") as f:
            q=csv.DictReader(f)
            for x in q:
                p=x["photo_id"]
                if p in selected or not p.isdecimal():
                    raise ValueError("Repeated or invalid original 400 photo ID")
                selected.add(p)
        if len(selected)!=(2 if test_only else 400):
            raise ValueError("Original 400-photo seasonal opportunity source changed")
    return receipt,selected


def audit(archive:Path,out:Path,*,400_source:Path|None=None,
          test_only=False)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Do not overwrite a SHA-frozen MOIS original photo temporal support audit")
    source_doc,selected=check_inputs(archive,400_source,test_only=test_only)
    base=Counter()
    monthly=defaultdict(Counter)
    observation={}
    selection_overlap=defaultdict(Counter)
    original_ids=set()
    fcp_taxa=set()
    with zipfile.ZipFile(archive) as z:
        member="mois_fcp_combined_original_247032_photo_hues.csv.gz"
        info=z.getinfo(member)
        with z.open(info) as stream:
            # Check actual uncompressed member byte SHA before looking at colour.
            raw=stream.read()
        if not test_only and hashlib.sha256(raw).hexdigest()!=PHOTO_SHA:
            raise ValueError("Genuine photo-level original measured union changed since authenticated CI")
        with gzip.open(io.BytesIO(raw),"rt",encoding="utf-8",newline="") as f:
            r=csv.DictReader(f)
            expected={"photo_id","inat_taxon_id","species","morph","old_mois_or_new_fcp",
                      "observation_id","observed_on","latitude","longitude"}
            if not expected.issubset(r.fieldnames or []):
                raise ValueError("Original measured photo source schema changed")
            for x in r:
                pid=x["photo_id"]
                if pid in original_ids or not pid.isascii() or not pid.isdecimal():
                    raise ValueError("Measured original photo ID duplicate or invalid")
                original_ids.add(pid)
                kind=x["old_mois_or_new_fcp"]
                morph=x["morph"]
                if kind not in ("historical","new_fcp") or morph not in (*COLOURS,"mixed_uncertain",""):
                    raise ValueError("Merged Mois/FCP original photo label not its fixed source")
                base[kind+"_photos"]+=1
                is_photo_colour=morph in COLOURS
                base[kind+"_fourstate_ROIs"]+=int(is_photo_colour)
                base[kind+"_unclassifiable_photo"]+=int(not is_photo_colour)
                if pid in selected:
                    selection_overlap[kind]["photo_ids"]+=1
                    selection_overlap[kind]["fourstate_labels"]+=int(is_photo_colour)
                    selection_overlap[kind]["white_labels"]+=int(morph=="white")
                if kind=="historical":
                    if x["observed_on"]:
                        base["historical_photos_with_historical_date_column"]+=1
                    continue
                # The already measured FCP cohort has original observation dates,
                # but the old shallow historical MOIS union does not.
                tid=x["inat_taxon_id"]
                day_text=x["observed_on"]
                try:
                    day=date.fromisoformat(day_text)
                except (TypeError,ValueError) as e:
                    raise ValueError("Measured FCP original photo has no valid frozen observation day") from e
                if not (date(1950,1,1)<=day<=SOURCE_YEAR and tid.isdecimal() and x["species"]):
                    raise ValueError("Source original observed-on beyond snapshot or source ID invalid")
                lat,lon=float(x["latitude"]),float(x["longitude"])
                if not (math.isfinite(lat) and math.isfinite(lon)):
                    raise ValueError("Original measured FCP photo has no reliable coordinate")
                grid(lat,lon,.25)
                ym=f"{day.year:04d}-{day.month:02d}"
                c=monthly[ym]
                c["photos"]+=1
                c["colour"]+=int(is_photo_colour)
                c["white"]+=int(morph=="white")
                c["chromatic"]+=int(is_photo_colour and morph!="white")
                fcp_taxa.add(tid)
                if is_photo_colour:
                    # Deduplicate images of the same photographed observation.
                    # Any conflicting photo hues from an observation are unknown
                    # instead of more independent flower observations.
                    obs=x["observation_id"]
                    if not obs:raise ValueError("Classified FCP photo missing independent source observation ID")
                    key=(tid,obs)
                    state=(lat,lon,day,morph)
                    if key in observation:
                        base["repeated_photo_from_same_observation"]+=1
                        old=observation[key]
                        if old[:3]!=(lat,lon,day) or old[3]!=morph:
                            base["conflicting_multi_photo_observation"]+=1
                            observation[key]=None
                    else:
                        observation[key]=state
    if base["historical_photos"]+base["new_fcp_photos"]!=len(original_ids):
        raise ValueError("Original measured source photo opportunities not conserved")
    if not test_only and (len(original_ids)!=UNION_IDS or base["historical_photos"]!=100543 or
        base["new_fcp_photos"]!=146489 or
        base["historical_fourstate_ROIs"]!=45953 or
        base["new_fcp_fourstate_ROIs"]!=74357 or
        len(fcp_taxa)!=1499 or
        base["historical_photos_with_historical_date_column"]!=0):
        raise ValueError("Actual frozen original 247032 source/photo-colour membership changed")
    observed_rows=[(tax,obs,*state) for (tax,obs),state in observation.items() if state is not None]
    # Freeze the eligible source-photo temporal support on independent
    # source observations and the local hemispheric calendar.
    grids=(.25,.5,1.0,5.0)
    samples={}
    for width in grids:
        site=defaultdict(lambda:defaultdict(Counter))
        for tax,obs,lat,lon,day,morph in observed_rows:
            season_name,season_year=season(day.year,day.month,lat)
            if season_name=="tropical_month_only":
                continue
            i,j=grid(lat,lon,width)
            site[(tax,i,j,season_year)][season_name]["obs"]+=1
            site[(tax,i,j,season_year)][season_name]["white"]+=int(morph=="white")
            site[(tax,i,j,season_year)][season_name]["nonwhite"]+=int(morph!="white")
        samples[width]=site
    summary=[]
    matching=[]
    for width,site in samples.items():
        for scope in ("all_source_years","climate_available_1958_to_2025"):
            filtered={k:v for k,v in site.items() if scope=="all_source_years" or
                      1958<=k[3]<=2025}
            for threshold in (1,2,3,5,10):
                qualified=[]
                for key,ss in filtered.items():
                    supported=[(s,v) for s,v in ss.items() if v["obs"]>=threshold]
                    if len(supported)>=2:
                        qualified.append((key,ss,supported))
                summary.append(dict(zip(SUMMARY_COLS,(
                    width,scope,threshold,len(qualified),
                    len({k[0] for k,_,_ in qualified}),
                    sum(sum(v["obs"] for v in ss.values()) for _,ss,_ in qualified),
                    False,False,False
                ))))
                if width==.25 and scope=="climate_available_1958_to_2025" and threshold==2:
                    for key,ss,supported in qualified:
                        tax,i,j,y=key
                        season_a,aa=supported[0]
                        season_b,bb=supported[1]
                        matching.append(dict(zip(MATCH_COLS,(
                            tax,i,j,y,season_a,season_b,aa["obs"],bb["obs"],
                            aa["white"],bb["white"],aa["nonwhite"],bb["nonwhite"],
                            width,False,False,True
                        ))))
    out.mkdir(parents=True,exist_ok=True)
    month_file=out/"mois_new_fcp_measured_photo_original_observed_month_ROI_opportunities.csv"
    with month_file.open("x",encoding="utf-8",newline="") as f:
        q=csv.DictWriter(f,fieldnames=MONTH_COLS,lineterminator="\n")
        q.writeheader()
        for ym,c in sorted(monthly.items()):
            n=c["colour"]
            q.writerow(dict(zip(MONTH_COLS,(
                ym,c["photos"],n,c["photos"]-n,c["white"],c["chromatic"],
                f"{c['white']/n:.8f}" if n else "",False,False
            ))))
    summary_file=out/"mois_measured_source_season_matched_observation_support_by_spatial_resolution.csv"
    with summary_file.open("x",encoding="utf-8",newline="") as f:
        q=csv.DictWriter(f,fieldnames=SUMMARY_COLS,lineterminator="\n")
        q.writeheader();q.writerows(summary)
    match_file=out/"mois_matched_0p25_original_observation_siteyear_two_season_ROI_support.csv"
    with match_file.open("x",encoding="utf-8",newline="") as f:
        q=csv.DictWriter(f,fieldnames=MATCH_COLS,lineterminator="\n")
        q.writeheader();q.writerows(sorted(matching,key=lambda r:(
            int(r["inat_taxon_id"]),int(r["local_meteorological_season_year"]),
            int(r["grid_lat_idx"]),int(r["grid_lon_idx"])))))
    horizons=Counter()
    for ym,c in monthly.items():
        year=int(ym[:4])
        k="pre_1958" if year<1958 else "actual_month_climate_to_2025" if year<=2025 else "2026_not_verified_climate"
        horizons[k+"_source_photos"]+=c["photos"]
        horizons[k+"_fourstate_ROIs"]+=c["colour"]
    original_overlap={key:dict(vals) for key,vals in selection_overlap.items()}
    result={
        "schema":"mois_genuine_measured_247032_photo_union_original_dated_ROI_season_support_v1",
        "status":"REAL_MEASURED_PHOTO_SEASONAL_SUPPORT_NOT_BOTANICAL_OPEN_FLOWER_SEASONAL_EFFECT",
        "source_union_archive_read_sha256":sha(archive),
        "original_photo_union_uncompressed_gzip_sha256":PHOTO_SHA,
        "original_20260927_400_photo_selection_SHA256":sha(400_source) if 400_source else None,
        "n_unique_original_photo_ids":len(original_ids),
        "source_photo_and_automatic_colour_funnel":dict(base),
        "n_source_taxa_new_fcp_with_any_original_photo":len(fcp_taxa),
        "n_classifiable_original_photos_with_source_observed_on":base["new_fcp_fourstate_ROIs"],
        "n_independent_original_colour_classified_observations":len(observed_rows),
        "source_monthly_climate_coverage_by_year":dict(horizons),
        "source_400_original_photo_ID_intersection":original_overlap,
        "minimum_independent_obs_each_season_thresholds":[1,2,3,5,10],
        "spatial_resolution_degrees":[.25,.5,1.0,5.0],
        "fine_grid_within_species_siteyear_climate_scope_min2_pair_count":len(matching),
        "source_photo_vs_true_flower_organ_anthesis_known":False,
        "source_photo_position_accuracy_not_provided_in_247k_union":True,
        "source_photographer_observer_disjointness_not_proven_in_union":True,
        "historical_mois_100543_observed_on_missing_in_frozen_union":True,
        "photo_white_rate_is_automatic_technical_conditional_on_roi_classifiability":True,
        "no_genetic_polymorphism_or_true_blooming_white_flower_rate_estimated":True,
        "no_monthly_weather_link_or_causal_eco_effect_performed":True,
        "dated_photo_month_ROI_opportunity_sha256":sha(month_file),
        "temporal_matching_sample_support_sha256":sha(summary_file),
        "fine_year_climate_matched_ROI_samples_sha256":sha(match_file),
    }
    (out/"mois_measured_original_247k_date_and_season_support_receipt.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-247032-measured-union-archive",type=Path,required=True)
    p.add_argument("--original-400-blind-photo-sample",type=Path)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    result=audit(a.original_247032_measured_union_archive,a.outdir,
                 a.original_400_blind_photo_sample)
    print("MOIS_REAL_247K_FOUR_HUE_SOURCE_SEASON_COMPARABILITY_NO_TRUE_BLOOM",
        json.dumps({k:result[k] for k in (
            "n_unique_original_photo_ids","n_classifiable_original_photos_with_source_observed_on",
            "n_independent_original_colour_classified_observations",
            "fine_grid_within_species_siteyear_climate_scope_min2_pair_count",
            "source_monthly_climate_coverage_by_year","source_400_original_photo_ID_intersection",
        )},sort_keys=True))


if __name__=="__main__":
    main()
