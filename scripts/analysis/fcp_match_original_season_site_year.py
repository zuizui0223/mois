#!/usr/bin/env python3
"""Compare ORIGINAL ROI-supported FCP white proportions across local seasons.

This post-ROI comparison is outcome-blind to unmatched geographic cells:
a within-cell season comparison contributes only when the SAME taxon,
same 5-degree cell, hemisphere, and meteorological season-year have
independently classifiable original observations in both seasons.

Five-degree cells are broad and DO NOT establish shared individual plants,
native populations, controlled sites, genetic polymorphism, adaptation, or
true weather anomalies. The effect is descriptive, without inferential
p-values or misleadingly precise photo-level confidence intervals.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import date
from pathlib import Path

SOURCE_STATUS = "DESCRIPTIVE_ORIGINAL_ROI_MONTHLY_PHENOTYPE_NOT_GLOBAL_COMPLETION"
TEST_STATUS = "SYNTHETIC_TEST_NO_BIOLOGICAL_FINDING"
SEASONS = ("local_winter", "local_spring", "local_summer", "local_autumn")
PAIRS = tuple((a,b) for i,a in enumerate(SEASONS) for b in SEASONS[i+1:])
INCOLS = (
    "inat_taxon_id","species","year_month","calendar_month","hemisphere",
    "season","lat5","lon5","dated_original_photo_opportunities",
    "new_ROI_photos_with_known_date","classifiable_flower_ROIs",
    "same_ROI_white_yellow_organ_mixed","unclassifiable_ROIs",
    "usable_flower_photo_labels","independent_usable_observations",
    "white_observations","nonwhite_observations",
    "contradictory_multiphoto_observations","distinct_usable_observers",
    "white_rate_per_observation_descriptive",
)
COUNT_COLS = (
    "dated_original_photo_opportunities","new_ROI_photos_with_known_date",
    "classifiable_flower_ROIs","same_ROI_white_yellow_organ_mixed",
    "unclassifiable_ROIs","usable_flower_photo_labels",
    "independent_usable_observations","white_observations",
    "nonwhite_observations","contradictory_multiphoto_observations",
)
ROWCOLS = (
    "inat_taxon_id","species","season_year","hemisphere","lat5","lon5",
    "local_season","n_calendar_months","source_date_photo_opportunities",
    "new_original_ROI_attempts","classifiable_flower_ROIs",
    "unclassifiable_ROIs","white_yellow_organ_mixed_ROIs",
    "independent_flower_observations","white_observations",
    "nonwhite_observations","white_fraction_observation_weighted",
    "roi_classifiability_among_attempted",
)
CONTRAST_COLS = (
    "inat_taxon_id","species","season_a","season_b","comparison_status",
    "matched_5degree_cell_seasonyear_strata",
    "matched_calendar_months_a","matched_calendar_months_b",
    "matched_independent_observations_a","matched_independent_observations_b",
    "balanced_min_observation_weight","balanced_white_fraction_a",
    "balanced_white_fraction_b","balanced_white_fraction_difference_a_minus_b",
    "unmatched_species_white_fraction_a","unmatched_species_white_fraction_b",
    "unmatched_species_white_fraction_difference_a_minus_b",
    "matched_roi_classifiability_a","matched_roi_classifiability_b",
    "geographic_cell_km_resolution_is_coarse",
    "source_observers_are_not_independent_across_seasons_verified",
    "causal_or_genetic_inference",
)


def sha(p:Path)->str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for part in iter(lambda:f.read(1024*1024),b""):
            h.update(part)
    return h.hexdigest()


def read_receipt(qc:Path,receipt:Path,*,test_only:bool=False)->dict:
    doc=json.loads(receipt.read_text(encoding="utf-8"))
    accepted=TEST_STATUS if test_only else SOURCE_STATUS
    if (doc.get("status")!=accepted or doc.get("source_snapshot_date")!="2026-09-27"
        or doc.get("monthly_site_QC_csv_sha256")!=sha(qc)
        or doc.get("original_source_is_not_confirmed_flowers") is not True
        or doc.get("roi_failed_is_not_called_nonflowering") is not True
        or doc.get("white_fraction_uses_independent_original_observations_not_photos") is not True
        or doc.get("calendar_date_was_never_used_to_select_or_classify_original_pixels") is not True
        or doc.get("all_93209351_photo_pixels_classified") is not False):
        raise ValueError("Original dated, observation-weighted, new ROI season QC not attested")
    return doc


def local_season_year(ym:str,hem:str,season:str)->int:
    if (len(ym)!=7 or ym[4]!="-"):
        raise ValueError("Not a source observed year-month")
    year,month=map(int,ym.split("-"))
    date(year,month,1)
    if hem not in ("north","south","tropical"):
        raise ValueError("Invalid hemisphere in source-date QC")
    if hem=="tropical":
        if season!="tropical_month_only":
            raise ValueError("Tropical locality assigned an unverified four-season label")
    else:
        northern=("local_winter" if month in (12,1,2) else
                  "local_spring" if month in (3,4,5) else
                  "local_summer" if month in (6,7,8) else "local_autumn")
        southern={"local_winter":"local_summer","local_spring":"local_autumn",
                  "local_summer":"local_winter","local_autumn":"local_spring"}
        expected=northern if hem=="north" else southern[northern]
        if season!=expected:
            raise ValueError("Hemisphere-month-local-season mismatch")
    # December is part of next calendar year's meteorological season (DJF)
    # in both hemispheres (north winter; south summer).
    return year+int(month==12) if hem!="tropical" else year


def _to_int(row:dict,key:str)->int:
    val=row[key]
    if not val.isascii() or not val.isdecimal():
        raise ValueError("Original month-site QC contains invalid count: "+key)
    return int(val)


def run(qc:Path,receipt:Path,out:Path,*,_test_only:bool=False)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Do not overwrite a sealed season comparison")
    source_doc=read_receipt(qc,receipt,test_only=_test_only)
    cells={}
    species={}
    temporal_rows=0
    with qc.open(newline="",encoding="utf-8") as h:
        reader=csv.DictReader(h)
        if tuple(reader.fieldnames or ())!=INCOLS:
            raise ValueError("Original dated source QC format changed")
        seen=set()
        for raw in reader:
            tid=raw["inat_taxon_id"]
            if not tid.isascii() or not tid.isdecimal() or not raw["species"]:
                raise ValueError("Missing original taxon/source species identity")
            name=raw["species"]
            if tid in species and species[tid]!=name:
                raise ValueError("One original taxon ID has inconsistent species names")
            species[tid]=name
            ym=raw["year_month"]
            hem=raw["hemisphere"]
            season=raw["season"]
            sy=local_season_year(ym,hem,season)
            if _to_int(raw,"calendar_month")!=int(ym[-2:]):
                raise ValueError("Observation date month changed in QC")
            lat5=_to_int(raw,"lat5")
            lon5=_to_int(raw,"lon5")
            if not (0<=lat5<=35 and 0<=lon5<=71):
                raise ValueError("5-degree original source site cell invalid")
            key=(tid,sy,hem,lat5,lon5,season)
            fingerprint=(tid,ym,hem,lat5,lon5)
            if fingerprint in seen:
                raise ValueError("Duplicate original observation month/cell outcome")
            seen.add(fingerprint)
            d={c:_to_int(raw,c) for c in COUNT_COLS}
            if (d["new_ROI_photos_with_known_date"] >
                d["dated_original_photo_opportunities"]
                or d["classifiable_flower_ROIs"] + d["unclassifiable_ROIs"] !=
                   d["new_ROI_photos_with_known_date"]
                or d["usable_flower_photo_labels"] +
                   d["same_ROI_white_yellow_organ_mixed"] !=
                   d["classifiable_flower_ROIs"]
                or d["white_observations"]+d["nonwhite_observations"] !=
                   d["independent_usable_observations"]
                or d["independent_usable_observations"]>
                   d["usable_flower_photo_labels"]
                or _to_int(raw,"distinct_usable_observers") >
                   d["independent_usable_observations"]):
                raise ValueError("Original flower, ROI and observation denominators conflict")
            n=d["independent_usable_observations"]
            expected=d["white_observations"]/n if n else None
            rate=raw["white_rate_per_observation_descriptive"].strip()
            if expected is None:
                if rate:
                    raise ValueError("White fraction computed from no flower observations")
            elif not rate or not math.isfinite(float(rate)) or abs(float(rate)-expected)>1e-7:
                raise ValueError("White fraction is not observation-weighted")
            block=cells.setdefault(key,{
                "species":name,"months":set(),"counts":defaultdict(int)})
            if block["species"]!=name:
                raise ValueError("Season block changed source species")
            block["months"].add(ym)
            for c in COUNT_COLS:
                block["counts"][c]+=d[c]
            temporal_rows+=1

    out.mkdir(parents=True,exist_ok=True)
    season_file=out/"observed_month_to_local_season_cell_year.csv"
    with season_file.open("x",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=ROWCOLS,lineterminator="\n")
        writer.writeheader()
        for (tax,year,hem,lat,lon,season),info in sorted(cells.items()):
            c=info["counts"]
            n=c["independent_usable_observations"]
            attempts=c["new_ROI_photos_with_known_date"]
            writer.writerow(dict(zip(ROWCOLS,(
                tax,info["species"],year,hem,lat,lon,season,
                len(info["months"]),c["dated_original_photo_opportunities"],
                attempts,c["classifiable_flower_ROIs"],
                c["unclassifiable_ROIs"],
                c["same_ROI_white_yellow_organ_mixed"],
                n,c["white_observations"],c["nonwhite_observations"],
                f"{c['white_observations']/n:.8f}" if n else "",
                f"{c['classifiable_flower_ROIs']/attempts:.8f}" if attempts else "",
            ))))

    # Counts pooled across ALL photographed geocells are only a comparator
    # indicating how much unbalanced geography can change the observed trend.
    totals=defaultdict(lambda:defaultdict(int))
    matched=defaultdict(dict)
    for key,info in cells.items():
        tax,year,hem,lat,lon,s=key
        if s not in SEASONS:
            continue # tropical months retained in season_file, NOT re-labelled
        n=info["counts"]["independent_usable_observations"]
        if n:
            totals[(tax,s)]["white"]+=info["counts"]["white_observations"]
            totals[(tax,s)]["obs"]+=n
            matched[(tax,year,hem,lat,lon)][s]=info

    pairs_file=out/"within_species_matched_coarse_cell_season_contrasts.csv"
    n_supported=0
    n_insufficient=0
    with pairs_file.open("x",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=CONTRAST_COLS,lineterminator="\n")
        writer.writeheader()
        for tax,name in sorted(species.items(),key=lambda v:int(v[0])):
            for a,b in PAIRS:
                strata=[]
                for key,blocks in matched.items():
                    if key[0]!=tax or a not in blocks or b not in blocks:
                        continue
                    strata.append((blocks[a],blocks[b]))
                weight=0
                weighted_a=weighted_b=0.0
                obs_a=obs_b=months_a=months_b=0
                attempted_a=attempted_b=good_a=good_b=0
                for A,B in strata:
                    x,y=A["counts"],B["counts"]
                    na,nb=(x["independent_usable_observations"],
                           y["independent_usable_observations"])
                    w=min(na,nb)
                    weight+=w
                    weighted_a+=w*x["white_observations"]/na
                    weighted_b+=w*y["white_observations"]/nb
                    obs_a+=na
                    obs_b+=nb
                    months_a+=len(A["months"])
                    months_b+=len(B["months"])
                    attempted_a+=x["new_ROI_photos_with_known_date"]
                    attempted_b+=y["new_ROI_photos_with_known_date"]
                    good_a+=x["classifiable_flower_ROIs"]
                    good_b+=y["classifiable_flower_ROIs"]
                naive_a=totals[(tax,a)]
                naive_b=totals[(tax,b)]
                pnaive_a=(naive_a["white"]/naive_a["obs"] if naive_a["obs"] else None)
                pnaive_b=(naive_b["white"]/naive_b["obs"] if naive_b["obs"] else None)
                out_a=weighted_a/weight if weight else None
                out_b=weighted_b/weight if weight else None
                n_supported+=bool(strata)
                n_insufficient+=not bool(strata)
                def fmt(value):
                    return f"{value:.8f}" if value is not None else ""
                writer.writerow(dict(zip(CONTRAST_COLS,(
                    tax,name,a,b,("MATCHED_BROAD_CELL_YEAR_DESCRIPTIVE_ONLY"
                                  if strata else "NO_MATCHED_CELL_YEAR_NO_TEMPORAL_IDENTIFICATION"),
                    len(strata),months_a,months_b,obs_a,obs_b,weight,
                    fmt(out_a),fmt(out_b),
                    fmt(out_a-out_b if strata else None),
                    fmt(pnaive_a),fmt(pnaive_b),
                    fmt(pnaive_a-pnaive_b if pnaive_a is not None and
                        pnaive_b is not None else None),
                    fmt(good_a/attempted_a if attempted_a else None),
                    fmt(good_b/attempted_b if attempted_b else None),
                    True,False,False,
                ))))

    report={
        "schema":"fcp_observation_season_year_geocell_overlap_v1",
        "status":("SYNTHETIC_ONLY_NOT_BIOLOGY" if _test_only else
                  "MATCHED_BROAD_CELL_YEAR_DESCRIPTIVE_NO_CAUSAL_EFFECT"),
        "source_snapshot_date":"2026-09-27",
        "input_original_month_qc_sha256":sha(qc),
        "input_original_month_qc_receipt_sha256":sha(receipt),
        "season_year_cell_csv_sha256":sha(season_file),
        "matched_season_contrasts_csv_sha256":sha(pairs_file),
        "taxa_with_at_least_one_valid_original_source_month":len(species),
        "n_input_original_month_cells":temporal_rows,
        "n_merged_season_year_cell_rows":len(cells),
        "n_species_season_pairs_with_matched_coarse_cell_year_support":n_supported,
        "n_species_season_pairs_with_no_matched_support":n_insufficient,
        "counted_5degree_cell_is_not_a_controlled_local_site":True,
        "no_observer_independence_across_seasons_inferable_from_aggregates":True,
        "no_climate_raster_join_or_seasonal_weather_anomaly_estimated":True,
        "no_uncertainty_or_causal_significance_estimated_from_clustered_photos":True,
        "no_genetic_FCP_or_pollinator_selection_inference":True,
        "input_original_roi_not_full_93m_pixel_census":True,
    }
    (out/"matched_season_overlap_receipt.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--month-site-qc",type=Path,required=True)
    p.add_argument("--month-site-receipt",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args(argv)
    print(json.dumps(run(a.month_site_qc,a.month_site_receipt,a.outdir),
                     sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
