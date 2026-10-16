#!/usr/bin/env python3
"""Predeclared season-by-flower-colour diagnostics from sealed Mois ROI results.

Northern and southern local spring/summer/autumn/winter are mirrored.
Tropical observations use calendar month/quarter ONLY: rainfall-defined
wet/dry seasons require independently pinned MONTHLY precipitation rasters
and must never be inferred from annual BIO12 or a hemisphere label.

Distinguishes seasonal floral-community TURNOVER from white-colour differences
among the same nominal species, including stricter species x original
geographic/precipitation-stratum pairs. Nothing identifies individual plants,
within-genotype plasticity, population flowering rates or climatic selection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from datetime import date
from itertools import combinations
from pathlib import Path

try:
    from scripts.acquisition.mois_stratified_photo_selection_20261016 import classify_site, sha256_file
    from scripts.analysis.mois_phenology_observer_bias_audit_20261016 import load_overlay
except ModuleNotFoundError:
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from scripts.acquisition.mois_stratified_photo_selection_20261016 import classify_site, sha256_file
    from scripts.analysis.mois_phenology_observer_bias_audit_20261016 import load_overlay

SCHEMA="mois_photo_season_opportunity_v1"
STATUS="COMPLETE_SYNTHETIC_OR_SOURCE_BOUND_SEASONAL_DESCRIPTIVE_AUDIT_NO_13ENV_FIT"
RECEIPT_SCHEMA="mois_full_selected_photo_roi_terminal_collection_v1"
RECEIPT_STATUS="COMPLETE_MOIS_SAMPLED_PHOTO_TRAITS_CONDITIONAL_CLASSIFIABILITY_NO_13ENV_FIT"
SOURCE_SNAPSHOT=date(2026,9,27)
REQUIRED={
    "inat_taxon_id","genus","measurement_id","observation_id","observer_id",
    "latitude","longitude","bio12_mm","region_stratum","precip_stratum",
    "observed_on","photo_design_weight","observation_inclusion_probability",
    "measurement_status","morph"
}
HUES={"white","yellow_orange","red_pink","blue_purple"}
LOCAL=("SPRING","SUMMER","AUTUMN","WINTER")
PAIRS=tuple(combinations(LOCAL,2))
DATE_PREFIX=re.compile(r"^(\d{4}-\d{2}-\d{2})(?:$|T[0-9]{2}:[0-9]{2}.*$)")
PRIVATE_FIELDS=[
    "measurement_id","inat_taxon_id","genus","observation_id","observer_id",
    "observed_on","year","month","domain","local_season","annotation_status",
    "region_stratum","precip_stratum","latitude","longitude","bio12_mm",
    "sample_observation_weight","sample_photo_weight","morph","measurement_status"
]

def observation_date(raw:str)->date|None:
    if not raw: return None
    found=DATE_PREFIX.fullmatch(raw)
    if not found:return None
    try:
        value=date.fromisoformat(found.group(1))
    except ValueError:
        return None
    if value>SOURCE_SNAPSHOT:
        raise ValueError("Observation date exceeds immutable 2026-09-27 photo-source snapshot")
    return value

def local_season(lat_raw:str, observed:str)->tuple[str,str,str,str]:
    when=observation_date(observed)
    if not lat_raw:
        return ("UNKNOWN_GEOGRAPHY","UNKNOWN_SEASON",
                when.strftime("%m") if when else "UNKNOWN_MONTH",
                str(when.year) if when else "UNKNOWN_YEAR")
    lat=float(lat_raw)
    if not math.isfinite(lat) or not -90<=lat<=90:
        raise ValueError("Original photographed-site latitude invalid")
    domain=("TROPICAL" if abs(lat)<23.5
            else "NORTH_EXTRATROPICAL" if lat>=23.5
            else "SOUTH_EXTRATROPICAL")
    if when is None:
        return domain,"UNKNOWN_SEASON","UNKNOWN_MONTH","UNKNOWN_YEAR"
    month=when.month
    if domain=="TROPICAL":
        # Q1..Q4 is a calendar grouping only, not a wet/dry label.
        season=f"CALENDAR_Q{(month-1)//3+1}"
    else:
        north={
            12:"WINTER",1:"WINTER",2:"WINTER",
            3:"SPRING",4:"SPRING",5:"SPRING",
            6:"SUMMER",7:"SUMMER",8:"SUMMER",
            9:"AUTUMN",10:"AUTUMN",11:"AUTUMN"}
        northern=north[month]
        season=northern if lat>=23.5 else {
            "WINTER":"SUMMER","SPRING":"AUTUMN",
            "SUMMER":"WINTER","AUTUMN":"SPRING"}[northern]
    return domain,season,f"{month:02}",str(when.year)

def wmean(num:float,den:float)->float|None:
    return num/den if den else None

def add(stats:dict,row:dict,obs_weight:float,photo_weight:float)->None:
    valid=row["measurement_status"]=="classified_four_state_morph"
    white=valid and row["morph"]=="white"
    stats["selected"]+=1
    stats["classifiable"]+=int(valid)
    stats["white"]+=int(white)
    stats["obs_total"]+=obs_weight
    stats["obs_class"]+=obs_weight*valid
    stats["obs_white"]+=obs_weight*white
    stats["photo_total"]+=photo_weight
    stats["photo_class"]+=photo_weight*valid
    stats["photo_white"]+=photo_weight*white

def newstat()->dict:
    return {x:0 for x in ("selected","classifiable","white",
                          "obs_total","obs_class","obs_white",
                          "photo_total","photo_class","photo_white")}

def describe(s:dict)->dict:
    return {
        "n_selected_observations":s["selected"],
        "n_classifiable_flower_photos":s["classifiable"],
        "n_white_among_classifiable":s["white"],
        "observation_weighted_classifiability":wmean(s["obs_class"],s["obs_total"]),
        "photo_weighted_classifiability":wmean(s["photo_class"],s["photo_total"]),
        "observation_weighted_white_given_classifiable":wmean(s["obs_white"],s["obs_class"]),
        "photo_weighted_white_given_classifiable":wmean(s["photo_white"],s["photo_class"]),
    }

def matched_season_comparisons(
    rows_by_key:dict, *, min_classifiable:int
)->list[dict]:
    """Equal-unit white contrasts ONLY among taxa/cells observed in both seasons."""
    if min_classifiable<1:raise ValueError("At least one classified observation per compared season required")
    units=defaultdict(dict)
    for key, s in rows_by_key.items():
        if s["classifiable"]>=min_classifiable and s["obs_class"]>0:
            units[key[:-1]][key[-1]]=s["obs_white"]/s["obs_class"]
    comparisons=[]
    groups=defaultdict(list)
    for identity, season_values in units.items():
        if identity[-1]=="TROPICAL": # no comparable spring/summer labels in tropics
            continue
        for first,second in PAIRS:
            if first in season_values and second in season_values:
                # Season 2 minus season 1, same taxon or taxon/site cell.
                groups[(identity[-1],first,second)].append(
                    season_values[second]-season_values[first])
    for (domain,first,second),differences in sorted(groups.items()):
        comparisons.append({
            "domain":domain,
            "first_season":first,"second_season":second,
            "n_matched_units":len(differences),
            "equal_unit_mean_white_change_second_minus_first":
                sum(differences)/len(differences),
            "n_positive":sum(x>0 for x in differences),
            "n_negative":sum(x<0 for x in differences),
            "n_zero":sum(x==0 for x in differences),
        })
    return comparisons

def audit(
    joined_csv:Path, source_receipt:Path,
    private_rows:Path, private_summary:Path, public_summary:Path,
    *, annotation_overlay:Path|None=None, annotation_sha256:str|None=None,
    min_classifiable_per_season:int=3, min_public_cell:int=10,
    synthetic:bool=False,
)->dict:
    if min_classifiable_per_season<1 or min_public_cell<1:
        raise ValueError("Seasonal matched support and public suppression thresholds must be positive")
    source=json.loads(source_receipt.read_text(encoding="utf-8"))
    if (source.get("schema")!=RECEIPT_SCHEMA
        or source.get("status")!=RECEIPT_STATUS
        or source.get("private_joined_photo_csv_sha256")!=sha256_file(joined_csv)
        or source.get("synthetic_only_test") is not synthetic
        or source.get("full_FCP_93m_photos_measured") is not False
        or source.get("all13_photo_site_predictors_joined") is not False):
        raise ValueError("Only a complete hash-bound Mois sample, never partial synthetic ecology")
    if (annotation_overlay is None)!=(annotation_sha256 is None):
        raise ValueError("Annotation overlay and its independent SHA256 must be supplied together")
    annotation,ann_provenance=({},None)
    if annotation_overlay is not None:
        annotation,ann_provenance=load_overlay(annotation_overlay,annotation_sha256)
    for path in (private_rows,private_summary,public_summary):
        if path.exists() or path.with_name(path.name+".partial").exists():
            raise ValueError("Refuse to overwrite sealed or interrupted seasonal evidence")
        path.parent.mkdir(parents=True,exist_ok=True)
    temporary=private_rows.with_name(private_rows.name+".partial")
    season=defaultdict(newstat)
    season_annotation=defaultdict(newstat)
    taxon_season=defaultdict(newstat)
    cell_season=defaultdict(newstat)
    month_counts=Counter()
    annotation_counts=Counter()
    observer_season=defaultdict(Counter)
    source_obs=set()
    all_observers=set()
    selected=classified=white=missing_dates=0
    try:
        with joined_csv.open(newline="",encoding="utf-8") as source_f, \
             temporary.open("w",newline="",encoding="utf-8") as output:
            reader=csv.DictReader(source_f)
            if not REQUIRED.issubset(reader.fieldnames or ()):
                raise ValueError("Season audit needs source photo coordinates, observed date and sealed ROI")
            writer=csv.DictWriter(output,fieldnames=PRIVATE_FIELDS)
            writer.writeheader()
            for row in reader:
                tid=row["inat_taxon_id"]
                obs=row["observation_id"]
                if not tid.isdecimal() or not row["genus"] or not obs or obs in source_obs:
                    raise ValueError("Duplicate/missing source observation or nominal taxon")
                source_obs.add(obs)
                if (annotation and obs not in annotation):
                    raise ValueError("Phenology overlay lacks a selected observation")
                ann=annotation.get(obs,"NOT_VERIFIED")
                region,precip=classify_site(row)
                if (row["region_stratum"],row["precip_stratum"])!=(region,precip):
                    raise ValueError("Observed photo geographic/BIO12 sampling stratum changed")
                domain,year_season,month,year=local_season(
                    row["latitude"].strip(),row["observed_on"].strip())
                if domain=="UNKNOWN_GEOGRAPHY" and row["longitude"].strip():
                    raise ValueError("Half-missing coordinates")
                try:
                    pi=float(row["observation_inclusion_probability"])
                    pweight=float(row["photo_design_weight"])
                except (TypeError,ValueError) as exc:
                    raise ValueError("Missing original observation/photo sampling weights") from exc
                if not (math.isfinite(pi) and 0<pi<=1 and math.isfinite(pweight)
                        and pweight>=1):
                    raise ValueError("Impossible original photo or independent-observation selection probability")
                obs_w=1/pi
                state=row["measurement_status"]
                morph=row["morph"]
                if (state=="classified_four_state_morph")!=(morph in HUES):
                    raise ValueError("Source photographic classifiable/morph consistency failed")
                if state!="classified_four_state_morph" and morph!="mixed_uncertain":
                    raise ValueError("Unclassified photo is not automatically coloured or white")
                attr=row["observer_id"].strip() or "UNKNOWN_OBSERVER"
                all_observers.add(attr)
                observer_season[(domain,year_season)][attr]+=1
                stage={"measurement_id":row["measurement_id"],
                       "inat_taxon_id":tid,"genus":row["genus"],
                       "observation_id":obs,"observer_id":attr,
                       "observed_on":row["observed_on"],"year":year,
                       "month":month,"domain":domain,
                       "local_season":year_season,"annotation_status":ann,
                       "region_stratum":region,"precip_stratum":precip,
                       "latitude":row["latitude"],"longitude":row["longitude"],
                       "bio12_mm":row["bio12_mm"],
                       "sample_observation_weight":f"{obs_w:.12g}",
                       "sample_photo_weight":f"{pweight:.12g}",
                       "morph":morph,"measurement_status":state}
                writer.writerow(stage)
                add(season[(domain,year_season)],row,obs_w,pweight)
                add(season_annotation[(domain,year_season,ann)],row,obs_w,pweight)
                add(taxon_season[(tid,row["genus"],domain,year_season)],row,obs_w,pweight)
                add(cell_season[(tid,row["genus"],region,precip,domain,year_season)],
                    row,obs_w,pweight)
                annotation_counts[ann]+=1
                month_counts[(domain,month)]+=1
                selected+=1
                classified+=int(state=="classified_four_state_morph")
                white+=int(morph=="white")
                missing_dates+=int(month=="UNKNOWN_MONTH")
        if annotation and set(annotation)!=source_obs:
            raise ValueError("Annotation overlay contains extra observations outside the selected sample")
        if selected!=source["original_photo_terminals_verified"] or classified!=source["original_colour_four_state_classifiable"]:
            raise ValueError("Photographic seasonal counts differ from original sealed ROI terminal source")
        season_records=[]
        for (domain,quarter),s in sorted(season.items()):
            matching=[
                item for (tid,gen,dom,sea),item in taxon_season.items()
                if (dom,sea)==(domain,quarter) and item["classifiable"]>0
            ]
            observers=observer_season[(domain,quarter)]
            flowering=season_annotation[(domain,quarter,"FLOWERS")]
            no_annot=season_annotation[(domain,quarter,"NO_ANNOTATION")]
            unverified=season_annotation[(domain,quarter,"NOT_VERIFIED")]
            season_records.append({
                "domain":domain,"local_season":quarter,
                **describe(s),
                "confirmed_flowers_annotated_observations":flowering["selected"],
                "confirmed_flowers_annotated_roi_classifiable":flowering["classifiable"],
                "confirmed_flowers_annotated_white_given_classifiable_observation_weighted":
                    wmean(flowering["obs_white"],flowering["obs_class"]),
                "checked_but_unannotated_observations":no_annot["selected"],
                "annotation_not_verified_observations":unverified["selected"],
                "n_taxa_with_classifiable_colour":len(matching),
                "n_taxa_with_at_least_3_classifiable_photos":
                    sum(k["classifiable"]>=3 for k in matching),
                "n_taxa_with_at_least_5_classifiable_photos":
                    sum(k["classifiable"]>=5 for k in matching),
                "n_distinct_observers":len(observers),
                "top_observer_photo_share":max(observers.values())/s["selected"],
            })
        # Same species but different observational seasons; this does not
        # itself identify within-individual colour change.
        species_pairs=matched_season_comparisons(taxon_season,
                        min_classifiable=min_classifiable_per_season)
        # Strict species x fixed geographic + original BIO12 stratum:
        # a much better support check for apparent time effects.
        cell_pairs=matched_season_comparisons(cell_season,
                        min_classifiable=min_classifiable_per_season)
        # Drop taxon and coarse-location identifiers from public output.
        def public_cells(records):
            result=[]
            for r in records:
                if r["n_selected_observations"]<min_public_cell:
                    result.append({"domain":r["domain"],"local_season":r["local_season"],
                                   "privacy_suppressed_count_below":min_public_cell})
                else:
                    item=r.copy()
                    if r["confirmed_flowers_annotated_observations"]<min_public_cell:
                        item["confirmed_flowers_annotated_observations"]=None
                        item["confirmed_flowers_annotated_roi_classifiable"]=None
                        item["confirmed_flowers_annotated_white_given_classifiable_observation_weighted"]=None
                        item["confirmed_flowers_subcell_privacy_suppressed"]=True
                    result.append(item)
            return result
        def public_pairs(records):
            return [
                r if r["n_matched_units"]>=min_public_cell
                  else {"domain":r["domain"],"first_season":r["first_season"],
                        "second_season":r["second_season"],
                        "privacy_suppressed_n_matched_units_below":min_public_cell}
                for r in records
            ]
        report={
            "schema":SCHEMA,"status":STATUS,
            "source_sample_roi_receipt_sha256":sha256_file(source_receipt),
            "source_sealed_roi_join_sha256":sha256_file(joined_csv),
            "selected_source_observations":selected,
            "source_classifiable_photographs":classified,
            "source_classified_white_photographs":white,
            "n_selected_with_unknown_observation_date":missing_dates,
            "n_source_unique_observers":len(all_observers),
            "annotation_overlay_attestation":ann_provenance,
            "annotation_status_selected_counts":dict(sorted(annotation_counts.items())),
            "n_matched_species_season_pairs_before_public_suppression":
                sum(p["n_matched_units"] for p in species_pairs),
            "n_stricter_matched_species_region_precip_season_pairs_before_public_suppression":
                sum(p["n_matched_units"] for p in cell_pairs),
            "seasonal_observation_photo_white_conditioned":public_cells(season_records),
            "same_species_season_matched":public_pairs(species_pairs),
            "same_species_same_geo_precip_cell_season_matched":public_pairs(cell_pairs),
            "public_minimum_cell_size":min_public_cell,
            "matched_min_classifiable_per_species_season":min_classifiable_per_season,
            "tropical_season_method":"calendar_month_or_quarter_only; NEVER wet/dry without pinned monthly rainfall climatology",
            "northern_southern_calendar_months_aligned_to_local_seasons":True,
            "old_genus_13_environment_estimates_updated":False,
            "all_13_photo_site_environment_predictors_available_here":False,
            "causal_seasonal_phenotype_or_wild_flora_incidence_identified":False,
            "original_image_colour_results_measured_in_this_audit":False,
            "synthetic_test_only":bool(synthetic),
            "note":"Season pairs are matched nominal species or coarsely matched cells, not repeated plant/genotype tests; observational source seasonality may reflect flowering, tourism and observer selection."
        }
        private={
            "schema":"mois_private_photo_seasonal_trait_support_v1",
            "source_sample_receipt_sha256":sha256_file(source_receipt),
            "source_join_sha256":sha256_file(joined_csv),
            "private_derived_photo_rows_sha256":sha256_file(temporary),
            "full_unsuppressed_seasonal_support":season_records,
            "matched_nominal_species_season_pairs":species_pairs,
            "matched_species_region_precip_cell_season_pairs":cell_pairs,
            "selected_obs_by_domain_month":[
                {"domain":d,"month":m,"n_selected":n}
                for (d,m),n in sorted(month_counts.items())],
            "contains_photo_observation_taxon_ids_and_coordinates":True,
            "sampling_bias_corrections_are_only_conditional_on_source":True,
        }
        tmp_priv=private_summary.with_name(private_summary.name+".partial")
        tmp_public=public_summary.with_name(public_summary.name+".partial")
        if tmp_priv.exists() or tmp_public.exists():
            raise ValueError("Interrupted seasonal report files must be inspected first")
        tmp_priv.write_text(json.dumps(private,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        tmp_public.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        os.replace(temporary,private_rows)
        os.replace(tmp_priv,private_summary)
        os.replace(tmp_public,public_summary)
        return report
    except Exception:
        # No failed partial analysis should be advertised as a complete
        # worldwide photographed season result.
        for candidate in (temporary,
                          private_summary.with_name(private_summary.name+".partial"),
                          public_summary.with_name(public_summary.name+".partial")):
            if candidate.exists():
                candidate.unlink()
        raise

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--private-joined-photo-colours",required=True,type=Path)
    p.add_argument("--sample-collection-receipt",required=True,type=Path)
    p.add_argument("--out-private-season-photo-table",required=True,type=Path)
    p.add_argument("--out-private-season-report",required=True,type=Path)
    p.add_argument("--out-public-aggregate",required=True,type=Path)
    p.add_argument("--annotation-overlay",type=Path)
    p.add_argument("--annotation-overlay-sha256")
    p.add_argument("--min-classified-per-season",type=int,default=3)
    p.add_argument("--min-public-cell",type=int,default=10)
    args=p.parse_args()
    r=audit(
        args.private_joined_photo_colours,args.sample_collection_receipt,
        args.out_private_season_photo_table,args.out_private_season_report,
        args.out_public_aggregate,annotation_overlay=args.annotation_overlay,
        annotation_sha256=args.annotation_overlay_sha256,
        min_classifiable_per_season=args.min_classified_per_season,
        min_public_cell=args.min_public_cell
    )
    print(json.dumps({k:r[k] for k in (
        "status","selected_source_observations","source_classifiable_photographs",
        "n_matched_species_season_pairs_before_public_suppression",
        "n_stricter_matched_species_region_precip_season_pairs_before_public_suppression"
    )},sort_keys=True))

if __name__=="__main__":
    main()
