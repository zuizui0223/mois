#!/usr/bin/env python3
"""Strictly provenance-bound MOIS audit: real 2025 seasonal ROI missingness,
photographic spatial support and TerraClimate observation-month anomalies.

This is an audit of *unreviewed original photo classification*, NOT estimates
of blooming flower colours, selection, local adaptation, or genotype frequencies.
A season-by-climate coefficient cannot be identified in one site-year.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

FCP_ZIP_CSV="postblind/first_original_source_matched_season_64_photo_roi_QC.csv"
FCP_ZIP_RECEIPT="postblind/first_original_source_matched_season_64_roi_receipt.json"
FCP_LOCK_SHA="1b16cdabd2936567dee87b8344019f40f315600c6510dfc1f25a03835024b6a0"
TC_ROOT="mois_physical_2025_original_temporal_climate/"
TC_CSV=TC_ROOT+"authentic_2025_asclepias_terraclimate_physical_month_grid.csv"
TC_RECEIPT=TC_ROOT+"verified_2025_terraclimate_source_month_receipt.json"
MONTH_COLS=("year_month","local_season","original_photo_attempts",
            "roi_four_state_automatic","roi_same_white_yellow_organs_mixed",
            "roi_unusable_unknown","automatic_white_label_UNREVIEWED",
            "automatic_chromatic_label_UNREVIEWED",
            "original_climate_gridcells_observed",
            "botanically_confirmed_open_flower_photos")
SEASON_COLS=("season","original_source_photos","independent_source_observations",
             "independent_source_observers","original_photo_sites_4km",
             "photo_attempts_in_7_shared_4km_sites",
             "automatic_four_state_roi_all_sites",
             "mixed_white_yellow_roi_all_sites",
             "automatic_classifiable_minus_mixed_all_sites",
             "automated_white_UNREVIEWED_all_sites",
             "automated_chromatic_UNREVIEWED_all_sites",
             "automatic_four_state_roi_shared_sites",
             "mixed_white_yellow_roi_shared_sites",
             "automatic_classifiable_minus_mixed_shared_sites",
             "automated_white_UNREVIEWED_shared_sites",
             "automated_chromatic_UNREVIEWED_shared_sites",
             "automatic_roi_classifiability_fraction",
             "year2025_tmean_C_weighted","year2025_srad_W_m2_weighted",
             "year2025_ppt_mm_weighted","year2025_def_mm_weighted",
             "year2025_vpd_kPa_weighted",
             "2025_tmean_anomaly_C_same_month",
             "2025_srad_anomaly_W_m2_same_month",
             "2025_ppt_anomaly_mm_same_month",
             "2025_def_anomaly_mm_same_month",
             "2025_vpd_anomaly_kPa_same_month",
             "phenology_stage_UNREVIEWED")
PHOTO_COLS=("original_photo_id","original_observation_id","original_observed_on",
            "source_photo_season","original_independent_observer_id",
            "original_terraclimate_lat_idx","original_terraclimate_lon_idx",
            "same_4km_cell_has_both_photo_seasons",
            "original_roi_four_state_automatic","roi_mixed_white_yellow",
            "original_roi_unusable_unknown","photo_automatic_roi_white_UNREVIEWED",
            "photo_automatic_roi_chromatic_UNREVIEWED",
            "original_photo_climate_year_month","tmean_C","srad_W_m2",
            "ppt_mm","def_mm","vpd_kPa",
            "tmean_C_anomaly_same_month","srad_W_m2_anomaly_same_month",
            "ppt_mm_anomaly_same_month","def_mm_anomaly_same_month",
            "vpd_kPa_anomaly_same_month",
            "open_flower_petal_tepal_verified")


def sha(b:bytes)->str:
    return hashlib.sha256(b).hexdigest()


def csv_rows(b:bytes)->list[dict]:
    with io.StringIO(b.decode("utf-8")) as f:
        q=csv.DictReader(f)
        if not q.fieldnames:
            raise ValueError("Expected SHA-sealed original CSV schema")
        return list(q)


def csv_write(file:Path,fields:tuple[str,...],rows:list[dict])->None:
    with file.open("x",newline="",encoding="utf-8") as f:
        q=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        q.writeheader()
        q.writerows(rows)


def coords(row:dict)->tuple[int,int]:
    lat=float(row["latitude"])
    lon=float(row["longitude"])
    acc=float(row["positional_accuracy_m"])
    if (not all(math.isfinite(v) for v in (lat,lon,acc)) or
        not -90<=lat<=90 or not -180<=lon<=180 or not 0<acc<=5000):
        raise ValueError("Genuine original date/photo source geolocation precision missing")
    return min(4319,math.floor((90-lat)*24)),math.floor(((lon+180)%360)*24)


def build(fcp_zip:Path,climate_zip:Path,out:Path,*,_test_only=False)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Do not overwrite genuine original source/CF archive results")
    with zipfile.ZipFile(fcp_zip) as z:
        fp=z.read(FCP_ZIP_CSV)
        pr=json.loads(z.read(FCP_ZIP_RECEIPT))
    with zipfile.ZipFile(climate_zip) as z:
        cc=z.read(TC_CSV)
        cr=json.loads(z.read(TC_RECEIPT))
    if (pr.get("status")!="FIRST_ORIGINAL_SEASONAL_PIXELS_MEASURED_BOTANICAL_ANTHESIS_NOT_VALIDATED" or
        pr.get("postblind_join_sha256")!=sha(fp) or
        (not _test_only and sha(fp)!=FCP_LOCK_SHA) or
        pr.get("original_photos_measured")!=64 or
        cr.get("status")!="AUTHENTIC_2025_MONTHLY_CLIMATE_VALUES_NOT_WHITE_FLOWER_ASSOCIATION" or
        cr.get("decoded_month_grid_climate_sha256")!=sha(cc) or
        cr.get("input_original_2025_photo_climate_request_sha256")!="4a181971255218aa2645313aaecd3d16739c1dd00ca7a365d6e1b4e1d974024d" or
        cr.get("original_source_photo_grid_months")!=40 or
        cr.get("unique_original_source_grid_cells")!=21 or
        cr.get("genuine_original_monthly_photo_observations")!=64 or
        cr.get("no_botanically_confirmed_open_flower_colour_data") is not True or
        cr.get("no_white_fraction_vs_climate_regression_or_causal_inference") is not True):
        raise ValueError("Mois genuine original ROI + actual CF-climate input SHA/receipt mismatch")
    photo=csv_rows(fp)
    climate=csv_rows(cc)
    if len(photo)!=64 or len(climate)!=40:
        raise ValueError("Missing photographed observations or source climate cells")
    climate_values={}
    for c in climate:
        ym=c["year_month"]
        ii,jj=int(c["tc_lat_idx"]),int(c["tc_lon_idx"])
        k=(ym,ii,jj)
        if k in climate_values or not ym.startswith("2025-"):
            raise ValueError("Climate source duplicates or borrowed another year")
        if (c["flowering_organ_is_botanically_unreviewed"]!="True" or
            c["climate_is_actual_2025_not_a_historical_normal"]!="True"):
            raise ValueError("Climatological normals or fake flowering stage mixed with real year")
        climate_values[k]=c
    count=Counter()
    observer=set()
    observations=set()
    photo_ids=set()
    grids_by_season=defaultdict(set)
    photo_rows=[]
    source=[]
    for p in photo:
        pid=p["photo_id"]
        obs=p["observation_id"]
        user=p["observer_id"]
        if pid in photo_ids or obs in observations or user in observer:
            raise ValueError("The original photo wave is not independent source observation and observer data")
        photo_ids.add(pid)
        observations.add(obs)
        observer.add(user)
        season=p["local_season"]
        ym=p["observed_year_month"]
        if season not in ("summer","autumn") or ym not in {
            f"2025-{n:02d}" for n in range(6,12)}:
            raise ValueError("Not the expected source-2025 two-temporal-season input")
        if (season=="summer") != (int(ym[5:])<=8):
            raise ValueError("Not exact source date/calendar season")
        if p["observed_on"][:7]!=ym:
            raise ValueError("Original observed_on disagrees with sampled month")
        if p["species"]!="Asclepias syriaca" or p["inat_taxon_id"]!="47911":
            raise ValueError("Original source sampled species changed")
        if p["flower_stage_verified_by_independent_botanists"]!="False" or p["flower_roi_organ_human_review"]!="UNREVIEWED":
            raise ValueError("External review not completed; no false anthesis assumptions")
        ii,jj=coords(p)
        c=climate_values.get((ym,ii,jj))
        if c is None:
            raise ValueError("Original source photo has no exact 2025 site/month climate cell")
        count[(ym,ii,jj)]+=1
        grids_by_season[season].add((ii,jj))
        four=p["measurement_status"]=="classified_four_state_morph"
        if four and p["morph"] not in {"white","yellow_orange","red_pink","blue_purple"}:
            raise ValueError("Unknown original four-state ROI")
        mix=p["same_roi_white_yellow_organ_mixed"]=="True"
        if p["same_roi_white_yellow_organ_mixed"] not in ("True","False") or mix and not four:
            raise ValueError("Original ROI floral-organ mixing or missingness schema drift")
        source.append((p,c,ii,jj,four,mix))
    if count!={(k):int(c["n_photo_observations"]) for k,c in climate_values.items()}:
        raise ValueError("Photographs and physical climate grid-month denominators do not match")
    shared=grids_by_season["summer"] & grids_by_season["autumn"]
    if not _test_only and (len(shared)!=7 or
        len(grids_by_season["summer"])!=16 or len(grids_by_season["autumn"])!=12):
        raise ValueError("Unexpected original seasonal geographic support")
    months=defaultdict(list)
    seasons=defaultdict(list)
    for p,c,ii,jj,four,mix in source:
        season=p["local_season"]
        ym=p["observed_year_month"]
        usable=four and not mix
        months[(ym,season)].append((p,c,ii,jj,four,mix))
        seasons[season].append((p,c,ii,jj,four,mix))
        def cl(key):
            return c[key]
        photo_rows.append(dict(zip(PHOTO_COLS,(
            p["photo_id"],p["observation_id"],p["observed_on"],season,p["observer_id"],
            ii,jj,(ii,jj) in shared,four,mix,not four,
            usable and p["morph"]=="white",
            usable and p["morph"]!="white",
            ym,cl("tmean_C"),cl("srad_W_m2"),cl("ppt_mm"),cl("def_mm"),cl("vpd_kPa"),
            cl("tmean_C_anomaly_1991_2020_same_month"),
            cl("srad_W_m2_anomaly_1991_2020_same_month"),
            cl("ppt_mm_anomaly_1991_2020_same_month"),
            cl("def_mm_anomaly_1991_2020_same_month"),
            cl("vpd_kPa_anomaly_1991_2020_same_month"),False))))
    monthrows=[]
    for (ym,season),rs in sorted(months.items()):
        monthrows.append(dict(zip(MONTH_COLS,(
            ym,season,len(rs),
            sum(four for _,_,_,_,four,_ in rs),
            sum(mix for _,_,_,_,_,mix in rs),
            sum(not four for _,_,_,_,four,_ in rs),
            sum(four and not mix and p["morph"]=="white" for p,_,_,_,four,mix in rs),
            sum(four and not mix and p["morph"]!="white" for p,_,_,_,four,mix in rs),
            len({(ii,jj) for _,_,ii,jj,_,_ in rs}),0))))
    by_season=[]
    for season in ("summer","autumn"):
        rs=seasons[season]
        common=[t for t in rs if (t[2],t[3]) in shared]
        if len(rs)!=32:raise ValueError("Missing one frozen original 32-photo summer/autumn fold")
        def count_good(rows):
            return (sum(four for _,_,_,_,four,_ in rows),
                    sum(mix for _,_,_,_,_,mix in rows),
                    sum(four and not mix for _,_,_,_,four,mix in rows),
                    sum(four and not mix and p["morph"]=="white" for p,_,_,_,four,mix in rows),
                    sum(four and not mix and p["morph"]!="white" for p,_,_,_,four,mix in rows))
        def mean_cl(col):
            return sum(float(c[col]) for _,c,_,_,_,_ in rs)/len(rs)
        yes,mi,usable,white,chrom=count_good(rs)
        yc,mc,uc,wc,cc=count_good(common)
        by_season.append(dict(zip(SEASON_COLS,(
            season,len(rs),len({p["observation_id"] for p,_,_,_,_,_ in rs}),
            len({p["observer_id"] for p,_,_,_,_,_ in rs}),
            len(grids_by_season[season]),len(common),
            yes,mi,usable,white,chrom,yc,mc,uc,wc,cc,
            f"{yes/len(rs):.6f}",
            *(f"{mean_cl(col):.6f}" for col in (
                "tmean_C","srad_W_m2","ppt_mm","def_mm","vpd_kPa",
                "tmean_C_anomaly_1991_2020_same_month",
                "srad_W_m2_anomaly_1991_2020_same_month",
                "ppt_mm_anomaly_1991_2020_same_month",
                "def_mm_anomaly_1991_2020_same_month",
                "vpd_kPa_anomaly_1991_2020_same_month")),
            True))))
    out.mkdir(parents=True,exist_ok=True)
    pfile=out/"2025_Asclepias_original_64_automated_ROI_climate_photo_join_UNREVIEWED.csv"
    mfile=out/"2025_Asclepias_original_month_ROI_classifiability.csv"
    sfile=out/"2025_Asclepias_seasonal_shared_grid_support_climate_missingness.csv"
    csv_write(pfile,PHOTO_COLS,sorted(photo_rows,key=lambda p:int(p["original_photo_id"])))
    csv_write(mfile,MONTH_COLS,monthrows)
    csv_write(sfile,SEASON_COLS,by_season)
    result={
        "schema":"mois_original_2025_summer_autumn_within_grid_roi_missingness_v1",
        "status":"REAL_ORIGINAL_PHOTOGRAPH_SEASON_CLIMATE_JOIN_NO_BOTANICALLY_CONFIRMED_COLOUR",
        "source_fcp_2025_photo_roi_sha256":sha(fp),
        "source_MOIS_2025_terraclimate_climate_sha256":sha(cc),
        "photo_attempts":64,
        "independent_original_observations":len(observations),
        "independent_original_observers":len(observer),
        "same_4km_terraclimate_gridcells_both_summer_and_autumn":len(shared),
        "summer_gridcells":len(grids_by_season["summer"]),
        "autumn_gridcells":len(grids_by_season["autumn"]),
        "all_original_source_gridcells":len(grids_by_season["summer"]|grids_by_season["autumn"]),
        "photo_month_grid_cell_counts_reconciled_to_full_climate_source":True,
        "no_photo_replacement_or_climate_conditioned_photo_selection":True,
        "n_independently_reviewed_visible_open_flower_petals":0,
        "all_64_botanical_stages_remain_UNREVIEWED":True,
        "raw_ROI_labels_are_not_white_flower_rate":True,
        "spatially_matched_photo_sample_is_not_climate_randomization":True,
        "observed_month_and_season_perfectly_confound_weather_and_temporal_stage":True,
        "no_seasonal_white_colour_climate_regression":True,
        "photo_climate_join_csv_sha256":sha(pfile.read_bytes()),
        "month_classifiability_csv_sha256":sha(mfile.read_bytes()),
        "shared_grid_season_QC_csv_sha256":sha(sfile.read_bytes()),
        "by_season":{r["season"]:r for r in by_season},
    }
    (out/"verified_2025_shared_sites_original_ROI_climate_QC_receipt.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fcp-original-2025-64-ROI-zip",type=Path,required=True)
    p.add_argument("--mois-original-2025-terraclimate-zip",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    r=build(a.fcp_original_2025_64_ROI_zip,
            a.mois_original_2025_terraclimate_zip,a.outdir)
    print(json.dumps({k:v for k,v in r.items() if k!="by_season"},sort_keys=True))
    for s,c in r["by_season"].items():
        print("MOIS_2025_SOURCE_SITE_SEASON_"+s.upper(),json.dumps(c,sort_keys=True))


if __name__=="__main__":
    main()
