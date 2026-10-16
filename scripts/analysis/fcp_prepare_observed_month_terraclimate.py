#!/usr/bin/env python3
"""Build a source-photo-bound observed-month TerraClimate request roster.

This is NOT a climate API query, weather dataset, ecological result, or
assumption that every iNaturalist original photo contains a flower.

Only after the real original monthly photo-ID source and the frozen original
ROI-v4 classifier both pass their independent SHA/identity gates may original
observed_on and source coordinate be joined for a local 1/24-degree (~4 km)
TerraClimate time-series request. Never expose date/site to image classifier.

Output: each independent original source observation represented once,
including unresolved/ROI-failed ones, and the deduplicated climate month-cell
requests. A 2026 observation stays in the dated source, but is not represented
as having 2025 or climatology substituted for 2026 weather.
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

SOURCE_DATE="2026-09-27"
REAL_ROI="FIRST_GENUINE_MONTHLY_ORIGINAL_GEO_MATCHED_PIXELS_MEASURED_NOT_1000_SPECIES"
REAL_MONTH="DESCRIPTIVE_ORIGINAL_ROI_MONTHLY_PHENOTYPE_NOT_GLOBAL_COMPLETION"
OK="classified_four_state_morph"
MORPHS={"white","yellow_orange","red_pink","blue_purple"}
REQ_COLS=(
    "inat_taxon_id","species","observation_id","observer_id",
    "observed_on","observed_year_month","latitude","longitude",
    "positional_accuracy_m","tc_lat_idx","tc_lon_idx",
    "tc_lat_center","tc_lon_center","n_new_roi_photo_attempts",
    "n_four_state_classifiable_roi_photos","n_unknown_roi_photos",
    "n_same_roi_white_yellow_organ_mixed","phenotype_status",
    "source_position_precision_eligible_5km",
    "exact_observed_month_supported_by_climate_release",
    "climate_join_status",
)
CLIM_COLS=(
    "tc_lat_idx","tc_lon_idx","year_month","tc_lat_center",
    "tc_lon_center","n_independent_original_observations",
    "n_photo_roi_attempts","n_photo_supported_white_observations",
    "n_photo_supported_chromatic_observations",
    "n_original_observations_without_confident_flower_colour",
    "source_month_is_not_a_climate_measurement",
)


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(2*1024*1024),b""):
            h.update(part)
    return h.hexdigest()


def signed(source:Path, roi_doc:Path, month_doc:Path, *,
           test_only=False)->tuple[dict,dict]:
    roi=json.loads(roi_doc.read_text(encoding="utf-8"))
    qc=json.loads(month_doc.read_text(encoding="utf-8"))
    accepted_roi=("SYNTHETIC_ORIGINAL_GEOGRAPHIC_PAIR_TEST_NOT_BIOLOGY"
                  if test_only else REAL_ROI)
    accepted_month=("SYNTHETIC_TEST_NO_BIOLOGICAL_FINDING"
                    if test_only else REAL_MONTH)
    if (roi.get("status")!=accepted_roi
        or roi.get("original_geographic_source_join_terminal_sha256")!=sha(source)
        or qc.get("status")!=accepted_month
        or qc.get("new_roi_terminal_sha256")!=sha(source)
        or qc.get("source_snapshot_date")!=SOURCE_DATE
        or qc.get("new_roi_terminal_photos")!=roi.get("new_original_photo_pixels_processed")
        or qc.get("original_source_is_not_confirmed_flowers") is not True
        or qc.get("white_fraction_uses_independent_original_observations_not_photos") is not True
        or qc.get("calendar_date_was_never_used_to_select_or_classify_original_pixels") is not True
        or qc.get("all_93209351_photo_pixels_classified") is not False
        or roi.get("old_measured_colours_imported") is not False):
        raise ValueError("Original monthly source/pixel identities not independently authenticated")
    return roi,qc


def observed_month(raw:str)->str:
    raw=str(raw or "").strip()
    if len(raw)!=10 or raw[4]!="-" or raw[7]!="-":
        return ""
    try:
        day=date.fromisoformat(raw)
        return day.strftime("%Y-%m") if day<=date.fromisoformat(SOURCE_DATE) else ""
    except ValueError:
        return ""


def grid(lat_raw:str, lon_raw:str)->tuple[int,int,str,str]:
    try:
        lat=float(lat_raw)
        lon=float(lon_raw)
    except (ValueError,TypeError) as e:
        raise ValueError("Invalid original photo source site") from e
    if (not math.isfinite(lat) or not math.isfinite(lon)
        or not -90<=lat<=90 or not -180<=lon<=180):
        raise ValueError("Original photo geographic coordinate outside globe")
    # 1/24-degree WGS84 cell centres are +89.979167.. and -179.979167..
    # The 180E/180W antimeridian is a single wrapped grid index.
    ii=min(4319,int(math.floor((90-lat)*24)))
    jj=int(math.floor(((lon+180)%360)*24))
    return ii,jj,f"{90-(ii+0.5)/24:.8f}",f"{-180+(jj+0.5)/24:.8f}"


def precise(raw:str)->bool:
    try:
        number=float(raw)
        return math.isfinite(number) and 0<number<=5000
    except (ValueError,TypeError):
        return False


def prepare(source:Path,roi_receipt:Path,monthly_receipt:Path,out:Path,*,
            last_verified_climate_year:int=2025,_synthetic_test=False)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Never overwrite a source-signed TerraClimate request roster")
    if not 1958<=last_verified_climate_year<=date.fromisoformat(SOURCE_DATE).year:
        raise ValueError("Invalid independently verified climate availability bound")
    roi,month=signed(source,roi_receipt,monthly_receipt,test_only=_synthetic_test)
    groups={}
    seen_photos=set()
    with source.open(newline="",encoding="utf-8") as h:
        reader=csv.DictReader(h)
        required={
            "inat_taxon_id","species","photo_id","observation_id",
            "observer_id","observed_on","latitude","longitude",
            "positional_accuracy_m","measurement_status","morph",
            "flower_fraction_white","flower_fraction_yellow",
        }
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("Post-ROI original photo source missing date/position/mask status")
        for row in reader:
            if (not row["photo_id"].isascii() or not row["photo_id"].isdecimal()
                or row["photo_id"] in seen_photos):
                raise ValueError("Original measured source photo ID duplicate or invalid")
            seen_photos.add(row["photo_id"])
            if (not row["inat_taxon_id"].isdecimal() or
                not row["observation_id"] or not row["observer_id"] or
                not row["species"]):
                raise ValueError("Source observation/taxon identity invalid")
            ym=observed_month(row["observed_on"])
            if ("source_observation_month" in row and
                row["source_observation_month"] and
                row["source_observation_month"]!=ym):
                raise ValueError("Original source day disagrees with previously recovered observation month")
            ii,jj,clat,clon=grid(row["latitude"],row["longitude"])
            key=(row["inat_taxon_id"],row["observation_id"])
            loc=(row["inat_taxon_id"],row["species"],row["observation_id"],
                 row["observer_id"],row["observed_on"],
                 row["latitude"],row["longitude"],
                 row["positional_accuracy_m"],ii,jj,clat,clon)
            if key not in groups:
                groups[key]={"loc":loc,"morphs":set(),"failed":0,"good":0,
                             "mixed":0,"photos":0}
            g=groups[key]
            if g["loc"]!=loc:
                raise ValueError("Multiple photographs of one original observation disagree on date/site/observer")
            good=row["measurement_status"]==OK
            morph=row["morph"]
            if good!=(morph in MORPHS) or (not good and morph!="mixed_uncertain"):
                raise ValueError("Original ROI-v4 four-state status/flower mask inconsistent")
            g["photos"]+=1
            if not good:
                g["failed"]+=1
                continue
            try:
                fw=float(row["flower_fraction_white"])
                fy=float(row["flower_fraction_yellow"])
            except ValueError as e:
                raise ValueError("Original nine-colour ROI fraction invalid") from e
            if not (math.isfinite(fw) and math.isfinite(fy) and
                    0<=fw<=1 and 0<=fy<=1):
                raise ValueError("Original flower ROI colour fraction invalid")
            g["good"]+=1
            if fw>=.1 and fy>=.1:
                g["mixed"]+=1
            else:
                g["morphs"].add(morph)
    if len(seen_photos)!=roi["new_original_photo_pixels_processed"]:
        raise ValueError("Original frozen pixel terminal changed number of photo IDs")

    out.mkdir(parents=True,exist_ok=True)
    obs_file=out/"real_original_roi_observation_terraclimate_month_keys.csv"
    request_file=out/"exact_observed_month_terraclimate_grid_requests.csv"
    requests=defaultdict(lambda:{"obs":0,"photos":0,"white":0,"chromatic":0,"unresolved":0})
    excluded=defaultdict(int)
    colour=defaultdict(int)
    with obs_file.open("x",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=REQ_COLS,lineterminator="\n")
        writer.writeheader()
        for _,g in sorted(groups.items(),key=lambda v:(int(v[0][0]),v[0][1])):
            (tax,name,obs,observer,day,lat,lon,accuracy,
             ii,jj,clat,clon)=g["loc"]
            ym=observed_month(day)
            precise_geo=precise(accuracy)
            if g["mixed"]:
                phenotype="same_flower_white_yellow_organ_mix_unresolved"
            elif len(g["morphs"])>1:
                phenotype="discordant_photo_morphs_unresolved"
            elif not g["morphs"]:
                phenotype="unclassifiable_or_not_flower_unknown"
            else:
                phenotype=("white" if g["morphs"]=={"white"} else "chromatic")
            colour[phenotype]+=1
            if not ym:
                join="missing_or_invalid_original_observation_day"
            elif not precise_geo:
                join="source_position_uncertainty_above_5km_or_unknown"
            elif int(ym[:4])<1958:
                join="source_observation_year_before_supported_climate_history"
            elif int(ym[:4])>last_verified_climate_year:
                join="year_not_in_independently_verified_monthly_climate_release"
            else:
                join="original_observation_month_grid_request_eligible_not_fetched"
            excluded[join]+=1
            if join=="original_observation_month_grid_request_eligible_not_fetched":
                req=requests[(ii,jj,ym)]
                req["obs"]+=1
                req["photos"]+=g["photos"]
                req["white"]+=int(phenotype=="white")
                req["chromatic"]+=int(phenotype=="chromatic")
                req["unresolved"]+=int(phenotype not in ("white","chromatic"))
            writer.writerow(dict(zip(REQ_COLS,(
                tax,name,obs,observer,day,ym,lat,lon,accuracy,ii,jj,
                clat,clon,g["photos"],g["good"],g["failed"],g["mixed"],
                phenotype,precise_geo,int(ym[:4])<=last_verified_climate_year
                if ym else False,join
            ))))
    with request_file.open("x",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=CLIM_COLS,lineterminator="\n")
        writer.writeheader()
        for (ii,jj,ym),counts in sorted(requests.items()):
            writer.writerow(dict(zip(CLIM_COLS,(
                ii,jj,ym,f"{90-(ii+0.5)/24:.8f}",
                f"{-180+(jj+0.5)/24:.8f}",
                counts["obs"],counts["photos"],counts["white"],
                counts["chromatic"],counts["unresolved"],True,
            ))))
    if sum(excluded.values())!=len(groups) or sum(g["photos"] for g in groups.values())!=len(seen_photos):
        raise ValueError("No valid observation denominators for source-ROI climate roster")
    result={
        "schema":"fcp_original_roi_postblind_terraclimate_observation_month_key_v1",
        "status":("SYNTHETIC_OBSERVED_MONTH_CLIMATE_REQUEST_NOT_BIOLOGY" if _synthetic_test
                  else "REAL_ORIGINAL_ROI_CLIMATE_REQUEST_KEYS_NOT_CLIMATE_MEASURED"),
        "snapshot_date":SOURCE_DATE,
        "original_pixel_terminal_sha256":sha(source),
        "original_roi_receipt_sha256":sha(roi_receipt),
        "original_monthly_site_qc_receipt_sha256":sha(monthly_receipt),
        "original_independent_observations":len(groups),
        "real_original_photo_pixels":len(seen_photos),
        "observations_by_photo_colour_phenotype":dict(colour),
        "observations_by_climate_request_eligibility":dict(excluded),
        "n_exact_source_month_raster_cells_to_request":len(requests),
        "climate_variables_requested":["tmin","tmax","ppt","srad","def","vpd"],
        "climate_expected_physical_units":{
            "tmin":"degC","tmax":"degC","ppt":"mm_month_total",
            "srad":"W_m2_month_mean","def":"mm_month_total","vpd":"kPa_month_mean"},
        "source_intended":"TerraClimate v1.1 (ERA5-era), 1/24 degree",
        "last_independently_verified_monthly_climate_year":last_verified_climate_year,
        "source_based_date_and_position_join_happens_after_blind_image_roi":True,
        "no_2026_climate_values_fabricated_from_2025_or_1991_2020_normals":True,
        "unknown_roi_not_treated_as_a_nonflower":True,
        "per_observation_colour_not_genetically_verified":True,
        "no_climate_raster_download_or_temp_radiation_water_effect_estimated":True,
        "all_93209351_original_source_photos_measured":False,
        "independent_observation_grid_key_csv_sha256":sha(obs_file),
        "unique_terra_cell_month_requests_csv_sha256":sha(request_file),
    }
    (out/"postblind_terraclimate_request_receipt.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
        encoding="utf-8")
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-roi-photo-terminal",type=Path,required=True)
    p.add_argument("--original-roi-receipt",type=Path,required=True)
    p.add_argument("--month-site-qc-receipt",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--last-verified-climate-year",type=int,default=2025)
    args=p.parse_args(argv)
    print(json.dumps(prepare(args.original_roi_photo_terminal,
                             args.original_roi_receipt,
                             args.month_site_qc_receipt,args.outdir,
                             last_verified_climate_year=args.last_verified_climate_year),
                     sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
