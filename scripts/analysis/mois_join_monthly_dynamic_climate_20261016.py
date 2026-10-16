#!/usr/bin/env python3
"""Exact original-site WorldClim 2.1 monthly dynamic T/S/P/vapour + VPD.

Use the SAME longitude/latitude pixel for all source months and variables.
Annual means are climatological means over ALL 12 months, not over highly
biased iNaturalist observation dates. Derived VPD uses saturation pressure at
mean monthly air temperature minus mean monthly vapour pressure; this is a
climate-normal APPROXIMATION, not an observed instantaneous daytime VPD.
Source belongs to the 1970--2000 WorldClim normal period, NOT photo year.

Only the 48-layer original source-admitted archive may enter production;
unclassified iNaturalist photo opportunities are never transformed to white.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts.acquisition.mois_stratified_photo_selection_20261016 import sha256_file
from scripts.analysis.mois_seasonal_flower_colour_opportunity_20261016 import observation_date
from scripts.analysis.mois_monthly_rainfall_flower_colour_20261016 import rainfall_label
from scripts.acquisition.mois_pin_worldclim_dynamic_monthly_climate_20261016 import (
    SCHEMA,VARIABLES,STATUS
)

REAL_MANIFEST_SHA="9d0cf37a64534542d7788debef28c13f106a15c3bf7bf91be0653f3ce106baae"
MAX_GRID_PIXELS=3_000_000
SOURCE_COLUMNS={
    "inat_taxon_id","genus","observation_id","measurement_id",
    "latitude","longitude","observed_on","morph","measurement_status"
}
COLOUR={"white","yellow_orange","red_pink","blue_purple"}
OUT_FIELDS=[
    "original_month", "climate_grid_row","climate_grid_col",
    "monthly_tavg_c","monthly_srad_kj_m2_day","monthly_prec_mm",
    "monthly_vapr_kpa","monthly_vpd_kpa",
    "delta_tavg_c","delta_srad_kj_m2_day","delta_prec_mm",
    "delta_vapr_kpa","delta_vpd_kpa",
    "site_mean_tavg_c","site_mean_srad_kj_m2_day","site_mean_prec_mm",
    "site_mean_vapr_kpa","site_mean_vpd_kpa",
    "source_site_annual_prec_mm","monthly_rainfall_phase",
    "vpd_raw_negative_clipped_at_observation_month","source_climate_join_status",
]


def saturation_pressure_kpa(tavg_c):
    """FAO-style saturation vapour pressure at MEAN MONTHLY air temperature."""
    t=np.asarray(tavg_c,dtype=np.float64)
    if np.any(~np.isfinite(t)) or np.any(t<=-237.3):
        raise ValueError("Invalid historical climatological Celsius temperature")
    return .6108*np.exp(17.27*t/(t+237.3))


class DynamicNormals:
    def __init__(self, root:Path, expected_manifest_sha:str, *, synthetic=False):
        manifest=root/"dynamic_monthly_manifest.json"
        if sha256_file(manifest)!=expected_manifest_sha:
            raise ValueError("Month-by-month 48-layer climate source manifest SHA256 mismatch")
        if not synthetic and expected_manifest_sha!=REAL_MANIFEST_SHA:
            raise ValueError("Original measured WorldClim 48-month source SHA changed")
        doc=json.loads(manifest.read_text(encoding="utf8"))
        if (doc.get("schema")!=SCHEMA or doc.get("status")!=STATUS
            or doc.get("synthetic_test_only") is not synthetic
            or doc.get("n_original_monthly_rasters")!=48
            or set(doc.get("months",{}))!=set(VARIABLES)
            or doc.get("climate_is_1970_2000_normal_not_photo_year_weather") is not True
            or doc.get("original_image_photo_ids_opened") is not False):
            raise ValueError("Original 48 climate-only month-source provenance invalid")
        self.arrays={}
        self.month_sha={}
        common=None
        for variable,rule in VARIABLES.items():
            months=doc["months"][variable]
            if set(months)!={f"{n:02}" for n in range(1,13)}:
                raise ValueError("A climate variable is missing original source month")
            band=[]
            for month in range(1,13):
                m=f"{month:02}"
                info=months[m]
                if info["original_variable_units"]!=rule["units"]:
                    raise ValueError("Original monthly climate source unit mismatch")
                path=root/info["path"]
                if sha256_file(path)!=info["sha256"]:
                    raise ValueError("Original TIFF source photo-climate byte changed")
                with rasterio.open(path) as tif:
                    transform=tif.transform
                    desc=(tif.width,tif.height,tuple(transform)[:6],
                          tif.crs.to_string() if tif.crs else None)
                    if (tif.count!=1 or tif.crs is None or tif.crs.to_epsg()!=4326
                        or transform.a<=0 or transform.e>=0
                        or abs(transform.b)>1e-12 or abs(transform.d)>1e-12
                        or tif.scales[0]!=1 or tif.offsets[0]!=0):
                        raise ValueError("Original source GeoTIFF uncalibrated or grid changed")
                    if desc[0]*desc[1]>MAX_GRID_PIXELS:
                        raise ValueError("Unexpected resampled/non-10m climate source")
                    if common is None:common=desc
                    elif desc!=common:
                        raise ValueError("T/S/P/vapour 48 original months do not share grid")
                    raw=tif.read(1,masked=True)
                    band.append(np.ma.filled(raw.astype(np.float32),np.nan))
                self.month_sha[variable+"_"+m]=info["sha256"]
            self.arrays[variable]=np.stack(band,axis=0)
        self.transform=rasterio.Affine(*common[2])
        self.height=common[1]
        self.width=common[0]
        self.manifest_sha=expected_manifest_sha
        self.cache={}
        self.synthetic=synthetic

    def cell(self, lat,lon):
        if lat=="" or lon=="" or lat is None or lon is None:
            if bool(lat)!=bool(lon):
                raise ValueError("Half-missing original photo-site geographic coordinates")
            return None
        la,lo=float(lat),float(lon)
        if not math.isfinite(la) or not math.isfinite(lo) or not -90<=la<=90 or not -180<=lo<=180:
            raise ValueError("Original observation location outside WGS84")
        x,y=(~self.transform)*(lo,la)
        row,col=math.floor(y),math.floor(x)
        if not 0<=row<self.height or not 0<=col<self.width:
            return None
        return row,col

    def profile(self, cell):
        if cell is None:return None
        if cell not in self.cache:
            row,col=cell
            p={}
            for key,b in self.arrays.items():
                values=np.asarray(b[:,row,col],dtype=np.float64)
                if np.any(~np.isfinite(values)) or np.any(values<VARIABLES[key]["min"]) \
                   or np.any(values>VARIABLES[key]["max"]):
                    self.cache[cell]=None
                    return None
                p[key]=values
            sat=saturation_pressure_kpa(p["tavg"])
            raw=sat-p["vapr"]
            p["vpd"]=np.maximum(0.,raw)
            p["vpd_months_clipped"]=int(np.sum(raw<0))
            self.cache[cell]=p
        return self.cache[cell]

    def site_month(self,lat,lon,observed_on)->dict:
        when=observation_date(observed_on)
        month=when.month if when else None
        cell=self.cell(lat,lon)
        p=self.profile(cell)
        if p is None:
            status="CLIMATE_UNAVAILABLE"
        elif month is None:
            status="OBSERVATION_DATE_UNKNOWN"
        else:
            status="VALID"
        out={name:"" for name in OUT_FIELDS}
        out["original_month"]=month if month is not None else ""
        out["climate_grid_row"]=cell[0] if cell is not None else ""
        out["climate_grid_col"]=cell[1] if cell is not None else ""
        out["source_climate_join_status"]=status
        if p is None:return out
        for source,col,unit in (
            ("tavg","tavg","c"),("srad","srad","kj_m2_day"),
            ("prec","prec","mm"),("vapr","vapr","kpa"),("vpd","vpd","kpa")
        ):
            avg=float(np.mean(p[source]))
            out["site_mean_"+col+"_"+unit]=avg
            if month is not None:
                val=float(p[source][month-1])
                out["monthly_"+col+"_"+unit]=val
                out["delta_"+col+"_"+unit]=val-avg
        out["source_site_annual_prec_mm"]=float(np.sum(p["prec"]))
        if month is not None:
            out["monthly_rainfall_phase"]=rainfall_label(p["prec"].tolist(),month)["label"]
            sat=saturation_pressure_kpa(p["tavg"][month-1])
            out["vpd_raw_negative_clipped_at_observation_month"]=bool(
                float(sat)<float(p["vapr"][month-1]))
        else:out["monthly_rainfall_phase"]="UNKNOWN_MONTH"
        return out


def join_original_sample(season_photo_csv:Path,season_private_receipt:Path,
                         climatology:DynamicNormals,out_private_csv:Path,
                         out_public_report:Path, *,
                         synthetic=False)->dict:
    """Strict private original photo table; public metadata-only support totals."""
    if out_private_csv.exists() or out_public_report.exists():
        raise ValueError("Source photo-climate join must never overwrite results")
    receipt=json.loads(season_private_receipt.read_text(encoding="utf8"))
    if (receipt.get("schema")!="mois_private_photo_seasonal_trait_support_v1"
        or receipt.get("private_derived_photo_rows_sha256")!=sha256_file(season_photo_csv)):
        raise ValueError("Not original source photo-season rows with SHA-sealed provenance")
    if climatology.synthetic is not synthetic:
        raise ValueError("Synthetic source cannot become an original data result")
    out_private_csv.parent.mkdir(parents=True,exist_ok=True)
    out_public_report.parent.mkdir(parents=True,exist_ok=True)
    staging=out_private_csv.with_name(out_private_csv.name+".INCOMPLETE")
    n=ok=classified=white=without_date=0
    seen=set()
    by_phase={}
    with season_photo_csv.open(newline="",encoding="utf8") as f, \
         staging.open("w",newline="",encoding="utf8") as target:
        reader=csv.DictReader(f)
        if not SOURCE_COLUMNS.issubset(reader.fieldnames or []):
            raise ValueError("Exact original observation source identity and true photo-ROI labels required")
        cols=reader.fieldnames+OUT_FIELDS
        writer=csv.DictWriter(target,fieldnames=cols)
        writer.writeheader()
        for row in reader:
            id_=row["measurement_id"]
            if not id_ or id_ in seen:
                raise ValueError("Original selected photograph identity missing or repeated")
            seen.add(id_)
            if (row["measurement_status"]=="classified_four_state_morph")!=(
                row["morph"] in COLOUR):
                raise ValueError("Original flower colour classification mismatch")
            result=climatology.site_month(row["latitude"],row["longitude"],row["observed_on"])
            writer.writerow({**row,**result})
            n+=1
            if result["source_climate_join_status"]=="VALID":
                ok+=1
                phase=result["monthly_rainfall_phase"]
                if phase not in by_phase:
                    by_phase[phase]={"n":0,"classifiable":0}
                by_phase[phase]["n"]+=1
                by_phase[phase]["classifiable"]+=int(
                    row["measurement_status"]=="classified_four_state_morph")
            if row["measurement_status"]=="classified_four_state_morph":
                classified+=1
                white+=int(row["morph"]=="white")
            without_date+=int(result["original_month"]=="")
    if not n or n<classified or classified<white:
        raise ValueError("No complete valid original photo source")
    result={
        "schema":"mois_original_photo_all_monthly_dynamic_climate_site_join_v1",
        "status":"COMPLETE_EXACT_ORIGINAL_SOURCE_PHOTO_MONTHLY_T_S_P_VAPR_JOIN_NOT_13ENV_REGRESSION",
        "original_sha_sealed_private_photo_csv_sha256":sha256_file(season_photo_csv),
        "original_sha_sealed_private_support_receipt_sha256":sha256_file(season_private_receipt),
        "original_48_raster_manifest_sha256":climatology.manifest_sha,
        "new_private_photo_monthly_climate_csv_sha256":sha256_file(staging),
        "selected_original_photographs":n,
        "source_four_colour_classifiable_photos":classified,
        "source_classifiable_white_photos":white,
        "observations_month_unknown":without_date,
        "n_successfully_joined_monthly_T_S_P_vapour":ok,
        "month_climate_wet_dry_support":by_phase,
        "VPD_monthly_approximate_from_mean_temperature_and_vapour_pressure":True,
        "seasonal_photo_white_environment_regression_fitted":False,
        "original_full13_photo_site_env_reconciled":False,
        "synthetic_test_only":bool(synthetic)
    }
    staging.replace(out_private_csv)
    out_public_report.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--verified-worldclim48-root",type=Path,required=True)
    p.add_argument("--verified-worldclim48-sha256",required=True)
    p.add_argument("--private-original-sample-season-photo",type=Path,required=True)
    p.add_argument("--private-season-receipt",type=Path,required=True)
    p.add_argument("--out-private-site-months",type=Path,required=True)
    p.add_argument("--out-public-support",type=Path,required=True)
    a=p.parse_args()
    weather=DynamicNormals(a.verified_worldclim48_root,a.verified_worldclim48_sha256)
    result=join_original_sample(a.private_original_sample_season_photo,
                 a.private_season_receipt,weather,
                 a.out_private_site_months,a.out_public_support)
    print(json.dumps({k:result[k] for k in (
        "status","selected_original_photographs",
        "n_successfully_joined_monthly_T_S_P_vapour")},sort_keys=True))


if __name__=="__main__":main()
