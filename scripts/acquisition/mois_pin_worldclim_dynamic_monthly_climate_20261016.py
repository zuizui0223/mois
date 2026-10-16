#!/usr/bin/env python3
"""Admit complete ORIGINAL WorldClim 2.1 10-minute monthly T/S/P/vapour rasters.

Four provider ZIPs tavg/srad/prec/vapr (1970--2000 climate normals) are
byte-hashed BEFORE any TIFF is used. All 48 GeoTIFFs must be present,
north-up, EPSG:4326 and on ONE identical global grid. No plant photos,
annotations or flower-classification data enter this source admission.

Expected original units, from WorldClim documentation:
tavg °C; srad kJ m^-2 day^-1; prec mm/month; vapr kPa.
No arbitrary source-value division by ten is allowed; the original raster
scales/offsets and plausible non-ocean source measurements are checked.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import re
import sys
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts.acquisition.mois_stratified_photo_selection_20261016 import sha256_file

SCHEMA="mois_worldclim21_original_monthly_tavg_srad_prec_vapr_10m_48bands_v1"
STATUS="VERIFIED_48_REAL_MONTHLY_GEOSPATIAL_CLIMATE_NORMALS_NO_PLANT_PHOTOS"
URL_ROOT="https://geodata.ucdavis.edu/climate/worldclim/2_1/base/"
CLIMATOLOGY="WorldClim 2.1 1970-2000 monthly 10m"
VARIABLES={
    "tavg":{"units":"degC","min":-95.,"max":65.},
    "srad":{"units":"kJ m-2 day-1","min":0.,"max":50000.},
    "prec":{"units":"mm/month","min":0.,"max":10000.},
    "vapr":{"units":"kPa","min":0.,"max":15.}
}
PREVIOUS_PREC_SHA="1090f578f6402672a2b6438e252156f6cf867a25002a8ad6ffcbfa19e516ddc6"
MAX_MEMBERS=13
MAX_FILE_BYTES=250_000_000
MAX_ARCHIVE_UNCOMPRESSED=1_500_000_000
HEADER="wc2.1_10m_"

def validate_originals(zips:dict[str,Path],known_digests:dict[str,str],
                       outroot:Path,*,synthetic:bool=False)->dict:
    stage=outroot.with_name(outroot.name+".INCOMPLETE")
    if outroot.exists() or stage.exists():
        raise ValueError("Cannot overwrite completed or interrupted source-climate admission")
    if set(zips)!=set(VARIABLES) or set(known_digests)!=set(VARIABLES):
        raise ValueError("Require all four WorldClim original monthly climate archives")
    import rasterio
    try:
        bounds_ref=None
        grid_ref=None
        bands={}
        archives={}
        readmes={}
        band_count=0
        stage.mkdir(parents=True)
        for kind in ("tavg","srad","prec","vapr"):
            zip_path=zips[kind]
            expected=known_digests[kind]
            if not synthetic and zip_path.name!=HEADER+kind+".zip":
                raise ValueError("Official archive basename must identify the original climate variable")
            if len(expected)!=64 or expected!=sha256_file(zip_path):
                raise ValueError(f"Provider {kind} original archive SHA256 mismatch")
            if not synthetic and kind=="prec" and expected!=PREVIOUS_PREC_SHA:
                raise ValueError("Original previously admitted precipitation archive changed")
            archives[kind]=expected
            with zipfile.ZipFile(zip_path) as z:
                members={}
                total=0
                readme=None
                for f in z.infolist():
                    if f.is_dir():continue
                    if f.filename=="readme.txt":
                        if readme is not None or f.file_size>65536 or (f.flag_bits&1):
                            raise ValueError("Duplicated, encrypted or large provider readme")
                        readme=hashlib.sha256(z.read(f)).hexdigest()
                        continue
                    pat=rf"(?:wc2\.1_10m_{kind}/)?(wc2\.1_10m_{kind}_(\d{{1,2}})\.tif)"
                    match=re.fullmatch(pat,f.filename)
                    if (match is None or f.flag_bits&1
                        or f.file_size<=0 or f.file_size>MAX_FILE_BYTES
                        or f.compress_size<=0):
                        raise ValueError("Unexpected, unsafe or non-monthly source GeoTIFF member")
                    month=int(match.group(2))
                    if month not in range(1,13) or month in members:
                        raise ValueError("Duplicate/non-calendar source monthly climate file")
                    members[month]=f
                    total+=f.file_size
                    if total>MAX_ARCHIVE_UNCOMPRESSED:
                        raise ValueError("Uncompressed source raster exceeds safety limit")
                if len(members)!=12 or len(z.infolist())>MAX_MEMBERS or z.testzip() is not None:
                    raise ValueError("Original monthly climate archive missing/malformed/extra members")
                readmes[kind]=readme
                entries={}
                d=stage/kind
                d.mkdir()
                for month in range(1,13):
                    member=members[month]
                    name=Path(member.filename).name
                    dst=d/name
                    with z.open(member) as src,dst.open("xb") as target:
                        for buf in iter(lambda:src.read(4*1024*1024),b""):
                            target.write(buf)
                    digest=sha256_file(dst)
                    with rasterio.open(dst) as tif:
                        t=tif.transform
                        grid=(tif.width,tif.height,tuple(t)[:6],tif.crs.to_string() if tif.crs else None)
                        if (tif.count!=1 or tif.crs is None or tif.crs.to_epsg()!=4326
                            or t.a<=0 or t.e>=0 or abs(t.b)>1e-12 or abs(t.d)>1e-12
                            or tif.scales[0]!=1 or tif.offsets[0]!=0
                            or (not synthetic and (
                                tif.width!=2160 or tif.height!=1080
                                or abs(t.a-1/6)>1e-8 or abs(t.e+1/6)>1e-8
                                or abs(t.c+180)>1e-6 or abs(t.f-90)>1e-6))):
                            raise ValueError("Original monthly raster units/grid/projection outside frozen 10m contract")
                        if grid_ref is None:grid_ref=grid
                        elif grid!=grid_ref:
                            raise ValueError("T/S/P/V original raster grids not identical")
                        value_samples=[]
                        # Read actual source grid at fixed sites for cheap
                        # additional plausibility, NOT for biological analysis.
                        for lat,lon in ((35.7,139.7),(-33.9,151.2),(0.0,28.0)):
                            x,y=(~t)*(lon,lat)
                            rr,cc=math.floor(y),math.floor(x)
                            if 0<=rr<tif.height and 0<=cc<tif.width:
                                v=tif.read(1,window=((rr,rr+1),(cc,cc+1)),masked=True)[0,0]
                                if not bool(getattr(v,"mask",False)):
                                    val=float(v)
                                    if not (math.isfinite(val)
                                      and VARIABLES[kind]["min"]<=val<=VARIABLES[kind]["max"]):
                                        raise ValueError("Official monthly variable source unit/value plausibility failed")
                                    value_samples.append(val)
                    entries[f"{month:02}"]={
                        "path":str(Path(kind)/name),"sha256":digest,
                        "original_variable_units":VARIABLES[kind]["units"]
                    }
                    band_count+=1
                bands[kind]=entries
        if band_count!=48:raise ValueError("Not forty-eight complete original monthly variables")
        manifest={
            "schema":SCHEMA,"status":STATUS,
            "source":CLIMATOLOGY,"units":{k:v["units"] for k,v in VARIABLES.items()},
            "source_public_base_url":URL_ROOT,
            "original_source_zip_sha256":archives,
            "source_original_readme_sha256":readmes,
            "months":bands,
            "grid":{"width":grid_ref[0],"height":grid_ref[1],
                    "affine6":grid_ref[2],"crs":grid_ref[3]},
            "n_original_monthly_rasters":band_count,
            "temperature_vapour_solar_precipitation_measured_at_photo_sites":False,
            "new_original_photo_colour_outcomes_measured":False,
            "original_image_photo_ids_opened":False,
            "temperature_derived_monthly_vpd_not_original_measurement":True,
            "climate_is_1970_2000_normal_not_photo_year_weather":True,
            "synthetic_test_only":bool(synthetic)
        }
        manifest_path=stage/"dynamic_monthly_manifest.json"
        manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf8")
        digest=sha256_file(manifest_path)
        receipt={
            "schema":"mois_source_attested_worldclim21_original_dynamic_monthlies_v1",
            "status":STATUS,
            "n_original_provider_zip_archives":4,
            "n_georeferenced_original_raster_layers":48,
            "original_provider_zip_sha256":archives,
            "manifest_sha256":digest,
            "climate_grid_shape":list(grid_ref[:2]),
            "provider_geotiff_data_units":manifest["units"],
            "biological_photo_pixels_measured":False,
            "iNaturalist_observation_date_joined":False,
            "global_white_fraction_estimated":False,
            "synthetic_test_only":bool(synthetic)
        }
        (stage/"dynamic_climate_source_receipt.json").write_text(
            json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf8")
        stage.replace(outroot)
        return receipt
    except Exception:
        # Leave incomplete stage visibly quarantined on interrupted source
        # operations. It is never reported as an admitted climate source.
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for kind in VARIABLES:
        p.add_argument(f"--original-{kind}-zip",type=Path,required=True)
        p.add_argument(f"--original-{kind}-sha256",required=True)
    p.add_argument("--out-private-source",type=Path,required=True)
    a=p.parse_args()
    orig={v:getattr(a,f"original_{v}_zip") for v in VARIABLES}
    sha={v:getattr(a,f"original_{v}_sha256") for v in VARIABLES}
    out=validate_originals(orig,sha,a.out_private_source)
    print("MOIS_ORIGINAL_48_MONTHLY_CLIMATE_BANDS_VERIFIED",
          json.dumps(out,sort_keys=True))


if __name__=="__main__":main()
