#!/usr/bin/env python3
"""Re-export ONLY scientific columns from the immutable 37,104-accepted
species MOIS environment/photographic nine-hue catalog into portable CSV.

This is a convenience adapter, not a source-data update. Explicit signed
source archive and inner Parquet checks; no iNaturalist photo pixels and no
third-party GIFT pollination species list included. That independent list
is supplied separately from the pre-existing original island H5C signed
qualification artifact at analysis time, not committed to the MOIS repo.

The environmental covariates describe each species' PRESENT GBIF
distribution footprint, not the photographed flowering organ or historical
native range.
"""
from __future__ import annotations
import argparse, gzip,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import pandas as pd
import pyarrow.parquet as pq
SRC_SHA="d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
MEMBER="BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
INNER_SHA="e772092266ebe603c26414e30a9f31f09990c2cfc608b590bb80a1e5c7ff306b"
COLS=[
 "key","species","genus","source_hue9","hue_status",
 "BIO01","BIO12","vpd_kPa_approx","MODIS_NDVI_2001_2021_mean",
 "ESA_AGB_woody_Mg_ha_2020","mean_abs_latitude_over_occupied_cells",
 "records_passed_spatial_QA","occupied_grid_cells",
 "srad_monthly_mean","wind_monthly_mean"
]
HUES={"white","yellow","orange","red","pink","magenta","purple","blue","bronze"}
def digest(x):return hashlib.sha256(x).hexdigest()
def require(ok,what):
    if not bool(ok):raise ValueError("MOIS_GIFT_INDEPENDENT_SOURCE: "+what)
def run(source,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwrite")
    require(digest(Path(source).read_bytes())==SRC_SHA,"Original source ZIP SHA changed")
    with ZipFile(source) as z:
        blob=z.read(MEMBER)
        meta=json.loads(z.read("receipt.json"))
    require(digest(blob)==INNER_SHA==meta["output_SHA256"][MEMBER],
            "Original Parquet SHA mismatch")
    table=pq.read_table(io.BytesIO(blob))
    require(set(COLS).issubset(table.column_names),"Required original ecological columns absent")
    d=table.select(COLS).to_pandas()
    d["key"]=d.key.astype(str)
    require(len(d)==37104 and d.key.is_unique,
            "COLXR global accepted species denominator changed")
    require(d.species.nunique()==len(d) and d.species.notna().all(),
            "Exact binomial species must be one per original source taxon")
    require(set(d.source_hue9)==HUES and
            d.genus.notna().all(),"Original nine photo hue states drifted")
    out.mkdir(parents=True)
    target=out/"MOIS_37104_accepted_species_photo9_and_current_environment.csv.gz"
    raw=d.to_csv(index=False,float_format="%.12g").encode("utf8")
    with target.open("wb") as stream:
        with gzip.GzipFile(fileobj=stream,mode="wb",mtime=0,
                           filename="") as zz:zz.write(raw)
    summary={
        "source_original_ZIP_SHA256":SRC_SHA,
        "source_parquet_inner_SHA256":INNER_SHA,
        "global_accepted_species":len(d),
        "global_photo_hue_counts":d.source_hue9.value_counts().sort_index().to_dict(),
        "global_genus_count":d.genus.nunique(),
        "export_n_columns":len(COLS),
        "not_independent_true_petal_colour":True,
        "GIFT_pollen_vector_mode_is_NOT_in_this_source_or_export":True,
        "source_CSV_SHA256":digest(target.read_bytes()),
        "original_GBIF_exposures_not_native_distribution":True
    }
    (out/"source_receipt.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("MOIS_INDEPENDENT_GIFT_CROSSWALK_SOURCE_READY",json.dumps(summary),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args();run(a.source,a.out)
