#!/usr/bin/env python3
"""MOIS: paired WorldClim 1.4 LGM-vs-current BIO1/BIO12 (3 GCMs) on true
GBIF 37,104 accepted-species x 100-km cells; genus-controlled species tests.

Current (1960-1990) and LGM (~22 ka) must share the WorldClim 1.4 baseline.
All source images are photographic 9-hue ROIs, NOT vetted petal pigmentation.
Current occupied GBIF cells are NOT evidence of LGM range/refugia.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

GBIF_SHA = "422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
ENV_SHA = "d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
GRID_FILE = "all_37104_species_x_EPSG6933_100km_cells.parquet"
ENV_FILE = "BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
GCMS = ("CCSM4", "MIROC-ESM", "MPI-ESM-P")
HUES = ("white", "yellow", "orange", "red", "pink",
        "magenta", "purple", "blue", "bronze")
OFFSETS = (-0.4, -0.2, 0.0, 0.2, 0.4)
SEED = 20261116

def require(value, message):
    if not bool(value):
        raise ValueError("MOIS_LGM: " + str(message))

def sha(b):
    if isinstance(b, (str, Path)):
        b = Path(b).read_bytes()
    return hashlib.sha256(b).hexdigest()

def archive_source(path, expected_sha):
    require(sha(path) == expected_sha, "signed existing source ZIP SHA differs")
    return zipfile.ZipFile(path)

def read_signed_sources(gbif_zip, env_zip):
    import pyarrow.parquet as pq
    with archive_source(gbif_zip, GBIF_SHA) as z:
        b = z.read(GRID_FILE)
        receipt = json.loads(z.read("receipt.json"))
    require(receipt["all_16_shards_source_verified"]
            and receipt["GBIF_occurrences_passed_all_filters"] == 396142825,
            "global GBIF verified 16-shard receipt absent")
    require(sha(b) == "e00b0408055984abfe4ec0ff0a24833fc148888afd1ea0e0c9b270a077d08b11",
            "GBIF species-grid source parquet SHA changed")
    grid = pq.read_table(io.BytesIO(b),
                         columns=["specieskey", "grid_x", "grid_y", "occurrence_count"]
                         ).to_pandas()
    require(len(grid) == 5493931 and int(grid.occurrence_count.sum()) == 396142825,
            "all original GBIF species-grid pairs not preserved")
    grid["specieskey"] = grid.specieskey.astype(str)

    with archive_source(env_zip, ENV_SHA) as z:
        b = z.read(ENV_FILE)
        receipt_env = json.loads(z.read("receipt.json"))
    require(sha(b) == receipt_env["output_SHA256"][ENV_FILE]
            and receipt_env["source_global_accepted_species"] == 37104,
            "signed complete climate+NDVI+woody AGB source changed")
    species = pq.read_table(io.BytesIO(b)).to_pandas()
    species["key"] = species.key.astype(str)
    require(len(species) == 37104 and species.key.is_unique
            and species.source_hue9.isin(HUES).all(), "photographic species table changed")
    return grid, species

def layer_number(name):
    """Parse WorldClim .bil/.tif, including bio1.bil or cclgmbi1.tif."""
    stem = Path(name).stem.lower()
    for pattern in (r"(?:bio|bi)(?:oclim)?[_-]?0?(1[0-9]|[1-9])(?:[_-]10m)?$",
                    r"(?:bio|bi)(?:oclim)?[_-]?0?(1[0-9]|[1-9])$"):
        m = re.search(pattern, stem)
        if m:
            return int(m.group(1))
    return None

def extract_two_bioclim(source_zip, tmpdir, label):
    """Extract only selected BIO1/12 plus BIL companions, not every 19 bands."""
    with zipfile.ZipFile(source_zip) as z:
        files = [x for x in z.namelist()
                 if Path(x).suffix.lower() in (".tif", ".tiff", ".bil")]
        matches = {}
        for n in files:
            k = layer_number(n)
            if k in (1, 12):
                matches.setdefault(k, []).append(n)
        require(set(matches) == {1, 12}
                and all(len(matches[i]) == 1 for i in (1, 12)),
                f"{label}: ambiguous/missing BIO1/12; examples={files[:14]}")
        out = {}
        for i in (1, 12):
            member = matches[i][0]
            base = Path(member).name
            target = Path(tmpdir) / f"{label}_{i}" / base
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(z.read(member))
            if target.suffix.lower() == ".bil":
                stem = str(Path(member).with_suffix(""))
                needed = [stem + ext for ext in (".hdr", ".prj")]
                have = {x.lower(): x for x in z.namelist()}
                for n in needed:
                    if n.lower() in have:
                        (target.parent / Path(n).name).write_bytes(z.read(have[n.lower()]))
                require(target.with_suffix(".hdr").exists(),
                        f"{label} BIO{i}: BIL header unavailable")
            out[i] = (target, sha(z.read(member)))
    return out

def read_band(path, index, label):
    import rasterio
    with rasterio.open(path) as ds:
        require(ds.count == 1 and ds.width >= 500 and ds.height >= 300,
                f"{label}: wrong 10m global raster dimensions")
        tr = ds.transform
        require(tr.a > 0 and tr.e < 0 and abs(tr.b) < 1e-12 and abs(tr.d) < 1e-12,
                f"{label}: raster not north-up")
        # Generic BIL may not advertise CRS; in that case require global WGS84 bounds.
        if ds.crs is not None:
            require(ds.crs.to_epsg() == 4326, f"{label}: not EPSG 4326")
        # WorldClim covers global land except Antarctica; bottom latitude is
        # typically -60, not -90. Do not incorrectly reject official rasters.
        require(abs(ds.bounds.left + 180) < 1.0 and abs(ds.bounds.right - 180) < 1.0
                and abs(ds.bounds.top - 90) < 1.0 and -90.1 <= ds.bounds.bottom <= -55,
                f"{label}: unexpected WorldClim longitude/latitude raster bounds")
        a = ds.read(1).astype(np.float32)
        nodata = ds.nodata
        if nodata is not None:
            a[a == nodata] = np.nan
        a[~np.isfinite(a)] = np.nan
        if index == 1:
            # WC 1.4 BIO1 is stored in tenths of a degree Celsius.
            a = a / 10.0
            a[(a < -90) | (a > 65)] = np.nan
        else:
            # BIO12: annual precipitation in mm.
            a[(a < 0) | (a > 30000)] = np.nan
        require(np.isfinite(a).sum() > 100000,
                f"{label} BIO{index}: unrealistic physical climate coverage")
        return a, (tr.c, tr.f, tr.a, tr.e)

def extract_at_points(values, coords, affine):
    """Numpy index into aligned EPSG4326 data; shape (25 x GBIF cells)."""
    lon, lat = coords
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    x0, y0, dx, dy = affine
    with np.errstate(invalid="ignore"):
        col = np.floor((lon - x0) / dx)
        row = np.floor((lat - y0) / dy)
    ok = (np.isfinite(col) & np.isfinite(row) & (col >= 0) & (row >= 0)
          & (col < values.shape[1]) & (row < values.shape[0]))
    out = np.full(lon.shape, np.nan, dtype=np.float32)
    ix = np.where(ok)
    out[ix] = values[row[ix].astype(np.int64), col[ix].astype(np.int64)]
    return out

def coords_in_grid(cells):
    from pyproj import Transformer
    t = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)
    x = cells.grid_x.to_numpy(dtype=float)
    y = cells.grid_y.to_numpy(dtype=float)
    coords = []
    for ox in OFFSETS:
        for oy in OFFSETS:
            lon, lat = t.transform((x + .5 + ox) * 100000.,
                                   (y + .5 + oy) * 100000.)
            coords.append((lon, lat))
    require(len(coords) == 25, "25 subsamples per equal-area 100km GBIF cell")
    lon = np.stack([pair[0] for pair in coords], axis=0)
    lat = np.stack([pair[1] for pair in coords], axis=0)
    return lon, lat

def paired_change(current, past):
    """Paired VALID same-land-position climate differences, per cell."""
    ct, cp = current
    pt, pp = past
    require(ct.shape == cp.shape == pt.shape == pp.shape,
            "paired climate arrays not aligned")
    mask = (np.isfinite(ct) & np.isfinite(cp) & np.isfinite(pt) & np.isfinite(pp))
    valid = np.sum(mask, axis=0)
    dt = np.where(mask, ct - pt, 0.)
    dp = np.where(mask, np.log1p(cp) - np.log1p(pp), 0.)
    tdiff = np.divide(np.sum(dt, axis=0), valid, out=np.full(valid.shape,np.nan),
                      where=valid > 0)
    pdiff = np.divide(np.sum(dp, axis=0), valid, out=np.full(valid.shape,np.nan),
                      where=valid > 0)
    return valid, tdiff, pdiff

def build_cell_changes(grid, current_zip, lgm_zips):
    cells = grid[["grid_x", "grid_y"]].drop_duplicates().reset_index(drop=True)
    require(len(cells) == 18577, "signed GBIF cell universe drift")
    coords = coords_in_grid(cells)
    receipt = {"input_archives_SHA256": {"WorldClim1p4_current":sha(current_zip)},
               "models": list(GCMS),
               "baseline_current":"1960-1990",
               "period_LGM":"about 22000 years before present",
               "points_per_cell":25}
    with tempfile.TemporaryDirectory() as td:
        paths = extract_two_bioclim(current_zip,td,"CURRENT")
        base = []
        for k in (1,12):
            values,affine = read_band(paths[k][0],k,"CURRENT")
            base.append(extract_at_points(values,coords,affine))
        receipt["current_BIO_layer_SHA256"] = {f"BIO{k}":paths[k][1] for k in (1,12)}
        for label, zipname in zip(GCMS,lgm_zips):
            p = extract_two_bioclim(zipname,td,label.replace("-","_"))
            lgm = []
            for k in (1,12):
                values,affine = read_band(p[k][0],k,label)
                lgm.append(extract_at_points(values,coords,affine))
            n, t, rain = paired_change(base, lgm)
            prefix=label.replace("-","_")
            cells[f"{prefix}_paired_points"] = n.astype("int16")
            cells[f"{prefix}_delta_temperature_C"] = t
            cells[f"{prefix}_delta_log_precip"] = rain
            receipt["input_archives_SHA256"][label]=sha(zipname)
            receipt[f"{label}_BIO_layer_SHA256"]={f"BIO{k}":p[k][1] for k in (1,12)}
            receipt[f"{label}_cells_valid"] = int(np.sum(n>0))
            print("MOIS_LGM_MODEL_PAIRED",label,"valid cells",int(np.sum(n>0)),
                  "median warming",float(np.nanmedian(t)),flush=True)
    cols_t=[f"{x.replace('-','_')}_delta_temperature_C" for x in GCMS]
    cols_p=[f"{x.replace('-','_')}_delta_log_precip" for x in GCMS]
    cells["threeGCM_delta_temperature_C_median"]=cells[cols_t].median(axis=1,skipna=False)
    cells["threeGCM_delta_log_precip_median"]=cells[cols_p].median(axis=1,skipna=False)
    cells["threeGCM_valid_all_models"]=cells[cols_t + cols_p].notna().all(axis=1)
    return cells,receipt

def attach_species(grid, species, cells):
    p=grid[["specieskey","grid_x","grid_y"]].merge(
        cells,on=["grid_x","grid_y"],how="left",validate="many_to_one")
    require(len(p)==5493931, "all GBIF species x cells must survive join")
    measure=[c for c in cells.columns
             if c not in ("grid_x","grid_y") and c!="threeGCM_valid_all_models"
             and not c.endswith("_paired_points")]
    d=p.groupby("specieskey", sort=False)[measure].mean()
    cover=p.assign(valid=p.threeGCM_valid_all_models.fillna(False).astype(int)).groupby(
        "specieskey",sort=False).agg(
            lgm_valid_occupied_cells=("valid","sum"),
            lgm_all_occupied_cells=("valid","size"))
    d=d.join(cover)
    d["lgm_climate_coverage"]=d.lgm_valid_occupied_cells/d.lgm_all_occupied_cells
    out=species.merge(d,left_on="key",right_index=True,how="left",
                      validate="one_to_one")
    require(len(out)==37104 and out.key.is_unique
            and (out.occupied_grid_cells.fillna(0).astype(int)
                 ==out.lgm_all_occupied_cells.fillna(0).astype(int)).all(),
            "join lost species or recorded GBIF occupation")
    return out

def genus_model(dat, historical_columns, equal_genus=False):
    """Species one-row within-genus model with genus-clustered uncertainty."""
    import statsmodels.api as sm
    y=dat.source_hue9.eq("white").astype(float)
    x=pd.DataFrame({
        "BIO01_C":dat.BIO01.astype(float),
        "log_BIO12":np.log1p(dat.BIO12.astype(float).clip(lower=0)),
        "NDVI":dat.MODIS_NDVI_2001_2021_mean.astype(float),
        "log_woody_AGB":np.log1p(dat.ESA_AGB_woody_Mg_ha_2020.astype(float).clip(lower=0)),
        "abs_lat":dat.mean_abs_latitude_over_occupied_cells.astype(float),
        "log_occ":np.log1p(dat.records_passed_spatial_QA.astype(float)),
    },index=dat.index)
    for c in historical_columns:
        x[c]=dat[c].astype(float)
    for status in sorted(dat.hue_status.astype(str).unique())[1:]:
        x["status_"+status]=dat.hue_status.astype(str).eq(status).astype(float)
    group=dat.genus.astype(str)
    numerical=[c for c in x if not c.startswith("status_")]
    for c in numerical:
        sd=x[c].std()
        require(np.isfinite(sd) and sd>1e-9,f"nonsensical model variable {c}")
        x[c]=(x[c]-x[c].mean())/sd
    xc=x-x.groupby(group).transform("mean")
    yc=y-y.groupby(group).transform("mean")
    xc=xc.loc[:,xc.std()>1e-9]
    require(set(historical_columns).issubset(xc.columns),"historic effect lost after genus adjustment")
    weights=(1/group.map(group.value_counts()).astype(float)
             if equal_genus else None)
    if equal_genus:
        model=sm.WLS(yc,xc,weights=weights).fit(
            cov_type="cluster",cov_kwds={"groups":group},use_t=False)
    else:
        model=sm.OLS(yc,xc).fit(
            cov_type="cluster",cov_kwds={"groups":group},use_t=False)
    estimates={}
    for c in historical_columns:
        lo,hi=model.conf_int().loc[c].tolist()
        estimates[c]={"beta_per_species_SD":float(model.params[c]),
                      "95CI":[float(lo),float(hi)],"p":float(model.pvalues[c])}
    R=np.zeros((len(historical_columns),len(xc.columns)))
    for i,k in enumerate(historical_columns):
        R[i,xc.columns.get_loc(k)]=1
    joint=model.wald_test(R,scalar=True)
    return {"n_species":len(dat),"n_genera":int(group.nunique()),
            "species_equal_vs_genus_equal":"genus_equal" if equal_genus else "species_equal",
            "r_squared_within_genus":float(model.rsquared),
            "p_joint_historic":float(joint.pvalue),"historic_effects":estimates}

def compare_species_models(species):
    from statsmodels.stats.multitest import multipletests
    cols = [f"{c.replace('-','_')}_delta_{metric}"
            for c in GCMS for metric in ("temperature_C","log_precip")]
    cols += ["threeGCM_delta_temperature_C_median","threeGCM_delta_log_precip_median"]
    other = ["BIO01","BIO12","MODIS_NDVI_2001_2021_mean",
             "ESA_AGB_woody_Mg_ha_2020","mean_abs_latitude_over_occupied_cells",
             "records_passed_spatial_QA","genus","source_hue9","hue_status"]
    cols_all=cols+other
    complete=(species[cols].notna().all(axis=1)
              &species[other[:6]].notna().all(axis=1)
              &(species.lgm_climate_coverage>=.75)
              &(species.BIO12>=0)
              &(species.ESA_AGB_woody_Mg_ha_2020>=0)
              &species.genus.notna()&species.hue_status.notna())
    d=species.loc[complete].copy()
    mixed=d.groupby("genus").source_hue9.transform(lambda a: len(set(a.eq("white")))==2)
    d=d.loc[mixed].copy()
    require(len(d)>10000 and d.genus.nunique()>1000,
            "not enough true species with within-genus hue contrast")
    models={}
    p_values=[]
    for name in list(GCMS)+["threeGCM"]:
        cols_here = ([f"{name.replace('-','_')}_delta_temperature_C",
                      f"{name.replace('-','_')}_delta_log_precip"]
                     if name in GCMS else
                     ["threeGCM_delta_temperature_C_median",
                      "threeGCM_delta_log_precip_median"])
        models[name]={}
        for lab,equal in (("species_equal",False),("genus_equal",True)):
            models[name][lab]=genus_model(d,cols_here,equal)
        p_values.append(models[name]["genus_equal"]["p_joint_historic"])
    q=multipletests(p_values,method="fdr_bh")[1].tolist()
    for (label,m),val in zip(models.items(),q):
        m["genus_equal"]["q_BH_four_historic_models"]=float(val)
    return {"n_complete_before_genus_filter":int(complete.sum()),
            "n_species_used":len(d),"n_genera":int(d.genus.nunique()),
            "models":models,
            "interpretation":"observational; current GBIF occupancy not ancient refugia",
            "source_photo_hue_is_vetted_local_flower_colour":False}

def run(gbif,env,current,lgms,out):
    out=Path(out)
    require(not out.exists(), "refuse overwrite output")
    grid,species=read_signed_sources(gbif,env)
    cells,receipt=build_cell_changes(grid,current,lgms)
    out_species=attach_species(grid,species,cells)
    model=compare_species_models(out_species)
    out.mkdir(parents=True)
    cfile=out/"WorldClim1p4_LGM_threeGCM_18577_cells.parquet"
    sfile=out/"GBIF37104_photo9_LGM_climate_species.parquet"
    cells.to_parquet(cfile,index=False)
    out_species.to_parquet(sfile,index=False)
    (out/"genus_model_results.json").write_text(json.dumps(model,indent=2,sort_keys=True)+"\n")
    audit={"schema":"mois_real_GBIF37104_threeGCM_LGM_climate_20261116",
           "gbif_accepted_species":37104,"species_x_grid_pairs":5493931,
           "source_valid_GBIF_occurrences":396142825,
           "n_current_occupied_global_cells":18577,
           "prior_climate_source_is_WorldClim1p4":True,
           "current_1960_1990_matches_LGM_calibration":True,
           "current_occupancy_is_not_refugium":True,
           "photo9_is_not_verified_true_petal_colour":True,
           "source_biomass_is_woody_AGB_not_NPP":True,
           "worldclim":receipt,"results":model,
           "output_SHA256":{f.name:sha(f) for f in out.iterdir() if f.is_file()}}
    (out/"receipt.json").write_text(json.dumps(audit,indent=2,sort_keys=True)+"\n")
    print("MOIS_LGM_FINISHED_REAL_SIGNED",json.dumps({
        "species":len(out_species),"cells":len(cells),
        "n_species_regression":model["n_species_used"],
        "n_genera":model["n_genera"],
        "joint_genus_equal_q":{k:v["genus_equal"]["q_BH_four_historic_models"]
                                for k,v in model["models"].items()}},
        sort_keys=True),flush=True)
    return audit

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for arg in ("gbif","env","current","ccsm4","miroc","mpi","output"):
        p.add_argument("--"+arg,required=True)
    a=p.parse_args()
    run(a.gbif,a.env,a.current,[a.ccsm4,a.miroc,a.mpi],a.output)
