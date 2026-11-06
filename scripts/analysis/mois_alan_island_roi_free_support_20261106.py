#!/usr/bin/env python3
"""Island-first ALAN feasibility on ORIGINAL MOIS observation IDs, without ROI selection.

This is a COUNTS/SUPPORT audit, not a fitted ALAN effect or native-flora census.
Input: SHA-pinned 1,032,882 canonical observations, accepted WCVP species,
independently documented flower colour, positive GloBI moth visit records,
and an independently sourced island polygon GPKG (prefer Island GSHHG 2.3.7).
No documented moth record is NOT evidence of no moth visit; positive visit
is NOT evidence of effective nocturnal pollination or dependency.

A review-approved island×accepted-species establishment table is optional;
WCVP broad TDWG3 native range must NEVER be silently used as island native.
The optional annual VIIRS source is distinct from the 2020/2024 FOUR-MONTH
proxy. Neither annual brightness nor local effective ALAN is invented.
No individual licensed species-level trait or image-level data is exported.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

SOURCE_SHA = {
    "canonical":"f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c",
    "names":"dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2",
    "accepted":"9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3",
    "direct":"7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577",
    "globi":"a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994",
}
COLOURS = {"white", "yellow_orange", "red_pink", "blue_purple", "green_brown_inconspicuous"}
YEARS = set(range(2012, 2025))

def require(test, message):
    if not bool(test):
        raise ValueError("MOIS_ISLAND_ROIFREE: "+message)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def cell005(lon, lat):
    return (np.floor((pd.to_numeric(lon)+180)/0.05).astype("int64").astype(str)
            + "_" + np.floor((pd.to_numeric(lat)+90)/0.05).astype("int64").astype(str))

def accepted_species_mapping(names, accepted):
    require({"inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"} <= set(names),
            "WCVP taxon table columns absent")
    require({"wcvp_accepted_plant_name_id","wcvp_accepted_name"} <= set(accepted),
            "WCVP accepted name table columns absent")
    a = accepted[["wcvp_accepted_plant_name_id","wcvp_accepted_name"]].drop_duplicates()
    require(a.wcvp_accepted_plant_name_id.is_unique, "Conflicting accepted WCVP ID")
    n = names.loc[names.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
                  ["inat_taxon_id","wcvp_accepted_plant_name_id"]].drop_duplicates()
    require(n.inat_taxon_id.is_unique, "One source taxon maps to multiple accepted IDs")
    return n.merge(a, on="wcvp_accepted_plant_name_id", how="left",
                   validate="many_to_one")[["inat_taxon_id","wcvp_accepted_name"]]

def independent_colour(direct):
    required={"accepted_species","trait_name","normalized_value","quality",
              "resolution_status","classification"}
    require(required <= set(direct), "Independent colour ledger columns absent")
    x=direct.loc[direct.trait_name.eq("flower_primary_color")
        & direct.quality.isin(["high","medium"])
        & direct.resolution_status.eq("resolved")
        & ~direct.classification.isin(["cultivar_contamination","partial_self_compatibility"]),
        ["accepted_species","normalized_value"]].copy()
    x=x.loc[x.normalized_value.isin(COLOURS)]
    require(x.accepted_species.is_unique, "Conflicting direct single-state flower colour")
    x["colour_state"]=np.where(x.normalized_value.eq("white"),"white","nonwhite")
    return x.rename(columns={"accepted_species":"wcvp_accepted_name"})[
        ["wcvp_accepted_name","colour_state"]]

def moth_evidence(globi):
    require({"plant_name","guild"} <= set(globi), "GloBI guild input missing")
    require(globi.guild.isin(["moth_family_visit","butterfly_family_visit"]).all(),
            "Not a verified GloBI positive-flower-visit input")
    x=globi[["plant_name","guild"]].drop_duplicates().assign(v=True).pivot(
        index="plant_name",columns="guild",values="v").fillna(False)
    for name in ("moth_family_visit","butterfly_family_visit"):
        if name not in x: x[name]=False
    return x.rename(columns={"moth_family_visit":"moth_positive",
           "butterfly_family_visit":"butterfly_positive"}).reset_index().rename(
               columns={"plant_name":"wcvp_accepted_name"})[
               ["wcvp_accepted_name","moth_positive","butterfly_positive"]]

def geolocate_islands(photos, polygons, minimum_km2=20):
    """Only genuine point-in-polygon; ambiguous boundary/overlap => unknown."""
    required={"island_id","area_km2","geometry"}
    require(required <= set(polygons), "Island polygon GPKG lacks required columns")
    require(polygons.crs is not None and polygons.crs.to_epsg()==4326,
            "Island polygon CRS must be EPSG:4326")
    p=polygons.loc[pd.to_numeric(polygons.area_km2,errors="coerce")>=minimum_km2,
                   ["island_id","area_km2","geometry"]].copy()
    require(p.island_id.notna().all() and p.island_id.is_unique,
            "Island IDs missing/repeated")
    require(p.geometry.is_valid.all() and (~p.geometry.is_empty).all(),
            "Invalid island polygon: refuse silent geometry repair")
    p["area_km2"]=pd.to_numeric(p.area_km2)
    candidates=photos.loc[photos.latitude.between(-90,90)
        & photos.longitude.between(-180,180),["photo_id","latitude","longitude"]]
    g=gpd.GeoDataFrame(candidates,geometry=gpd.points_from_xy(
        candidates.longitude,candidates.latitude),crs="EPSG:4326")
    match=gpd.sjoin(g,p,how="inner",predicate="intersects")[["photo_id","island_id"]]
    ambiguity=match.photo_id.value_counts()
    bad=set(ambiguity[ambiguity>1].index)
    clean=match.loc[~match.photo_id.isin(bad)].drop_duplicates("photo_id")
    require(clean.photo_id.is_unique, "Multiple island IDs per source image")
    return clean, {"island_polygon_units_area_filtered":len(p),
       "photos_with_valid_coordinates":len(candidates),
       "photos_in_multiple_island_polygons_excluded":len(bad),
       "photos_unambiguously_inside_islands":len(clean)}

def island_native_status(x, establishment=None):
    """Unknown is NOT introduced; no broad TDWG native status substituted."""
    if establishment is None:
        x["establishment"]="unknown"
        return x
    required={"island_id","wcvp_accepted_name","establishment_status","review_status"}
    require(required<=set(establishment),"Audited island-native source columns absent")
    e=establishment.loc[establishment.review_status.eq("accepted"),
        ["island_id","wcvp_accepted_name","establishment_status"]].copy()
    require(e.establishment_status.isin(
        ["native","introduced","cultivated","transient","unknown"]).all(),
        "Invalid establishment source status")
    require(not e.duplicated(["island_id","wcvp_accepted_name"]).any(),
        "Conflicting island-specific establishment")
    e=e.rename(columns={"establishment_status":"establishment"})
    x=x.merge(e,on=["island_id","wcvp_accepted_name"],how="left",validate="many_to_one")
    x["establishment"]=x.establishment.fillna("unknown")
    return x

def legacy_viirs_support(v20=None,v24=None):
    frames=[]
    if v20 is not None:
        require({"cell_id","source_complete","viirs_quarterly_median"} <= set(v20),
                "2020 quarter-month proxy columns absent")
        require(v20.cell_id.is_unique,"2020 light cell duplicated")
        a=v20.loc[v20.source_complete.eq(True),
                  ["cell_id","viirs_quarterly_median"]].copy()
        a.columns=["cell_id","radiance"]
        a["year"]=2020; a["source"]="quarter_months_NOT_annual"
        frames.append(a)
    if v24 is not None:
        require({"cell_id","VIIRS_2024_quarter_median","VIIRS_2024_n_valid_months"}<=set(v24),
                "2024 quarter-month proxy columns absent")
        require(v24.cell_id.is_unique,"2024 light cell duplicated")
        a=v24.loc[pd.to_numeric(v24.VIIRS_2024_n_valid_months)>=3,
                  ["cell_id","VIIRS_2024_quarter_median"]].copy()
        a.columns=["cell_id","radiance"]
        a["year"]=2024; a["source"]="quarter_months_NOT_annual"
        frames.append(a)
    return pd.concat(frames,ignore_index=True) if frames else pd.DataFrame({
        "cell_id":pd.Series(dtype="string"),"radiance":pd.Series(dtype="float64"),
        "year":pd.Series(dtype="int64"),"source":pd.Series(dtype="string")})

def signed_annual_viirs(annual):
    required={"cell_id","year","radiance","cf_cvg","product_version"}
    require(required<=set(annual),"Annual VIIRS sourced cells/coverage absent")
    a=annual.copy()
    a["year"]=pd.to_numeric(a.year,errors="coerce")
    require(a.year.isin(YEARS).all(),"Annual VIIRS years outside declared frame")
    require(a.product_version.isin(["VNL_V2.1","VNL_V2.2"]).all(),
            "Untracked annual VIIRS product version")
    require(not a.duplicated(["year","cell_id"]).any(),"Annual VIIRS cells not unique")
    a["radiance"]=pd.to_numeric(a.radiance,errors="coerce")
    a["cf_cvg"]=pd.to_numeric(a.cf_cvg,errors="coerce")
    a=a.loc[np.isfinite(a.radiance)&a.cf_cvg.gt(0)].copy()
    a["source"]="annual_versioned_VNL"
    return a[["cell_id","radiance","year","source"]]

def prepare_photo_frame(canonical, mapping, colour, moth, assignment, establishment=None):
    need={"photo_id","observation_id","inat_taxon_id","observed_on",
          "latitude","longitude","source_subarm","strict9"}
    require(need<=set(canonical),"Canonical observation columns absent")
    require(canonical.photo_id.is_unique and canonical.observation_id.is_unique,
            "Original photo/observation IDs not unique")
    x=canonical.loc[canonical.source_subarm.eq("historical_main")].copy()
    x["year"]=pd.to_numeric(x.observed_on.astype(str).str[:4],errors="coerce")
    x=x.loc[x.year.isin(YEARS)].copy()
    x=x.merge(assignment,on="photo_id",how="inner",validate="one_to_one")
    x=x.merge(mapping,on="inat_taxon_id",how="left",validate="many_to_one")
    x=x.merge(colour,on="wcvp_accepted_name",how="left",validate="many_to_one")
    x=x.merge(moth,on="wcvp_accepted_name",how="left",validate="many_to_one")
    x["moth_positive"]=x.moth_positive.eq(True)
    x["butterfly_positive"]=x.butterfly_positive.eq(True)
    x=island_native_status(x,establishment)
    require(x.photo_id.is_unique and x.observation_id.is_unique,
            "Taxon or establishment join multiplied original photos")
    x["cell_id"]=cell005(x.longitude,x.latitude)
    return x

def summarise(x, light, *, privacy_min_species=5):
    """Per-island/year PHOTO support, distinct species and actual light ranges."""
    require(not light.duplicated(["year","cell_id"]).any(), "Duplicate site-year light")
    x=x.merge(light,on=["year","cell_id"],how="left",validate="many_to_one")
    out=[]
    for (island,year),s in x.groupby(["island_id","year"],sort=True):
        known=s.loc[s.colour_state.isin(["white","nonwhite"])]
        moth=s.loc[s.moth_positive]
        native=s.loc[s.establishment.eq("native")]
        documented_species=known.wcvp_accepted_name.nunique()
        m=s[["cell_id","radiance"]].dropna().drop_duplicates("cell_id")
        n_lit=len(m)
        # Do not substitute missing annual radiance with zero or with legacy months.
        val=pd.to_numeric(m.radiance,errors="coerce")
        val=val[np.isfinite(val)]
        logval=np.log1p(np.maximum(val,0))
        out.append({
           "island_id":island,"year":int(year),
           "photos_all_original_no_roi_gate":len(s),
           "observations":s.observation_id.nunique(),
           "photos_strict9_success_diagnostic_only":int(s.strict9.eq(True).sum()),
           "accepted_species":int(s.wcvp_accepted_name.nunique()),
           "photos_with_independent_single_colour":len(known),
           "documented_single_colour_species":int(documented_species),
           "photos_with_positive_moth_visit_record":len(moth),
           "moth_visit_record_species":int(moth.wcvp_accepted_name.nunique()),
           "photos_colour_and_moth_record":int((s.colour_state.notna()&s.moth_positive).sum()),
           "species_colour_and_moth_record":int(known.loc[known.moth_positive,
                                                "wcvp_accepted_name"].nunique()),
           "species_both_moth_and_butterfly_records":int(s.loc[
             s.moth_positive&s.butterfly_positive,"wcvp_accepted_name"].nunique()),
           "accepted_island_native_status_photos":len(native),
           "accepted_island_native_status_species":int(native.wcvp_accepted_name.nunique()),
           "accepted_island_introduced_status_photos":int(s.establishment.eq("introduced").sum()),
           "unknown_island_establishment_photos":int(s.establishment.eq("unknown").sum()),
           "same_year_viirs_matched_photos":int(s.radiance.notna().sum()),
           "same_year_viirs_distinct_cells":n_lit,
           "viirs_min":float(val.min()) if len(val) else None,
           "viirs_max":float(val.max()) if len(val) else None,
           "viirs_log1p_q90_minus_q10":float(np.quantile(logval,.9)-np.quantile(logval,.1))
                   if len(logval)>=3 else None,
           "light_source":str(s.source.dropna().iloc[0]) if s.source.notna().any()
                          else "NO_same_year_light",
           # Release only support aggregates; never species-direct restricted values.
           "known_white_species":int(known.loc[known.colour_state.eq("white"),
                                      "wcvp_accepted_name"].nunique())
                   if documented_species>=privacy_min_species else None,
           "direct_colour_counts_suppressed":bool(documented_species<privacy_min_species)
        })
    return pd.DataFrame(out)

def run(args):
    inputs={key:Path(getattr(args,key)) for key in SOURCE_SHA}
    for key,path in inputs.items():
        require(sha(path)==SOURCE_SHA[key],"Frozen original SHA changed: "+key)
    target=Path(args.output)
    require(not target.exists(),"Refuse to overwrite existing audit")
    c=pd.read_csv(inputs["canonical"],compression="gzip",low_memory=False,
      dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
      usecols=["photo_id","observation_id","inat_taxon_id","observed_on",
               "latitude","longitude","source_subarm","strict9"])
    require(len(c)==1_032_882,"Canonical photographed observation frame changed")
    historical=c.loc[c.source_subarm.eq("historical_main")].copy()
    historical["year"]=pd.to_numeric(historical.observed_on.astype(str).str[:4],errors="coerce")
    historical=historical.loc[historical.year.isin(YEARS)]
    islands=gpd.read_file(args.islands,layer=args.islands_layer)
    assignment,geo_audit=geolocate_islands(historical,islands,args.min_area_km2)
    n=pd.read_csv(inputs["names"],compression="gzip",dtype=str,
                  usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    a=pd.read_csv(inputs["accepted"],compression="gzip",dtype=str,
                  usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    d=pd.read_csv(inputs["direct"],compression="gzip",dtype=str,
                  usecols=["accepted_species","trait_name","normalized_value","quality",
                           "resolution_status","classification"])
    g=pd.read_csv(inputs["globi"],compression="gzip",usecols=["plant_name","guild"])
    native=pd.read_csv(args.island_establishment,dtype=str) if args.island_establishment else None
    x=prepare_photo_frame(historical,accepted_species_mapping(n,a),
       independent_colour(d),moth_evidence(g),assignment,native)
    legacy=legacy_viirs_support(
        pd.read_csv(args.viirs2020,dtype={"cell_id":str}) if args.viirs2020 else None,
        pd.read_csv(args.viirs2024,dtype={"cell_id":str}) if args.viirs2024 else None)
    annual=(signed_annual_viirs(pd.read_csv(args.viirs_annual,dtype={"cell_id":str}))
            if args.viirs_annual else pd.DataFrame(columns=legacy.columns))
    # An annual record, if available, takes priority. The quarterly proxy is a
    # separate diagnostic; never splice it into an "annual" time series.
    source=annual if args.viirs_annual else legacy
    table=summarise(x,source)
    target.mkdir(parents=True)
    table.to_csv(target/"island_year_ROI_free_source_support.csv",index=False)
    year=x.groupby("year").agg(
        attempted_observations=("observation_id","nunique"),
        island_units=("island_id","nunique"),
        accepted_species=("wcvp_accepted_name","nunique")).reset_index()
    year.to_csv(target/"yearly_photo_source_support.csv",index=False)
    receipt={
      "schema":"mois_ALAN_island_ROIfree_source_support_v1",
      "input_sha256":SOURCE_SHA,
      "island_geometry_sha256":sha(args.islands),
      "island_geometry_min_area_km2":args.min_area_km2,
      "geo_audit":geo_audit,
      "historical_source_2012_2024_observations_before_polygon":int(len(historical)),
      "actual_island_original_photo_observations":int(len(x)),
      "islands_with_original_photographs":int(x.island_id.nunique()),
      "photographic_years":sorted(map(int,x.year.dropna().unique().tolist())),
      "island_photo_independent_direct_colour":int(x.colour_state.notna().sum()),
      "island_photo_moth_visit_positive":int(x.moth_positive.sum()),
      "island_photo_direct_colour_and_moth_positive":int(
          (x.colour_state.notna()&x.moth_positive).sum()),
      "island_photo_strict9_success_not_used_for_inclusion":int(x.strict9.eq(True).sum()),
      "island_specific_establishment_review_source_supplied":args.island_establishment is not None,
      "annual_VNL_product_source_supplied":args.viirs_annual is not None,
      "only_four_month_VIIRS_proxy_available_if_no_annual":args.viirs_annual is None,
      "annual_albedo_or_ground_illuminance_NOT_measured":True,
      "no_moth_record_is_not_moth_absence":True,
      "moth_visit_is_not_independent_effective_nocturnal_pollination":True,
      "photograph_not_observed_is_not_flora_absence":True,
      "no_original_RGB_download_no_ROI_reinference":True,
      "individual_colour_or_native_restricted_rows_published":False,
      "output_sha256":{
         p.name:sha(p) for p in sorted(target.glob("*.csv"))}}
    (target/"source_receipt.json").write_text(json.dumps(
        receipt,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("MOIS_ALAN_ISLAND_REAL_ROI_FREE_SUPPORT",json.dumps({
        k:receipt[k] for k in ("historical_source_2012_2024_observations_before_polygon",
        "actual_island_original_photo_observations","islands_with_original_photographs",
        "island_photo_independent_direct_colour","island_photo_moth_visit_positive",
        "island_photo_direct_colour_and_moth_positive",
        "annual_VNL_product_source_supplied")},sort_keys=True),flush=True)
    return table,receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in SOURCE_SHA:p.add_argument("--"+key,required=True)
    p.add_argument("--islands",required=True)
    p.add_argument("--islands-layer",default="islands")
    p.add_argument("--min-area-km2",type=float,default=20)
    p.add_argument("--viirs2020")
    p.add_argument("--viirs2024")
    p.add_argument("--viirs-annual",help="CSV: year,cell_id,radiance,cf_cvg,product_version")
    p.add_argument("--island-establishment",help="review-accepted island×species native status only")
    p.add_argument("--output",required=True)
    run(p.parse_args())

if __name__=="__main__":
    main()
