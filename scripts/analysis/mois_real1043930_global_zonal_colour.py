#!/usr/bin/env python3
"""Real source-backed 1,043,930 photo MOIS global white-first colour atlas.

Global gradients: 0/23.5/55/66.56 degree latitude in BOTH hemispheres,
meteorological seasons by hemisphere, WorldClim annual humidity/rain,
elevation, and GSHHG island/mainland geography. EVERY photo ID remains in
the sampling/ROI success denominator. Only 390,541 original admissible nine
palette ROIs enter primary colour data, 36,224 ambiguous as sensitivity only.

Main responses: continuous white palette area, dominant flower-ROI colour and
fraction, runner-up fraction, within-ROI nine-colour entropy, bronze/yellow/
pink etc percentages. The nine bins lack GREEN, therefore the presence of
green floral petals is NOT estimable from this source; actual sRGB CIELAB on
the 256 K1 photo pilot is a separately audited technical layer. CIELAB L*,
C*, hue must never be invented from the stored nine palette bins.

These are photographed *floral union ROI* appearances, not native wild flora
frequency, phenotypic plasticity, causal selection, guaranteed petals,
anthocyanin content or photo-count-independent population colours.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd

PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
COLUMNS=tuple("continuous9_"+s for s in PALETTE)
N=1_043_930
STRICT=390_541
AMBIGUOUS=36_224
FOUR=459_235
PHOTO_SHA="b20bbefe15020c21c3b5b0f7edf3f8e08148a4ee7dda8bf772a89bf6b126f447"
CLIMATE_SHA="ec39b49649acc2ca4a468c00c86d05800e33ea2db1f91a3b119076be38432735"
GEO_SHA=""; # Retrieved from independently frozen geography receipt, not a guess.
SCHEMA="mois_genuine1043930_global_zonal9_maincolour_white_first_v1"

def check(ok,msg):
    if not ok:raise ValueError("MOIS_GLOBAL_COLOUR_ZONES: "+msg)

def sha(file):
    h=hashlib.sha256()
    with Path(file).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def latitude_zone(lat):
    """Preserve actual geographic hemispheres, including polar photo support."""
    if not np.isfinite(lat) or abs(lat)>90:
        return "unknown"
    hemi="N" if lat>=0 else "S"
    absolute=abs(lat)
    if absolute>=66.5622:band="polar_66.56_90"
    elif absolute>=55:band="subpolar_55_66.56"
    elif absolute>=23.5:band="temperate_23.5_55"
    else:band="tropical_0_23.5"
    return hemi+"_"+band

def season_label(lat,month):
    if not np.isfinite(lat) or abs(lat)>90 or pd.isna(month) or not 1<=month<=12:
        return "unknown"
    m=int(month)
    north=lat>=0
    if m in (12,1,2):
        return "winter" if north else "summer"
    if m in (3,4,5):
        return "spring" if north else "autumn"
    if m in (6,7,8):
        return "summer" if north else "winter"
    return "autumn" if north else "spring"

def elev_band(v):
    if not np.isfinite(v):return "unknown"
    if v<0:return "below_sealevel"
    if v<500:return "low_0_500m"
    if v<1500:return "mid_500_1500m"
    if v<2500:return "high_1500_2500m"
    return "alpine_2500m_plus"

def q5(values):
    a=np.asarray(values,dtype=float)
    finite=np.isfinite(a)
    result=np.full(len(a),"unknown",dtype=object)
    if finite.any():
        cut=np.quantile(a[finite],[.2,.4,.6,.8])
        idx=np.searchsorted(cut,a[finite],side="right")+1
        result[finite]=np.array(["Q"+str(x) for x in idx],dtype=object)
        return result,[float(x) for x in cut]
    return result,[]

def main_colour_fraction_table(photos,source,climate,geography,geo_receipt):
    photos,source,climate,geography,geo_receipt=map(
        Path,(photos,source,climate,geography,geo_receipt))
    rec=json.loads(source.read_text(encoding="utf-8"))
    coast=json.loads(geo_receipt.read_text(encoding="utf-8"))
    check(rec.get("schema")=="mois_real1043930_original9_with_FCP447500_independent_source_reuse_v1"
          and rec.get("latest_real_photo_ID_attempts")==N
          and rec.get("strict_original9_vectors_now")==STRICT
          and rec.get("ambiguous_original9_vectors_now")==AMBIGUOUS
          and rec.get("output_photo_lake_SHA256")==PHOTO_SHA
          and sha(photos)==PHOTO_SHA,
          "Original genuine nine-palette photo/receipt SHA mismatch")
    check(sha(climate)==CLIMATE_SHA
          and coast.get("real_original_source_photo_attempts")==N
          and coast["output_geography_photo_ID_parquet_SHA256"]==sha(geography),
          "Not original source-aligned annual climatic or high-res GSHHG location")
    photo_cols=["photo_id","inat_taxon_id","species","arm","date",
                "latitude","longitude","morph","original_four_hue",
                "continuous9_strict","continuous9_exploratory_ambiguous",
                "continuous9_entropy",*COLUMNS]
    p=pd.read_parquet(photos,columns=photo_cols)
    climate_cols=["photo_id","wc2_bio12_10m","wc2_annual_vapr_10m",
                  "wc2_elev_2p5m","wc2_bio05_10m","annual_solar_Wm2"]
    c=pd.read_parquet(climate,columns=climate_cols)
    g=pd.read_parquet(geography,columns=["photo_id","land_type","landmass","area_km2"])
    check(len(p)==len(c)==len(g)==N
          and p.photo_id.astype(str).tolist()==c.photo_id.astype(str).tolist()
          ==g.photo_id.astype(str).tolist(),
          "Latitude/WorldClim/island mixed photo IDs or row ordering")
    check(p.photo_id.astype(str).is_unique
          and int(p.continuous9_strict.sum())==STRICT
          and int(p.continuous9_exploratory_ambiguous.sum())==AMBIGUOUS,
          "True photo ID original nine-comp outcome drift")
    check(int(p.morph.isin(["white","yellow_orange","red_pink","blue_purple"]).sum())==FOUR,
          "Original source four-hue photo count altered")
    P=p[list(COLUMNS)].to_numpy(dtype=float)
    admission=p.continuous9_strict.to_numpy(bool)
    ambiguous=p.continuous9_exploratory_ambiguous.to_numpy(bool)
    check(not (admission&ambiguous).any()
          and np.isfinite(P[admission|ambiguous]).all()
          and np.allclose(P[admission|ambiguous].sum(axis=1),1,atol=2e-5),
          "Photo palettes are not actual normalized valid ROI measurements")
    check(np.isnan(P[~(admission|ambiguous)]).all(),
          "A non-evaluable photo was silently assigned a palette fraction")
    largest=P.argmax(axis=1) # nan for unmeasured excluded below
    largest[~(admission|ambiguous)]=-1
    ordering=np.sort(np.where(np.isfinite(P),P,0),axis=1)
    p["main_fraction"]=np.where(admission,ordering[:,-1],np.nan)
    p["second_fraction"]=np.where(admission,ordering[:,-2],np.nan)
    names=np.array(PALETTE,dtype=object)
    p["dominant_colour"]=np.where(admission,names[np.clip(largest,0,8)],"unknown")
    p["strict"]=admission
    p["ambiguous"]=ambiguous
    p["classified4"]=p.morph.isin(["white","yellow_orange","red_pink","blue_purple"])
    p["white_fraction"]=np.where(admission,P[:,0],np.nan)
    p["bronze_fraction"]=np.where(admission,P[:,8],np.nan)
    p["pink_magenta_red_fraction"]=np.where(admission,np.sum(P[:,3:6],axis=1),np.nan)
    p["yellow_orange_fraction"]=np.where(admission,np.sum(P[:,1:3],axis=1),np.nan)
    p["blue_purple_fraction"]=np.where(admission,np.sum(P[:,6:8],axis=1),np.nan)
    p["dominant_white"]=np.where(admission,(p.dominant_colour=="white").to_numpy(float),np.nan)
    p["dominant_bronze"]=np.where(admission,(p.dominant_colour=="bronze").to_numpy(float),np.nan)
    p["mixed_top_below70"]=np.where(admission,(p.main_fraction<.7).to_numpy(float),np.nan)
    p["entropy9"]=np.where(admission,
        pd.to_numeric(p.continuous9_entropy,errors="coerce").to_numpy(float),np.nan)
    p["latitude"]=pd.to_numeric(p.latitude,errors="coerce")
    p["longitude"]=pd.to_numeric(p.longitude,errors="coerce")
    p["clim_vapour"]=pd.to_numeric(c.wc2_annual_vapr_10m,errors="coerce")
    p["clim_rain"]=pd.to_numeric(c.wc2_bio12_10m,errors="coerce")
    p["clim_elevation"]=pd.to_numeric(c.wc2_elev_2p5m,errors="coerce")
    p["clim_heat"]=pd.to_numeric(c.wc2_bio05_10m,errors="coerce")
    p["clim_solar"]=pd.to_numeric(c.annual_solar_Wm2,errors="coerce")
    p["land_type"]=g.land_type.astype(str).to_numpy()
    p["landmass"]=g.landmass.astype(str).to_numpy()
    p["area_km2"]=pd.to_numeric(g.area_km2,errors="coerce")
    p["month"]=pd.to_datetime(p.date,errors="coerce").dt.month
    p["hemisphere"]=np.select([p.latitude>=0,p.latitude<0],["north","south"],default="unknown")
    p["latitude_zone"]=[latitude_zone(v) for v in p.latitude.to_numpy(float)]
    p["zone_polar_or_other"]=np.where(p.latitude.abs()>=66.5622,"polar","nonpolar")
    p["season"]=[season_label(lat,m) for lat,m in zip(p.latitude.to_numpy(float),p.month.to_numpy(float))]
    p["altitude_band"]=[elev_band(v) for v in p.clim_elevation.to_numpy(float)]
    p["vapour_quintile"],vap_cut=q5(p.clim_vapour.to_numpy(float))
    p["rainfall_quintile"],rain_cut=q5(p.clim_rain.to_numpy(float))
    p["source_cohort"]=np.where(p.arm.astype(str).str.startswith("FCP5_"),"FCP5",p.arm)
    p["taxon"]=p.inat_taxon_id.fillna("").astype(str)
    p["genus"]=p.species.fillna("").astype(str).str.split().str[0]
    check(not p.land_type.isna().any(),"GSHHG may be unassigned, never NaN")
    return p,{"vapour_source_all_photo_quintiles":vap_cut,
              "rainfall_source_all_photo_quintiles":rain_cut,
              "input_9roi_sha256":sha(photos),"input_climate_sha256":sha(climate),
              "input_landmass_sha256":sha(geography)}

OUTCOMES=("white_fraction","entropy9","dominant_white","dominant_bronze",
          "main_fraction","second_fraction","mixed_top_below70",
          "bronze_fraction","yellow_orange_fraction",
          "pink_magenta_red_fraction","blue_purple_fraction")

def describe_stratum(d,key,value):
    original=len(d)
    x=d[d.strict]
    known=int(d.classified4.sum())
    result={"axis":key,"stratum":str(value),
            "source_photo_attempts":original,
            "known_4hue_ROIs":known,
            "strict_original9_ROIs":int(len(x)),
            "ambiguity_only_original9_ROIs":int(d.ambiguous.sum()),
            "strict9_success_pct_of_all_attempts":100*float(len(x)/original) if original else None,
            "four_hue_success_pct_of_all_attempts":100*known/original if original else None,
            "n_source_photo_taxa_strict9":int(x.taxon.nunique()),
            "n_source_photo_genera_strict9":int(x.genus.nunique()),
            "strict9_photo_based_not_wild_flora_truth":True}
    for outcome in OUTCOMES:
        vals=pd.to_numeric(x[outcome],errors="coerce").dropna()
        result[outcome+"_photo_mean_pct"]=100*float(vals.mean()) if len(vals) else None
        if not x.empty:
            taxa=x.groupby("taxon",observed=True)[outcome].mean()
            result[outcome+"_taxon_equal_pct"]=100*float(taxa.mean())
        else:
            result[outcome+"_taxon_equal_pct"]=None
    if len(x):
        counts=x.dominant_colour.value_counts()
        for colour in PALETTE:
            result["n_dominant_"+colour]=int(counts.get(colour,0))
            result["share_dominant_"+colour+"_pct"]=100*int(counts.get(colour,0))/len(x)
    else:
        for colour in PALETTE:
            result["n_dominant_"+colour]=0
            result["share_dominant_"+colour+"_pct"]=None
    return result

def profile(p):
    keys=[
      ("hemisphere",["hemisphere"]),
      ("latitude_zone",["latitude_zone"]),
      ("season_by_hemisphere",["hemisphere","season"]),
      ("season_by_latitude_zone",["latitude_zone","season"]),
      ("latitude_island",["latitude_zone","land_type"]),
      ("island_mainland",["land_type"]),
      ("moisture_vapour_quintiles",["vapour_quintile"]),
      ("rainfall_quintiles",["rainfall_quintile"]),
      ("elevation",["altitude_band"]),
      ("hemisphere_moisture",["hemisphere","vapour_quintile"]),
      ("hemisphere_elevation",["hemisphere","altitude_band"]),
      ("island_moisture",["land_type","vapour_quintile"]),
      ("island_elevation",["land_type","altitude_band"]),
      ("hemisphere_island",["hemisphere","land_type"]),
      ("season_island",["land_type","season"]),
    ]
    rows=[]
    for display,group in keys:
        for k,v in p.groupby(group,dropna=False,observed=True,sort=True):
            if not isinstance(k,tuple):k=(k,)
            name="|".join(map(str,k))
            rows.append(describe_stratum(v,display,name))
    return rows

def taxon_local_season(p):
    """Retain how frequently the exact source taxon is photographed in each season.

    Paired 5° geocell same-named-taxon summer/winter is still cross-sectional
    photo sampling, NOT genetic evolution or same-plant plasticity.
    """
    x=p.loc[p.strict & p.hemisphere.isin(["north","south"]) &
         p.season.isin(["summer","winter"]) &
         p.latitude.between(-90,90)&p.longitude.between(-180,180)].copy()
    x["lat5"]=np.floor((x.latitude+90)/5).astype("int16")
    x["lon5"]=np.floor((x.longitude+180)/5).astype("int16")
    group=["hemisphere","source_cohort","taxon","genus","lat5","lon5","season"]
    d=x.groupby(group,observed=True).agg(
        n_photos=("photo_id","size"),
        white_mean=("white_fraction","mean"),
        entropy_mean=("entropy9","mean"),
        dominant_white=("dominant_white","mean")).reset_index()
    pairs=[]
    paired={}
    for min_photo in (1,2):
        piv=d.pivot_table(index=group[:-1],columns="season",
                          values=["n_photos","white_mean","entropy_mean","dominant_white"],
                          aggfunc="first")
        good=(piv[("n_photos","summer")].ge(min_photo)&
              piv[("n_photos","winter")].ge(min_photo))
        q=piv.loc[good].copy()
        for hemi in ("north","south"):
            ix=q.index.get_level_values("hemisphere")==hemi
            v=q.loc[ix]
            response={}
            for trait in ("white_mean","entropy_mean","dominant_white"):
                if len(v):
                    difference=100*(v[(trait,"summer")]-v[(trait,"winter")])
                    response[trait+"_mean_summer_minus_winter_pp"]=float(difference.mean())
                    response[trait+"_median_pp"]=float(difference.median())
                else:
                    response[trait+"_mean_summer_minus_winter_pp"]=None
                    response[trait+"_median_pp"]=None
            paired[f"{hemi}_minimum_photo_each_season_{min_photo}"]={
                "matched_taxon_by_5deg_by_source_cohort_cells":len(v),
                "different_source_taxa":int(v.index.get_level_values("taxon").nunique()),
                "different_genera":int(v.index.get_level_values("genus").nunique()),
                **response,
                "same_named_taxon_region_NOT_individual_or_same_population":True}
    return paired

def cohort_axis(p):
    return {
        "per_source_arm":{str(k):{
            "attempted":int(len(v)),"strict":int(v.strict.sum()),
            "ambiguous":int(v.ambiguous.sum()),
            "four_hue":int(v.classified4.sum())}
            for k,v in p.groupby("source_cohort",observed=True)},
        "photographic_unknowns_not_imputed":True,
    }

def run(photo,receipt,climate,geography,geo_receipt,out):
    out=Path(out)
    check(not out.exists(),"Immutable global colour zones cannot be overwritten")
    p,cuts=main_colour_fraction_table(photo,receipt,climate,geography,geo_receipt)
    rows=profile(p)
    season=taxon_local_season(p)
    arctic=p.loc[p.latitude.ge(66.5622)]
    antarctic=p.loc[p.latitude.le(-66.5622)]
    report={
      "schema":SCHEMA,
      "source_original_photo_IDs_attempted":N,
      "genuine_four_hue_ROIs":FOUR,
      "original_strict_nine_colour_ROIs":STRICT,
      "ambiguity_nine_colour_sensitivity_only":AMBIGUOUS,
      "no_source_colour_photo_IDs_excluded_from_missing_denominator":True,
      "strict_nine_ROI_colour_proxy_coverage":STRICT/N,
      "outcome_columns":[*OUTCOMES],
      "source_pins":cuts,
      "all_crossed_geographic_climate_zone_profiles":rows,
      "matched_source_taxon_5degree_summer_winter":season,
      "source_cohort_counts":cohort_axis(p),
      "northern_arctic_original_photo_attempts":int(len(arctic)),
      "northern_arctic_original_strict9":int(arctic.strict.sum()),
      "southern_antarctic_original_photo_attempts":int(len(antarctic)),
      "southern_antarctic_original_strict9":int(antarctic.strict.sum()),
      "unavailable_from_saved_source":{
          "major_flower_roi_CIELAB_L_star":True,
          "major_flower_roi_CIELAB_C_star":True,
          "major_flower_roi_CIELAB_hue":True,
          "green_colour_palette_component":True,
          "botanically_verified_major_petal_colour":True,
      },
      "frozen_srgb_device_lab_pilot_is_not_global1043930_Lab":True,
      "landmass_antarctica_not_automatically_island":True,
      "white_vs_dark_or_green_not_one_drab_bin":True,
      "source_photo_image_downloads":0,
      "original_ROI_model_reinferences":0,
    }
    out.mkdir(parents=True)
    outfile=out/"mois1043930_global_zonal_white_dominant_palette_profiles.csv"
    pd.DataFrame(rows).to_csv(outfile,index=False)
    (out/"mois1043930_global_zonal_white_colour_result.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    with (out/"mois1043930_global_zonal_run_receipt.json").open("w") as f:
        f.write(json.dumps({
          "schema":SCHEMA,
          "full_photo_attempts":N,"strict_colour":STRICT,"ambiguous":AMBIGUOUS,
          "zonal_profiles":len(rows),
          "northern_arctic_strict_photo9":int(arctic.strict.sum()),
          "southern_antarctic_strict_photo9":int(antarctic.strict.sum()),
          "zonal_csv_sha256":sha(outfile),
          "result_json_sha256":sha(out/"mois1043930_global_zonal_white_colour_result.json"),
        },indent=2,sort_keys=True)+"\n")
    wanted={
      "source":N,"strict":STRICT,"profiles":len(rows),
      "north_arctic":describe_stratum(arctic,"polar_photo","Arctic_geographic"),
      "south_antarctic":describe_stratum(antarctic,"polar_photo","Antarctic_geographic"),
      "hemispheres":[x for x in rows if x["axis"]=="hemisphere"],
      "island_mainland":[x for x in rows if x["axis"]=="island_mainland"],
      "season_by_hemisphere":[x for x in rows if x["axis"]=="season_by_hemisphere"],
      "source_taxon_paired_seasons":season,
    }
    for name,items in wanted.items():
        if isinstance(items,list):
            wanted[name]=[{k:v for k,v in a.items() if k in (
              "stratum","source_photo_attempts","strict_original9_ROIs",
              "strict9_success_pct_of_all_attempts",
              "white_fraction_photo_mean_pct","white_fraction_taxon_equal_pct",
              "entropy9_photo_mean_pct","main_fraction_photo_mean_pct",
              "share_dominant_white_pct","share_dominant_bronze_pct",
              "share_dominant_pink_pct","share_dominant_yellow_pct"
            )} for a in items]
        elif name in ("north_arctic","south_antarctic"):
            wanted[name]={k:v for k,v in items.items() if k in (
              "stratum","source_photo_attempts","strict_original9_ROIs",
              "white_fraction_photo_mean_pct","white_fraction_taxon_equal_pct",
              "entropy9_photo_mean_pct","n_source_photo_taxa_strict9",
              "share_dominant_white_pct","share_dominant_bronze_pct")}
    print("MOIS_LATEST1043930_GLOBAL_WHITECOLOR_LAT_POLAR_SEASON_ISLAND_REAL_RESULTS",
        json.dumps(wanted,sort_keys=True,allow_nan=False),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("photo","receipt","climate","geography","geo-receipt","output"):
        p.add_argument("--"+key,required=True,type=Path)
    a=p.parse_args()
    run(a.photo,a.receipt,a.climate,a.geography,a.geo_receipt,a.output)

if __name__=="__main__":main()
