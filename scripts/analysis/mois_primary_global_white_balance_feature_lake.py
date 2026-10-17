#!/usr/bin/env python3
"""Cache genuine MOIS ROI hues and environmental covariates for global ecology.

Near-half WHITE versus ALL coloured flowers is descriptive, not a 50:50
evolutionary equilibrium or evidence of pigment cost/fitness. Keep unknown ROIs
unknown, keep season-enriched source out of the primary global denominator.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd

HUES = ("white","yellow_orange","red_pink","blue_purple")
N, CLASS, DATED, PARTIAL = 309439,149662,149653,248312
SHA_PHOTO="546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"
SHA_CLIMATE="e8cf9ab1c935b906731b6c783a45193c9c13cca8cb226fe46c4335d0e1c75cee"
SHA_STATIC="43dcd2608f5a25e0a5d45d35810604a62ab659fb2cb9fb1a4e4257ef756d56c3"
P_COLS=("photo_id","inat_taxon_id","species","morph","observation_id",
        "qualified_date","old_mois_or_new_fcp","latitude","longitude",
        "wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay",
        "soil_available_water_proxy")
ENV=("wc2_bio05_10m","wc2_bio12_10m","wc2_bio14_10m","wc2_bio15_10m",
     "wc2_annual_srad_normal_10m_Wm2")
E_COLS=("photo_id","inat_taxon_id","species","morph","latitude","longitude",*ENV)
TERRAIN=("wc2_elev_2p5m","terrain_slope_2p5m_degree","terrain_TRI_2p5m_m",
         "terrain_TPI_2p5m_m","terrain_roughness_2p5m_m")
S_COLS=("photo_id","inat_taxon_id","species","morph",*TERRAIN)
SOIL=("wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay",
      "soil_available_water_proxy")
P_OUT="primary309439_photo_hue_source.parquet"
E_OUT="primary309439_annual_climate_source.parquet"
S_OUT="primary309439_partial_terrain_soil_source.parquet"
STATS="primary309439_white_chromatic_balance.json"
DECILES="primary309439_abiotic_decile_balance.csv"
MAP="primary309439_latitude_white_balance.csv"
MANIFEST="primary309439_immutable_photo_feature_lake_manifest.json"


def require(ok, why):
    if not ok: raise ValueError(why)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()


def source(p_path,e_path,s_path):
    paths=(p_path,e_path,s_path)
    for path, expected in zip(paths,(SHA_PHOTO,SHA_CLIMATE,SHA_STATIC)):
        require(Path(path).is_file() and sha(path)==expected,
                "Source actual ROI photo or original environment SHA drift")
    p=pd.read_csv(p_path,compression="gzip",dtype=str,keep_default_na=False,
                  usecols=list(P_COLS),low_memory=False)
    e=pd.read_csv(e_path,compression="gzip",dtype=str,keep_default_na=False,
                  usecols=list(E_COLS),low_memory=False)
    s=pd.read_csv(s_path,compression="gzip",dtype=str,keep_default_na=False,
                  usecols=list(S_COLS),low_memory=False)
    require(len(p)==len(e)==N and len(s)==PARTIAL and
            all(x.photo_id.is_unique for x in (p,e,s)),
            "Original source photo counts or photo IDs differ")
    known=p.morph.isin(HUES)
    require(int(known.sum())==CLASS and
            int((known&p.qualified_date.str.fullmatch(
                r"\d{4}-\d{2}-\d{2}")).sum())==DATED,
            "Original classified photo or recorded observed_on differs")
    for title,x in (("climate",e),("static",s)):
        matched=p[["photo_id","inat_taxon_id","species","morph"]].merge(
            x[["photo_id","inat_taxon_id","species","morph"]],on="photo_id",
            suffixes=("_a","_b"),validate="one_to_one")
        require(len(matched)==len(x),
                "Unmatched site-feature original photo IDs: "+title)
        for field in ("inat_taxon_id","species","morph"):
            require(matched[field+"_a"].eq(matched[field+"_b"]).all(),
                    "ROI hue or named taxon changed across feature join")
    m=p.merge(e.drop(columns=["inat_taxon_id","species","morph"]),
              on="photo_id",validate="one_to_one",suffixes=("","_env"))
    for fld in ("latitude","longitude"):
        aa=pd.to_numeric(m[fld],errors="coerce").to_numpy(float)
        bb=pd.to_numeric(m[fld+"_env"],errors="coerce").to_numpy(float)
        overlap=np.isfinite(aa)&np.isfinite(bb)
        require(np.all(np.abs(aa[overlap]-bb[overlap])<=0.001),
                "Old photo site position differs from annual WorldClim site")
    m=m.drop(columns=["latitude_env","longitude_env"]).merge(
        s.drop(columns=["inat_taxon_id","species","morph"]),
        on="photo_id",how="left",validate="one_to_one")
    require(len(m)==N and m.photo_id.is_unique,
            "Environmental photo join changed original measured denominator")
    for k in (*ENV,*TERRAIN,*SOIL,"latitude","longitude"):
        m[k]=pd.to_numeric(m[k],errors="coerce")
    # The pinned historical 309439 annual srad column is original
    # WorldClim kJ m-2 day-1 despite its earlier Wm2 name.
    m["annual_solar_Wm2"]=m["wc2_annual_srad_normal_10m_Wm2"]/86.4
    m=m.drop(columns=["wc2_annual_srad_normal_10m_Wm2"])
    m["classified"]=m.morph.isin(HUES)
    m["white"]=m.morph.eq("white")
    return m


def fractions(d):
    n=len(d);c=int(d.classified.sum());w=int(d.white.sum())
    return {"photo_attempts":n,"classifiable":c,"white":w,
            "chromatic":c-w,"unknown":n-c,
            "white_fraction":w/c if c else None,
            "classification_yield":c/n if n else None}


def genus(name):
    if not isinstance(name,str) or not name.strip(): return None
    t=name.split()[0]
    return t if re.fullmatch("[A-Z][A-Za-z]+",t) else None


def balances(d,*,draws=999):
    names=d[d.classified & d.inat_taxon_id.ne("")].copy()
    names["genus"]=names.species.map(genus)
    species=names.groupby("inat_taxon_id").agg(
        p=("white","mean"),n=("photo_id","size"),genus=("genus","first"))
    genera=species.dropna(subset=["genus"]).groupby("genus").p.mean()
    require(len(species)>0 and len(genera)>0,"Missing photographed taxa/genera")
    g=genera.to_numpy(float)
    rng=np.random.default_rng(20261017)
    b=np.array([g[rng.integers(0,len(g),len(g))].mean()
                for _ in range(draws)])
    rates=fractions(d)
    nearhalf={ "photos":rates["white_fraction"],
               "species_equal":float(species.p.mean()),
               "genus_equal":float(g.mean())}
    return {
        "photo_rates":rates,
        "near_half_rather_than_exact_50_50_fraction":nearhalf,
        "distance_from_50_percent_in_percentage_points":
            {k:100*(v-.5) for k,v in nearhalf.items()},
        "within_5_percentage_points_of_half_descriptive":
            {k:bool(abs(v-.5)<=.05) for k,v in nearhalf.items()},
        "genus_bootstrap_percentile95_fraction":
            [float(x) for x in np.percentile(b,[2.5,97.5])],
        "named_taxa_with_classified_photos":int(len(species)),
        "nominal_genera_with_classified_taxa":int(len(genera)),
        "taxa_with_only_one_classified_photo":int(species.n.eq(1).sum()),
        "unknown_photos_identification_bounds_total_original_attempts":
            [rates["white"]/len(d),(rates["white"]+rates["unknown"])/len(d)],
        "balance_is_not_maintained_50_50_equilibrium":True,
        "reproductive_fitness_pigment_cost_not_observed":True,
    }



def balance_sampling_sensitivity(d):
    """Separate photo-depth and genus-richness effects on apparent near-half."""
    classified=d[d.classified & d.inat_taxon_id.ne("")].copy()
    classified["genus"]=classified.species.map(genus)
    taxa=classified.groupby("inat_taxon_id").agg(
        white_fraction=("white","mean"),photos=("photo_id","size"),
        genus=("genus","first"))
    summary=[]
    for min_photo in (1,2,5):
        eligible=taxa[(taxa.photos>=min_photo)&taxa.genus.notna()]
        genera=eligible.groupby("genus").white_fraction.agg(["mean","size"])
        for min_taxa in (1,2,5):
            g=genera[genera["size"]>=min_taxa]
            selected=eligible[eligible.genus.isin(g.index)]
            summary.append({
                "min_classified_original_photos_per_taxon":min_photo,
                "min_named_taxa_per_genus":min_taxa,
                "classifiable_named_taxa":int(len(selected)),
                "eligible_genera":int(len(g)),
                "taxon_equal_white_fraction":float(
                    selected.white_fraction.mean()) if len(selected) else None,
                "genus_equal_white_fraction":
                    float(g["mean"].mean()) if len(g) else None,
                "genus_equal_distance_from_half_pp":
                    100*(float(g["mean"].mean())-.5) if len(g) else None,
            })
    return summary


def geographic_cell_sensitivity(d, *, degrees=10, min_classified=30):
    """Spatially equal cells are a bias diagnostic, not area/plant weighting."""
    x=d[d.classified & d.latitude.between(-90,90) &
        d.longitude.between(-180,180)].copy()
    x["latcell"]=np.floor((x.latitude+90)/degrees).astype(int)
    x["loncell"]=np.floor((x.longitude+180)/degrees).astype(int)
    cells=x.groupby(["latcell","loncell"]).white.agg(["size","mean"])
    cells=cells.loc[cells["size"]>=min_classified]
    return {
        "degrees_per_geographic_cell":degrees,
        "minimum_classified_photos_per_cell":min_classified,
        "eligible_cells":int(len(cells)),
        "photos_in_eligible_cells":int(cells["size"].sum()),
        "equal_cell_mean_white_fraction":
            float(cells["mean"].mean()) if len(cells) else None,
        "cell_white_rate_5th_50th_95th_percentile":
            [float(x) for x in np.percentile(cells["mean"],[5,50,95])]
            if len(cells) else None,
        "not_area_weighted_or_representative_biological_population":True,
    }


def deciles(d):
    variables=(*ENV[:-1],"annual_solar_Wm2")
    results=[]
    for v in variables:
        q=d.dropna(subset=[v]).copy()
        if len(q)<20 or q[v].nunique()<10:continue
        q["decile"]=pd.qcut(q[v],10,duplicates="drop",labels=False)
        for k,group in q.groupby("decile",observed=True):
            cls=group[group.classified]
            taxa=cls.groupby("inat_taxon_id").white.mean()
            results.append({
                "variable":v,"decile":int(k),
                "environment_median":float(group[v].median()),
                **fractions(group),
                "taxon_equal_photo_white_fraction":
                    float(taxa.mean()) if len(taxa) else None,
                "named_taxa_with_known_hue":len(taxa),
                "annual_climatic_1970_2000_not_observation_weather":True,
            })
    return pd.DataFrame(results)


def latitude_table(d):
    x=d[d.latitude.between(-90,90)].copy()
    x["hemisphere"]=np.where(x.latitude.ge(0),"north","south")
    x["abs_lat_band"]=pd.cut(x.latitude.abs(),[0,30,60,90.001],
        right=False,labels=["0-30","30-60","60-90"])
    result=[]
    for (band,hemi),g in x.groupby(["abs_lat_band","hemisphere"],observed=True):
        result.append({"band":str(band),"hemisphere":hemi,**fractions(g)})
    return pd.DataFrame(result)


def run(photo,climate,static,out):
    require(not out.exists() or not any(out.iterdir()),
            "Existing immutable feature lake must never be overwritten")
    d=source(photo,climate,static)
    summary=balances(d)
    summary["taxon_and_genus_depth_sensitivity"]=balance_sampling_sensitivity(d)
    summary["geographic_equal_10degree_cell_sensitivity"]=geographic_cell_sensitivity(d)
    summary.update({
        "schema":"mois_309439_world_photographic_white_chromatic_nearhalf_v1",
        "known_hue_counts":{h:int(d.morph.eq(h).sum()) for h in HUES},
        "photo_RGB_hue_not_true_flower_organ_or_UV_pigment":True,
        "missing_photo_hues_are_not_nonwhite":True,
        "separate_targeted256_season_photos_not_in_global_balance":True,
        "static_soil_terrain_source_only_248312_partially_missing":True,
        "annual_climate_not_observed_year_month":True,
        "fitness_cost_benefit_are_hypotheses_not_estimated":True,
    })
    out.mkdir(parents=True,exist_ok=True)
    outputs={
        P_OUT:d[["photo_id","inat_taxon_id","species","morph",
                 "observation_id","qualified_date","old_mois_or_new_fcp",
                 "latitude","longitude","classified","white"]],
        E_OUT:d[["photo_id",*ENV[:-1],"annual_solar_Wm2"]],
        S_OUT:d[["photo_id",*SOIL,*TERRAIN]],
    }
    for name,table in outputs.items():
        table.to_parquet(out/name,index=False,engine="pyarrow",compression="zstd")
        check=pd.read_parquet(out/name,columns=["photo_id"])
        require(check.photo_id.equals(d.photo_id.reset_index(drop=True)),
                "Saved Parquet no longer keyed in original photo ID order")
    (out/STATS).write_text(json.dumps(summary,sort_keys=True,indent=2,
        ensure_ascii=False,allow_nan=False)+"\n")
    deciles(d).to_csv(out/DECILES,index=False,lineterminator="\n")
    latitude_table(d).to_csv(out/MAP,index=False,lineterminator="\n")
    files=(P_OUT,E_OUT,S_OUT,STATS,DECILES,MAP)
    receipt={
        "schema":"mois_immutable_309439_photo_RGB_environment_Parquet_v1",
        "original_photo_SHA":SHA_PHOTO,
        "original_annual_climate_SHA":SHA_CLIMATE,
        "original_partial_248312_terrain_SHA":SHA_STATIC,
        "original_measured_RGB_photo_IDs":N,
        "known_four_hue_original_ROIs":CLASS,
        "original_observation_day_and_known_hue":DATED,
        "all_photo_attempts_retained_in_parquet":True,
        "soil_terrain_not_extrapolated_to_unobserved_sites":True,
        "raw_licensed_image_pixels_not_published":True,
        "annual_solar_source_unit_corrected_from_kJ_to_Wm2":True,
        "targeted_season_sampling_excluded_from_primary_world":True,
        "global_50_50_adaptive_equilibrium_test_not_run":True,
        "files_SHA256":{file:sha(out/file) for file in files},
    }
    (out/MANIFEST).write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_309439_WORLD_WHITE_CHROMATIC_BALANCE",
          json.dumps({"photo_count":N,"label_count":CLASS,
                      "hues":summary["known_hue_counts"],
                      "fractions":summary[
                          "near_half_rather_than_exact_50_50_fraction"],
                      "taxa":summary["named_taxa_with_classified_photos"],
                      "genera":summary["nominal_genera_with_classified_taxa"],
                      "genus_bootstrap95":summary["genus_bootstrap_percentile95_fraction"],
                      "depth_gate_sensitivity":summary["taxon_and_genus_depth_sensitivity"],
                      "equal_10degree_cells":summary["geographic_equal_10degree_cell_sensitivity"]},
                     sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-photo",type=Path,required=True)
    ap.add_argument("--original-annual-climate",type=Path,required=True)
    ap.add_argument("--partial-soil-terrain",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    run(a.original_photo,a.original_annual_climate,
        a.partial_soil_terrain,a.output)
