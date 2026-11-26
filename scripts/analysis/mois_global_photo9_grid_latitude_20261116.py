#!/usr/bin/env python3
"""Species- and genus-controlled MOIS photographic hues vs REAL GBIF global grid.

Sources: exact signed November 2026 anonymized GBIF 10,346-file S3 full snapshot
and signed COL XR 37,104 accepted species with nine measured photo ROI fractions.

No GBIF coordinate is treated as a local petal hue, a native-only range or a
biological pigment trait. GBIF records count is only a proxy for observation
effort confounded with geographic range.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

GLOBAL_ZIP_SHA="422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
SOURCE_ZIP_SHA="abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a"
GLOBAL_GRID_SHA="e00b0408055984abfe4ec0ff0a24833fc148888afd1ea0e0c9b270a077d08b11"
GLOBAL_COVERAGE_SHA="6ce4ac5289dccfdfce3dae946de3b1617bbc646eba76de75086c59795e5a4603"
SOURCE_TRAITS_SHA="b06d8980860ab11d1e9cff1af7dd1d92b3bfa9de4f2d227bc7e58809f931663a"
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
BANDS=tuple(range(-90,90,10))
GRID_METRES=100000

def require(ok,msg):
    if not bool(ok): raise ValueError("MOIS_GBIF_GLOBAL_ECOGEOGRAPHY: "+str(msg))

def sha(x):
    if isinstance(x,(str,Path)):x=Path(x).read_bytes()
    return hashlib.sha256(x).hexdigest()

def lat_band_floor(value):
    if not np.isfinite(value):return float("nan")
    return max(-90,min(80,math.floor(float(value)/10)*10))

def within_group(frame,groups):
    return frame-frame.groupby(groups).transform("mean")

def photo_mean_top(values):
    v=np.asarray(values,dtype=float)
    require(v.shape==(9,) and np.all(np.isfinite(v))
            and (v>=0).all() and abs(v.sum()-1)<1e-3,
            "invalid photo9 numeric mean vector")
    return HUES[int(v.argmax())]

def signed_data(grid_source,accepted_source):
    import pyarrow.parquet as pq
    require(sha(grid_source)==GLOBAL_ZIP_SHA and sha(accepted_source)==SOURCE_ZIP_SHA,
            "signed GBIF source artifact changed")
    with ZipFile(grid_source) as z:
        rawcov=z.read("all_37104_accepted_species_grid_coverage_and_photo9.csv.gz")
        rawgrid=z.read("all_37104_species_x_EPSG6933_100km_cells.parquet")
        source_receipt=json.loads(z.read("receipt.json"))
    with ZipFile(accepted_source) as z:
        rawtraits=z.read("current_exact_accepted_COLXR_numeric9_traits.csv.gz")
    require(sha(rawcov)==GLOBAL_COVERAGE_SHA and sha(rawgrid)==GLOBAL_GRID_SHA
            and sha(rawtraits)==SOURCE_TRAITS_SHA,"internal signed sources changed")
    require(source_receipt["all_16_shards_source_verified"]
            and source_receipt["complete_public_S3_object_count"]==10346
            and source_receipt["GBIF_occurrences_passed_all_filters"]==396142825,
            "source S3 global completeness not certified")
    c=pd.read_csv(io.BytesIO(rawcov),compression="gzip",
                  dtype={"key":str,"inat_taxon_id":str,"specieskey":str})
    t=pd.read_csv(io.BytesIO(rawtraits),compression="gzip",
                  dtype={"gbif_colxr_accepted_species_key":str,"inat_taxon_id":str})
    require(len(c)==len(t)==37104 and c.key.is_unique
            and t.gbif_colxr_accepted_species_key.is_unique,
            "accepted species roster incomplete")
    keep=["gbif_colxr_accepted_species_key","inat_taxon_id","photo_source_species_name",
          "photo_evidence_status","photo_modal_palette9",
          "independent_roi9_observations","after_strict9",
          *("observed_mean_"+p for p in HUES)]
    u=t[keep].rename(columns={"inat_taxon_id":"photo_taxon_id"})
    d=c.merge(u,left_on="key",right_on="gbif_colxr_accepted_species_key",
              how="left",validate="one_to_one",indicator=True)
    require(d._merge.eq("both").all()
        and d.inat_taxon_id.eq(d.photo_taxon_id).all()
        and d.source_hue9.eq(d.photo_modal_palette9).all()
        and d.hue_status.eq(d.photo_evidence_status).all(),
        "one-to-one GBIF photo trait join not exact")
    mean=d[["observed_mean_"+p for p in HUES]].to_numpy(float)
    require(np.isfinite(mean).all() and (mean>=0).all()
            and np.allclose(mean.sum(axis=1),1,atol=5e-5),
            "9 photo ROI fractions nonfinite or not compositional")
    d["mean9_top"]=np.asarray(HUES)[mean.argmax(axis=1)]
    d["modal_mean9_disagree"]=d["mean9_top"].ne(d.source_hue9)
    require(int(d.modal_mean9_disagree.sum())==2466,
            "signed 37,104 modality-vs-mean mismatch changed")
    require(d.occupied_grid_cells.gt(0).sum()==36758
            and d.records_passed_spatial_QA.sum()==396142825,
            "full 37,104 species actual GBIF occupancy mismatch")
    d["genus"]=d.species.str.split().str[0]
    d["white_photo"]=(d.source_hue9=="white").astype(float)
    d["log_occ"]=np.log1p(d.records_passed_spatial_QA)
    d["log_cells"]=np.log1p(d.occupied_grid_cells)
    d=d.drop(columns="_merge")
    table=pq.read_table(io.BytesIO(rawgrid),
             columns=["specieskey","grid_x","grid_y","occurrence_count"]).to_pandas()
    require(len(table)==5493931
            and table.occurrence_count.sum()==396142825
            and table.groupby(["specieskey","grid_x","grid_y"]).size().max()==1,
            "5,493,931 distinct GBIF species-cell pairs or source records changed")
    require(set(table.specieskey)==set(d.loc[d.occupied_grid_cells.gt(0),"key"]),
            "source species × grid keys differ from accepted species coverage")
    group=table.groupby("specieskey",sort=False).agg(
        cells=("grid_y","size"),occurrence_count=("occurrence_count","sum"))
    dpos=d.set_index("key").loc[group.index]
    require(np.array_equal(dpos.occupied_grid_cells.to_numpy(),group.cells.to_numpy())
            and np.array_equal(dpos.records_passed_spatial_QA.to_numpy(),
                               group.occurrence_count.to_numpy()),
            "species-cell totals are not consistent with 37,104 source atlas")
    return d,table

def model(d,yname,within=False,genus_equal=False,nine=False,mean_white=False):
    y=d[yname].astype(float).copy()
    x=pd.DataFrame(index=d.index)
    x["log_occ"]=d.log_occ.to_numpy(float)
    if mean_white:
        x["photo_mean_white_fraction"]=d.observed_mean_white.to_numpy(float)
    elif nine:
        for hue in HUES[1:]:
            x["vs_white_"+hue]=d.source_hue9.eq(hue).astype(float)
    else:x["white_photo"]=d.white_photo.to_numpy(float)
    states=sorted(d.hue_status.unique())
    for status in states[1:]:
        x["source_status_"+status]=d.hue_status.eq(status).astype(float)
    if within:
        groups=d.genus
        x=within_group(x,groups)
        y=within_group(y.to_frame("response"),groups)["response"]
        if genus_equal:
            weights=1/d.groupby("genus")["genus"].transform("size").to_numpy(float)
            m=sm.WLS(y,x,weights=weights).fit(
                cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
        else:
            m=sm.OLS(y,x).fit(cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
    else:
        x=sm.add_constant(x)
        m=sm.OLS(y,x).fit(cov_type="cluster",
                           cov_kwds={"groups":d.genus},use_t=False)
    return m

def coef(m,n):
    a=m.conf_int().loc[n]
    return {"beta":float(m.params[n]),"CI95":[float(a.iloc[0]),float(a.iloc[1])],
            "p":float(m.pvalues[n])}

def run(*,grid_source,accepted_source,output):
    from pyproj import Transformer
    out=Path(output);require(not out.exists(),"do not overwrite published audit")
    d,grid=signed_data(grid_source,accepted_source)
    positive=d.loc[d.occupied_grid_cells.gt(0)].copy()
    require(len(positive)==36758,"positive species denominator changed")
    grid["grid_y"]=grid["grid_y"].astype("int32")
    unique_y=sorted(grid.grid_y.unique())
    _,lat=Transformer.from_crs("EPSG:6933","EPSG:4326",
        always_xy=True).transform([0.0]*len(unique_y),[(y+.5)*GRID_METRES for y in unique_y])
    lat=dict(zip(unique_y,lat))
    grid["approx_cell_center_latitude"]=grid.grid_y.map(lat).astype(float)
    # A 100km raster cell's centre need not lie inside projection limits at
    # the extreme polar edge. Never convert a missing latitude to 0°.
    ok=np.isfinite(grid.approx_cell_center_latitude)
    invalid_lat=int((~ok).sum())
    geo=grid.loc[ok].copy()
    geo["abs_lat"]=geo.approx_cell_center_latitude.abs()
    geo["latitude_10degree_lower"]=np.floor(
        geo.approx_cell_center_latitude/10).clip(-9,8).astype(int)*10
    geo["latitude_10degree_lower"]=geo.latitude_10degree_lower.clip(-90,80)
    latmean=geo.groupby("specieskey").agg(
        mean_abs_latitude_over_occupied_cells=("abs_lat","mean"),
        mean_signed_latitude_over_occupied_cells=("approx_cell_center_latitude","mean"),
        occupied_cells_with_valid_lat=("grid_y","size")).reset_index()
    d=d.merge(latmean,left_on="key",right_on="specieskey",
              how="left",validate="one_to_one",suffixes=("","_lat"))
    d["lat_model_valid"]=d.occupied_cells_with_valid_lat.gt(0)
    # Each species contributes once per occupied 10-degree band. This
    # is a photographic species-hue label copied across its GBIF band presences,
    # not local colour measurement.
    b=(geo[["specieskey","latitude_10degree_lower"]].drop_duplicates()
       .merge(d[["key","source_hue9","genus"]],left_on="specieskey",
              right_on="key",how="left",validate="many_to_one"))
    require(b.source_hue9.notna().all(),"grid band source species trait missing")
    byband=b.groupby(["latitude_10degree_lower","source_hue9"]).size().unstack(fill_value=0)
    byband=byband.reindex(index=BANDS,columns=HUES,fill_value=0)
    band_report=byband.reset_index()
    band_report["unique_species_band_total"]=byband.sum(axis=1).to_numpy()
    band_report["white_photo_fraction_of_observed_species"]=(
        byband.white/byband.sum(axis=1).replace(0,np.nan)).to_numpy(float)
    genus_band=b.groupby(["latitude_10degree_lower","genus"]).agg(
        white_photo_fraction_in_genus=("source_hue9",lambda x:float(x.eq("white").mean())))
    equal=genus_band.groupby("latitude_10degree_lower").white_photo_fraction_in_genus.mean()
    band_report["genus_equal_white_photo_fraction"]=band_report.latitude_10degree_lower.map(equal)
    band_report["photo_hue_source_is_local_at_GBIF_grid"]=False

    # Exclude species with no usable georeferenced GBIF record from
    # conditional range models. Their absence is NOT a biological zero range.
    mixed_genus=positive.groupby("genus").white_photo.nunique().eq(2)
    mixed=positive.loc[positive.genus.map(mixed_genus)].copy()
    stable=mixed.loc[~mixed.modal_mean9_disagree].copy()
    stable=stable.loc[stable.genus.map(stable.groupby("genus").white_photo.nunique().eq(2))]
    white_cluster=model(positive,"log_cells")
    white_genus=model(positive,"log_cells",within=True)
    white_equal=model(mixed,"log_cells",within=True,genus_equal=True)
    white_stable=model(stable,"log_cells",within=True,genus_equal=True)
    fractions=model(positive,"log_cells",within=True,mean_white=True)
    nine=model(positive,"log_cells",within=True,nine=True)
    names=[n for n in nine.params.index if n.startswith("vs_white_")]
    R=np.zeros((len(names),len(nine.params)))
    for i,n in enumerate(names):R[i,list(nine.params.index).index(n)]=1
    wald=nine.wald_test(R,scalar=True)
    q=multipletests([nine.pvalues[n] for n in names],method="fdr_bh")[1]
    nine_details=[{"hue":n.removeprefix("vs_white_"),
                   **coef(nine,n),"q_BH_eight_contrasts":float(fdr)}
                  for n,fdr in zip(names,q)]
    latitude=d.loc[d.lat_model_valid].copy()
    lati=model(latitude,"mean_abs_latitude_over_occupied_cells",within=True)
    lati_equal=model(latitude.loc[latitude.genus.map(
        latitude.groupby("genus").white_photo.nunique().eq(2))],
        "mean_abs_latitude_over_occupied_cells",within=True,genus_equal=True)

    hues=d.groupby("source_hue9").agg(accepted_species=("key","size"),
        with_recorded_grid_cells=("occupied_grid_cells",lambda x:int((x>0).sum())),
        median_species_cells_all=("occupied_grid_cells","median"),
        total_GBIF_qualifying_occurrences=("records_passed_spatial_QA","sum")).reset_index()
    descript=positive.groupby("source_hue9").agg(
        median_species_100km_cells=("occupied_grid_cells","median"),
        median_GBIF_record_count=("records_passed_spatial_QA","median"),
        median_recorded_countries=("observed_countries","median")).reset_index()
    hues=hues.merge(descript,on="source_hue9",validate="one_to_one")
    receipt={
      "schema":"mois_GBIF_37104_9photo_hues_global_range_latitude_genus_FE_v1",
      "source_full_10346_S3_Parquets_and_16_shard_validated":True,
      "GBIF_source_snapshot":"2026-11-01",
      "source_global_GBIF_ZIP_SHA256":GLOBAL_ZIP_SHA,
      "source_COLXR_37104_traits_ZIP_SHA256":SOURCE_ZIP_SHA,
      "species_accepted_total":37104,
      "species_with_observed_valid_GBIF_range":len(positive),
      "species_without_valid_GBIF_record":int(d.occupied_grid_cells.eq(0).sum()),
      "source_species_grid_cells":len(grid),
      "source_occurrences_after_GBIF_QA":396142825,
      "source_photo_palette9_modal_vs_mean_disagreements":int(d.modal_mean9_disagree.sum()),
      "genera_in_positive_species":int(positive.genus.nunique()),
      "genera_with_both_white_and_nonwhite":int(mixed.genus.nunique()),
      "species_in_mixed_hue_genera":len(mixed),
      "white_photo_range_median_grid_cells":float(positive.loc[
          positive.white_photo.eq(1),"occupied_grid_cells"].median()),
      "nonwhite_photo_range_median_grid_cells":float(positive.loc[
          positive.white_photo.eq(0),"occupied_grid_cells"].median()),
      "statistical_models":{
          "white_country_unadjusted_not_analysed":True,
          "range_log_cells_photo_white_effort_status_genus_cluster":
              coef(white_cluster,"white_photo"),
          "range_within_genus_white_effect_genus_cluster":
              coef(white_genus,"white_photo"),
          "range_within_mixed_genus_equal_weight_white_effect":
              coef(white_equal,"white_photo"),
          "range_within_mixed_genus_equal_exclude_modal_mean_disagree":
              coef(white_stable,"white_photo"),
          "range_within_genus_mean9_white_fraction":
              coef(fractions,"photo_mean_white_fraction"),
          "range_9photo_hues_within_genus_joint_Wald":
              {"chi2":float(wald.statistic),"df":len(names),"p":float(wald.pvalue)},
          "range_9photo_hues_vs_white_within_genus":nine_details,
          "mean_absolute_grid_latitude_white_within_genus_degrees":
              coef(lati,"white_photo"),
          "mean_absolute_grid_latitude_white_genus_equal_degrees":
              coef(lati_equal,"white_photo")
      },
      "latitude_derived_from_EPSG6933_cell_center_not_GBIF_point":True,
      "latitude_cell_centers_invalid_outside_projection":invalid_lat,
      "individual_species_cells_in_18_latitude_bands":
          len(b),
      "global_37104_photo_hues_botanical_petal_confirmed":False,
      "local_GBIF_occurrence_hue_measured":False,
      "native_wild_range_verified":False,
      "species_genetic_fixed_hue_verified":False,
      "genus_FE_is_not_full_phylogenetic_correction":True,
      "GBIF_occurrence_count_is_effort_and_range_proxy":True,
      "no_climate_soil_elevation_fitted":True,
      "grid_type":"EPSG:6933 100km square equal area, not GBIF ISEA3H r8"
    }
    out.mkdir(parents=True)
    d[["key","inat_taxon_id","species","genus","source_hue9","hue_status",
       "mean9_top","modal_mean9_disagree","observed_mean_white",
       "occupied_grid_cells","records_passed_spatial_QA","observed_countries",
       "mean_abs_latitude_over_occupied_cells","mean_signed_latitude_over_occupied_cells",
       "occupied_cells_with_valid_lat"]].to_csv(
        out/"all_37104_accepted_species_global_GBIF_photo9_georange_latitudes.csv.gz",
        index=False,compression={"method":"gzip","mtime":0})
    hues.to_csv(out/"global_photo9_hue_100km_occupancy_summary.csv",index=False)
    band_report.to_csv(out/"global_10degree_latitude_photo_hue_species_presence.csv",index=False)
    pd.DataFrame(nine_details).to_csv(out/"within_genus_9hue_vs_white_FDR.csv",index=False)
    receipt["output_SHA256"]={p.name:sha(p) for p in out.iterdir()}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_GLOBAL_37104_9PHOTO_HUE_GBIF_RANGE_LATITUDE_MODELS",
      json.dumps({k:receipt[k] for k in (
        "species_accepted_total","species_with_observed_valid_GBIF_range",
        "genera_with_both_white_and_nonwhite",
        "latitude_cell_centers_invalid_outside_projection",
        "source_species_grid_cells","statistical_models")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--GBIF-global-ZIP",required=True)
    p.add_argument("--COLXR-photo-ZIP",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    run(grid_source=a.GBIF_global_ZIP,accepted_source=a.COLXR_photo_ZIP,
        output=a.output)
