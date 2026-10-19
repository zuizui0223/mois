#!/usr/bin/env python3
"""Join actual MOIS photo sites to TDWG3 polygons and Kew WCVP *regional* range status.

A native_at_TDWG3 record is a statement about an accepted taxon's broad
botanical country, NOT proof that a photographed individual is wild/native,
that a tiny island population is native, or that flower colour is botanically
verified. Preserve geographic ambiguity, source photo uncertainty, unmatched
names, native/introduced uncertainty and original ROI unknowns.

Inputs are individually SHA-authenticated releases. Spatial boundaries are
TDWG/wgsrpd level3.geojson at an immutable Git commit, not guessed country names.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import make_valid, union_all
from shapely.geometry import GeometryCollection, Polygon, MultiPolygon

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, assemble, require, sha256,
)

PHOTO_OUT="mois350399_original_photo_TDWG3_WCVP_region_source_status.parquet"
CELL_OUT="mois350399_TDWG3_native_introduced_source_photo_composition.csv"
AUDIT_OUT="mois350399_TDWG3_WCVP_geography_source_receipt.json"
N_PHOTOS=350399
N_CLASSIFIED=164297
N_KNOWN_SOURCE_TAXA=74821
WCVP_MATCH="mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"
WCVP_DIST="mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz"
WCVP_RECEIPT="mois350399_wcvp_TDWG3_source_only_taxon_native_receipt.json"
VALID_STATUSES={
    "native_at_TDWG3","introduced_at_TDWG3","doubtful","extinct",
    "unknown_flag","contradictory_TDWG3_records"
}
TDWG_SOURCE_COMMIT="52da7828aba9d461dd133c27b3bd7a4407161f54"
TDWG_GIT_BLOB="91104e5159e31f88154833a51d3b0d1c9271083f"


def _polygonal_only(geometry):
    """Drop linear artefacts produced while repairing historical TDWG rings."""
    if isinstance(geometry, (Polygon, MultiPolygon)):
        return geometry
    if isinstance(geometry, GeometryCollection):
        pieces = [_polygonal_only(g) for g in geometry.geoms]
        polygonal = [g for g in pieces if g is not None and not g.is_empty]
        if not polygonal:
            return None
        return union_all(polygonal)
    return None


def polygons_read(path, *, minimum_codes=300):
    frame = gpd.read_file(path)
    # The historical GIS source spells this field Level3_cod, but GeoJSON
    # exports can change its case. Require exactly one canonical field.
    column = [k for k in frame.columns if k.casefold() == "level3_cod"]
    require(len(column) == 1 and "geometry" in frame,
            "Original TDWG3 code or polygon field schema changed")
    require(frame.crs is not None and frame.crs.to_epsg() == 4326,
            "TDWG3 original source geometry must declare EPSG:4326")
    frame = frame[[column[0], "geometry"]].rename(
        columns={column[0]: "LEVEL3_COD"}).copy()
    frame["LEVEL3_COD"] = frame.LEVEL3_COD.astype(str).str.strip()
    require(frame.LEVEL3_COD.str.fullmatch(r"[A-Z]{3}").all(),
            "Original TDWG Level-3 country codes are not valid three-letter codes")
    missing = int(frame.geometry.isna().sum())
    empty = int(frame.geometry.is_empty.sum())
    require(missing == 0 and empty == 0,
            f"Original TDWG3 has {missing} missing and {empty} empty shapes; "
            "do not classify points to those undefined country boundaries")
    # Real, pinned TDWG/WGSRPD 2011 coast/botanical-country shapes contain
    # invalid ring self-crossings. Record AND repair them deterministically
    # (GEOS make_valid), rather than assuming invalid shapely predicates
    # reflect genuine native locations. Never modify original GeoJSON bytes.
    invalid_mask = ~frame.geometry.is_valid
    invalid_count = int(invalid_mask.sum())
    collection_count = 0
    repaired = 0
    for idx in frame.index[invalid_mask]:
        original = frame.at[idx, "geometry"]
        corrected = make_valid(original)
        if isinstance(corrected, GeometryCollection):
            collection_count += 1
        corrected = _polygonal_only(corrected)
        require(corrected is not None and not corrected.is_empty,
                f"Invalid historical TDWG3 polygon {idx} has no repairable polygon")
        require(corrected.is_valid,
                f"Historical TDWG3 polygon {idx} remains invalid after make_valid")
        frame.at[idx, "geometry"] = corrected
        repaired += 1
    require(frame.geometry.is_valid.all() and not frame.geometry.is_empty.any(),
            "Not all original TDWG3 polygons can be repaired for geographic overlay")
    require(frame.LEVEL3_COD.nunique() >= minimum_codes,
            "Original TDWG3 polygon table covers too few Level-3 regions")
    frame.attrs["source_geometry_repair_audit"] = {
        "original_feature_count": len(frame),
        "original_unique_TDWG3_codes": int(frame.LEVEL3_COD.nunique()),
        "original_missing_shape_count": missing,
        "original_empty_shape_count": empty,
        "original_invalid_polygon_count": invalid_count,
        "repair_attempted_count": repaired,
        "original_invalid_repair_geometry_collections": collection_count,
        "remaining_invalid_after_polygon_only_make_valid": 0,
        "geometry_repair": "GEOS make_valid_polygon_components_only",
        "original_input_file_modified": False,
        "coastline_or_island_specific_native_status_inferred": False,
    }
    print("MOIS_ORIGINAL_TDWG3_GEOMETRY_REPAIR_AUDIT",
          json.dumps(frame.attrs["source_geometry_repair_audit"],
                     sort_keys=True), flush=True)
    return frame

def assign_photo_regions(photos, geo):
    """Unambiguous point-within-polygon TDWG3, never snapping ocean/coasts."""
    require(photos.photo_id.is_unique and "latitude" in photos
            and "longitude" in photos, "Invalid original photo coordinates")
    original=photos[["photo_id","latitude","longitude"]].copy().reset_index(drop=True)
    lat=pd.to_numeric(original.latitude,errors="coerce")
    lon=pd.to_numeric(original.longitude,errors="coerce")
    valid=lat.between(-90,90)&lon.between(-180,180)
    result=pd.DataFrame({
        "photo_id":original.photo_id.astype(str).to_numpy(),
        "tdwg_level3_code":[""]*len(original),
        "tdwg3_geometry_assignment_state":np.where(
            valid,"no_polygonic_interior_match_or_boundary",
            "original_coordinates_missing_or_invalid"),
    })
    result["source_original_latitude"]=lat
    result["source_original_longitude"]=lon
    # Positional uncertainty is preserved; the dominant legacy cohort has no
    # authenticated accuracy field, so most entries must remain unavailable.
    accuracy=(pd.to_numeric(photos["positional_accuracy"],errors="coerce")
              if "positional_accuracy" in photos else
              pd.Series(np.nan,index=photos.index))
    acc=accuracy.reset_index(drop=True)
    result["original_reported_accuracy_m"]=acc
    result["positional_accuracy_status"]=np.select(
        [acc.notna()&acc.ge(0)&acc.le(5000),
         acc.notna()&acc.gt(5000)],
        ["source_reported_le_5000m","source_reported_gt_5000m"],
        default="not_source_reported")
    if not valid.any():
        return result
    idx=np.flatnonzero(valid.to_numpy())
    coordinates=gpd.GeoDataFrame(
        {"source_row_id":idx},
        geometry=gpd.points_from_xy(lon.iloc[idx],lat.iloc[idx]),
        crs="EPSG:4326")
    hits=gpd.sjoin(coordinates,
                  geo[["LEVEL3_COD","geometry"]],
                  how="left",predicate="within")
    known=hits[["source_row_id","LEVEL3_COD"]].dropna().drop_duplicates()
    if len(known):
        counts=known.groupby("source_row_id",observed=True).LEVEL3_COD.agg(
            ["size","first"])
        unique=counts.index[counts["size"].eq(1)].to_numpy(int)
        ambiguous=counts.index[counts["size"].gt(1)].to_numpy(int)
        result.loc[unique,"tdwg_level3_code"]=counts.loc[unique,"first"].to_numpy()
        result.loc[unique,"tdwg3_geometry_assignment_state"]="unique_TDWG3_polygon"
        result.loc[ambiguous,"tdwg3_geometry_assignment_state"]="multiple_TDWG3_polygons"
    require(len(result)==len(photos) and result.photo_id.is_unique,
            "Original photo to TDWG3 map lost source identity")
    return result


def verified_wcvp_tables(directory, original_photo_files, photos):
    directory=Path(directory)
    manifest=json.loads((directory/WCVP_RECEIPT).read_text())
    require(manifest["schema"]==
            "mois350399_kew_wcvp_official_accepted_species_TDWG3_native_source_v1"
            and manifest["source_photo_attempts"]==N_PHOTOS
            and manifest["source_four_hue_classified"]==N_CLASSIFIED
            and manifest["wcvp_distributions_resolution"]=="TDWG_level3_only"
            and manifest["island_native_status_inferred"] is False
            and manifest["photo_coordinate_overlaid_TDWG3_polygons"] is False
            and manifest["H3_native_island_hypothesis_tested"] is False,
            "Not the authenticated real350399 original Kew TDWG3 source-only Release")
    for file in (WCVP_MATCH,WCVP_DIST):
        require(sha256(directory/file)==manifest["data_files_sha256"][file],
                "Published Kew names/distribution SHA changed: "+file)
    for arm,file in original_photo_files.items():
        require(sha256(file)==manifest["source_parquet_sha256"][arm],
                "Source MOIS photo pixel cohort changed after WCVP taxonomy was measured")
    taxa=pd.read_csv(directory/WCVP_MATCH,dtype=str,keep_default_na=False)
    distribution=pd.read_csv(directory/WCVP_DIST,dtype=str,keep_default_na=False)
    require(len(taxa)==manifest["taxa_with_any_source_photos"]==N_KNOWN_SOURCE_TAXA
            and taxa.inat_taxon_id.is_unique
            and taxa.wcvp_accepted_plant_name_id.notna().all(),
            "WCVP accepted-source taxon table is incomplete or duplicated")
    sums=photos.groupby("inat_taxon_id",observed=True).agg(
        total=("photo_id","size"),hues=("classified","sum"))
    matched=taxa.set_index("inat_taxon_id").loc[sums.index]
    require((pd.to_numeric(matched.source_original_photo_attempts).to_numpy()
             ==sums.total.to_numpy()).all()
            and (pd.to_numeric(matched.source_four_hue_classified_photos).to_numpy()
                 ==sums.hues.to_numpy()).all(),
            "WCVP accepted species taxon-to-photo source identity count drift")
    require({"wcvp_accepted_plant_name_id","tdwg_level3_code",
             "native_introduced_state"}<=set(distribution),
            "WCVP source Level3 distribution fields missing")
    require(distribution.native_introduced_state.isin(VALID_STATUSES).all(),
            "Unexpected source distribution native or introduced flags")
    require(not distribution.duplicated([
        "wcvp_accepted_plant_name_id","tdwg_level3_code"]).any(),
        "Contradictory WCVP source rows not reconciled by original names importer")
    return taxa,distribution,manifest


def attach_region_status(photos, geo_assign, taxa, distribution):
    """The result is a photo's TAXON'S regional WCVP status, NOT that photo's nativeness."""
    require(geo_assign.photo_id.tolist()==photos.photo_id.astype(str).tolist(),
            "TDWG geography not aligned with immutable original photo IDs")
    result=photos[["photo_id","inat_taxon_id","species","arm","observation_id",
                   "morph","classified","white"]].copy().reset_index(drop=True)
    result=result.merge(geo_assign,on="photo_id",how="left",
                        sort=False,validate="one_to_one")
    require(result.photo_id.tolist()==photos.photo_id.tolist(),
            "Geographic status join reordered original source photo roster")
    result=result.merge(
        taxa[["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"]],
        how="left",on="inat_taxon_id",sort=False,validate="many_to_one")
    native=distribution[["wcvp_accepted_plant_name_id","tdwg_level3_code",
                         "native_introduced_state"]]
    result=result.merge(native,on=["wcvp_accepted_plant_name_id","tdwg_level3_code"],
                        how="left",sort=False,validate="many_to_one")
    match=result.wcvp_match_state.eq("unique_exact_name_to_accepted_id")
    result.loc[~match,"native_introduced_state"]=pd.NA
    result["taxon_TDWG3_range_annotation"]=np.select([
        result.tdwg3_geometry_assignment_state.ne("unique_TDWG3_polygon"),
        ~match,
        result.native_introduced_state.isna(),
        result.native_introduced_state.notna()
    ],[
        "TDWG3_photo_location_unassigned",
        "accepted_species_name_unresolved",
        "WCVP_range_not_listed_for_this_TDWG3",
        "WCVP_"+result.native_introduced_state.fillna("unknown").astype(str)
    ], default="other_uncertainty")
    result["native_on_specific_island"]="not_established_by_WCVP_TDWG3"
    result["photographed_individual_nativeness"]="not_established_from_area_level_range"
    require(len(result)==len(photos) and result.photo_id.is_unique,
            "Attaching broad TDWG3 native metadata duplicated original photo pixels")
    require(result.classified.eq(result.morph.isin(HUES)).all(),
            "Original photographic RGB ROI classifiability was changed")
    return result


def summarise_source_status(d):
    keys=["arm","taxon_TDWG3_range_annotation"]
    def counts(s):
        c=s.loc[s.classified]
        return pd.Series({
            "photo_attempts":len(s),
            "four_hue_ROI_count":len(c),
            "white_label_ROIs":int(c.white.sum()),
            "yellow_orange_label_ROIs":int(c.morph.eq("yellow_orange").sum()),
            "red_pink_label_ROIs":int(c.morph.eq("red_pink").sum()),
            "blue_purple_label_ROIs":int(c.morph.eq("blue_purple").sum()),
            "unclassifiable_original_ROIs":int((~s.classified).sum()),
            "photo_white_among_classifiable":(float(c.white.mean()) if len(c) else None),
            "nominal_taxa_with_original_photos":int(s.inat_taxon_id.nunique()),
            "source_distinct_observations":int(s.observation_id.nunique())
        })
    rows=[]
    for group,g in d.groupby(keys,observed=True,dropna=False):
        vals=counts(g)
        rows.append({**dict(zip(keys,group)),**vals.to_dict()})
    out=pd.DataFrame(rows)
    require(int(out.photo_attempts.sum())==len(d)
            and int(out.four_hue_ROI_count.sum())==int(d.classified.sum()),
            "TDWG3 regional status denominator not equal all actually attempted photos")
    return out


def run(core,broad,wave2,cached8192,tail16384, wcvp_dir,
        tdwg_geojson, output):
    output=Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "Refuse overwriting previous signed original-photo TDWG3 archive")
    photos,arms=assemble(core,broad,wave2,cached8192,tail16384,strict=True)
    paths=dict(core=core,broad=broad,wave2=wave2,
               cached8192=cached8192,tail16384=tail16384)
    names,dist,source=verified_wcvp_tables(wcvp_dir,paths,photos)
    geo=polygons_read(tdwg_geojson)
    mapped=assign_photo_regions(photos,geo)
    combined=attach_region_status(photos,mapped,names,dist)
    table=summarise_source_status(combined)
    output.mkdir(parents=True,exist_ok=True)
    combined.to_parquet(output/PHOTO_OUT,index=False,compression="zstd")
    table.to_csv(output/CELL_OUT,index=False,lineterminator="\n")
    counts=combined.taxon_TDWG3_range_annotation.value_counts(dropna=False)
    receipt={
        "schema":"mois_original350399_photo_taxon_WCVP_regional_TDWG3_geoposition_v1",
        "original_actual_photo_attempts":len(photos),
        "four_hue_classifiable_original_ROIs":int(photos.classified.sum()),
        "unique_photo_ids":True,
        "WCVP_source_receipt_sha256":sha256(Path(wcvp_dir)/WCVP_RECEIPT),
        "WCVP_input_data_shas":source["data_files_sha256"],
        "WCVP_original_archive_sha256":source["wcvp_archive_sha256"],
        "TDWG_GeoJSON_source_repository":"https://github.com/tdwg/wgsrpd",
        "TDWG_source_commit":TDWG_SOURCE_COMMIT,
        "TDWG_original_Git_blob_SHA1":TDWG_GIT_BLOB,
        "TDWG_original_GeoJSON_SHA256":sha256(tdwg_geojson),
        "TDWG_original_unique_level3_codes":int(geo.LEVEL3_COD.nunique()),
        "TDWG_original_geometry_repair_audit":geo.attrs.get(
            "source_geometry_repair_audit",{}),
        "WCVP_native_distributions_are_only_TDWG_level3":True,
        "photo_geographic_status_counts":{
            str(k):int(v) for k,v in
            mapped.tdwg3_geometry_assignment_state.value_counts().items()},
        "photo_regional_taxon_annotation_counts":{
            str(k):int(v) for k,v in counts.items()},
        "source_colour_technical_unknowns_retained":int((~combined.classified).sum()),
        "native_on_specific_island_determined":False,
        "inferred_wild_vs_cultivated_individual_status":False,
        "H3_native_island_adaptation_tested":False,
        "TDWG_polygon_coastline_scale_matches_exact_island":False,
        "positional_accuracy_metadata_retained":True,
        "no_unmatched_WCVP_name_or_missing_range_forced_to_introduced":True,
        "no_photo_colour_reclassification_or_new_image_download":True,
        "files_sha256":{
            key:sha256(output/key) for key in (PHOTO_OUT,CELL_OUT)}
    }
    (output/AUDIT_OUT).write_text(json.dumps(
        receipt,sort_keys=True,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("MOIS_REAL350399_WCVP_TDWG3_PHOTO_COORDINATE_SOURCE_ONLY_COMPLETE",
          json.dumps({k:receipt[k] for k in (
              "original_actual_photo_attempts",
              "TDWG_original_unique_level3_codes",
              "photo_geographic_status_counts",
              "photo_regional_taxon_annotation_counts")},
              sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("core","broad","wave2","cached8192","tail16384",
                 "wcvp-dir","tdwg-geojson","output"):
        p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args()
    run(args.core,args.broad,args.wave2,args.cached8192,args.tail16384,
        args.wcvp_dir,args.tdwg_geojson,args.output)


if __name__=="__main__":
    main()
