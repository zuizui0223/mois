"""Synthetic TDWG3 boundaries test original-photo name/range status, without island assertions."""
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from scripts.analysis.mois_real350399_wcvp_geographic_native_status_20261019 import (
    assign_photo_regions, attach_region_status, polygons_read, summarise_source_status,
)


def polygons():
    return gpd.GeoDataFrame(
        {"LEVEL3_COD":["AAA","BBB"]},
        geometry=[box(10,10,20,20),box(20,10,30,20)],
        crs="EPSG:4326"
    )


def photo(pid, tax, lon, lat, hue, arm="core", acc=None):
    return dict(photo_id=str(pid),inat_taxon_id=str(tax),species="Campanula punctata",
                observation_id="o"+str(pid),arm=arm,
                latitude=lat,longitude=lon,
                positional_accuracy=acc,morph=hue,
                classified=hue in ("white","red_pink","blue_purple","yellow_orange"),
                white=hue=="white")


def names():
    return pd.DataFrame([
        {"inat_taxon_id":"1","wcvp_match_state":"unique_exact_name_to_accepted_id",
         "wcvp_accepted_plant_name_id":"kew-001"},
        {"inat_taxon_id":"2","wcvp_match_state":"multiple_accepted_ids_ambiguous",
         "wcvp_accepted_plant_name_id":""},
        {"inat_taxon_id":"3","wcvp_match_state":"unique_exact_name_to_accepted_id",
         "wcvp_accepted_plant_name_id":"kew-003"}
    ])


def ranges():
    return pd.DataFrame([
        {"wcvp_accepted_plant_name_id":"kew-001",
         "tdwg_level3_code":"AAA",
         "native_introduced_state":"native_at_TDWG3"},
        {"wcvp_accepted_plant_name_id":"kew-001",
         "tdwg_level3_code":"BBB",
         "native_introduced_state":"introduced_at_TDWG3"}
    ])


def test_read_tdwg3_geojson_geometry_and_code(tmp_path):
    path=tmp_path/"source.geojson"
    polygons().to_file(path,driver="GeoJSON")
    out=polygons_read(path,minimum_codes=2)
    assert out.crs.to_epsg()==4326
    assert sorted(out.LEVEL3_COD)==["AAA","BBB"]


def test_geographical_uncertainty_is_never_assigned_as_mainland(tmp_path):
    photos=pd.DataFrame([
        photo(1,1,15,15,"white",acc=100),
        photo(2,1,25,15,"red_pink",acc=10000),
        photo(3,2,20,15,"white"),   # exact boundary: within neither
        photo(4,1,-40,0,"blue_purple"),
        photo(5,1,np.nan,5,"white"),
        photo(6,1,15,15,"mixed_uncertain")
    ])
    geo=assign_photo_regions(photos,polygons())
    assert geo.tdwg_level3_code.tolist()==["AAA","BBB","","","","AAA"]
    assert geo.loc[2,"tdwg3_geometry_assignment_state"]=="no_polygonic_interior_match_or_boundary"
    assert geo.loc[4,"tdwg3_geometry_assignment_state"]=="original_coordinates_missing_or_invalid"
    assert geo.loc[0,"positional_accuracy_status"]=="source_reported_le_5000m"
    assert geo.loc[1,"positional_accuracy_status"]=="source_reported_gt_5000m"
    assert geo.loc[5,"positional_accuracy_status"]=="not_source_reported"


def test_photo_region_join_not_native_individual_or_islet():
    photos=pd.DataFrame([
        photo(1,1,15,15,"white"),
        photo(2,1,25,15,"red_pink",arm="tail16384"),
        photo(3,2,15,15,"white"),
        photo(4,3,15,15,"white"),
        photo(5,1,20,15,"mixed_uncertain")
    ])
    geo=assign_photo_regions(photos,polygons())
    out=attach_region_status(photos,geo,names(),ranges())
    assert len(out)==5
    assert out.taxon_TDWG3_range_annotation.tolist()==[
        "WCVP_native_at_TDWG3","WCVP_introduced_at_TDWG3",
        "accepted_species_name_unresolved",
        "WCVP_range_not_listed_for_this_TDWG3",
        "TDWG3_photo_location_unassigned"
    ]
    assert out.native_on_specific_island.eq("not_established_by_WCVP_TDWG3").all()
    assert out.photographed_individual_nativeness.eq(
        "not_established_from_area_level_range").all()
    assert out.classified.sum()==4
    assert out.white.sum()==3
    table=summarise_source_status(out)
    assert table.photo_attempts.sum()==5
    assert table.four_hue_ROI_count.sum()==4
    assert table.unclassifiable_original_ROIs.sum()==1


def test_multiple_polygons_different_codes_kept_ambiguous():
    overlap=gpd.GeoDataFrame(
        {"LEVEL3_COD":["AAA","BBB"]},
        geometry=[box(10,10,20,20),box(15,15,25,25)],
        crs="EPSG:4326")
    data=pd.DataFrame([photo(1,1,17,17,"white")])
    out=assign_photo_regions(data,overlap)
    assert out.loc[0,"tdwg3_geometry_assignment_state"]=="multiple_TDWG3_polygons"
    assert out.loc[0,"tdwg_level3_code"]==""


def test_several_photos_same_observation_are_not_pseudoreplicated():
    photos=pd.DataFrame([photo(1,1,15,15,"white"),
                         photo(2,1,15,15,"red_pink")])
    photos.loc[1,"observation_id"]="o1"
    status=attach_region_status(photos,assign_photo_regions(photos,polygons()),
                                names(),ranges())
    table=summarise_source_status(status)
    assert table.photo_attempts.sum()==2
    assert table.source_distinct_observations.sum()==1
