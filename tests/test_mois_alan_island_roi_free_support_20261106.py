"""Exact-island ROI-free source-support contracts; synthetic only, no external API."""
import unittest

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

from scripts.analysis import mois_alan_island_roi_free_support_20261106 as m

class IslandROIFreeTests(unittest.TestCase):
    def frames(self):
        canonical=pd.DataFrame({
           "photo_id":["10","11","12","13","14"],
           "observation_id":["o10","o11","o12","o13","o14"],
           "inat_taxon_id":["1","1","2","3","4"],
           "observed_on":["2020-06-10","2020-06-11","2024-06-15","2024-07-01","2020-01-01"],
           "latitude":[1,1.2,1.3,1.4,2],
           "longitude":[140,140.2,140.3,140.4,141],
           "source_subarm":["historical_main"]*4+["adaptive_refill"],
           "strict9":[True,False,False,False,True]})
        names=pd.DataFrame({"inat_taxon_id":["1","2","3","4"],
             "wcvp_match_state":["unique_exact_name_to_accepted_id"]*3+["ambiguous"],
             "wcvp_accepted_plant_name_id":["100","200","300","400"]})
        accepted=pd.DataFrame({"wcvp_accepted_plant_name_id":["100","200","300"],
             "wcvp_accepted_name":["Plant a","Plant b","Plant c"]})
        direct=pd.DataFrame({"accepted_species":["Plant a","Plant b","Plant c"],
            "trait_name":["flower_primary_color"]*3,
            "normalized_value":["white","blue_purple","multistate"],
            "quality":["high","medium","medium"],
            "resolution_status":["resolved"]*3,
            "classification":["accepted"]*3})
        interactions=pd.DataFrame({"plant_name":["Plant a","Plant a","Plant b"],
            "guild":["moth_family_visit","butterfly_family_visit","butterfly_family_visit"]})
        polygons=gpd.GeoDataFrame({"island_id":["islandA"],"area_km2":[101.0]},
              geometry=[Polygon([(139,0),(142,0),(142,3),(139,3)])],crs="EPSG:4326")
        return canonical,names,accepted,direct,interactions,polygons

    def test_photo_original_roi_unclassified_still_retained(self):
        c,n,a,d,g,p=self.frames()
        assigned,counts=m.geolocate_islands(c.loc[c.source_subarm.eq("historical_main")],p)
        x=m.prepare_photo_frame(c,m.accepted_species_mapping(n,a),
           m.independent_colour(d),m.moth_evidence(g),assigned)
        self.assertEqual(counts["photos_unambiguously_inside_islands"],4)
        self.assertEqual(len(x),4)
        self.assertEqual(x.colour_state.notna().sum(),3)
        self.assertEqual(x.strict9.eq(True).sum(),1)
        self.assertEqual(int(x.moth_positive.sum()),2)
        self.assertEqual(int((x.moth_positive & x.butterfly_positive).sum()),2)
        self.assertEqual(x.establishment.unique().tolist(),["unknown"])

    def test_one_accepted_species_counted_once_despite_multiple_photos(self):
        c,n,a,d,g,p=self.frames()
        assigned,_=m.geolocate_islands(c.iloc[:4],p)
        x=m.prepare_photo_frame(c,m.accepted_species_mapping(n,a),
            m.independent_colour(d),m.moth_evidence(g),assigned)
        out=m.summarise(x,m.legacy_viirs_support(),privacy_min_species=1)
        row=out.loc[out.year.eq(2020)].iloc[0]
        self.assertEqual(int(row.photos_all_original_no_roi_gate),2)
        self.assertEqual(int(row.accepted_species),1)
        self.assertEqual(int(row.photos_with_independent_single_colour),2)
        self.assertEqual(int(row.moth_visit_record_species),1)
        self.assertTrue(row.direct_colour_counts_suppressed==False)

    def test_overlap_is_unknown_not_double_island(self):
        c,*_=self.frames()
        p=gpd.GeoDataFrame({"island_id":["A","B"],"area_km2":[101,102]},
            geometry=[Polygon([(139,0),(141,0),(141,3),(139,3)]),
                      Polygon([(139.5,0),(142,0),(142,3),(139.5,3)])],
            crs="EPSG:4326")
        matched,s=m.geolocate_islands(c.iloc[:4],p)
        self.assertEqual(s["photos_in_multiple_island_polygons_excluded"],4)
        self.assertTrue(matched.empty)

    def test_source_mapping_and_unresolved_traits_are_not_nonwhite(self):
        _,n,a,d,g,_=self.frames()
        x=m.accepted_species_mapping(n,a)
        self.assertNotIn("4",set(x.inat_taxon_id))
        labels=m.independent_colour(d)
        self.assertEqual(set(labels.colour_state),{"white","nonwhite"})
        self.assertNotIn("Plant c",set(labels.wcvp_accepted_name))

    def test_native_requires_accepted_island_species_records(self):
        c,n,a,d,g,p=self.frames()
        assignments,_=m.geolocate_islands(c.iloc[:4],p)
        native=pd.DataFrame({
            "island_id":["islandA","islandA","islandA"],
            "wcvp_accepted_name":["Plant a","Plant b","Plant c"],
            "establishment_status":["native","introduced","native"],
            "review_status":["accepted","accepted","pending"]})
        x=m.prepare_photo_frame(c,m.accepted_species_mapping(n,a),
            m.independent_colour(d),m.moth_evidence(g),assignments,native)
        self.assertEqual(x.loc[x.inat_taxon_id.eq("1"),"establishment"].tolist(),
                         ["native","native"])
        self.assertEqual(x.loc[x.inat_taxon_id.eq("2"),"establishment"].tolist(),
                         ["introduced"])
        self.assertEqual(x.loc[x.inat_taxon_id.eq("3"),"establishment"].tolist(),
                         ["unknown"])

    def test_monthly_quarters_are_not_annual_time_series(self):
        v20=pd.DataFrame({"cell_id":["1_1","1_2"],
            "source_complete":[True,False],"viirs_quarterly_median":[2.3,100]})
        t=m.legacy_viirs_support(v20,None)
        self.assertEqual(len(t),1)
        self.assertEqual(t.iloc[0]["year"],2020)
        self.assertEqual(t.iloc[0]["source"],"quarter_months_NOT_annual")
        with self.assertRaisesRegex(ValueError,"Annual VIIRS"):
            m.signed_annual_viirs(pd.DataFrame({"cell_id":["x"],"year":[2020]}))

    def test_annual_viirs_needs_real_coverage_and_unique_ids(self):
        annual=pd.DataFrame({"cell_id":["x","y"],"year":[2012,2024],
            "radiance":[.2,4.1],"cf_cvg":[0,15],
            "product_version":["VNL_V2.1","VNL_V2.2"]})
        good=m.signed_annual_viirs(annual)
        self.assertEqual(len(good),1)
        self.assertEqual(good.year.iloc[0],2024)
        with self.assertRaisesRegex(ValueError,"not unique"):
            m.signed_annual_viirs(pd.concat([annual,annual],ignore_index=True))

    def test_invalid_island_id_or_geometry_fails_closed(self):
        c,*_,p=self.frames()
        p=pd.concat([p,p],ignore_index=True)
        with self.assertRaisesRegex(ValueError,"Island IDs"):
            m.geolocate_islands(c,p)

    def test_ambiguous_same_species_native_status_fails_closed(self):
        c,n,a,d,g,p=self.frames()
        assignment,_=m.geolocate_islands(c.iloc[:4],p)
        status=pd.DataFrame({"island_id":["islandA","islandA"],
         "wcvp_accepted_name":["Plant a","Plant a"],
         "establishment_status":["native","introduced"],
         "review_status":["accepted","accepted"]})
        with self.assertRaisesRegex(ValueError,"Conflicting island-specific"):
            m.prepare_photo_frame(c,m.accepted_species_mapping(n,a),
                m.independent_colour(d),m.moth_evidence(g),assignment,status)

if __name__ == "__main__":
    unittest.main()
