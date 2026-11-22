"""Failure-mode and scientific-inference unit tests for the v4.0 chronology."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_firstrecords_v40_interregional_firstrecords_photo9_20261122 import (
    region_curves, model_sample, strict_links, test_hue_association as evaluate_photo_hue,
    MAIN_YEAR_MIN, YEAR_ENDPOINTS
)

class FirstRecordsV40Tests(unittest.TestCase):
    def setUp(self):
        self.flora=pd.DataFrame({
            "key":["k1","k2","k3","k4"],
            "species":["Acacia alpha","Acacia beta","Betula gamma","Betula delta"],
            "genus":["Acacia","Acacia","Betula","Betula"],
            "source_hue9":["white","yellow","white","blue"],
            "hue_status":["provisional_single"]*4,
            "occupied_grid_cells":[10,20,31,300]})
    def test_chronology_same_year_multiple_regions_not_a_spread_event(self):
        r=pd.DataFrame([
            ("k1","loc1",1900,"Europe"),("k1","loc2",1900,"Asia"),
            ("k1","loc3",1940,"Europe"),("k1","loc4",2001,"Asia"),
            ("k2","loc5",1950,"Africa"),("k2","loc6",2002,"Africa"),
            ("k3","loc7",2000,"Europe"),
            ("k4","loc8",1870,"Europe"),("k4","loc9",1899,"Europe"),
        ],columns=["key","locationID","year","continent"])
        curves=region_curves(r,self.flora).set_index("key")
        self.assertEqual(curves.loc["k1","same_year_regions"],2)
        self.assertEqual(curves.loc["k1","new_regions_50yr"],1)
        self.assertEqual(curves.loc["k1","new_regions_100yr"],1)
        self.assertEqual(curves.loc["k1","first_record_continent"],"MULTIPLE_CONTINENTS")
        self.assertEqual(curves.loc["k2","any_new_regions_50yr"],0)
        self.assertEqual(curves.loc["k2","any_new_regions_100yr"],1)
        self.assertEqual(curves.loc["k3","same_year_regions"],1)
        self.assertNotIn("k3",model_sample(curves.reset_index(),50).key.tolist())
        self.assertIn("k1",model_sample(curves.reset_index(),50).key.tolist())
        self.assertNotIn("k1",model_sample(curves.reset_index(),50,single_start=True).key.tolist())
        self.assertEqual(YEAR_ENDPOINTS[50],1973)
    def test_outcome_matches_exact_source_confidence_establishment_and_region(self):
        loc=pd.DataFrame([
            ("L1","Europe mainland","Europe"),
            ("L1","Europe mainland","Europe"),
            ("L2","Island independent region","Asia"),
        ],columns=["locationID","location","continent"])
        tax=pd.DataFrame([
            ("1","Acacia alpha","Vascular plants","ACCEPTED","EXACT","SPECIES"),
            ("2","Acacia beta","Vascular plants","SYNONYM","EXACT","SPECIES"),
        ],columns=["taxonID","species","taxaGroup","GBIFstatus","GBIFmatchType","GBIFtaxonRank"])
        base=["1","Acacia alpha"]
        r=pd.DataFrame([
            ("L1","1","Acacia alpha","1900","1900","high confidence","present","introduced","established","Free"),
            ("L2","1","Acacia alpha","1940","1940","high confidence","present","introduced","established","Free"),
            ("L2","1","Acacia alpha","1942","1942","high confidence","present","introduced","established","Free"),
            ("L1","1","Acacia alpha","1950","1940 - 1950","high confidence","present","introduced","established","Free"),
            ("L1","1","Acacia alpha","1980","1980","high confidence","present","introduced","casual","Free"),
            ("L2","2","Acacia beta","1910","1910","high confidence","present","introduced","established","Free"),
        ],columns=["locationID","taxonID","taxon","firstRecordEvent","verbatimFirstRecordEvent",
                   "confidenceFirstRecordEvent","occurrenceStatus","establishmentMeans",
                   "degreeOfEstablishment","accessRights"])
        r["location"]=r.locationID.map({"L1":"Europe mainland","L2":"Island independent region"})
        matched,_,audit=strict_links(r,tax,loc,self.flora)
        self.assertEqual(len(matched),2)
        self.assertEqual(int(matched.year.min()),1900)
        self.assertEqual(int(matched.year.max()),1940)
        self.assertEqual(matched.locationID.nunique(),2)
        self.assertEqual(audit["species_location_duplicate_rows_removed"],1)
        self.assertEqual(audit["source_location_IDs_after_validated_dedup"],2)
    def test_insufficient_exchangeability_does_not_create_fake_p(self):
        curves=pd.DataFrame({
            "key":["k1","k2"],"first_year":[1900,1901],
            "source_hue9":["white","yellow"],"genus":["Acacia"]*2,
            "first_record_continent":["Europe"]*2,
            "hue_status":["provisional_single"]*2,
            "same_year_regions":[1,1],
            "new_regions_25yr":[0,1],"new_regions_50yr":[1,1],"new_regions_100yr":[1,1],
            "any_new_regions_25yr":[0,1],"any_new_regions_50yr":[1,1],"any_new_regions_100yr":[1,1]})
        x=model_sample(curves,50)
        result,draws=evaluate_photo_hue(x,self.flora,50,"genus_breadth",999,20261122)
        self.assertEqual(result["status"],"insufficient primary exchangeability or sample")
        self.assertIsNone(draws)
        self.assertEqual(len(x),2)

if __name__=="__main__":
    unittest.main()
