"""No-network tests for source-backed COL XR GBIF equal-area occurrence grid."""
import unittest
import pandas as pd
from pyproj import Transformer
from scripts.analysis.mois_gbif_colxr_occurrence_grid_pilot import (
    allowed_coordinates,grid_100km,choose_species
)


class GbifGridPilotTests(unittest.TestCase):
    def test_coordinates_with_uncertainty_filter(self):
        self.assertTrue(allowed_coordinates({
            "decimalLatitude":35.5,"decimalLongitude":139.5,
            "coordinateUncertaintyInMeters":1000}))
        self.assertFalse(allowed_coordinates({
            "decimalLatitude":0,"decimalLongitude":0}))
        self.assertFalse(allowed_coordinates({
            "decimalLatitude":35.5,"decimalLongitude":139.5,
            "coordinateUncertaintyInMeters":85000}))
        self.assertFalse(allowed_coordinates({
            "decimalLatitude":99,"decimalLongitude":139.5}))

    def test_equal_area_projection_cells_not_raw_latitude_grid(self):
        trans=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
        a=grid_100km(139.0,35.0,trans)
        self.assertTrue(all(isinstance(x,int) for x in a))
        self.assertNotEqual(a,(139,35))

    def test_only_actual_colxr_species_accepted(self):
        d=pd.DataFrame({"status":["colxr_exact_accepted_species",
                                  "synonym_or_doubtful_review",
                                  "colxr_exact_accepted_species"],
                        "key":["AAA","OLD","BBB"],
                        "inat_taxon_id":["1","2","3"],
                        "source_species":["Abies firma","Abies sp","Abies koreana"],
                        "dominant_hue9":["white","white","purple"]})
        got=choose_species(d)
        self.assertEqual(set(got.key),{"AAA","BBB"})


if __name__=="__main__":
    unittest.main()
