"""No temporal-change model: test source and spatial floral ALAN correspondence."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis import mois_viirs_spatial_levels_20261105 as a

class TestSpatialLevels(unittest.TestCase):
    def template(self):
        rows=[]
        for b in range(42):
            for j in range(6):
                night=(j+.3*b)/6
                rows.append(dict(cell_id=str(b)+"-"+str(j),block_id=str(b),
                    white_taxon_equal=.40 + .01*night,
                    viirs_quarterly_median=np.expm1(night),
                    population_density=np.expm1(.2*j+.02*b),
                    n_attempted=5+j,bio5=20+j/4,bio12=900+j*70,elev=100+j*10))
        return pd.DataFrame(rows)
    def test_response_is_visible_ROI_not_species_loss(self):
        x=a.get_frame(self.template(),"white_taxon_equal","viirs_quarterly_median","n_attempted")
        self.assertTrue(x.white_pp.between(0,100).all())
        self.assertFalse("species_extinct" in x)
        self.assertEqual(len(x),42*6)
    def test_outcome_year_change_is_not_a_regressor(self):
        self.assertNotIn("delta",a.MODEL)
        self.assertIn("C(block_id)",a.MODEL)
        self.assertIn("population",a.MODEL)
    def test_missing_population_is_not_zero_and_fails(self):
        x=self.template();x.loc[0,"population_density"]=np.nan
        with self.assertRaisesRegex(ValueError,"Nonfinite"):
            a.get_frame(x,"white_taxon_equal","viirs_quarterly_median","n_attempted")
    def test_source_hash_drift_fails_before_read(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"old.csv";p.write_text("cell_id,light\n1,10")
            with self.assertRaisesRegex(ValueError,"source SHA drift"):
                a.source(p,p)
    def test_spatial_cluster_model_requires_independent_blocks(self):
        x=self.template().iloc[:60].copy()
        model=a.get_frame(x,"white_taxon_equal","viirs_quarterly_median","n_attempted")
        with self.assertRaisesRegex(ValueError,"support"):
            a.fit(model,"test",2020,"n_attempted")
    def test_new_photo_or_moth_trait_not_estimated_by_model(self):
        self.assertNotIn("moth",a.MODEL.lower())
        self.assertNotIn("new_photo",a.MODEL.lower())
        self.assertNotIn("delta_white",a.MODEL)
        self.assertTrue(a.SHA_2020.startswith("3562db"))
        self.assertTrue(a.SHA_2024.startswith("1c5f03"))

if __name__=="__main__":
    unittest.main()
