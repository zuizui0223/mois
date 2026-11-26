"""Synthetic regressions for source-only third-wave breadth strategy."""
import unittest
import pandas as pd

from scripts.analysis.mois_species_first_real1280_cumulative import (
    choose_third_wave, catalogue_row_ids,
)


class Wave3BreadthFirstTests(unittest.TestCase):
    def test_breadth_before_retried_taxa_and_no_observation_reuse(self):
        atlas=pd.DataFrame({
            "inat_taxon_id":["10","11","12","13","14","15"],
            "hue_status":["provisional_single","no_strict9","mixed_image_unresolved",
                          "no_strict9","no_strict9","no_strict9"],
            "dominant_hue9":["white","","","","",""],
            "after_strict9":[1,0,1,0,0,0],
            "after_photo_attempts":[1,1,1,0,0,0]})
        previous=pd.DataFrame({"photo_id":["100","110","120"],
                               "observation_id":["a","b","c"],
                               "inat_taxon_id":["10","11","12"]})
        candidates=pd.DataFrame({
            "photo_id":["101","111","121","130","131","140","150"],
            "observation_id":["p","b","q","r","s","t","u"],
            "inat_taxon_id":["10","11","12","13","13","14","15"],
            "source_observation_rank":[2,2,2,1,2,1,1]})
        selected=choose_third_wave(atlas,previous,candidates,quota=2)
        self.assertTrue(set(selected.inat_taxon_id).issubset({"13","14","15"}))
        self.assertTrue(selected.acquisition_priority.eq("first_unmeasured").all())
        self.assertTrue(selected.photo_id.is_unique and selected.observation_id.is_unique)
        self.assertFalse(set(selected.photo_id)&set(previous.photo_id))

    def test_older_species_nine_pixel_ids_deduplicated(self):
        old=pd.DataFrame({"inat_taxon_id":["1","2","2"],
                          "strict9":[True,True,False]})
        first=pd.DataFrame({"inat_taxon_id":["2","3"],
                            "strict9":[False,True]})
        self.assertEqual(set(catalogue_row_ids(first,old)),{"1","2","3"})

    def test_rejects_duplicate_candidate_photos(self):
        atlas=pd.DataFrame({"inat_taxon_id":["1","2"],
                            "hue_status":["no_strict9","no_strict9"],
                            "dominant_hue9":["",""],
                            "after_strict9":[0,0],
                            "after_photo_attempts":[0,0]})
        prior=pd.DataFrame({"photo_id":["9"],"observation_id":["z"],
                            "inat_taxon_id":["9"]})
        candidates=pd.DataFrame({
            "inat_taxon_id":["1","2"],"photo_id":["5","5"],
            "observation_id":["a","b"],"source_observation_rank":[1,1]})
        with self.assertRaisesRegex(ValueError,"duplicate or malformed"):
            choose_third_wave(atlas,prior,candidates,quota=1)


if __name__=="__main__":
    unittest.main()
