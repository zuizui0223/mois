"""Offline synthetic regression for wave4 original nine-ROI species-first append."""
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_wave4_incremental import add_and_queue
from scripts.analysis.mois_species_first_palette9 import PALETTE


class FourthWaveIncrementTests(unittest.TestCase):
    def fixture(self):
        df=pd.DataFrame({
            "inat_taxon_id":[str(i) for i in range(1,42112)],
            "hue_status":["no_strict9"]*42111,
            "dominant_hue9":[""]*42111,
            "observed_hue9_set":[""]*42111,
            "after_strict9":[0]*42111,
            "after_photo_attempts":[0]*42111,
            "independent_roi9_observations":[0]*42111,
            "distinct_palette_votes":[0]*42111,
            "ambiguous_observation_rois":[0]*42111,
            "next_measure_priority":[0]*42111
        })
        for c in PALETTE:df["observed_mean_"+c]=np.nan
        delta=pd.DataFrame({
            "inat_taxon_id":["1","2"],"photo_id":["101","201"],
            "observation_id":["a","b"],"strict9":[True,False],
            "dominant_label9":["white",""]
        })
        for c in PALETTE:
            delta["strict9_"+c]=[1. if c=="white" else 0.,np.nan]
        old=pd.DataFrame({
            "inat_taxon_id":["99","100"],
            "photo_id":["901","902"],
            "observation_id":["old1","old2"]
        })
        can=pd.DataFrame({
            "inat_taxon_id":["2","3","4","5"],
            "photo_id":["202","301","401","501"],
            "observation_id":["c","d","e","f"],
            "source_observation_rank":[2,1,1,1]
        })
        return df,delta,old,can

    def test_new_nine_colour_and_missing_photo_are_separate(self):
        atlas,delta,old,can=self.fixture()
        out,q=add_and_queue(atlas,delta,old,can,quota=2,expected_previous=2)
        x=out.set_index("inat_taxon_id")
        self.assertEqual(x.loc["1","after_strict9"],1)
        self.assertEqual(x.loc["1","dominant_hue9"],"white")
        self.assertEqual(x.loc["2","after_strict9"],0)
        self.assertEqual(x.loc["2","hue_status"],"no_strict9")
        self.assertTrue(set(q.inat_taxon_id).issubset({"3","4","5"}))
        self.assertTrue(q.acquisition_priority.eq("first_unmeasured").all())
        self.assertTrue(q.photo_id.is_unique and q.observation_id.is_unique)

    def test_repeat_photo_or_source_observation_cannot_gain_species(self):
        atlas,delta,old,can=self.fixture()
        x=delta.copy();x.loc[0,"photo_id"]="901"
        with self.assertRaisesRegex(ValueError,"previously measured original photo"):
            add_and_queue(atlas,x,old,can,quota=2,expected_previous=2)
        x=delta.copy();x.loc[0,"observation_id"]="old1"
        with self.assertRaisesRegex(ValueError,"already measured previously"):
            add_and_queue(atlas,x,old,can,quota=2,expected_previous=2)

    def test_prior_numeric_colour_cannot_be_overwritten(self):
        atlas,delta,old,can=self.fixture()
        atlas.loc[atlas.inat_taxon_id=="1","after_strict9"]=1
        atlas.loc[atlas.inat_taxon_id=="1","hue_status"]="provisional_single"
        with self.assertRaisesRegex(ValueError,"already had a strict9"):
            add_and_queue(atlas,delta,old,can,quota=2,expected_previous=2)


if __name__=="__main__":
    unittest.main()
