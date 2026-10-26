"""Regression: actual nine-component original source photos are not 4-hue guesses."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd

from scripts.analysis.mois_real1043930_continuous9_source_lake import (
    read_nine, PALETTE, FRACTIONS, PROVENANCE,
)

class SourceNinePaletteTests(unittest.TestCase):
    def setUp(self):
        self.t=TemporaryDirectory()
        self.addCleanup(self.t.cleanup)
        self.file=Path(self.t.name)/"roi.parquet"

    def fixture(self,with_fractions=True):
        statuses=[
            "classified_four_state_morph",
            "classified_four_state_morph",
            "not_evaluable_ambiguous_palette_composition",
            "not_evaluable_roi_or_flip_gate",
            "image_acquisition_failed"
        ]
        data=pd.DataFrame({
            "photo_id":["10","11","12","13","14"],
            "morph":["white","red_pink","mixed_uncertain","mixed_uncertain","mixed_uncertain"],
            "measurement_status":statuses,
        })
        if with_fractions:
            for f in FRACTIONS:data[f]=np.zeros(len(data),dtype=float)
            data.loc[0,FRACTIONS[0]]=0.7
            data.loc[0,FRACTIONS[1]]=0.3
            data.loc[1,FRACTIONS[3]]=1.
            data.loc[2,FRACTIONS[0]]=0.4
            data.loc[2,FRACTIONS[4]]=0.6
            data.loc[3,FRACTIONS[0]]=0.99
            # Other quality gate failure with positive palette is NOT rescuable.
        data.to_parquet(self.file,index=False)
        return data

    def test_original_known_and_ambiguous_only_have_continuous_data(self):
        self.fixture()
        r,a=read_nine(self.file,source_arm="TEST",original_four_hues=2)
        self.assertEqual(a["strict_original_4hue_ROIs_with_all_9_fractions"],2)
        self.assertEqual(a["ambiguous_explicit_palette_usable_exploratory"],1)
        self.assertEqual(a["original_photos_with_positive_all9_numerics"],4)
        self.assertEqual(r.continuous9_status.iloc[2],PROVENANCE["ambiguous"])
        self.assertEqual(r.continuous9_status.iloc[3],PROVENANCE["other"])
        self.assertTrue(np.isnan(r.continuous9_white.iloc[3]))
        self.assertTrue(np.isnan(r.continuous9_white.iloc[4]))
        self.assertAlmostEqual(r.continuous9_white.iloc[0],.7,places=6)
        self.assertGreater(r.continuous9_entropy.iloc[0],0)
        self.assertTrue(r.continuous9_strict.iloc[0])
        self.assertFalse(r.continuous9_strict.iloc[2])

    def test_missing_fraction_columns_does_not_map_four_hues_to_palette(self):
        self.fixture(with_fractions=False)
        r,a=read_nine(self.file,source_arm="OLD",original_four_hues=2)
        self.assertEqual(a["strict_original_4hue_ROIs_with_all_9_fractions"],0)
        self.assertFalse(a["all_nine_original_fraction_columns_present"])
        self.assertEqual(int(r.continuous9_strict.sum()),0)
        self.assertTrue(r.continuous9_white.isna().all())

    def test_partial_fraction_source_not_interpreted_as_complete_vector(self):
        x=self.fixture()
        x.drop(columns=[FRACTIONS[-1]]).to_parquet(self.file,index=False)
        r,a=read_nine(self.file,source_arm="PARTIAL",original_four_hues=2)
        self.assertEqual(a["strict_original_4hue_ROIs_with_all_9_fractions"],0)
        self.assertEqual(a["available_fraction_columns"],list(FRACTIONS[:-1]))

    def test_source_status_disagreement_is_error_not_silent_recoding(self):
        x=self.fixture()
        x.loc[0,"measurement_status"]="not_evaluable_roi_or_flip_gate"
        x.to_parquet(self.file,index=False)
        with self.assertRaisesRegex(ValueError,"status contradicts 4-hue"):
            read_nine(self.file,source_arm="BAD",original_four_hues=2)

if __name__=="__main__":
    unittest.main()
