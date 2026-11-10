"""Source-only synthetic tests for authentic first256 nine-palette delta and wave2."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_incremental_first256 import (
    read_first_wave, next_wave, parse_boolean, FIRST_ATTEMPTS, HUES,
)
from scripts.analysis.mois_species_first_palette9 import PALETTE, sha256


class First256IncrementalTests(unittest.TestCase):
    def test_boolean_does_not_convert_string_false_to_true(self):
        a = parse_boolean(pd.Series(["True", "False", "false", "true"]))
        self.assertEqual(a.tolist(), [True, False, False, True])
        with self.assertRaisesRegex(ValueError, "Boolean"):
            parse_boolean(pd.Series(["False", "maybe"]))

    def test_actual_observation_level_new_photo_source_pairing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = pd.DataFrame([{
                "photo_id": str(1000+i), "inat_taxon_id": str(2000+i),
                "observation_id": "independent-"+str(i)} for i in range(FIRST_ATTEMPTS)])
            photos = []
            for i in range(FIRST_ATTEMPTS):
                row = dict(src.iloc[i])
                good = i < 91
                row.update(strict9_original_photo_numeric="True" if good else "False",
                           morph="white" if good else "mixed_uncertain")
                for color in PALETTE:
                    row["flower_fraction_" + color] = (
                        1.0 if color=="white" else 0.0) if good else ""
                photos.append(row)
            measured = pd.DataFrame(photos)
            sf = root/"source.csv"; mf = root/"photos.csv"; rf = root/"receipt.json"
            src.to_csv(sf, index=False)
            measured.to_csv(mf, index=False)
            rf.write_text(json.dumps({
                "attempted": 256, "downloaded": 256, "photo_colours_measured": True,
                "source_sha256": sha256(sf), "photos_sha256": sha256(mf)}))
            got = read_first_wave(sf, mf, rf, sha_check=False)
            self.assertEqual(int(got.strict9.sum()),91)
            self.assertEqual(int(got.measured_palette_dominant_label.eq("white").sum()),91)
            self.assertEqual(got.observation_id.nunique(),256)
            modified = measured.copy(); modified.loc[0,"observation_id"] = "WRONG"
            modified.to_csv(mf,index=False)
            rf.write_text(json.dumps({
                "attempted": 256, "downloaded": 256, "photo_colours_measured": True,
                "source_sha256": sha256(sf), "photos_sha256": sha256(mf)}))
            with self.assertRaisesRegex(ValueError, "taxon/observation conflicted"):
                read_first_wave(sf,mf,rf,sha_check=False)

    def test_wave2_retries_unknown_then_different_new_taxa(self):
        first = pd.DataFrame({
            "inat_taxon_id":["1","2","3"], "photo_id":["101","201","301"],
            "observation_id":["a","b","c"], "strict9":[True,False,True],
            "measured_palette_dominant_label":["white","",""]})
        combined = first[["inat_taxon_id","photo_id","observation_id","strict9"]].copy()
        atlas = pd.DataFrame({"inat_taxon_id":["1","2","3","4","5"],
                              "hue_status":["provisional_single","no_strict9",
                                            "mixed_image_unresolved","no_strict9","no_strict9"],
                              "after_strict9":[1,0,1,0,0]})
        candidates = pd.DataFrame({
            "inat_taxon_id":["1","2","2","3","4","4","5"],
            "photo_id":["102","201","202","302","401","402","501"],
            "observation_id":["d","b","e","f","g","h","i"],
            "source_observation_rank":[2,1,2,2,1,2,1]})
        result = next_wave(atlas,combined,first,candidates,quota=3)
        self.assertEqual(set(result.inat_taxon_id.tolist()[:2]), {"2","3"})
        self.assertEqual(len(set(result.inat_taxon_id)-{"2","3"}),1)
        self.assertTrue((set(result.inat_taxon_id)-{"2","3"}) <= {"4","5"})
        self.assertNotIn("1",set(result.inat_taxon_id))
        self.assertEqual(result.photo_id.nunique(),len(result))
        self.assertEqual(result.observation_id.nunique(),len(result))


if __name__=="__main__":
    unittest.main()
