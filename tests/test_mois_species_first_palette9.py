"""Synthetic tests: no photo acquisition; no GBIF remote calls."""
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_palette9 import (
    build_catalogue, plan_next_photos, PALETTE
)


def row(pid, taxon, obs, *, hue="white", strict=True):
    o = dict(photo_id=str(pid), inat_taxon_id=str(taxon), observation_id=str(obs), strict9=strict)
    for p in PALETTE:
        o["strict9_" + p] = (1. if p == hue else 0.) if strict else np.nan
    return o


class SpeciesFirstPalette9Tests(unittest.TestCase):
    def catalogue(self):
        photos = pd.DataFrame([
            row(1, 10, "a"), row(2, 10, "a"), row(3, 10, "b"),
            row(4, 10, "c"), row(5, 11, "d", hue="yellow"),
            row(6, 12, "e", hue="red"), row(7, 12, "f", hue="red"),
            row(8, 12, "g", hue="blue"), row(9, 12, "h", hue="blue"),
            row(10, 13, "i", strict=False), row(11, 14, "", hue="magenta")
        ])
        taxa = pd.DataFrame({"inat_taxon_id": ["10", "11", "12", "13", "14", "15"],
                             "species": ["Plant " + str(i) for i in range(6)],
                             "full_93m_original_observation_count": [10] * 6,
                             "full_93m_original_photo_count": [20] * 6})
        coverage = pd.DataFrame({"inat_taxon_id": ["10", "11", "12", "13", "14", "15"],
                                 "after_photo_attempts": [3, 1, 4, 1, 1, 0],
                                 "after_strict9": [3, 1, 4, 0, 1, 0]})
        return photos, taxa, coverage

    def test_two_photos_one_observation_do_not_fake_three_independent(self):
        p, n, c = self.catalogue()
        atlas = build_catalogue(p, n, c).set_index("inat_taxon_id")
        self.assertEqual(atlas.loc["10", "hue_status"], "sampling_consensus_single")
        self.assertEqual(atlas.loc["10", "independent_roi9_observations"], 3)
        self.assertEqual(atlas.loc["11", "hue_status"], "provisional_single")
        self.assertEqual(atlas.loc["12", "hue_status"], "sampled_multicolour")
        self.assertEqual(atlas.loc["12", "observed_hue9_set"], "blue|red")
        self.assertEqual(atlas.loc["13", "hue_status"], "no_strict9")
        self.assertEqual(atlas.loc["15", "hue_status"], "no_strict9")
        self.assertEqual(atlas.loc["14", "hue_status"], "provisional_single")
        self.assertEqual(atlas.loc["14", "independent_roi9_observations"], 0)

    def test_queue_prioritizes_missing_and_does_not_repeat_photo_observation_or_taxon(self):
        p, n, c = self.catalogue()
        atlas = build_catalogue(p, n, c)
        candidate = pd.DataFrame({"photo_id": ["1", "20", "21", "22", "23", "24", "25", "26"],
                                  "observation_id": ["z", "a", "p", "q", "r", "s", "t", "u"],
                                  "inat_taxon_id": ["11", "11", "13", "13", "15", "10", "11", "12"]})
        queue = plan_next_photos(atlas, p, candidate, quota=10)
        self.assertEqual(set(queue.inat_taxon_id), {"13", "15", "11"})
        self.assertEqual(len(queue), 3)
        self.assertEqual(queue.iloc[0].hue_status, "no_strict9")
        self.assertEqual(queue.iloc[1].hue_status, "no_strict9")
        self.assertEqual(len(set(queue.photo_id)), len(queue))

    def test_missing_numeric_vectors_and_photo_duplicates_are_not_silently_accepted(self):
        p, n, c = self.catalogue()
        p.loc[p.photo_id == "1", "strict9_white"] = np.nan
        with self.assertRaisesRegex(ValueError, "numeric"):
            build_catalogue(p, n, c)
        p, n, c = self.catalogue()
        p.loc[p.photo_id == "2", "photo_id"] = "1"
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            build_catalogue(p, n, c)


if __name__ == "__main__":
    unittest.main()
