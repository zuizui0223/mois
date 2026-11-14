"""Source-only depth256 tests; no original image pixels are fetched."""
import unittest
from scripts.acquisition.mois_species_first_unlabelled_depth256 import choose_unlabelled
from scripts.acquisition.mois_species_first_zero9_source import select_rows, output_rows
from scripts.acquisition.mois_onepass_full93m_target_tiers import obs_draw
import pandas as pd


class SignedDepth256SourceTests(unittest.TestCase):
    def test_only_unlabelled_taxa_in_current_signed_disjoint_group(self):
        a=pd.DataFrame({
            "inat_taxon_id":[str(x) for x in range(1,42112)],
            "dominant_hue9":["white"]*42111,
            "hue_status":["provisional_single"]*42111
        })
        a.loc[a.inat_taxon_id=="16",["dominant_hue9","hue_status"]]=["","no_strict9"]
        a.loc[a.inat_taxon_id=="17",["dominant_hue9","hue_status"]]=["","mixed_image_unresolved"]
        a.loc[a.inat_taxon_id=="23",["dominant_hue9","hue_status"]]=["","no_strict9"]
        a.loc[a.inat_taxon_id=="24",["dominant_hue9","hue_status"]]=["","multicolour_candidate"]
        a.loc[a.inat_taxon_id=="18",["dominant_hue9","hue_status"]]=["blue","no_strict9"]
        self.assertEqual(choose_unlabelled(a,2),{"16","17","23"})
        self.assertEqual(choose_unlabelled(a,3),set())
        a.loc[0,"inat_taxon_id"]="16"
        with self.assertRaisesRegex(ValueError,"incomplete"):
            choose_unlabelled(a,2)

    def test_256_independent_orig_photos_but_no_observation_repeats_or_unlicensed(self):
        rows=[]
        for i in range(281):
            uuid=f"raw-original-obs-{i}"
            rows.append({
                "inat_taxon_id":"16","photo_id":str(500000+i),
                "observation_uuid":uuid,
                "observation_inclusion_score":str(obs_draw(uuid)),
                "frozen_loader_supported":"true",
                "photo_license":"cc-by",
                "image_extension":"jpg",
                "observer_id":"100","observed_on":"2023-07-01",
                "latitude":"35","longitude":"135","tier":"0"
            })
        # Repeat a photo from one previously seen original observation.
        seen_obs={rows[1]["observation_uuid"]}
        seen_photo={rows[2]["photo_id"]}
        rows[3]["photo_license"]="all-rights-reserved"
        pool,stats=select_rows(rows,{"16"},seen_photo,seen_obs,limit=256)
        self.assertEqual(len(pool["16"]),256)
        selected=output_rows(pool)
        self.assertEqual(len(selected),256)
        self.assertTrue(selected.photo_id.is_unique)
        self.assertTrue(selected.observation_id.is_unique)
        self.assertTrue(selected.photo_url_large.str.startswith(
            "https://inaturalist-open-data.s3.amazonaws.com/photos/").all())
        self.assertFalse(set(selected.photo_id)&seen_photo)
        self.assertFalse(set(selected.observation_id)&seen_obs)
        self.assertNotIn(rows[3]["photo_id"],set(selected.photo_id))
        self.assertEqual(stats["already_measured_or_observed"],2)
        self.assertEqual(stats["missing_nonCC_photo_licence"],1)

if __name__=="__main__":
    unittest.main()
