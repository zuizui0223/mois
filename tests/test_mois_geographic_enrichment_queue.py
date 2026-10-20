"""Source-only region enrichment contracts; no fixtures are measured flowers."""
import unittest
from unittest.mock import patch
import json
from pathlib import Path

from scripts.acquisition import mois_queue_from_full93m_target_tiers as q
from test_mois_queue_from_full93m_target_tiers import TestFromTargetTier


class TestGeographicEnrichment(unittest.TestCase):
    def test_predeclared_rates_and_latitude_boundaries(self):
        cases = [(-90, "south_extratropical", 0.25),
                 (-23.441, "south_extratropical", 0.25),
                 (-23.44, "tropical", 0.5),
                 (0, "tropical", 0.5),
                 (23.44, "tropical", 0.5),
                 (23.441, "north_extratropical", 1/64),
                 (90, "north_extratropical", 1/64),
                 ("", "unknown", 1/64)]
        for latitude, band, rate in cases:
            self.assertEqual(q.geographic_enrichment_band(latitude), band)
            self.assertEqual(q.GEO_ENRICHMENT_RATES[band],rate)

    def test_out_of_bounds_or_nonnumeric_latitude_is_unknown_not_native_filter(self):
        for v in ("nan", "inf", "-91", "91", "NA", ""):
            self.assertEqual(q.geographic_enrichment_band(v), "unknown")

    def test_independent_observation_hash_and_joint_photo_selection(self):
        d=q.geographic_enrichment_draw("uuid-photo-A")
        self.assertEqual(d,q.geographic_enrichment_draw("uuid-photo-A"))
        self.assertGreaterEqual(d,0);self.assertLess(d,1)
        self.assertNotEqual(d,q.geographic_enrichment_draw("uuid-photo-B"))

    def test_country_and_botanical_nativeness_never_used_in_rate(self):
        self.assertNotIn("native", q.GEO_ENRICHMENT_SEED.lower())
        self.assertEqual(set(q.GEO_ENRICHMENT_RATES),{
            "north_extratropical","south_extratropical","tropical","unknown"})

    def test_geographic_selection_omits_screen_negatives_only_by_metadata_hash(self):
        t=TestFromTargetTier("test_signed_source_tier_only_unseen_supported_original_images")
        t.setUp()
        try:
            root=t.root
            # Existing fixture has obs A with photo1+2 at latitude=20,
            # obs B at 15; photo1 was historically attempted.
            # The scientific source upper bound stays explicit and stable.
            with patch.object(q,'EXPECTED_PHOTOS',3),patch.object(q,'EXPECTED_TAXA',2), \
                 patch.object(q,'obs_draw',lambda _:0.0001), \
                 patch.object(q,'geographic_enrichment_draw',lambda _:0.1):
                r=q.queue_tier(t.registry,[root/'tiers'/'receipt.json'],
                    root/'tiers'/'targets',t.counts,[t.prior],
                    root/'geo-source',tier=0,compute_bucket=0,
                    compute_bucket_count=1,check_taxa_sha=False,
                    expected_taxa=2,expected_shards=2,
                    geographic_enrichment=True)
            self.assertEqual(r['source_tier_total_photo_candidates'],3)
            self.assertEqual(r['new_photo_attempts_planned'],1)
            self.assertEqual(r['previously_measured_photo_ids_skipped'],1)
            self.assertEqual(r['selected_new_photos_unsupported_by_frozen_loader'],1)
            self.assertEqual(r['geographic_selection_excluded_photo_candidates'],0)
            self.assertEqual(r['selection_probability'],None)
            self.assertEqual(r['source_tier_probability_lower_bound'],0)
            self.assertEqual(r['geographic_selection_rates']['tropical'],0.5)
            output=(root/'geo-source'/'shards'/'shard-000.csv').read_text()
            self.assertIn("global93m_geographic_enrichment_tier_00",output)
            self.assertIn("0.0001",output) # tier span 0.0002 * 0.5
            self.assertFalse(r['native_photo_status_not_inferred_or_hard_filtered']==False)
        finally:
            t.tearDown()

    def test_still_processes_no_photo_if_observation_not_drawn_and_records_exclusions(self):
        t=TestFromTargetTier("test_signed_source_tier_only_unseen_supported_original_images")
        t.setUp()
        try:
            root=t.root
            with patch.object(q,'EXPECTED_PHOTOS',3),patch.object(q,'EXPECTED_TAXA',2), \
                 patch.object(q,'obs_draw',lambda _:0.0001), \
                 patch.object(q,'geographic_enrichment_draw',lambda _:0.9):
                r=q.queue_tier(t.registry,[root/'tiers'/'receipt.json'],
                    root/'tiers'/'targets',t.counts,[t.prior],
                    root/'notdrawn',tier=0,compute_bucket=0,
                    compute_bucket_count=1,check_taxa_sha=False,
                    expected_taxa=2,expected_shards=2,
                    geographic_enrichment=True)
            self.assertEqual(r['new_photo_attempts_planned'],0)
            self.assertEqual(r['geographic_selection_excluded_photo_candidates'],3)
            self.assertEqual(r['work_units'],0)
        finally:
            t.tearDown()

if __name__=="__main__":
    unittest.main()
