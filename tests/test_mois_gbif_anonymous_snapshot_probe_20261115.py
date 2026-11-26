import unittest
from scripts.analysis.mois_gbif_anonymous_snapshot_probe_20261115 import (
    snapshot_name, choose_recent_snapshot, key_overlap, canonical_pair, classify_join
)
class TestGBIFAnonymousSnapshot(unittest.TestCase):
    def test_date_only_s3_snapshot(self):
        self.assertEqual(snapshot_name("occurrence/2026-11-01/"),"2026-11-01")
        self.assertIsNone(snapshot_name("occurrence/latest/"))
    def test_choose_latest_not_in_future(self):
        dates=["occurrence/2025-08-01/","occurrence/2026-10-01/","occurrence/2026-12-01/"]
        self.assertEqual(choose_recent_snapshot(dates,asof="2026-11-15"),"2026-10-01")
    def test_canonical_species_not_authorship_or_higher_taxon(self):
        self.assertEqual(canonical_pair("Solanum polyacanthos Lam."),"Solanum polyacanthos")
        self.assertIsNone(canonical_pair("Solanum"))
        self.assertIsNone(canonical_pair(""))
        self.assertIsNone(canonical_pair("Solanum sp."))
    def test_key_overlap_is_distinct_from_canonical_match(self):
        obs=["12345","67890","12345"]
        col=["4Y3TB","3WDFP"]
        self.assertEqual(key_overlap(obs,col),0)
        self.assertEqual(classify_join(0,10),"taxonomy_mismatch_key_bridge_required")
        self.assertEqual(classify_join(5,10),"actual_colxr_key_overlap")
        self.assertEqual(classify_join(0,0),"no_matched_species_in_sample")
if __name__=="__main__":unittest.main()
