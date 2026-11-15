import unittest
from scripts.analysis.mois_gbif_37104_one_sql_validate import cube_sql,payload,COUNTS

class TestOneGBIFSQL(unittest.TestCase):
    def test_one_query_preserves_all_unique_taxon_keys(self):
        keys=['4Y3TB','ABC9','222A']
        s=cube_sql(keys)
        self.assertIn("specieskey IN ('222A', '4Y3TB', 'ABC9')",s)
        self.assertIn("GBIF_ISEA3HCode(8",s)
        self.assertIn("occurrenceStatus = 'PRESENT'",s)
        self.assertEqual(payload(keys)['format'],'SQL_TSV_ZIP')
    def test_duplicate_taxon_key_refused(self):
        with self.assertRaisesRegex(ValueError,"unsafe or duplicated"):
            cube_sql(['4Y3TB','4Y3TB'])
    def test_injected_species_key_refused(self):
        with self.assertRaisesRegex(ValueError,"unsafe or duplicated"):
            cube_sql(["4Y3TB') OR TRUE --"])
    def test_full_query_length_well_defined(self):
        keys=[f"A{n:04d}" for n in range(10000)]
        s=cube_sql(keys)
        self.assertEqual(s.count("'"),len(keys)*2+2)
        self.assertGreater(len(s),80000)
    def test_probing_full_before_smaller_fallback(self):
        self.assertEqual(COUNTS[0],37104)
        self.assertEqual(COUNTS[-1],250)
        self.assertTrue(all(a>b for a,b in zip(COUNTS,COUNTS[1:])))

if __name__=='__main__':unittest.main()
