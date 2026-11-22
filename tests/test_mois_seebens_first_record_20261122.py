"""Contract tests for Seebens initial time-axis source integration."""
import unittest
import pandas as pd
from scripts.analysis.mois_seebens_v31_first_record_link_20261122 import countries,link,dataset
class FirstRecords(unittest.TestCase):
    def test_not_fuzzy_countries(self):
        got=countries(pd.Series(["Germany","United States","California","Isle of Wight"]))
        self.assertEqual(got.tolist()[:2],["DE","US"])
        self.assertTrue(pd.isna(got.iloc[2]))
        self.assertTrue(pd.isna(got.iloc[3]))
    def test_original_date_interval_not_exact(self):
        s=pd.DataFrame({"key":["a","b","c"],"species":["Acacia alpha","Acacia beta","Acacia gamma"],
                        "source_hue9":["white","yellow","blue"],
                        "hue_status":["provisional_single"]*3,"occupied_grid_cells":[1,3,9]})
        f=pd.DataFrame({"Region":["Germany"]*3,"LifeForm":["Vascular plants"]*3,
                        "TaxonName":s.species,"FirstRecord":["1800","1905","2000"],
                        "FirstRecord_orig":["1800","1900-1910","2000"],
                        "PresentStatus":["established"]*3,"Phylum":["Tracheophyta"]*3,
                        "Source":["x"]*3,"DataUsage":["Free"]*3})
        _,x,audit=link(s,f,min_matches=1)
        self.assertEqual(int(x.original_exact.sum()),2)
        self.assertEqual(audit["MOIS_matched_country_original_exact_year"],2)
        self.assertEqual(dataset(x)[0].shape[0],2)
if __name__=="__main__":unittest.main()
