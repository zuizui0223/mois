"""GloBI pagination sentinel test: 1001 rows for limit=1000 means MORE."""
import unittest
import pandas as pd
from scripts.analysis import mois_alan_globi_stage2_uncap_20261105 as s

class SourceStage2Tests(unittest.TestCase):
    def sample(self):
        return pd.DataFrame([
            dict(family="Sphingidae",relation="visitsFlowersOf",offset=2000,
                 limit=1000,n_rows=893,failed=False),
            dict(family="Papilionidae",relation="visitsFlowersOf",offset=12000,
                 limit=1000,n_rows=1001,failed=False),
            dict(family="Papilionidae",relation="visitsFlowersOf",offset=13000,
                 limit=1000,n_rows=1001,failed=False),
            dict(family="Noctuidae",relation="visitsFlowersOf",offset=2000,
                 limit=1000,n_rows=1001,failed=True)])
    def test_limit_plus_one_is_truncated(self):
        d=s.find_capped(self.sample())
        self.assertEqual(len(d),1)
        self.assertEqual(d[0]["family"],"Papilionidae")
        self.assertEqual(d[0]["offset"],13000)
    def test_short_final_page_not_truncated(self):
        d=self.sample()
        d.loc[d.family.eq("Papilionidae")&d.offset.eq(13000),"n_rows"]=980
        self.assertEqual(s.find_capped(d),[])
    def test_source_SHA_pinned(self):
        self.assertEqual(len(s.SOURCE),64)
        self.assertEqual(len(s.PAGE),64)
    def test_extending_requires_real_capped_queries(self):
        old=pd.DataFrame({"plant_name":["Silene latifolia"],
             "guild":["moth_family_visit"],"source_family":["Sphingidae"],
             "relation":["visitsFlowersOf"],"insect_name":["Moth sp."],"study_title":["A"]})
        with self.assertRaisesRegex(ValueError,"capped provider"):
            s.add_pages(old,self.sample().loc[self.sample().family.eq("Sphingidae")])

if __name__=="__main__":unittest.main()
