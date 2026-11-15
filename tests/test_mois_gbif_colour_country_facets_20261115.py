import unittest
import pandas as pd
from scripts.analysis.mois_gbif_colour_country_facets_20261115 import extract,choose,SEED,HUES

class TestGBIFCountryFacets(unittest.TestCase):
    def test_colxr_facets_are_georeferenced_record_counts_not_native_range(self):
        sample={"count":17,"results":[],"facets":[
            {"field":"COUNTRY","counts":[{"name":"JP","count":11},{"name":"US","count":5}]},
            {"field":"CONTINENT","counts":[{"name":"ASIA","count":11},{"name":"NORTH_AMERICA","count":5}]}]}
        found=extract(sample)
        self.assertEqual([x["name"] for x in found["country"]],["JP","US"])
        self.assertEqual(sum(x["count"] for x in found["country"]),16)
    def test_missing_facet_is_not_assumed_empty_distribution(self):
        with self.assertRaisesRegex(ValueError,"missing requested facet"):
            extract({"count":1,"results":[],"facets":[{"field":"country","counts":[{"name":"JP","count":1}]}]})
    def test_facet_cap_refuses_truncated_range(self):
        with self.assertRaisesRegex(ValueError,"potentially capped"):
            extract({"count":200,"results":[],"facets":[
                {"field":"country","counts":[{"name":str(i),"count":1} for i in range(200)]},
                {"field":"continent","counts":[]}]})
    def test_facet_sum_cannot_exceed_occurrence_count(self):
        with self.assertRaisesRegex(ValueError,"exceeds"):
            extract({"count":1,"results":[],"facets":[
                {"field":"country","counts":[{"name":"JP","count":2}]},
                {"field":"continent","counts":[]}]})
    def test_hash_blind_equal_hue_selection(self):
        rows=[]
        for i,h in enumerate(HUES):
            for j in range(4):
                rows.append({"source_hue9":h,"key":f"K{i}_{j}","inat_taxon_id":str(100*i+j),
                             "species":f"sp{i}_{j}"})
        d=pd.DataFrame(rows)
        a=choose(d,2);b=choose(d.sample(frac=1,random_state=2),2)
        self.assertEqual(len(a),18)
        self.assertEqual(a.key.tolist(),b.key.tolist())
        self.assertTrue(a.groupby("source_hue9").size().eq(2).all())

if __name__=="__main__":unittest.main()
