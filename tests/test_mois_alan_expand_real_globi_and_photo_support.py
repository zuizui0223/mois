"""GloBI incremental source coverage tests: only existing visits, no nocturnal imputation."""
import unittest
from io import StringIO
import pandas as pd
from scripts.analysis import mois_alan_expand_real_globi_and_photo_support as s

class MoreGloBISourceTests(unittest.TestCase):
    def minimal(self):
        return pd.DataFrame([
           {"plant_name":"Silene latifolia","guild":"moth_family_visit"},
           {"plant_name":"Silene vulgaris","guild":"moth_family_visit"},
           {"plant_name":"Daucus carota","guild":"butterfly_family_visit"}])
    def test_recorded_plant_with_both_visitors_is_not_moth_only(self):
        d=pd.concat([self.minimal(),pd.DataFrame(
             [{"plant_name":"Silene latifolia","guild":"butterfly_family_visit"}])])
        q=s.classify_groups(d)
        self.assertEqual(q["Silene latifolia"],"both")
        self.assertEqual(q["Silene vulgaris"],"moth_only")
        self.assertEqual(q["Daucus carota"],"butterfly_only")
    def test_no_GloBI_record_is_not_assigned_daytime(self):
        q=s.classify_groups(self.minimal())
        self.assertNotIn("Nicotiana tabacum",q.index)
    def test_unknown_guild_source_fails_closed(self):
        d=self.minimal().copy();d.loc[0,"guild"]="nocturnal_exclusive"
        with self.assertRaisesRegex(ValueError,"unverified"):
             s.classify_groups(d)
    def test_only_species_direct_known_colour_is_counted(self):
        d=pd.DataFrame({"photo_id":["1","2","3"],"inat_taxon_id":["x","y","z"],
          "observation_id":["o1","o2","o3"],"source_subarm":["historical_main"]*3,
          "strict9":[True]*3,"latitude":[36.1]*3,"longitude":[140.1]*3,
          "observed_on":["2020-04-01"]*3,"wcvp_accepted_name":["Silene latifolia","Daucus carota","Silene vulgaris"],
          "strict9_white":[.5,.2,.4]})
        mapping=d[["inat_taxon_id","wcvp_accepted_name"]]
        d=d.drop(columns=["wcvp_accepted_name"])
        colours=pd.DataFrame([
            {"accepted_species":"Silene latifolia","trait_name":"flower_primary_color","quality":"high","resolution_status":"resolved","classification":"direct","normalized_value":"white"},
            {"accepted_species":"Daucus carota","trait_name":"flower_primary_color","quality":"high","resolution_status":"resolved","classification":"direct","normalized_value":"yellow_orange"},
            {"accepted_species":"Silene vulgaris","trait_name":"flower_primary_color","quality":"high","resolution_status":"resolved","classification":"direct","normalized_value":"mixed_or_variable"}])
        x,rec=s.assess(self.minimal(),mapping,colours,d)
        self.assertTrue(rec["unknown_species_colour_not_imputed"])
        q=x.loc[(x.period=="2020")&(x.grid_degrees==.25)].iloc[0]
        self.assertEqual(int(q["cells_with_both_visit_groups_ge1"]),1)
        self.assertEqual(int(q["cells_both_groups_independently_documented_white_ge1"]),0)
        self.assertEqual(int(q["cells_four_direct_colour_and_visit_groups_ge1"]),0)
    def test_no_photo_status_promoted_from_duplicate_guild_rows(self):
        c=pd.concat([self.minimal(),self.minimal()],ignore_index=True)
        q=s.classify_groups(c)
        self.assertEqual(len(q),3)
    def test_bad_source_PINS_refuse_unverified_data(self):
        self.assertEqual(len(s.SHA["status"]),64)
        self.assertEqual(len(s.SHA["canonical"]),64)
        self.assertEqual(len(s.SHA["globi"]),64)
    def test_photo_daynight_is_not_inferred_from_moth_visit(self):
        self.assertNotIn("nocturnal",s.classify_groups(self.minimal()).unique().tolist())
    def test_expansion_parameters_guard(self):
        with self.assertRaisesRegex(ValueError,"courteous"):
            s.extend(self.minimal(),pd.DataFrame({"family":["Sphingidae"],"relation":["visitsFlowersOf"],
               "page_cap_possible_truncation":[True]}),workers=0)

if __name__=="__main__":
    unittest.main()
