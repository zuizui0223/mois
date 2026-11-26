"""Synthetic no-network tests for incremental GBIF COL XR accepted taxon identity."""
import unittest
from unittest.mock import patch

import pandas as pd

from scripts.analysis import mois_gbif_colxr_new_hue_incremental as m


class IncrementalGBIFAcceptedTaxonomy(unittest.TestCase):
    def source(self):
        d=pd.DataFrame([
            ("1","Plantus alba","white","provisional_single"),
            ("2","Plantus purpurea","purple","sampled_multicolour"),
            ("3","Plantus flava","yellow","provisional_single"),
            ("4","Plantus rosea","pink","provisional_single"),
            ("5","Plantus rubra","red","provisional_single")],
            columns=["inat_taxon_id","species","dominant_hue9","hue_status"])
        old=pd.DataFrame([
            {"inat_taxon_id":"1","species":"Plantus alba",
             "source_hue9":"white","hue_status":"provisional_single",
             "status":"colxr_exact_accepted_species","key":"ABCD",
             "COL_XR_checklistKey":m.COL_XR},
            {"inat_taxon_id":"2","species":"Plantus purpurea",
             "source_hue9":"purple","hue_status":"sampled_multicolour",
             "status":"synonym_or_doubtful_review","key":"",
             "COL_XR_checklistKey":m.COL_XR}])
        return d,old

    def test_new_only_cannot_relabel_prior_colour_taxa(self):
        atlas, old=self.source()
        with patch.object(m,"PRIOR_NAMED",2),patch.object(m,"PRIOR_EXACT",1),patch.object(m,"PRIOR_SOURCE_HUES",2),patch.object(m,"UPDATED_SOURCE_HUES",5):
            new=m.select_new_named(atlas,old)
        self.assertEqual(set(new.inat_taxon_id),{"3","4","5"})
        self.assertEqual(len(new),3)
        wrong=atlas.copy();wrong.loc[0,"dominant_hue9"]="orange"
        with patch.object(m,"PRIOR_NAMED",2),patch.object(m,"PRIOR_EXACT",1),patch.object(m,"PRIOR_SOURCE_HUES",2),patch.object(m,"UPDATED_SOURCE_HUES",5):
            with self.assertRaisesRegex(ValueError,"Known name/hue"):
                m.select_new_named(wrong,old)

    def test_collisions_review_not_double_count_accepted_species(self):
        _,old=self.source()
        add=pd.DataFrame([
            {"inat_taxon_id":"3","species":"Plantus flava",
             "source_hue9":"yellow","hue_status":"provisional_single",
             "status":"colxr_exact_accepted_species","key":"EFGH",
             "COL_XR_checklistKey":m.COL_XR},
            {"inat_taxon_id":"4","species":"Plantus rosea",
             "source_hue9":"pink","hue_status":"provisional_single",
             "status":"colxr_exact_accepted_species","key":"ABCD",
             "COL_XR_checklistKey":m.COL_XR},
            {"inat_taxon_id":"5","species":"Plantus rubra",
             "source_hue9":"red","hue_status":"provisional_single",
             "status":"no_match_or_unrecognized","key":"",
             "COL_XR_checklistKey":m.COL_XR}])
        full,good=m.combine_names(old,add)
        self.assertEqual(len(full),5)
        self.assertEqual(len(good),1)
        self.assertEqual(good.key.item(),"EFGH")
        self.assertEqual(int(full.status.eq("duplicate_accepted_species_key_review").sum()),2)
        self.assertEqual(int(full.status.eq("no_match_or_unrecognized").sum()),1)

    def test_binomial_source_colour_strata_and_taxon_identity(self):
        atlas,_=self.source()
        with patch.object(m,"EXPECTED_SOURCE_TAXA",5),patch.object(m,"UPDATED_SOURCE_HUES",5):
            x=m.exact_name_source(atlas)
            self.assertEqual(len(x),5)
            bad=atlas.copy();bad.loc[0,"inat_taxon_id"]="2"
            with self.assertRaisesRegex(ValueError,"duplicated"):
                m.exact_name_source(bad)


if __name__=="__main__":
    unittest.main()
