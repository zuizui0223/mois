"""Synthetic no-API tests for COL XR incremental names after 11 signed new original photos."""
import unittest
from unittest.mock import patch

import pandas as pd
from scripts.analysis import mois_postsnapshot_11_original_new2_COLXR as m


class RealPostSnapshotNewNames(unittest.TestCase):
    def fixtures(self):
        a=pd.DataFrame([
            {"inat_taxon_id":"1","species":"Plantus alba","dominant_hue9":"white",
             "hue_status":"provisional_single"},
            {"inat_taxon_id":"2","species":"Plantus rubra","dominant_hue9":"red",
             "hue_status":"provisional_single"},
            {"inat_taxon_id":"3","species":"Plantus flava","dominant_hue9":"yellow",
             "hue_status":"provisional_single"},
            {"inat_taxon_id":"4","species":"Plantus noflower","dominant_hue9":"",
             "hue_status":"no_strict9"}])
        p=pd.DataFrame([{
            "inat_taxon_id":"1","species":"Plantus alba",
            "source_hue9":"white","hue_status":"provisional_single",
            "COL_XR_checklistKey":m.COL_XR,
            "status":"colxr_exact_accepted_species","key":"g1"}])
        return a,p

    def test_match_only_new_clear_original_photo_species(self):
        a,p=self.fixtures()
        with patch.object(m,"TAXA",4),patch.object(m,"AFTER_HUES",3),             patch.object(m,"PRIOR_HUES",1),patch.object(m,"PRIOR_ACCEPTED",1):
            extra=m.select_new(a,p)
        self.assertEqual(set(extra.inat_taxon_id),{"2","3"})
        self.assertTrue(extra.dominant_hue9.isin(m.PALETTE).all())

    def test_prior_hue_cannot_be_rewritten(self):
        a,p=self.fixtures()
        a.loc[0,"dominant_hue9"]="purple"
        with patch.object(m,"TAXA",4),patch.object(m,"AFTER_HUES",3),             patch.object(m,"PRIOR_HUES",1),patch.object(m,"PRIOR_ACCEPTED",1):
            with self.assertRaisesRegex(ValueError,"drifted"):
                m.select_new(a,p)

    def test_outside_accepted_binomial_not_promoted(self):
        a,p=self.fixtures()
        a.loc[2,"species"]="unresolved taxon sp."
        with patch.object(m,"TAXA",4),patch.object(m,"AFTER_HUES",3),             patch.object(m,"PRIOR_HUES",1),patch.object(m,"PRIOR_ACCEPTED",1):
            selected=m.select_new(a,p)
        self.assertEqual(len(selected),1)
        self.assertEqual(selected.inat_taxon_id.iloc[0],"2")


if __name__=="__main__":
    unittest.main()
