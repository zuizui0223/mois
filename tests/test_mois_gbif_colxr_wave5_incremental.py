"""No-network tests of wave5 incremental GBIF name/colour identity."""
import unittest
from unittest.mock import patch
import pandas as pd

from scripts.analysis import mois_gbif_colxr_wave5_incremental as m


class Wave5GBIFIncrementTests(unittest.TestCase):
    def fixture(self):
        atlas=pd.DataFrame([
            {"inat_taxon_id":"1","species":"Plantus alba",
             "dominant_hue9":"white","hue_status":"provisional_single"},
            {"inat_taxon_id":"2","species":"Plantus flava",
             "dominant_hue9":"yellow","hue_status":"sampled_multicolour"},
            {"inat_taxon_id":"3","species":"Plantus rubra",
             "dominant_hue9":"red","hue_status":"provisional_single"},
            {"inat_taxon_id":"4","species":"Plantus rosea",
             "dominant_hue9":"pink","hue_status":"provisional_single"},
            {"inat_taxon_id":"5","species":"Plantus nodata",
             "dominant_hue9":"","hue_status":"no_strict9"}
        ])
        old=pd.DataFrame([
            {"inat_taxon_id":"1","species":"Plantus alba","source_hue9":"white",
             "hue_status":"provisional_single","status":"colxr_exact_accepted_species",
             "key":"AAAA","COL_XR_checklistKey":m.COL_XR},
            {"inat_taxon_id":"2","species":"Plantus flava","source_hue9":"yellow",
             "hue_status":"sampled_multicolour","status":"synonym_or_doubtful_review",
             "key":"","COL_XR_checklistKey":m.COL_XR}
        ])
        return atlas,old

    def test_preserves_existing_and_matches_only_new_hue_taxa(self):
        atlas,old=self.fixture()
        with patch.object(m,"SOURCE_TAXA",5),patch.object(m,"AFTER_WAVE5_HUE_LABELS",4),patch.object(m,"PRIOR_HUE_LABELS",2),patch.object(m,"PRIOR_ACCEPTED",1):
            new=m.select_incremental_hue_names(atlas,old)
        self.assertEqual(set(new.inat_taxon_id),{"3","4"})

    def test_prior_hue_change_is_error(self):
        atlas,old=self.fixture()
        atlas.loc[0,"dominant_hue9"]="blue"
        with patch.object(m,"SOURCE_TAXA",5),patch.object(m,"AFTER_WAVE5_HUE_LABELS",4),patch.object(m,"PRIOR_HUE_LABELS",2),patch.object(m,"PRIOR_ACCEPTED",1):
            with self.assertRaisesRegex(ValueError,"drifted"):
                m.select_incremental_hue_names(atlas,old)

    def test_duplicate_prior_colxr_taxa_are_rejected(self):
        atlas,old=self.fixture()
        old.loc[1,"status"]="colxr_exact_accepted_species"
        old.loc[1,"key"]="AAAA"
        with patch.object(m,"SOURCE_TAXA",5),patch.object(m,"AFTER_WAVE5_HUE_LABELS",4),patch.object(m,"PRIOR_HUE_LABELS",2),patch.object(m,"PRIOR_ACCEPTED",2):
            with self.assertRaisesRegex(ValueError,"duplicate"):
                m.select_incremental_hue_names(atlas,old)


if __name__=="__main__":
    unittest.main()
