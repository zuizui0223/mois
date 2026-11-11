"""Offline checks for GBIF accepted-species nine-fraction trait preservation."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_gbif_colxr_nine_colour_traits_bridge import (
    build_traits,PALETTE
)


class NinePaletteGbifBridgeTests(unittest.TestCase):
    def fixture(self):
        rows=[]
        states=[("sampled_multicolour","yellow","yellow|white"),
                ("provisional_single","white","white")]
        for i,(status,hue,multi) in enumerate(states):
            d={"inat_taxon_id":str(i+1),"species":"Plantus "+chr(97+i),
               "hue_status":status,"dominant_hue9":hue,
               "observed_hue9_set":multi,"after_strict9":2,
               "after_photo_attempts":5,"independent_roi9_observations":2,
               "distinct_palette_votes":2 if i==0 else 1,
               "ambiguous_observation_rois":0}
            d.update({"observed_mean_"+h:(1.0 if h==hue else 0.) for h in PALETTE})
            rows.append(d)
        # Authentic atlas shape: all other 42,109 source taxa have no numerical ROI.
        for j in range(3,42112):
            r={"inat_taxon_id":str(j),"species":"Plantus unnamed",
               "hue_status":"no_strict9","dominant_hue9":"",
               "observed_hue9_set":"","after_strict9":0,
               "after_photo_attempts":0,"independent_roi9_observations":0,
               "distinct_palette_votes":0,"ambiguous_observation_rois":0}
            r.update({"observed_mean_"+h:np.nan for h in PALETTE})
            rows.append(r)
        atlas=pd.DataFrame(rows)
        matches=pd.DataFrame({"inat_taxon_id":["1","2"],
            "species":["Plantus a","Plantus b"],
            "source_hue9":["yellow","white"],
            "hue_status":["sampled_multicolour","provisional_single"],
            "key":["ABCD","EFGH"],"COL_XR_checklistKey":[
                "7ddf754f-d193-4cc9-b351-99906754a03b"]*2,
            "status":["colxr_exact_accepted_species"]*2})
        return atlas,matches

    def test_never_force_multicolour_taxon_into_single_hue(self):
        atlas,matches=self.fixture()
        d=build_traits(atlas,matches,verify_counts=False)
        self.assertEqual(len(d),2)
        self.assertEqual(d.hue_evidence_class.value_counts().to_dict(),
                         {"sampled_multicolour":1,"sampled_single_hue":1})
        self.assertEqual(d.loc[d.inat_taxon_id=="1","photo_observed_palette9_states"].iloc[0],
                         "yellow|white")
        self.assertFalse(d.all_range_hue_fixed_proven.any())
        self.assertFalse(d.hue_measured_at_gbif_occurrence.any())
        self.assertEqual(len([c for c in d.columns if c.startswith("observed_mean_")]),9)

    def test_mismatched_source_species_colour_or_colliding_colxr_key_refused(self):
        atlas,matches=self.fixture()
        bad=matches.copy();bad.loc[0,"source_hue9"]="red"
        with self.assertRaisesRegex(ValueError,"hue changed"):
            build_traits(atlas,bad,verify_counts=False)
        bad=matches.copy();bad.loc[0,"key"]="EFGH"
        with self.assertRaisesRegex(ValueError,"not exact accepted"):
            build_traits(atlas,bad,verify_counts=False)


if __name__=="__main__":
    unittest.main()
