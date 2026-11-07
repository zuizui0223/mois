"""Source identity, ROI technical unknown, photo-level vs medium500 evidence."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_bulk128_integrate_original_roi9_and_medium_lab import (
    integrate_frames,strict_nine,summarize_strata,FRAC,REFERENCE,TARGETED
)

class Bulk128PixelJoinTests(unittest.TestCase):
    def originals(self):
        records=[]
        for i in range(4):
            p={
                "photo_id":str(70000+i),
                "observation_id":f"uuid-{i}",
                "inat_taxon_id":str(200+i),
                "source_arm":REFERENCE if i==0 else TARGETED,
                "morph":["white","red_pink","mixed_uncertain","mixed_uncertain"][i],
                "measurement_status":[
                    "classified_four_state_morph","classified_four_state_morph",
                    "not_evaluable_ambiguous_palette_composition","image_acquisition_failed"][i],
                "image_sha256":"a"*64 if i<3 else "",
            }
            for k in FRAC:p[k]="0"
            p[FRAC[0]]="1" if i==0 else ("0.2" if i==1 else ("1" if i==2 else ""))
            p[FRAC[4]]="0.8" if i==1 else ("0" if i!=3 else "")
            records.append(p)
        return pd.DataFrame(records)
    def medium(self):
        r=[]
        for i in range(4):
            r.append(dict(photo_id=str(70000+i),observation_id=f"uuid-{i}",
                inat_taxon_id=str(200+i),design_cohort="reference" if i==0 else "targeted",
                original_priority="P5_existing_10plus_strict9" if i==0 else
                    "P2_source_has_photos_no_colour_10plus",
                latitude_region="S_tropical",local_season="winter",
                screen_status=["negative","candidate","candidate","error"][i],
                medium_decoded="True" if i!=3 else "False",
                floral_box_main_L_star="61.4" if i in (1,2) else "",
                floral_box_main_C_star="20" if i in (1,2) else "",
                floral_box_main_rgb_hex="#CC0055" if i in (1,2) else "",
                unverified_floral_box_not_petal_ROI="True"))
        return pd.DataFrame(r)

    def test_medium_negative_cannot_delete_strict_full_roi9(self):
        d=integrate_frames(self.medium(),self.originals())
        self.assertEqual(len(d),4)
        self.assertEqual(d.strict_original_roi9.tolist(),[True,True,False,False])
        self.assertEqual(int(d.medium_screen_negative_with_strict_original_roi9.sum()),1)
        self.assertEqual(int(d.medium_screen_candidate_but_no_strict_original_roi9.sum()),1)
        self.assertEqual(int(d.newly_colour_supported_taxon_in_frozen_P1P2_targeted.sum()),1)
        self.assertTrue(pd.isna(d.loc[2,"roi9_white"]))
        self.assertTrue(pd.isna(d.loc[3,"roi9_red"]))
        self.assertAlmostEqual(d.loc[1,"roi9_white"],0.2)
        out=summarize_strata(d,["design_cohort"])
        by={x["design_cohort"]:x for x in out}
        self.assertEqual(by["reference"]["medium_screen_negative_but_strict_roi9"],1)
        self.assertEqual(by["targeted"]["unique_P1P2_prior_no_colour_taxa_newly_strict_roi9"],1)

    def test_conflicting_source_taxon_or_observation_id_fails_closed(self):
        d=self.medium();x=self.originals()
        x.loc[1,"inat_taxon_id"]="999"
        with self.assertRaisesRegex(ValueError,"disagreeing observation or plant taxon"):
            integrate_frames(d,x)
        x=self.originals()
        x.loc[0,"observation_id"]="some-other-observation"
        with self.assertRaisesRegex(ValueError,"disagreeing observation or plant taxon"):
            integrate_frames(d,x)

    def test_no_technical_photo_failure_can_be_strict_roi9(self):
        x=self.originals()
        x.loc[3,"morph"]="white"
        x.loc[3,FRAC[0]]="1"
        with self.assertRaisesRegex(ValueError,"acquisition failure coloured"):
            strict_nine(x)

    def test_unverified_box_and_source_design_fail_closed(self):
        medium=self.medium()
        medium.loc[1,"design_cohort"]="reference"
        with self.assertRaisesRegex(ValueError,"role differs"):
            integrate_frames(medium,self.originals())
        medium=self.medium()
        medium.loc[1,"unverified_floral_box_not_petal_ROI"]="False"
        with self.assertRaisesRegex(ValueError,"not source-verifiable"):
            integrate_frames(medium,self.originals())

if __name__=="__main__":
    unittest.main()
