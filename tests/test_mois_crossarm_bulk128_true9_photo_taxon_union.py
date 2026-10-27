"""Cross-arm actual original photo/observation IDs and nine-palette missingness.

All fixtures are synthetic local records: no S3, image decoding, or model inference.
"""
import unittest
import pandas as pd
import numpy as np

from scripts.analysis.mois_crossarm_bulk128_true9_photo_taxon_union import (
    merge_arms,supplement,validated_arm,PALETTE,FINAL_COLOURS,SOURCE
)

class CrossArmOriginalPhotoTaxonTests(unittest.TestCase):
    @staticmethod
    def canonical(pids,taxa,obs,valid,*,source_design="original_1043930"):
        data=[]
        for pid,taxon,oid,strict in zip(pids,taxa,obs,valid):
            d=dict(photo_id=str(pid),inat_taxon_id=str(taxon),
                observation_id=oid,source_design=source_design,
                source_subarm="legacy" if source_design==SOURCE else "targeted",
                observed_on="",latitude=np.nan,longitude=np.nan,
                strict9=bool(strict),
                original_roi_four_hue="white" if strict else "mixed_uncertain")
            for k in PALETTE:
                d["strict9_"+k]=1.0 if strict and k=="white" else np.nan
            data.append(d)
        return pd.DataFrame(data)

    def register(self):
        return pd.DataFrame(dict(
            inat_taxon_id=["10","11","12","13","14"],
            n_unique_source_photo_ids=[25,50,75,100,200],
            n_strict_roi9=[1,0,0,0,0],
            priority=[
                "P3_thin_colour_1to4",
                "P2_source_has_photos_no_colour_10plus",
                "P2_source_has_photos_no_colour_10plus",
                "P2_source_has_photos_no_colour_10plus",
                "P2_source_has_photos_no_colour_10plus",
            ]
        ))

    def arms(self):
        baseline=self.canonical([1100,1101,1102],["10","10","11"],
                                ["old","old","old-2"],[True,False,False])
        geo=self.canonical([2100,2101],["11","12"],["n1","n2"],
                           [True,False],source_design="bias_refill_2048")
        p2=self.canonical([3100,3101],["13","11"],["n3","n4"],
                          [True,True],source_design="p2_first_4096")
        bulk=self.canonical([4100,4101],["12","14"],["n5","n6"],
                            [True,True],source_design="bulk128_131072")
        return baseline,geo,p2,bulk

    def test_photo_level_union_not_sum_of_newly_covered_taxa_per_arm(self):
        joined,taxa,r=merge_arms(*self.arms(),self.register(),
                                 verify_source_counts=False)
        self.assertEqual(len(joined),9)
        self.assertEqual(joined.photo_id.nunique(),9)
        self.assertEqual(int(joined.strict9.sum()),5)
        self.assertEqual(r["original_source_nominal_taxa_with_strict9_before"],1)
        self.assertEqual(r["original_source_nominal_taxa_with_strict9_after"],5)
        self.assertEqual(r["net_newly_measured_nominal_source_taxa"],4)
        self.assertEqual(r["nominal_source_taxa_without_any_strict9_after"],0)
        self.assertEqual(r["baseline_legacy_photo_rows_reusing_same_observation"],1)
        self.assertEqual(sum(taxa.newly_colour_supported),4)
        self.assertEqual(taxa.loc[taxa.inat_taxon_id=="11","after_strict9"].item(),2)
        self.assertEqual(len(joined.loc[~joined.strict9]),4)
        self.assertTrue(joined.loc[~joined.strict9,list(FINAL_COLOURS)]
                        .isna().all().all())

    def test_both_photo_overlap_and_new_observation_reuse_are_blocked(self):
        base,geo,p2,bulk=self.arms()
        p2.loc[0,"photo_id"]=geo.loc[0,"photo_id"]
        with self.assertRaisesRegex(ValueError,"Exact original photo ID appears"):
            merge_arms(base,geo,p2,bulk,self.register(),verify_source_counts=False)
        base,geo,p2,bulk=self.arms()
        bulk.loc[0,"observation_id"]=geo.loc[0,"observation_id"]
        with self.assertRaisesRegex(ValueError,"previously attempted another new source observation"):
            merge_arms(base,geo,p2,bulk,self.register(),verify_source_counts=False)

    def test_source_false9_and_malformed_numeric_vector_are_not_strict(self):
        bad=self.canonical([9999],["15"],["obs"],[False])
        bad.loc[0,"strict9_white"]=1.0
        clean=validated_arm(bad,SOURCE)
        self.assertTrue(clean[list(FINAL_COLOURS)].isna().all().all())
        bad.loc[0,"strict9"]=True
        with self.assertRaisesRegex(ValueError,"Strict9 measured in a four-hue-unknown"):
            validated_arm(bad,SOURCE)
        bad.loc[0,"original_roi_four_hue"]="white"
        bad.loc[0,"strict9_white"]=np.nan
        with self.assertRaisesRegex(ValueError,"numeric vector invalid"):
            validated_arm(bad,SOURCE)

    def test_supplement_must_preserve_actual_boolean_strict9(self):
        row={"photo_id":"101","observation_id":"uuid-101","inat_taxon_id":"404",
             "continuous9_strict":True,"original_four_hue":"white",
             "sampling_design_cohort":"reference"}
        for name in PALETTE:
            row["continuous9_"+name]=1.0 if name=="white" else 0.
        got=supplement(pd.DataFrame([row]),"bias_refill_2048")
        self.assertEqual(got.source_design.iloc[0],"bias_refill_2048")
        self.assertEqual(got.source_subarm.iloc[0],"reference")
        self.assertTrue(got.strict9.iloc[0])
        with self.assertRaisesRegex(ValueError,"genuine measured boolean"):
            bad=pd.DataFrame([{**row,"continuous9_strict":"False"}])
            validated_arm(supplement(bad,"bias_refill_2048"),"bias_refill_2048")

if __name__=="__main__":
    unittest.main()
