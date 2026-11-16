import unittest
import pandas as pd

from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
    match_gbif_to_wcvp, summarize_photo_species_by_tdwg
)

class TestWCVPNativeRegion(unittest.TestCase):
    def test_exact_name_only_no_ambiguous_wcvp_id(self):
        photo = pd.DataFrame({
            "key":["1","2","3","4"],
            "species":["Campanula punctata","Gentiana scabra",
                       "Viola mandshurica","Notfound species"],
            "source_hue9":["white","blue","yellow","pink"],
            "hue_status":["sampling_consensus_single"]*4
        })
        names=pd.DataFrame({
            "source_species_exact_binomial":[
                "Campanula punctata","Gentiana scabra","Viola mandshurica",
                "Viola mandshurica","Notfound species"],
            "wcvp_match_state":["unique_exact_name_to_accepted_id"]*4+
                               ["no_accepted_species_match"],
            "wcvp_accepted_plant_name_id":["10","11","12","13",""]
        })
        d=match_gbif_to_wcvp(photo,names,expected_species=4)
        self.assertEqual(d.wcvp_match_exact_binomial.tolist(),
                         [True,True,False,False])
        self.assertEqual(d.wcvp_accepted_plant_name_id.iloc[0],"10")
        self.assertTrue(pd.isna(d.wcvp_accepted_plant_name_id.iloc[2]))

    def test_conservative_native_introduced_species_counts(self):
        p=pd.DataFrame({
            "specieskey":["1","2","3","1","2","3"],
            "grid_x":[0,0,0,1,1,1],
            "grid_y":[0,0,0,0,0,0],
            "occurrence_count":[1,2,3,4,5,6]
        })
        sp=pd.DataFrame({
            "key":["1","2","3"],
            "source_hue9":["white","yellow","white"],
            "wcvp_accepted_plant_name_id":["10","11",pd.NA],
            "wcvp_match_exact_binomial":[True,True,False]})
        region=pd.DataFrame({
            "grid_x":[0,1],"grid_y":[0,0],
            "code":["AAA","BBB"],
            "region_strict":[True,True],"region_moderate":[True,True]})
        dist=pd.DataFrame({
            "wcvp_accepted_plant_name_id":["10","10","11","11"],
            "tdwg_level3_code":["AAA","BBB","AAA","BBB"],
            "native_introduced_state":["native_at_TDWG3","introduced_at_TDWG3",
                                       "introduced_at_TDWG3","native_at_TDWG3"]})
        d,a=summarize_photo_species_by_tdwg(p,sp,region,dist,"strict")
        rows=d.set_index("grid_x")
        self.assertEqual(rows.loc[0,"native_species_n"],1)
        self.assertEqual(rows.loc[0,"native_white_n"],1)
        self.assertEqual(rows.loc[1,"native_species_n"],1)
        self.assertEqual(rows.loc[1,"native_white_n"],0)
        self.assertEqual(a["uncertain_pairs"],2)
        self.assertEqual(d.GBIF_occurrence_count.sum(),21)

if __name__=="__main__":
    unittest.main()
