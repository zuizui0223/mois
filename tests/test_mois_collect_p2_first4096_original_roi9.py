"""No new colours from unknown original ROI status; every P2 taxon adds breadth."""
import tempfile
import unittest
from pathlib import Path
import json

import numpy as np
import pandas as pd

from scripts.analysis.mois_collect_p2_first4096_original_roi9 import (
    classify_nine,summarize,attach_provenance,PALETTE,sha
)
from scripts.acquisition.mois_p2_first_4096_source import ARM

class StrictSourceP2OriginalROITests(unittest.TestCase):
    def test_continuous9_strict_preserves_technical_and_ambiguous_unknown(self):
        d=pd.DataFrame({
            "photo_id":["11","12","13","14"],"inat_taxon_id":["1","2","3","4"],
            "observation_id":["a","b","c","d"],"source_arm":[ARM]*4,
            "morph":["white","red_pink","mixed_uncertain","mixed_uncertain"],
            "measurement_status":["classified_four_state_morph",
                                  "classified_four_state_morph",
                                  "not_evaluable_ambiguous_palette_composition",
                                  "image_acquisition_failed"],
        })
        for x in PALETTE:d["flower_fraction_"+x]=0.0
        d.loc[0,"flower_fraction_white"]=1.0
        d.loc[1,"flower_fraction_pink"]=.8
        d.loc[1,"flower_fraction_white"]=.2
        d.loc[2,"flower_fraction_bronze"]=1.0
        for x in PALETTE:d.loc[3,"flower_fraction_"+x]=np.nan
        c=classify_nine(d)
        self.assertEqual(c.strict_original9.tolist(),[True,True,False,False])
        self.assertEqual(c.technical_four_hue_unknown.tolist(),
                         [False,False,True,True])
        self.assertAlmostEqual(c.loc[1,"strict9_white_fraction"],.2)
        self.assertTrue(pd.isna(c.loc[2,"strict9_bronze"]))
        self.assertTrue(pd.isna(c.loc[3,"strict9_white"]))
        self.assertEqual(int(c.photo_acquisition_failed.sum()),1)
        c["latitude_region"]=["N_tropical"]*4
        c["recorded_local_season"]=["summer"]*4
        groups=summarize(c)
        self.assertEqual(len(groups),1)
        self.assertEqual(groups[0]["new_original_photos"],4)
        self.assertEqual(groups[0]["P2_taxa_newly_strict9_coloured"],2)
        self.assertEqual(groups[0]["new_four_hue_unknown"],2)

    def test_source_provenance_matches_photo_and_taxon_without_guessing(self):
        with tempfile.TemporaryDirectory() as temporary:
            src=Path(temporary)
            p=pd.DataFrame({
                "photo_id":["101","102"],
                "original_observation_id":["ob1","ob2"],
                "inat_taxon_id":["21","22"],
                "source_gap_priority":["P2","P2"],
                "prior_strict9":[0,0],
                "latitude_region":["S_tropical","N_tropical"],
                "recorded_local_season":["winter","summer"],
                "prior_attempt_effort_bin":["E0_1","E2_4"],
            })
            f=src/"selected_P2_photo_source_provenance.csv"
            p.to_csv(f,index=False)
            (src/"manifest.json").write_text(json.dumps({
                "new_photo_attempts_planned":2,
                "photo_source_provenance_SHA256":sha(f),
            }))
            d=[{"photo_id":"102","observation_id":"ob2",
                "inat_taxon_id":"22","source_arm":ARM},
               {"photo_id":"101","observation_id":"ob1",
                "inat_taxon_id":"21","source_arm":ARM}]
            a=attach_provenance(d,src)
            self.assertEqual(len(a),2)
            self.assertEqual(a.inat_taxon_id.tolist(),["22","21"])
            self.assertEqual(a.latitude_region.tolist(),["N_tropical","S_tropical"])
            self.assertTrue(a.prior_strict9.eq(0).all())
            d[0]["inat_taxon_id"]="23"
            with self.assertRaisesRegex(ValueError,"differ"):
                attach_provenance(d,src)

if __name__=="__main__":
    unittest.main()
