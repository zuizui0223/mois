"""Frozen ROI9 observation identity is outcome blind; no manual review needed."""
import unittest
import pandas as pd
from scripts.analysis.mois_roi9_observation_frame_20261029 import (
    observation_frame, coarse_zone, canonical_rank,
)

class NoHumanObservationIdentityTests(unittest.TestCase):
    def data(self):
        rows=[]
        for i in range(20):
            rows.append({
                "photo_id":str(5000+i),
                "observation_id":"obs-"+str(i//2),
                "inat_taxon_id":str(100+i//2),
                "source_design":"historical_K1_FCP_source_heterogeneous",
                "observed_on":"2025-07-20",
                "latitude":-13.0,"longitude":121.0,
                "strict9":bool(i%3==0),
            })
        return pd.DataFrame(rows)
    def test_canonical_same_across_order_or_result_changes(self):
        d=self.data()
        a,s=observation_frame(d)
        self.assertEqual(len(a),10)
        self.assertEqual(s["additional_photos_of_already_represented_observations"],10)
        original=dict(zip(a.observation_id,a.photo_id))
        b,_=observation_frame(d.sample(frac=1,random_state=123))
        self.assertEqual(dict(zip(b.observation_id,b.photo_id)),original)
        d["strict9"]=~d.strict9
        c,_=observation_frame(d)
        self.assertEqual(dict(zip(c.observation_id,c.photo_id)),original)
        self.assertEqual(coarse_zone("S_tropical"),"S_tropical")
        self.assertNotEqual(canonical_rank("obs-a","1"),canonical_rank("obs-a","2"))
    def test_missing_ids_never_replaced_by_fabricated_observations(self):
        d=self.data()
        d.loc[0,"observation_id"]=""
        a,s=observation_frame(d)
        self.assertEqual(s["photos_without_observation_uuid"],1)
        self.assertEqual(len(a),10)
        d.loc[0,"photo_id"]=d.loc[1,"photo_id"]
        with self.assertRaisesRegex(ValueError,"unique"):
            observation_frame(d)

if __name__=="__main__":
    unittest.main()
