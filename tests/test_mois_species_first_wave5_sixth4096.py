"""Synthetic original-photo tests for fifth→sixth species-first 4096 target pipeline."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_wave5_sixth4096 import (
    read_fifth, append_fifth_and_choose_sixth
)
from scripts.analysis.mois_species_first_palette9 import PALETTE,sha256


class FifthWaveAndSixthSpeciesBreadth(unittest.TestCase):
    def test_four_signed_256_original_photo_terminals_and_unknowns(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            src=root/"source";result=root/"measured"
            (src/"shards").mkdir(parents=True)
            sha={}
            for n in range(4):
                unit=pd.DataFrame({
                    "photo_id":[str(100000+n*256+i) for i in range(256)],
                    "inat_taxon_id":[str(30000+n*256+i) for i in range(256)],
                    "observation_id":[f"obs-{n}-{i}" for i in range(256)],
                })
                unit["shard"]=str(n)
                original=unit.copy()
                original["morph"]=["white"]+["mixed_uncertain"]*255
                for p in PALETTE:
                    original["flower_fraction_"+p]=[
                        "1.0" if p=="white" else "0.0"]+[""]*255
                path=src/"shards"/f"shard-{n:03d}.csv"
                path.parent.mkdir(parents=True,exist_ok=True)
                unit.to_csv(path,index=False)
                sha[str(n)]=sha256(path)
                folder=result/f"shard_{n:03d}";folder.mkdir(parents=True)
                out=folder/"photos.csv";original.to_csv(out,index=False)
                (folder/"receipt.json").write_text(json.dumps({
                    "schema":"mois_full_original93m_real_RGB_photo_unit_v1",
                    "shard":n,"attempted":256,
                    "source_sha256":sha256(path),
                    "photos_sha256":sha256(out),
                    "photo_colours_measured":True
                }))
            (src/"manifest.json").write_text(json.dumps({
                "schema":"mois_full_source_observation_sampling_cohort_v1",
                "source_arm":"MOIS_SPECIES_FIRST_FIFTH_NEW_ZERO9_1024",
                "new_photo_attempts_planned":1024,
                "source_previous_fourth_photo_attempts_plus_baseline":1184474,
                "source_prior_no_reused_original_photo_or_observation":True,
                "worker_ids":[0,1,2,3],
                "photo_shard_sha256":sha}))
            out,proof=read_fifth(src,result)
            self.assertEqual(len(out),1024)
            self.assertEqual(out.photo_id.nunique(),1024)
            self.assertEqual(int(out.strict9.sum()),4)
            self.assertEqual(int(out.dominant_label9.eq("white").sum()),4)
            self.assertEqual(int((~out.strict9).sum()),1020)
            self.assertEqual(len(proof),4)

    def test_only_new_taxa_are_reached_in_following_source(self):
        taxa=[str(i) for i in range(1,42112)]
        count=31235
        d=pd.DataFrame({
            "inat_taxon_id":taxa,
            "hue_status":["provisional_single"]*count+["no_strict9"]*(42111-count),
            "dominant_hue9":["white"]*count+[""]*(42111-count),
            "observed_hue9_set":["white"]*count+[""]*(42111-count),
            "after_strict9":[1]*count+[0]*(42111-count),
            "after_photo_attempts":[1]*count+[0]*(42111-count),
            "independent_roi9_observations":[1]*count+[0]*(42111-count),
            "distinct_palette_votes":[1]*count+[0]*(42111-count),
            "ambiguous_observation_rois":[0]*42111,
            "next_measure_priority":[2]*count+[0]*(42111-count)
        })
        for p in PALETTE:d["observed_mean_"+p]=[1.0 if p=="white" else 0.0]*count+[np.nan]*(42111-count)
        newtaxa=taxa[count:count+1024]
        new=pd.DataFrame({
            "inat_taxon_id":newtaxa,
            "photo_id":[str(10000+i) for i in range(1024)],
            "observation_id":[f"newobs-{i}" for i in range(1024)],
            "strict9":[i<360 for i in range(1024)],
            "dominant_label9":["yellow" if i<360 else "" for i in range(1024)]
        })
        for p in PALETTE:new["strict9_"+p]=[1.0 if p=="yellow" else 0.0]*360+[np.nan]*664
        old=pd.DataFrame({"photo_id":["1","2"],"observation_id":["old-1","old-2"],
                          "inat_taxon_id":["1","2"]})
        opt=pd.DataFrame({
            "inat_taxon_id":taxa[count+1024:count+1024+8],
            "photo_id":[str(20000+i) for i in range(8)],
            "observation_id":[f"later-{i}" for i in range(8)],
            "source_observation_rank":[1]*8
        })
        atlas,queue=append_fifth_and_choose_sixth(d,new,old,opt,quota=4,
                                                  expected_previous=2)
        self.assertEqual(int(atlas.after_strict9.gt(0).sum()),count+360)
        self.assertEqual(len(queue),4)
        self.assertTrue(queue.acquisition_priority.eq("first_unmeasured").all())
        self.assertTrue(queue.inat_taxon_id.is_unique)
        self.assertFalse(set(queue.inat_taxon_id)&set(newtaxa))
        self.assertFalse(set(queue.photo_id)&set(new.photo_id))

    def test_previously_measured_original_photo_not_reused(self):
        from scripts.analysis.mois_species_first_wave5_sixth4096 import require
        with self.assertRaisesRegex(ValueError,"duplicate"):
            require(False,"duplicate original photo")


if __name__=="__main__":
    unittest.main()
