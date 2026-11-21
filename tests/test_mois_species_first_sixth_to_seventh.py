"""Signed sixth ROI-v4 to seventh species breadth regression tests."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd

from scripts.analysis import mois_species_first_sixth_to_seventh as w
from scripts.analysis.mois_species_first_palette9 import PALETTE,sha256

class SignedSixthRealPhoto(unittest.TestCase):
    def make_fixture(self,root):
        source=root/"source";(source/"shards").mkdir(parents=True)
        worker=root/"workers"
        raw=[];signatures={}
        for i in range(16):
            photo=[str(1000000+i*256+j) for j in range(256)]
            taxon=[str(50000+i*256+j) for j in range(256)]
            obs=[f"obs-sixth-{i}-{j}" for j in range(256)]
            original=pd.DataFrame({"photo_id":photo,"inat_taxon_id":taxon,
                                   "observation_id":obs,
                                   "morph":["white" if i*256+j<1406 else "unknown" for j in range(256)],
                                   "measurement_status":["measured" if i*256+j<1406 else "not_detected" for j in range(256)],
                                   "image_sha256":["test-image"]*256})
            for hue in PALETTE:
                values=[]
                for j in range(256):
                    k=i*256+j
                    values.append(("1" if hue=="white" else "0") if k<1087
                                  else (("0.5" if hue in ("white","yellow") else "0") if k<1406 else ""))
                original["flower_fraction_"+hue]=values
            planned=original[["photo_id","inat_taxon_id","observation_id"]].copy()
            planned["shard"]=str(i)
            src=source/"shards"/f"shard-{i:03d}.csv";planned.to_csv(src,index=False)
            f=worker/f"mois-sixth4096-genuine-original-worker-{i//4}"/f"shard_{i:03d}"
            f.mkdir(parents=True)
            p=f/"photos.csv";original.to_csv(p,index=False)
            r=f/"receipt.json"
            r.write_text(json.dumps({
                "schema":"mois_full_original93m_real_RGB_photo_unit_v1",
                "shard":i,"attempted":256,"source_sha256":sha256(src),
                "photos_sha256":sha256(p),"photo_colours_measured":True}))
            signatures[str(i)]=sha256(r);raw.append(original)
        (source/"manifest.json").write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "source_arm":"MOIS_SPECIES_FIRST_SIXTH_NEW4096",
            "colour_measured":False,"work_units":16,
            "worker_ids":list(range(16)),"new_photo_attempts_planned":4096,
            "source_previous_attempts":w.PRIOR_SOURCE,
            "photo_shard_sha256":{str(i):sha256(source/"shards"/f"shard-{i:03d}.csv") for i in range(16)}}))
        derivative=pd.concat(raw,ignore_index=True)
        derivative["strict9_original_photo_numeric"]=[str(i<1406) for i in range(4096)]
        derivative["dominant_hue9_label"]=["white" if i<1087 else "" for i in range(4096)]
        derivative.loc[1406:,["flower_fraction_"+p for p in PALETTE]]=""
        audited=root/"audited";audited.mkdir()
        derivative.to_csv(audited/"mois_sixth4096_actual_original_nine_roi.csv.gz",
                          index=False,compression="gzip")
        (audited/"receipt.json").write_text(json.dumps({
            "schema":w.AUDIT_SCHEMA,"actual_independent_original_photo_attempts":4096,
            "actual_genuine_strict9_numeric_success":1406,
            "actual_clear_nine_palette_hues":1087,"actual_technical_unknowns":2690,
            "photo_observation_source_taxon_ID_collisions":0,
            "sixteen_original_photo_worker_receipt_sha256":signatures}))
        return source,worker,audited

    def test_sixteen_authentic_worker_receipts_numeric_and_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);src,worker,audited=self.make_fixture(root)
            frame,proof=w.read_sixth(src,worker,audited)
            self.assertEqual(len(frame),4096)
            self.assertEqual(int(frame.strict9.sum()),1406)
            self.assertEqual(int(frame.dominant_label9.ne("").sum()),1087)
            self.assertTrue(frame.loc[~frame.strict9,["strict9_"+h for h in PALETTE]].isna().all().all())
            self.assertEqual(len(proof["sixteen_original_photo_worker_receipt_sha256"]),16)
            file=worker/"mois-sixth4096-genuine-original-worker-0"/"shard_000"/"receipt.json"
            file.write_text(file.read_text()+" ")
            with self.assertRaisesRegex(ValueError,"SHA changed"):
                w.read_sixth(src,worker,audited)

    def test_append_only_zero_hue_unseen_priority(self):
        taxa=[str(i) for i in range(1,42112)]
        n=31618;label=29168
        d=pd.DataFrame({
            "inat_taxon_id":taxa,
            "hue_status":["provisional_single"]*label+["mixed_image_unresolved"]*(n-label)+["no_strict9"]*(42111-n),
            "dominant_hue9":["white"]*label+[""]*(42111-label),
            "observed_hue9_set":["white"]*label+[""]*(42111-label),
            "after_strict9":[1]*n+[0]*(42111-n),
            "after_photo_attempts":[1]*n+[0]*(42111-n),
            "independent_roi9_observations":[1]*n+[0]*(42111-n),
            "distinct_palette_votes":[1]*label+[0]*(42111-label),
            "ambiguous_observation_rois":[0]*42111,
            "next_measure_priority":[2]*label+[1]*(n-label)+[0]*(42111-n)
        })
        for p in PALETTE:d["observed_mean_"+p]=[1.0 if p=="white" else 0.0]*n+[np.nan]*(42111-n)
        fresh=taxa[n:n+4096]
        sixth=pd.DataFrame({
            "inat_taxon_id":fresh,"photo_id":[str(800000+i) for i in range(4096)],
            "observation_id":[f"6obs-{i}" for i in range(4096)],
            "strict9":[i<1406 for i in range(4096)],
            "dominant_label9":["yellow" if i<1087 else "" for i in range(4096)]
        })
        for p in PALETTE:
            sixth["strict9_"+p]=[1.0 if p=="yellow" else 0.0]*1087 + [0.5 if p in ("white","yellow") else 0.0]*319 + [np.nan]*2690
        old=pd.DataFrame({"photo_id":["1","2"],"observation_id":["oldobs-1","oldobs-2"],"inat_taxon_id":["1","2"]})
        others=taxa[n+4096:n+4096+16]
        options=pd.DataFrame({"inat_taxon_id":others,
                              "photo_id":[str(900000+i) for i in range(16)],
                              "observation_id":[f"7obs-{i}" for i in range(16)],
                              "source_observation_rank":[1]*16})
        with patch.object(w,"PRIOR_SOURCE",2):
            updated,queue=w.append_and_select(d,sixth,old,options,quota=8)
        self.assertEqual(int(updated.after_strict9.gt(0).sum()),33024)
        self.assertEqual(int(updated.dominant_hue9.ne("").sum()),30255)
        self.assertTrue(queue.acquisition_priority.eq("first_unmeasured").all())
        self.assertFalse(set(queue.inat_taxon_id)&set(fresh))

if __name__=="__main__":
    unittest.main()
