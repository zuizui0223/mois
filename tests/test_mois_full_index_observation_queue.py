"""Actual FCP global source selection: full frame, nested inclusion, no photo quotas."""
from pathlib import Path
from unittest.mock import patch
import csv,gzip,hashlib,json,tempfile,unittest
from scripts.acquisition import mois_full_index_observation_queue as q

class TestFullSourceQueue(unittest.TestCase):
    def setUp(self):
        self.t=tempfile.TemporaryDirectory();self.root=Path(self.t.name)
        (self.root/"index"/"photos").mkdir(parents=True)
        (self.root/"index"/"observations").mkdir()
        original=[
            ("101","OBS-A","1","cc-by","jpg","user1"),
            ("102","OBS-A","1","cc-by","jpg","user1"),
            ("103","OBS-B","2","cc-by","jpg","user2"),
            ("104","OBS-C","2","cc-by","jpg","user3"),
        ]
        observations=[
            ("OBS-A","1","user1","2025-07-01","30","135","100"),
            ("OBS-B","2","user2","2025-05-02","43","142","200"),
            ("OBS-C","2","user3","2025-08-13","38","134","10")
        ]
        shards=[]
        for i in range(2):
            photo=self.root/"index"/"photos"/f"part-{i:03d}.csv.gz"
            obs=self.root/"index"/"observations"/f"part-{i:03d}.csv.gz"
            p=[r for r in original if int(r[2])%2==i]
            o=[r for r in observations if int(r[1])%2==i]
            for fn,header,data in (
                (photo,["photo_id","observation_uuid","inat_taxon_id","photo_license","image_extension","photographer_id"],p),
                (obs,["observation_uuid","inat_taxon_id","observer_id","observed_on","latitude","longitude","positional_accuracy"],o)):
                with gzip.open(fn,"wt",newline="") as f:
                    w=csv.writer(f);w.writerow(header);w.writerows(data)
            shards.append({"part":i,"photo_count":len(p),"observation_count":len(o),
                 "photos_sha256":q.sha(photo),"observations_sha256":q.sha(obs)})
        (self.root/"index"/"manifest.json").write_text(json.dumps({
            "schema":q.SCHEMA,"photos_emitted":4,
            "taxa_in_frame":2,"photo_caps":None,"photo_colours_measured":False,
            "shards":shards
        }))
        self.count=self.root/"counts.csv"
        self.count.write_text("inat_taxon_id,species\n1,Plant alpha\n2,Plant beta\n")
        self.prior=self.root/"prior.csv"
        self.prior.write_text("photo_id\n103\n")
    def tearDown(self):self.t.cleanup()
    def run_queue(self,p,output):
        with patch.object(q,"N_PHOTOS",4),patch.object(q,"N_TAXA",2),patch.object(q,"N_SHARDS",2):
            return q.prepare(self.root/"index",self.count,[self.prior],self.root/output,p=p,batch_size=3)
    def test_observation_level_p1_retains_both_photos_and_previous_status(self):
        doc=self.run_queue(1.0,"done")
        self.assertEqual((doc["full_frame_photo_ids"],doc["sampled_photo_ids_before_prior_dedup"],
                          doc["previously_measured_photo_ids_skipped"],doc["new_photo_attempts_planned"]),
                         (4,4,1,3))
        self.assertTrue(doc["no_per_taxon_or_observation_photo_count_cap"])
        with (self.root/"done"/"shards"/"shard-000.csv").open(newline="") as f:
            data=list(csv.DictReader(f))
        self.assertEqual({r["photo_id"] for r in data},{"101","102","104"})
        self.assertEqual(sum(r["observation_id"]=="OBS-A" for r in data),2)
        self.assertTrue(all(r["selection_probability"]=="1.0" for r in data))
        self.assertFalse(doc["colour_measured"])
    def test_nested_probability_milestones_do_not_drop_prior_included_observations(self):
        ids=["OBS-A","OBS-B","OBS-C","OBS-D","OBS-E"]
        low={x for x in ids if q.obs_draw(x)<0.1}
        high={x for x in ids if q.obs_draw(x)<0.5}
        self.assertTrue(low.issubset(high))
    def test_source_hash_drift_fails_closed(self):
        file=self.root/"index"/"photos"/"part-001.csv.gz"
        file.write_bytes(file.read_bytes()+b"altered")
        with self.assertRaisesRegex(ValueError,"SHA"):
            self.run_queue(1.0,"bad")
        self.assertFalse((self.root/"bad"/"manifest.json").exists())
    def test_incomplete_count_only_frame_not_accepted(self):
        doc=self.root/"index"/"manifest.json"
        v=json.loads(doc.read_text());v["photos_emitted"]=93209351
        doc.write_text(json.dumps(v))
        with self.assertRaisesRegex(ValueError,"complete"):
            self.run_queue(1.0,"bad")

if __name__=="__main__":unittest.main()
