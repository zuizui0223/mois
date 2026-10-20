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

    def test_full_frame_selected_unsupported_formats_are_source_only_not_inferred_hues(self):
        # Simulate exact upstream source records of GIF, empty extension,
        # and case/dot variations of a genuinely supported JPEG.
        photos=self.root/"index"/"photos"/"part-001.csv.gz"
        with gzip.open(photos,"rt",newline="") as stream:
            reader=csv.DictReader(stream)
            header=reader.fieldnames
            contents=list(reader)
        for pid,extension in (("105","gif"),("106",".JPEG "),("107","")):
            contents.append(dict(zip(header,(pid,"OBS-A","1","cc-by",extension,"user1"))))
        with gzip.open(photos,"wt",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=header)
            writer.writeheader();writer.writerows(contents)
        manifest=self.root/"index"/"manifest.json"
        m=json.loads(manifest.read_text())
        m["photos_emitted"]=7
        m["shards"][1]["photo_count"]+=3
        m["shards"][1]["photos_sha256"]=q.sha(photos)
        manifest.write_text(json.dumps(m))
        with patch.object(q,"N_PHOTOS",7),patch.object(q,"N_TAXA",2),patch.object(q,"N_SHARDS",2):
            result=q.prepare(self.root/"index",self.count,[self.prior],
                             self.root/"unusual",p=1.0,batch_size=3)
        self.assertEqual(result["sampled_photo_ids_before_prior_dedup"],7)
        self.assertEqual(result["previously_measured_photo_ids_skipped"],1)
        self.assertEqual(result["new_photo_attempts_planned"],4)
        self.assertEqual(result["selected_new_photos_unsupported_by_frozen_loader"],2)
        self.assertEqual(result["unsupported_original_source_photo_ext_counts"],
                         {"":1,"gif":1})
        self.assertTrue(result["no_per_taxon_or_observation_photo_count_cap"])
        unsupported=self.root/"unusual"/"unsupported_selected_original_photo_formats.jsonl"
        self.assertEqual(q.sha(unsupported),
                         result["selected_new_photo_formats_unsupported_sha256"])
        rows=[json.loads(s) for s in unsupported.read_text().splitlines()]
        self.assertEqual({r["photo_id"] for r in rows},{"105","107"})
        self.assertTrue(all("morph" not in r and "white" not in r for r in rows))
        selected=[]
        for path in sorted((self.root/"unusual"/"shards").glob("*.csv")):
            with path.open(newline="") as stream:
                selected.extend(csv.DictReader(stream))
        self.assertEqual({r["photo_id"] for r in selected},{"101","102","104","106"})
        jpeg=next(r for r in selected if r["photo_id"]=="106")
        self.assertTrue(jpeg["photo_url_large"].endswith("/106/large.jpeg"))
        self.assertTrue(all("white" not in r and "morph" not in r for r in selected))

    def test_only_original_frozen_loader_extensions_are_normalized(self):
        self.assertEqual(q.frozen_original_extension(" .JPG "),"jpg")
        self.assertEqual(q.frozen_original_extension(".JPEG"),"jpeg")
        self.assertEqual(q.frozen_original_extension("PNG"),"png")
        self.assertIsNone(q.frozen_original_extension("gif"))
        self.assertIsNone(q.frozen_original_extension("webp"))
        self.assertIsNone(q.frozen_original_extension(""))
        self.assertIsNone(q.frozen_original_extension("../large.jpg?format=PNG"))


    def test_second_p0004_source_draw_adds_only_new_observations_not_first_pixels(self):
        # Frozen observation draw is nested: A in both p0002/p0004;
        # B enters only p0004; C in neither. Previous original photo IDs
        # must not be acquired again even when their observation is selected.
        chance={"OBS-A":0.0001,"OBS-B":0.0003,"OBS-C":0.9}
        with patch.object(q,"N_PHOTOS",4),patch.object(q,"N_TAXA",2), \
             patch.object(q,"N_SHARDS",2), \
             patch.object(q,"obs_draw",side_effect=lambda uuid:chance[uuid]):
            first=q.prepare(self.root/"index",self.count,
                [self.prior],self.root/"first_p0002",p=0.0002,batch_size=3)
            self.assertEqual(first["new_photo_attempts_planned"],2)
            first_roster=[]
            for source_path in (self.root/"first_p0002"/"shards").glob("*.csv"):
                with source_path.open(newline="") as stream:
                    first_roster.extend(csv.DictReader(stream))
            self.assertEqual({r["photo_id"] for r in first_roster},{"101","102"})
            # Existing original sampled observations are already genuine
            # measurements; include one additional historical photo 104.
            previous=self.root/"all_finished_prior.csv"
            previous.write_text("photo_id\n101\n102\n104\n")
            second=q.prepare(self.root/"index",self.count,
                [previous],self.root/"second_p0004",p=0.0004,batch_size=3)
        self.assertEqual(second["sampled_photo_ids_before_prior_dedup"],3)
        self.assertEqual(second["previously_measured_photo_ids_skipped"],2)
        self.assertEqual(second["new_photo_attempts_planned"],1)
        self.assertEqual(second["selected_new_photos_unsupported_by_frozen_loader"],0)
        with (self.root/"second_p0004"/"shards"/"shard-000.csv").open(newline="") as stream:
            selected=list(csv.DictReader(stream))
        self.assertEqual([r["photo_id"] for r in selected],["103"])
        self.assertEqual(selected[0]["observation_id"],"OBS-B")
        self.assertEqual(selected[0]["selection_probability"],"0.0004")
        self.assertFalse(second["colour_measured"])
        self.assertFalse(second["ecological_analysis_run"])
        self.assertTrue(second["all_unselected_source_photos_eligible_at_higher_probability"])


    def test_third_p0008_nested_observation_queue_excludes_all_prior_original_pixels(self):
        # A in p0004 and p0008, B newly admitted ONLY at p0008, C
        # excluded from both. Already processed A photos must never be
        # fetched again, even though A passes the larger hash threshold.
        draw={"OBS-A":0.0003,"OBS-B":0.0006,"OBS-C":0.9}
        previous=self.root/"complete_prior_RBG.csv"
        previous.write_text("photo_id\n101\n102\n104\n")
        with patch.object(q,"N_PHOTOS",4),patch.object(q,"N_TAXA",2), \
             patch.object(q,"N_SHARDS",2), \
             patch.object(q,"obs_draw",side_effect=lambda uuid:draw[uuid]):
            source=q.prepare(self.root/"index",self.count,
                [previous],self.root/"third_p0008",p=0.0008,batch_size=3)
        self.assertEqual(source["full_frame_photo_ids"],4)
        self.assertEqual(source["selection_probability"],0.0008)
        self.assertEqual(source["sampled_photo_ids_before_prior_dedup"],3)
        self.assertEqual(source["previously_measured_photo_ids_skipped"],2)
        self.assertEqual(source["new_photo_attempts_planned"],1)
        self.assertTrue(source["no_per_taxon_or_observation_photo_count_cap"])
        self.assertFalse(source["colour_measured"])
        with (self.root/"third_p0008"/"shards"/"shard-000.csv").open(newline="") as f:
            selected=list(csv.DictReader(f))
        self.assertEqual([z["photo_id"] for z in selected],["103"])
        self.assertEqual(selected[0]["observation_id"],"OBS-B")
        self.assertEqual(selected[0]["selection_probability"],"0.0008")

if __name__=="__main__":unittest.main()
