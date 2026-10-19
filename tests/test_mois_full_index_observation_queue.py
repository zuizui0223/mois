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
    def test_unfamiliar_source_extensions_are_retained_in_precolour_photo_queue(self):
        # Original 93m has heterogeneous image_extension metadata. Previously
        # an unsupported selected photo aborted the entire source-frame draw.
        updates={"101":"WEBP","102":"GIF","104":"JPG"}
        for i in range(2):
            filename=self.root/"index"/"photos"/f"part-{i:03d}.csv.gz"
            with gzip.open(filename,"rt",newline="") as stream:
                reader=csv.DictReader(stream)
                header=reader.fieldnames
                rows=list(reader)
            for row in rows:
                row["image_extension"]=updates.get(row["photo_id"],row["image_extension"])
            with gzip.open(filename,"wt",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=header)
                writer.writeheader();writer.writerows(rows)
            manifest_file=self.root/"index"/"manifest.json"
            manifest=json.loads(manifest_file.read_text())
            manifest["shards"][i]["photos_sha256"]=q.sha(filename)
            manifest_file.write_text(json.dumps(manifest))
        receipt=self.run_queue(1.0,"nonjpeg")
        self.assertEqual(receipt["new_photo_attempts_planned"],3)
        self.assertTrue(receipt["all_selected_unsupported_extension_original_photo_IDs_retained"])
        self.assertEqual(receipt["selected_original_image_format_status"],
                         {"frozen_original_roi_supported":1,
                          "outside_frozen_original_roi_loader":2})
        with (self.root/"nonjpeg"/"shards"/"shard-000.csv").open(newline="") as stream:
            rows={r["photo_id"]:r for r in csv.DictReader(stream)}
        self.assertEqual(set(rows),{"101","102","104"})
        self.assertEqual(rows["101"]["source_image_extension_original"],"WEBP")
        self.assertEqual(rows["101"]["source_image_extension_status"],
                         "outside_frozen_original_roi_loader")
        self.assertTrue(rows["101"]["photo_url_large"].endswith("/large.webp"))
        self.assertTrue(rows["104"]["photo_url_large"].endswith("/large.jpg"))

    def test_malformed_extension_has_safe_deterministic_nonfetchable_url(self):
        ext, status=q.extension_for_frozen_original_loader(
            "jpg?redirect=https://evil.example/test")
        self.assertEqual(ext,"unsupported")
        self.assertEqual(status,"malformed_or_blank_original_extension")
        self.assertEqual(q.extension_for_frozen_original_loader("  .JPEG "),
                         ("jpeg","frozen_original_roi_supported"))
        self.assertEqual(q.extension_for_frozen_original_loader(""),
                         ("unsupported","malformed_or_blank_original_extension"))

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
