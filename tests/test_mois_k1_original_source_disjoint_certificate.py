"""Signed historic 929138 versus future K1 all-source original photo ID exclusion."""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_k1_original_source_disjoint_certificate as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class TestEntireK1PrecolourOriginalIDExclusion(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.source=self.root/"source"
        (self.source/"shards").mkdir(parents=True)
        self.manifest=self.source/"manifest.json"
        self.shard=self.source/"shards/shard-000.csv"
        self.out=self.root/"output"
        self.old={"10":("main","OLD_OBS","white","classified_four_state_morph","a"*64)}
        self.old_obs={"OLD_OBS":("main","10")}
        self.rows=[
            dict(photo_id="20",observation_id="NEW_OBS_1",inat_taxon_id="1",
                 measurement_id="NEWPHOTO_1",selection_probability="0.1",
                 source_arm=m.SOURCE_ARM,shard="0"),
            dict(photo_id="21",observation_id="OLD_OBS",inat_taxon_id="2",
                 measurement_id="NEWPHOTO_2",selection_probability="0.004",
                 source_arm=m.SOURCE_ARM,shard="0")
        ]
        self.write_source()

    def write_source(self):
        with self.shard.open("w",newline="",encoding="utf-8") as stream:
            w=csv.DictWriter(stream,fieldnames=list(self.rows[0]),lineterminator="\n")
            w.writeheader();w.writerows(self.rows)
        doc=dict(schema="mois_full_source_observation_sampling_cohort_v1",
            full_frame_photo_ids=93209351,
            source_arm=m.SOURCE_ARM,
            source_precolour_balanced_manifest_SHA256=m.PARENT_SOURCE_SHA,
            new_photo_attempts_planned=len(self.rows),
            all_balanced_selected_original_observations=m.SELECTED_SOURCE_OBSERVATIONS,
            reused_original_photo_ids=m.REUSED_OLD_ORIGINALS,
            work_units=1,worker_ids=[0],group_ids=[0],
            units_per_group=4,original_model_batch_size=256,
            photo_shard_sha256={"0":sha(self.shard)},
            known_previous_full_photo_and_observation_IDs_source_SHA_verified=True,
            previously_counted_photo_ID_duplicates_are_NOT_new_inferences=True,
            original_photo_source_probabilities_unchanged_by_prior_pixel_reuse=True,
            technical_unknown_not_other_colour=True,
            images_downloaded_by_this_source_bridge=0,
            new_colour_ROIs_inferred_by_this_source_bridge=0)
        self.manifest.write_text(json.dumps(doc),encoding="utf-8")

    def certify(self,where=None,**kwargs):
        with unittest.mock.patch.object(m,"GROUPS",1),unittest.mock.patch.object(m,"UNITS",1):
            return m.certify_source(self.source,self.old,self.old_obs,[],
                where or self.out,expected_old=1,expected_hues=1,
                expected_new=2,expected_sha=sha(self.manifest),**kwargs)

    def test_sha_pin_every_original_photo_AND_historical_ID_noncollision(self):
        result=self.certify()
        self.assertEqual(result["K1_all_new_precolour_source_original_photo_IDs"],2)
        self.assertEqual(result["K1_exact_historical_photo_ID_collisions"],0)
        self.assertEqual(result["K1_historical_original_observation_ID_overlaps_are_not_independence"],1)
        self.assertEqual(result["historical_verified_four_hue_original_photo_ROIs"],1)
        self.assertTrue(result["original_photo_joint_selection_probability_retained"])
        self.assertFalse(result["global_white_flora_fraction_unweighted_estimated"])
        gz=self.out/result["K1_new_source_photo_roster"]
        self.assertEqual(sha(gz),result["K1_new_source_photo_roster_sha256"])
        with gzip.open(gz,"rt",newline="") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual([x["photo_id"] for x in rows],["20","21"])
        self.assertEqual([x["selection_probability"] for x in rows],["0.1","0.004"])
        self.assertEqual([x["original_worker_group"] for x in rows],["0","0"])

    def test_real_duplicate_historical_photo_never_reclassified(self):
        self.rows[0]["photo_id"]="10"
        self.write_source()
        with self.assertRaisesRegex(ValueError,"photo-ID/observation collision"):
            self.certify()
        self.assertFalse(self.out.exists())

    def test_original_precolour_source_shard_byte_change_rejected(self):
        self.shard.write_text(self.shard.read_text()+"\nchanged")
        with self.assertRaisesRegex(ValueError,"shard SHA mismatch"):
            self.certify()
        self.assertFalse(self.out.exists())

    def test_nonpositive_probability_never_measured_or_rewritten(self):
        self.rows[1]["selection_probability"]="0"
        self.write_source()
        with self.assertRaisesRegex(ValueError,"inclusion probability"):
            self.certify()
        self.assertFalse(self.out.exists())

    def test_duplicate_photo_or_observation_source_rejected(self):
        self.rows[1]["observation_id"]=self.rows[0]["observation_id"]
        self.write_source()
        with self.assertRaisesRegex(ValueError,"photo-ID/observation collision"):
            self.certify()
        self.assertFalse(self.out.exists())

    def test_immutable_gzip_is_byte_identical_for_same_signed_source(self):
        self.certify()
        second=self.certify(self.root/"another")
        first=json.loads((self.out/"K1_original_source_disjointness_certificate.json").read_text())
        self.assertEqual(first["K1_new_source_photo_roster_sha256"],
                         second["K1_new_source_photo_roster_sha256"])
        with self.assertRaisesRegex(ValueError,"cannot overwrite"):
            self.certify()


if __name__=="__main__":
    unittest.main()
