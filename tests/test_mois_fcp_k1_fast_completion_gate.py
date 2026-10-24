"""K1 scheduled group archiver: archived-only fast metadata and full fallback."""
from __future__ import annotations
import hashlib
import json
import unittest
from scripts.acquisition import mois_fcp_k1_fast_completion_gate as m


def archive(n, groups=None):
    groups=list(range(n)) if groups is None else sorted(groups)
    per=[{"original_worker_group":g,
          "original_group_receipt_SHA256":format(g+1,"064x"),
          "new_original_photo_attempts":min(1024,m.SOURCE_TARGET-1024*g),
          "four_hue_classified":400} for g in groups]
    digest=hashlib.sha256("\n".join(
        f'{r["original_worker_group"]}:{r["original_group_receipt_SHA256"]}'
        for r in per).encode()).hexdigest()
    tag=f"mois-K1-genuine-completed-original-ROI-{n}groups-{digest[:16]}-v1"
    summary={
        "schema":m.CUM_SCHEMA,
        "source_run":m.SOURCE_RUN,
        "total_scheduled_original_worker_groups":m.MAX_GROUPS,
        "new_original_photo_target_total":m.SOURCE_TARGET,
        "previously_measured_K1_photo_IDs_reused_elsewhere_NOT_NEW":m.REUSE,
        "source_precolour_probability_photo_sampling_preserved":True,
        "within_K1_original_source_group_units_disjoint":True,
        "new_original_images_downloaded_or_colour_inferred_in_THIS_ledger":0,
        "this_is_NOT_a_crossarm_photo_ID_deduplicated_global_total":True,
        "global_wild_flora_or_blooming_petal_status_not_independently_verified":True,
        "completed_and_SHA_sealed_original_worker_groups":n,
        "completed_original_worker_group_IDs":groups,
        "per_group":per,
        "source_original_group_receipt_aggregate_SHA256":digest,
    }
    return tag,summary


def inventory(success=(0,1,2,3)):
    jobs=[]
    art=[]
    for g in range(75):
        status="success" if g in success else None
        jobs.append({"name":f"warm-original-ROI ({g})","conclusion":status})
        if status=="success":
            art.append({"name":f"mois-balanced-onephoto-K1-real-RGB-group-{g}",
                        "digest":"sha256:"+"a"*64,"size_in_bytes":1234,
                        "expired":False})
    return jobs,art


class FastOriginalPhotoGroupArchiverTests(unittest.TestCase):
    def test_one_latest_signed_cumulative_subset_skips_all_old_worker_lookups(self):
        old,snapshot=archive(3)
        tags=["unrelated-release",old,"mois-K1-genuine-completed-original-ROI-1groups-"+("0"*16)+"-v1"]
        self.assertEqual(m.newest_snapshot_tag(tags),old)
        self.assertEqual(m.verified_archived_groups(snapshot,old),{0,1,2})
        jobs,art=inventory()
        actual=m.successful_finished_groups(jobs,art)
        looked=[]
        fresh=m.groups_needing_archive(actual,{0,1,2},lambda x: looked.append(x) or False)
        self.assertEqual(fresh,[3])
        self.assertEqual(looked,[3])

    def test_already_signed_missing_from_cumulative_group_is_not_republished(self):
        jobs,art=inventory()
        done=m.successful_finished_groups(jobs,art)
        looked=[]
        pending=m.groups_needing_archive(done,{0,1},lambda g: looked.append(g) or (g==2))
        self.assertEqual(pending,[3])
        self.assertEqual(looked,[2,3])

    def test_no_snapshot_fallback_checks_complete_tag_inventory(self):
        j,a=inventory()
        done=m.successful_finished_groups(j,a)
        complete={m.SEALED_TAG.format(g) for g in (0,2,3)}
        self.assertEqual(m.missing_from_all_releases(done,complete),[1])
        self.assertIsNone(m.newest_snapshot_tag(["new-nonK1-release"]))

    def test_changed_signed_cumulative_sha_does_not_skip_real_pixels(self):
        tag,snapshot=archive(2)
        snapshot["per_group"][0]["original_group_receipt_SHA256"]="9"*64
        with self.assertRaisesRegex(ValueError,"SHA/tag mismatch"):
            m.verified_archived_groups(snapshot,tag)

    def test_cumulative_old_pixel_source_and_group_are_strict(self):
        tag,snapshot=archive(2)
        snapshot["source_run"]=99
        with self.assertRaisesRegex(ValueError,"source receipts"):
            m.verified_archived_groups(snapshot,tag)
        tag,snapshot=archive(2)
        snapshot["completed_original_worker_group_IDs"]=[0,0]
        with self.assertRaisesRegex(ValueError,"group source identity"):
            m.verified_archived_groups(snapshot,tag)
        tag,snapshot=archive(2)
        snapshot["per_group"][0]["new_original_photo_attempts"]=9
        with self.assertRaisesRegex(ValueError,"worker terminal"):
            m.verified_archived_groups(snapshot,tag)

    def test_missing_original_zip_or_wrong_SHA_rejected(self):
        j,a=inventory()
        a.pop()
        with self.assertRaisesRegex(ValueError,"lack their saved ZIPs"):
            m.successful_finished_groups(j,a)
        j,a=inventory()
        a[0]["digest"]="no-real-original-image-proof"
        with self.assertRaisesRegex(ValueError,"authenticated SHA"):
            m.successful_finished_groups(j,a)
        j,a=inventory()
        a[0]["expired"]=True
        with self.assertRaisesRegex(ValueError,"expired"):
            m.successful_finished_groups(j,a)

    def test_all75_groups_and_last_partial977_images_valid(self):
        tag,snapshot=archive(75)
        groups=m.verified_archived_groups(snapshot,tag)
        self.assertEqual(len(groups),75)
        self.assertEqual(snapshot["per_group"][-1]["new_original_photo_attempts"],977)
        j,a=inventory(tuple(range(75)))
        self.assertEqual(m.successful_finished_groups(j,a),groups)
        self.assertEqual(m.groups_needing_archive(groups,groups,lambda _: self.fail("Already archived")),[])

    def test_no_false_archive_if_artifact_job_not_successful(self):
        j,a=inventory()
        a[0]["name"]="mois-balanced-onephoto-K1-real-RGB-group-70"
        with self.assertRaisesRegex(ValueError,"without matching successful"):
            m.successful_finished_groups(j,a)
        with self.assertRaisesRegex(ValueError,"not among successful"):
            m.groups_needing_archive({1,2},{0},lambda _: False)

    def test_conflicting_same_count_immutable_snapshot_refused(self):
        t,_=archive(2)
        wrong=t.replace(t.split("-")[-2],"a"*16)
        with self.assertRaisesRegex(ValueError,"Conflicting equally sized"):
            m.newest_snapshot_tag([t,wrong])


if __name__=="__main__":
    unittest.main()
