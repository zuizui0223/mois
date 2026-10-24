"""Photo-ID-immutable 93m balanced sample reuse and never-downloaded source bridge."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_93m_balanced_prior_ROI_reuse_bridge as bridge
from scripts.acquisition.mois_completed_pixel_checkpoint import sha
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS


class BalancingExistingOriginalPixelsTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/"fullsource.csv"
        self.manifest=self.root/"source.json"
        self.technical=self.root/"technical.csv"
        self.output=self.root/"new-photo-shards"
        self.priors=[]
        self.prior_counts={}
        specs=[
            ("main",[("1","obsA","white"),("2","obsM","mixed_uncertain")]),
            ("prior_group00",[("3","obsC","white")]),
            ("prior_R1",[("4","obsD","yellow_orange")]),
            ("prior_R2",[("5","obsE","mixed_uncertain")])
        ]
        for arm,items in specs:
            frame=pd.DataFrame([{
                "photo_id":pid, "observation_id":obs, "morph":hue,
                "measurement_status":"classified_four_state_morph" if hue!="mixed_uncertain"
                      else "not_evaluable_roi_or_flip_gate",
                "image_sha256":"d"*64
            } for pid,obs,hue in items])
            p=self.root/(arm+".parquet")
            rec=self.root/(arm+".json")
            frame.to_parquet(p,index=False)
            count=len(frame);hue=sum(z!="mixed_uncertain" for _,_,z in items)
            if arm=="main":
                rec.write_text(json.dumps({
                    "schema":"mois_sealed_actual_pixel_checkpoint_v1",
                    "combined_unique_photo_attempts":count,
                    "combined_four_hue_classified":hue,
                    "output_sha256":{p.name:sha(p)},
                }))
            else:
                rec.write_text(json.dumps({
                    "schema":"mois_global93m_observation_baseline_real_RGB_photo_results_v1",
                    "attempted":count,"classifiable":hue,
                    "measured_photo_parquet_sha256":sha(p),
                }))
            self.priors.append((arm,p,rec))
            self.prior_counts[arm]=(count,hue)
        self.technical.write_text(
            "photo_id,observation_id,acquired,status\n8,obsT,true,technical_cohort\n"
        )
        self.rows=[]
        fixtures=[
            ("1","obsA","supported",2,3),
            ("6","obsB","supported",4,2),
            ("7","obsM","supported",1,1),  # observation already measured, distinct photo -> STILL selected
            ("8","obsT","supported",1,1),  # technical benchmark -> no re-download
            ("9","obsI","unsupported_original_extension",1,1),
        ]
        for pid,obs,status,m,n in fixtures:
            self.rows.append({
                "photo_id":pid,"observation_id":obs,
                "inat_taxon_id":"1","species":"Campanula test",
                "observer_id":"photographer","observed_on":"2025-06-01",
                "latitude":"35","longitude":"140",
                "positional_accuracy":"5","photo_license":"cc-by",
                "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                "photo_source_disposition":status,
                "source_arm":bridge.SOURCE_ARM,
                "photos_in_original_observation":str(m),
                "observations_in_original_source_stratum":str(n),
                "joint_photo_selection_probability":format(.0064/m/n,".17g"),
            })
        with self.source.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=self.rows[0].keys(),lineterminator="\n")
            writer.writeheader();writer.writerows(self.rows)
        self.manifest.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "sampling_schema":bridge.SOURCE_SCHEMA,
            "source_arm":bridge.SOURCE_ARM,
            "all_upstream_original_photo_IDs":93_209_351,
            "original_tier06_source_photo_IDs_accounted":597070,
            "sampled_distinct_original_observations":5,
            "sampled_per_stratum_K":1,
            "source_sample_decided_without_photo_pixels":True,
            "selection_probability_is_per_photo_in_signed_shards":True,
            "original_source_sample_full_table_SHA256":sha(self.source),
        }))

    def bridge(self,**kwargs):
        return bridge.build(
            self.manifest,self.source,self.priors,self.technical,self.output,
            expected_sample=5,
            expected_sample_sha=sha(self.manifest),
            expected_selected_sha=sha(self.source),
            expected_prior_counts=self.prior_counts,
            expected_technical_sha=sha(self.technical),**kwargs)

    def test_every_selected_photo_disposition_conserved_and_old_pixels_not_refetched(self):
        r=self.bridge()
        self.assertEqual(r["all_balanced_selected_original_observations"],5)
        self.assertEqual(r["new_photo_attempts_planned"],2)
        self.assertEqual(r["reused_original_photo_ids"],1)
        self.assertEqual(r["technical_holdout_images_not_redownloaded"],1)
        self.assertEqual(r["unsupported_original_formats_technical_unknown"],1)
        self.assertEqual(r["previously_measured_original_observations_with_different_sample_photo"],1)
        self.assertEqual(r["work_units"],1)
        self.assertEqual(r["group_ids"],[0])
        self.assertFalse(r["colour_measured"])
        self.assertTrue(r["original_photo_source_probabilities_unchanged_by_prior_pixel_reuse"])
        self.assertEqual(len(r["prior_original_pixel_arms"]),4)
        self.assertEqual(r["prior_original_crossarm_photo_overlap"]["unique_old_photo_IDs"],5)
        self.assertEqual(r["prior_original_crossarm_photo_overlap"]["signed_arm_photo_rows"],5)
        with (self.output/"shards/shard-000.csv").open(newline="") as f:
            kept=list(csv.DictReader(f))
        self.assertEqual(tuple(kept[0].keys()),WORKER_COLUMNS)
        self.assertEqual({v["photo_id"] for v in kept},{"6","7"})
        self.assertEqual(set(v["observation_id"] for v in kept),{"obsB","obsM"})
        self.assertEqual(
            {v["photo_id"]:float(v["selection_probability"]) for v in kept},
            {"6":.0064/8,"7":.0064})
        self.assertEqual(sha(self.output/"shards/shard-000.csv"),
                         r["photo_shard_sha256"]["0"])
        with (self.output/"already_measured_original_photo_reuse.csv").open(newline="") as f:
            old=list(csv.DictReader(f))
        self.assertEqual(len(old),1)
        self.assertEqual(old[0]["photo_id"],"1")
        self.assertEqual(old[0]["prior_automatic_ROI_morph"],"white")

    def change_prior_to_shared_photo(self, *, new_obs="obsA", new_hue="white",
                                     image_sha="d"*64):
        path=self.priors[1][1]
        table=pd.read_parquet(path)
        table.loc[0,"photo_id"]="1"
        table.loc[0,"observation_id"]=new_obs
        table.loc[0,"morph"]=new_hue
        table.loc[0,"measurement_status"]="classified_four_state_morph"
        table.loc[0,"image_sha256"]=image_sha
        table.to_parquet(path,index=False)
        receipt=self.priors[1][2]
        d=json.loads(receipt.read_text())
        d["measured_photo_parquet_sha256"]=sha(path)
        receipt.write_text(json.dumps(d))

    def test_historical_exact_same_photo_and_observation_is_reused_once(self):
        self.change_prior_to_shared_photo()
        d=self.bridge()
        self.assertEqual(d["prior_original_crossarm_photo_overlap"][
            "old_cross_arm_exact_duplicate_photo_rows"],1)
        self.assertEqual(d["prior_original_crossarm_photo_overlap"][
            "old_cross_arm_same_photo_evidence_conflicts"] if
            "old_cross_arm_same_photo_evidence_conflicts" in
            d["prior_original_crossarm_photo_overlap"] else 0,0)
        self.assertEqual(d["reused_original_photo_ids"],1)
        self.assertEqual(d["new_photo_attempts_planned"],2)
        file=self.output/"prior_original_photo_crossarm_overlap.csv"
        self.assertEqual(d["prior_original_crossarm_overlap_rows_SHA256"],sha(file))
        with file.open(newline="") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(rows[0]["photo_id"],"1")
        self.assertEqual(rows[0]["first_source_arm"],"main")
        self.assertEqual(rows[0]["additional_source_arm"],"prior_group00")
        self.assertEqual(rows[0]["evidence_conflict"],"false")

    def test_same_original_photo_with_different_observation_fails_closed(self):
        self.change_prior_to_shared_photo(new_obs="obsC")
        with self.assertRaisesRegex(ValueError,
                                    "Conflicting original observation UUID"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def test_contradictory_prior_photo_colour_not_reused_as_white_flower(self):
        self.change_prior_to_shared_photo(new_hue="red_pink")
        data=self.bridge()
        self.assertEqual(data["prior_original_crossarm_photo_overlap"][
            "old_cross_arm_same_photo_evidence_conflicts"],1)
        self.assertEqual(data["selected_reused_photos_with_conflicting_prior_ROI_status"],1)
        conflict=self.output/"prior_original_photo_ROI_conflicts.csv"
        self.assertEqual(data["prior_original_crossarm_conflict_rows_SHA256"],sha(conflict))
        with (self.output/"already_measured_original_photo_reuse.csv").open(newline="") as f:
            reused=list(csv.DictReader(f))
        self.assertEqual(reused[0]["prior_automatic_ROI_morph"],"mixed_uncertain")
        self.assertEqual(reused[0]["prior_ROI_status"],"conflicting_prior_ROI_evidence")
        self.assertEqual(reused[0]["prior_photo_ROI_provenance_conflict"],"true")

    def set_legacy_blank_prior_ROI(self, *, status=""):
        # Historical projection has a terminal original photo ID but no
        # eligible ROI-v4 four-hue inference. It must remain NONCLASSIFIABLE.
        p=self.priors[0][1]
        old=pd.read_parquet(p)
        old.loc[1,"morph"]=""
        old.loc[1,"measurement_status"]=status
        old.loc[1,"image_sha256"]=""
        old.to_parquet(p,index=False)
        receipt=self.priors[0][2]
        prior=json.loads(receipt.read_text())
        prior["output_sha256"][p.name]=sha(p)
        receipt.write_text(json.dumps(prior))

    def test_legacy_original_checkpoint_blank_morph_is_unknown_not_white(self):
        self.set_legacy_blank_prior_ROI()
        result=self.bridge()
        self.assertEqual(result["prior_original_pixel_arms"][0]["four_hue_classified"],1)
        self.assertEqual(
            result["prior_original_pixel_arms"][0][
                "historical_legacy_blank_ROI_technical_unknown"],1)
        self.assertEqual(
            result["prior_original_crossarm_photo_overlap"][
                "old_source_technical_unknown_without_four_hue"],1)
        self.assertEqual(result["new_photo_attempts_planned"],2)
        self.assertFalse(result["colour_measured"])
        # A blank source projection cannot silently turn into an original
        # image classified as a white or coloured flowering petal.
        self.assertEqual(
            bridge.canonical_prior_original_ROI_hue("", "")[0],
            "mixed_uncertain")

    def test_classified_morph_may_not_be_blank_in_existing_pixel_receipts(self):
        self.set_legacy_blank_prior_ROI(status="classified_four_state_morph")
        with self.assertRaisesRegex(ValueError,
                                    "may not be interpreted as a classified petal"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def test_unrecognized_historical_photo_colour_fails_closed(self):
        p=self.priors[0][1]
        frame=pd.read_parquet(p)
        frame.loc[1,"morph"]="green"
        frame.to_parquet(p,index=False)
        prior=self.priors[0][2]
        report=json.loads(prior.read_text())
        report["output_sha256"][p.name]=sha(p)
        prior.write_text(json.dumps(report))
        with self.assertRaisesRegex(ValueError,
                                    "Unknown prior image colour"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def change_reused_original_photo_observation_to_legacy_missing(self, unknown="nan"):
        p=self.priors[0][1]
        frame=pd.read_parquet(p)
        frame.loc[0,"observation_id"]=unknown
        frame.to_parquet(p,index=False)
        receipt=self.priors[0][2]
        record=json.loads(receipt.read_text())
        record["output_sha256"][p.name]=sha(p)
        receipt.write_text(json.dumps(record))

    def test_signed_full93m_photo_UUID_backfills_unknown_old_roi_source_only(self):
        self.change_reused_original_photo_observation_to_legacy_missing("nan")
        data=self.bridge()
        self.assertEqual(data["reused_original_photo_ids"],1)
        self.assertEqual(data["reused_photo_legacy_observation_unresolved"],1)
        self.assertTrue(data["missing_legacy_observation_backfilled_only_from_signed_93m_source"])
        with (self.output/"already_measured_original_photo_reuse.csv").open(newline="") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]["photo_id"],"1")
        self.assertEqual(rows[0]["observation_id"],"obsA")
        self.assertEqual(rows[0]["prior_observation_id_missing_in_legacy_source"],"true")
        self.assertEqual(rows[0]["prior_automatic_ROI_morph"],"white")
        self.assertEqual(data["new_photo_attempts_planned"],2)
        self.assertFalse(data["colour_measured"])

    def test_blank_and_null_tokens_are_not_trusted_as_real_original_UUIDs(self):
        for missing in ("", "nan", "None", "NULL", "null", "na", "N/A"):
            self.assertTrue(bridge.prior_observation_unknown(missing))
        for real in ("obsA", "1234abcd", "actual-photo-observation"):
            self.assertFalse(bridge.prior_observation_unknown(real))

    def test_genuinely_conflicting_nonmissing_photo_observation_remains_refused(self):
        self.change_reused_original_photo_observation_to_legacy_missing("a-real-but-different-UUID")
        with self.assertRaisesRegex(ValueError,"conflicting original observation UUID"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def rewrite_balanced_original_source_photo_UUID(self, value):
        source=list(self.rows)
        source[0]=dict(source[0],observation_id=value)
        with self.source.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(source[0]),lineterminator="\n")
            writer.writeheader()
            writer.writerows(source)
        manifest=json.loads(self.manifest.read_text())
        manifest["original_source_sample_full_table_SHA256"]=sha(self.source)
        self.manifest.write_text(json.dumps(manifest))

    def test_exact_photo_ID_reused_across_numeric_iNat_and_signed_UUID_representations(self):
        uuid="15dbeb6f-493a-4a3a-91bf-a832eddf1001"
        self.rewrite_balanced_original_source_photo_UUID(uuid)
        self.change_reused_original_photo_observation_to_legacy_missing("1234567")
        report=self.bridge()
        self.assertEqual(report["reused_original_photo_ids"],1)
        self.assertEqual(report["reused_photo_prior_numeric_observation_UUID_alias"],1)
        self.assertEqual(report["reused_photo_legacy_observation_unresolved"],0)
        self.assertTrue(report[
            "same_photo_ID_prior_numeric_observation_alias_not_a_validated_global_UUID_crosswalk"
        ])
        self.assertEqual(report["new_photo_attempts_planned"],2)
        with (self.output/"already_measured_original_photo_reuse.csv").open(newline="") as f:
            old=list(csv.DictReader(f))
        self.assertEqual(len(old),1)
        self.assertEqual(old[0]["photo_id"],"1")
        self.assertEqual(old[0]["observation_id"],uuid)
        self.assertEqual(old[0]["prior_old_observation_ID_reference"],"1234567")
        self.assertEqual(old[0]["observation_ID_crosswalk_relation"],
                         "old_numeric_ID_vs_signed_UUID_same_photo_ID_unverified_crosswalk")
        self.assertEqual(old[0]["same_photo_ID_does_not_prove_independent_observation_crosswalk"],
                         "true")
        self.assertEqual(old[0]["prior_automatic_ROI_morph"],"white")
        self.assertFalse(report["colour_measured"])

    def test_numeric_ID_cannot_be_accepted_if_signed_source_uses_another_text_scheme(self):
        self.change_reused_original_photo_observation_to_legacy_missing("1234567")
        with self.assertRaisesRegex(ValueError,"conflicting original observation UUID"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def test_two_different_nonmissing_signed_UUIDs_are_never_equated(self):
        uuid="15dbeb6f-493a-4a3a-91bf-a832eddf1001"
        self.rewrite_balanced_original_source_photo_UUID(uuid)
        self.change_reused_original_photo_observation_to_legacy_missing(
            "15dbeb6f-493a-4a3a-91bf-a832eddf1002")
        with self.assertRaisesRegex(ValueError,"conflicting original observation UUID"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def test_unchanged_selected_photo_probability_is_mandatory(self):
        d=list(self.rows)
        d[1]=dict(d[1],joint_photo_selection_probability="0.5")
        with self.source.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=d[0],lineterminator="\n")
            writer.writeheader();writer.writerows(d)
        manifest=json.loads(self.manifest.read_text())
        manifest["original_source_sample_full_table_SHA256"]=sha(self.source)
        self.manifest.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,"probability changed"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def test_byte_tampering_of_original_93m_sample_is_refused(self):
        with self.source.open("ab") as f:
            f.write(b"fake-extra-source")
        with self.assertRaisesRegex(ValueError,"(manifest or selected roster SHA changed|Not the exact SHA-signed)"):
            self.bridge()

    def test_sample_is_immutable_after_success(self):
        self.bridge()
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.bridge()


if __name__=="__main__":
    unittest.main()
