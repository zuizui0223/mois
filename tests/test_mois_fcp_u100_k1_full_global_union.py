"""Real photo-ID/terminal/SHA regression for signed 75 K1 + FCP U100 3141."""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_fcp_u100_k1_full_global_union as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class TestFullK1AndNeverMeasuredFCPGlobalPhotoUnion(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.output=self.root/"output"
        self.k1_global=self.root/"k1.json"
        self.fcp_receipt=self.root/"fcp.json"
        self.fcp_parquet=self.root/"real_fcp.parquet"
        self.preflight=self.root/"preflight.json"
        self.disp=self.root/"dispositions.csv"
        self.prefix=self.root/"prefix.csv"
        self.k1_roster=self.root/"k1_source.csv.gz"
        self.rows=[]
        for taxon,pid,action in [
            ("1","100","prior_original_photo_ROI_reuse"),
            ("2","200","reserved_by_inflight_K1_original_photo"),
            ("3","300","new_original_photo_ROI_queued"),
        ]:
            n=100+int(taxon)
            p=1/n/2
            self.rows.append(dict(
                photo_id=pid,observation_id="OBS-"+taxon,
                inat_taxon_id=taxon,
                full_93m_distinct_original_observations_of_taxon=str(n),
                photos_in_original_observation="2",
                depth1_joint_original_photo_inclusion_probability=format(p,".17g"),
            ))
        self.write_all()

    @staticmethod
    def write_csv(file,rows):
        with file.open("w",encoding="utf-8",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def write_all(self):
        self.write_csv(self.prefix,self.rows)
        self.disposition=[]
        for i,row in enumerate(self.rows):
            self.disposition.append(dict(
                photo_id=row["photo_id"],observation_id=row["observation_id"],
                inat_taxon_id=row["inat_taxon_id"],source_arm=m.PREFIX_ARM,
                reuse_or_new_photo_outcome=[
                    "prior_original_photo_ROI_reuse",
                    "reserved_by_inflight_K1_original_photo",
                    "new_original_photo_ROI_queued",
                ][i],
                depth1_joint_original_photo_inclusion_probability=
                    row["depth1_joint_original_photo_inclusion_probability"],
            ))
        self.write_csv(self.disp,self.disposition)
        ids=[dict(photo_id="200"),dict(photo_id="201")]
        with gzip.open(self.k1_roster,"wt",newline="",encoding="utf-8") as stream:
            w=csv.DictWriter(stream,fieldnames=["photo_id"],lineterminator="\n")
            w.writeheader();w.writerows(ids)
        self.pixel={
            "photo_id":"300", "observation_id":"OBS-3", "inat_taxon_id":"3",
            "morph":"white", "measurement_status":"classified_four_state_morph",
            "image_sha256":"a"*64,
            "selection_probability":self.rows[2]["depth1_joint_original_photo_inclusion_probability"],
            "source_arm":m.PREFIX_ARM,
        }
        pd.DataFrame([self.pixel]).to_parquet(self.fcp_parquet,index=False)
        self.k1_global.write_text(json.dumps({
            "schema":m.K1_GLOBAL_SCHEMA,
            "signed_historical_true_distinct_original_photo_IDs":1,
            "K1_original_new_photo_attempts":2,
            "genuine_K1_completed_and_SHA_sealed_groups":75,
            "source_K1_all_precolour_original_photo_IDs_disjoint_from_historical":2,
            "zero_all_source_K1_vs_historical_photo_ID_overlap":True,
            "counts_reproduce_independent_full_929138_ID_union":True,
            "exact_worldwide_historical_plus_completed_K1_distinct_original_photo_IDs":3,
            "exact_worldwide_historical_plus_completed_K1_four_hue_ROIs":1,
            "additional_raw_image_downloads_or_ROI_inference_during_this_update":0,
        }))
        self.fcp_receipt.write_text(json.dumps({
            "schema":m.FCP_ROI_SCHEMA,
            "original_93m_photo_ID_source_population":m.UPSTREAM,
            "eligible_global_first_FCP_nominal_taxon_IDs":3,
            "first_original_photo_candidates_retained":3,
            "previously_measured_original_image_ROIs_reused_elsewhere":1,
            "inflight_K1_original_photo_IDs_reserved_elsewhere":1,
            "newly_true_image_measured_original_photo_attempts":1,
            "individual_photo_design_probabilities_1_over_N_over_M_preserved":True,
            "photo_colour_inferred_by_unchanged_frozen_FCP_ROI_v4":True,
            "technical_unknowns_kept_unclassified":True,
            "individual_botanical_native_wild_anthesis_verified":False,
            "globally_unbiased_white_flora_fraction_estimated":False,
            "actual_new_image_parquet_SHA256":sha(self.fcp_parquet),
            "new_images_downloaded_OR_ROI_inferred_by_this_data_only_collector":0,
            "previous_real_photo_pixel_Parquets_rewritten":False,
            "newly_true_image_FCP_four_hue_classified":1,
            "newly_true_image_acquisition_failed":0,
            "newly_true_image_unclassifiable_acquired":0,
            "new_original_photo_hue_counts":{
                "white":1,"yellow_orange":0,"red_pink":0,
                "blue_purple":0,"mixed_uncertain":0
            }
        }))
        self.preflight.write_text(json.dumps({
            "schema":m.FCP_PREFLIGHT_SCHEMA,
            "eligible_nominal_FCP_U100_source_taxa":3,
            "already_measured_old_photo_IDs_to_reuse_NOT_RECOMPUTED":1,
            "already_running_K1_new_photo_IDs_to_reuse_NOT_RECOMPUTED":1,
            "never_before_measured_FCP_U100_original_photo_ROI_target_IDs":1,
            "new_original_photos_overlap_with_K1":0,
            "source_original_observations_one_per_taxon":3,
            "all_photo_inclusion_probabilities_individually_source_checked":True,
            "new_original_image_bytes_downloaded_here":0,
            "flower_colour_models_run_here":0,
            "ready_for_inference_ONLY_AFTER_K1_75_WORKERS_AND_75_SIGNED_GROUP_RELEASES":True,
        }))

    def run_union(self,**kw):
        args=dict(
            k1_global_sha=sha(self.k1_global),
            fcp_receipt_sha=sha(self.fcp_receipt),
            fcp_parquet_sha=sha(self.fcp_parquet),
            preflight_sha=sha(self.preflight),
            disposition_sha=sha(self.disp),
            prefix_sha=sha(self.prefix),
            k1_roster_sha=sha(self.k1_roster),
            old_n=1,k1_n=2,new_n=1,eligible=3,
            old_reuse=1,k1_reuse=1,k1_hues_expected=1,
            fcp_hues_expected=1,k1_groups_expected=75,
        )
        args.update(kw)
        return m.union(self.k1_global,self.fcp_receipt,self.fcp_parquet,
            self.preflight,self.disp,self.prefix,self.k1_roster,
            self.output,**args)

    def test_actual_distinct_photo_ID_total_and_hues(self):
        result=self.run_union()
        self.assertEqual(result["exact_global_unique_original_photo_ID_attempts"],4)
        self.assertEqual(result["exact_global_automatic_four_hue_photo_ROIs"],2)
        self.assertEqual(result["FCP_old_exact_original_photo_IDs_reused_without_new_inference"],1)
        self.assertEqual(result["FCP_K1_exact_photo_IDs_reused_without_new_inference"],1)
        self.assertEqual(result["new_FCP_images_overlapping_historical_or_K1"],0)
        self.assertFalse(result["globally_unbiased_white_flora_rate_not_estimated"] is False)
        self.assertEqual(result["previous_raw_image_downloads_or_image_inference_repeated_for_this_union"],0)
        self.assertTrue((self.output/"mois_K1_FCP_U100_source_audited_global_original_photo_union.json").is_file())

    def test_repeated_K1_photo_must_not_be_labelled_as_new_FCP(self):
        pixel=dict(self.pixel,photo_id="200")
        pd.DataFrame([pixel]).to_parquet(self.fcp_parquet,index=False)
        with self.assertRaisesRegex(ValueError,"photo IDs|new measured"):
            self.run_union()
        self.assertFalse(self.output.exists())

    def test_changed_status_cannot_turn_failed_photo_to_classified_white(self):
        pixel=dict(self.pixel,measurement_status="image_acquisition_failed")
        pd.DataFrame([pixel]).to_parquet(self.fcp_parquet,index=False)
        with self.assertRaisesRegex(ValueError,"terminal image SHA|technical nonresponse"):
            self.run_union()
        self.assertFalse(self.output.exists())

    def test_wrong_photo_probability_does_not_survive_signed_source(self):
        pixel=dict(self.pixel,selection_probability="0.4")
        pd.DataFrame([pixel]).to_parquet(self.fcp_parquet,index=False)
        with self.assertRaisesRegex(ValueError,"selection probability"):
            self.run_union()
        self.assertFalse(self.output.exists())

    def test_exact_FCP_source_disposition_required_not_candidate_only(self):
        disp=[dict(x) for x in self.disposition]
        disp[2]["reuse_or_new_photo_outcome"]="reserved_by_inflight_K1_original_photo"
        self.write_csv(self.disp,disp)
        with self.assertRaisesRegex(ValueError,"K1 already-selected|first-photo"):
            self.run_union()
        self.assertFalse(self.output.exists())

    def test_forged_source_SHA_must_fail_even_with_plausible_counts(self):
        self.fcp_parquet.write_bytes(self.fcp_parquet.read_bytes()+b"tampered")
        with self.assertRaisesRegex(ValueError,"Published original photo SHA"):
            self.run_union(fcp_parquet_sha="0"*64)
        self.assertFalse(self.output.exists())

    def test_old_K1_global_completion_must_be_75_proven_groups(self):
        doc=json.loads(self.k1_global.read_text())
        doc["genuine_K1_completed_and_SHA_sealed_groups"]=74
        self.k1_global.write_text(json.dumps(doc))
        with self.assertRaisesRegex(ValueError,"all 75 actual"):
            self.run_union()
        self.assertFalse(self.output.exists())

    def test_immutable_reconciliation_does_not_overwrite(self):
        self.run_union()
        with self.assertRaisesRegex(ValueError,"Refuse overwrite"):
            self.run_union()


if __name__=="__main__":
    unittest.main()
