"""Data-only recovery of 170 cached and 342 previous-photo medium box statuses."""
import csv,hashlib,json,tempfile,unittest
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_512_prior_roi9_unknown_technical_screen_integrate import (
    combine,sha,N_CACHE,N_NEW
)
from scripts.analysis.mois_512_prior_unknown_technical_screen_plan import (
    SCHEMA as SOURCE_SCHEMA,ROI_REVIEW_SHA,BULK_PARQUET_SHA
)
from scripts.acquisition.mois_existing_roi_unknown_342_medium_flower_diagnostic import (
    SCHEMA as WORKER_SCHEMA,MODEL_SHA
)

class SignedOld512DiagnosticMergeTests(unittest.TestCase):
    def source(self,root):
        root=Path(root)
        src=root/"source";jobs=root/"workers"
        (src/"workers").mkdir(parents=True)
        jobs.mkdir()
        header=["photo_id","observation_id","inat_taxon_id","source_design",
                "latitude_zone","local_season","postoutcome_source_gap",
                "previous_original_photo_attempts","source_photo_opportunities",
                "reason"]
        def row(i):
            return dict(photo_id=str(100000+i),
                observation_id=f"obs-of-already-measured-{i}",
                inat_taxon_id=str(20000+i),
                source_design="bulk128_131072" if i<N_CACHE else "p2_first_4096",
                latitude_zone="S_tropical" if i%2 else "N_temperate",
                local_season="winter" if i%2 else "spring",
                postoutcome_source_gap=
                    "P2_source_has_photos_no_colour_10plus",
                previous_original_photo_attempts=7,
                source_photo_opportunities=500,
                reason="FOUR_HUE_UNKNOWN_AND_NO_STRICT9")
        originals=[row(i) for i in range(512)]
        def write(path,rows,columns):
            with path.open("w",newline="") as file:
                w=csv.DictWriter(file,fieldnames=columns)
                w.writeheader();w.writerows(rows)
        cache=[{**r,"screen_status":"candidate" if i%3 else "negative",
                "medium_actual_RGB_decoded":True,
                "floral_box_main_L_star":70.2 if i%3 else ""}
               for i,r in enumerate(originals[:N_CACHE])]
        q=originals[N_CACHE:]
        cache_file=src/"already_screened_500px_170_no_refetch.csv"
        q_file=src/"needs_bounded_medium500_screen_old_photoIDs.csv"
        write(cache_file,cache,header+["screen_status",
             "medium_actual_RGB_decoded","floral_box_main_L_star"])
        write(q_file,q,header)
        sizes=[43,43,43,43,43,43,42,42]
        shard_hashes={}
        for i,size in enumerate(sizes):
            batch=q[sum(sizes[:i]):sum(sizes[:i+1])]
            f=src/"workers"/f"shard_{i:02d}.csv"
            write(f,batch,header)
            shard_hashes[str(i)]=sha(f)
        m={"schema":SOURCE_SCHEMA,"source_review_photo_rows":512,
            "existing_500px_photo_screen_results_reused":170,
            "existing_unknown_photos_needing_small_bounded_500px_recheck":342,
            "new_photo_IDs_selected_for_world_census":0,
            "input_review_csv_SHA256":ROI_REVIEW_SHA,
            "existing_bulk128_RGB_pixel_parquet_SHA256":BULK_PARQUET_SHA,
            "cached_medium_result_SHA256":sha(cache_file),
            "new_diagnostic_photo_queue_SHA256":sha(q_file),
            "eight_worker_photo_counts":sizes,
            "worker_source_SHA256":shard_hashes}
        (src/"manifest.json").write_text(json.dumps(m))
        for i,size in enumerate(sizes):
            folder=jobs/f"mois-oldROI9-medium-diagnostic-{i}"
            folder.mkdir()
            batch=q[sum(sizes[:i]):sum(sizes[:i+1])]
            records=[{**row,"screen_status":"candidate" if j%2 else "negative",
                      "technical_diagnostic_class":
                        "detectable_large_flower_box_ROI9_failed" if j%2 else
                        "detector_negative_requires_botanical_review",
                      "500px_pixel_decoded":True,
                      "largest_box_area_fraction":.2 if j%2 else 0.,
                      "flower_box_count_002":1 if j%2 else 0,
                      "strict_original_roi9_remains_unknown":True
                      } for j,row in enumerate(batch)]
            photo_file=folder/f"shard_{i:02d}_original_medium_diagnostic.csv"
            write(photo_file,records,list(records[0]))
            receipt={
                "schema":WORKER_SCHEMA,"worker":i,
                "model_SHA256":MODEL_SHA,
                "original_source_manifest_SHA256":sha(src/"manifest.json"),
                "individual_original_source_worker_SHA256":shard_hashes[str(i)],
                "measurement_CSV_SHA256":sha(photo_file),
                "no_raw_image_bytes_archived":True,
                "prior_ROI9_unknowns_never_imputed_as_white_or_pigmented":True,
                "original_old_photo_ids_attempted_for_first_medium_diagnostic":size,
                "decoded_medium500":size}
            (folder/f"shard_{i:02d}_source_receipt.json").write_text(json.dumps(receipt))
        return src,jobs

    def test_every_original_photo_once_no_strict9_colour_imputation(self):
        with tempfile.TemporaryDirectory() as temp:
            src,workers=self.source(temp)
            output=Path(temp)/"integrated"
            r=combine(src,workers,output)
            self.assertEqual(r["old_original_photo_IDs_with_all_strictROI9_unknown"],512)
            self.assertEqual(r["distinct_source_taxa_reviewed_technically_not_botanically"],512)
            self.assertEqual(r["old_photo_ID_500px_detector_outputs_reused_without_fetch"],170)
            self.assertEqual(r["old_photo_IDs_newly_diagnosed_500px_but_already_original_fullROI_failed"],342)
            self.assertEqual(r["genuine_new_strict9_ROI_labels_generated"],0)
            self.assertEqual(sum(r["detector_statuses_512"].values()),512)
            combined=pd.read_csv(output/"all512_oldROI9_unknown_photoID_detector_status.csv")
            self.assertEqual(combined.photo_id.nunique(),512)
            self.assertEqual(combined.inat_taxon_id.nunique(),512)
            self.assertTrue(combined.strict_original9_pigment_fraction.isna().all())
            self.assertTrue(combined.old_ROI9_missing.all())

    def test_new_diagnostic_original_photo_ID_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            src,workers=self.source(temp)
            target=workers/"mois-oldROI9-medium-diagnostic-0"/"shard_00_original_medium_diagnostic.csv"
            rows=pd.read_csv(target)
            rows.loc[0,"photo_id"]=99999999
            rows.to_csv(target,index=False)
            with self.assertRaisesRegex(ValueError,"changed"):
                combine(src,workers,Path(temp)/"rejected")

if __name__=="__main__":
    unittest.main()
