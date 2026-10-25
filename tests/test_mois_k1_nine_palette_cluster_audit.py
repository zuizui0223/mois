"""Original signed MOIS continuous nine-palette diagnostic/cluster safety tests."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd

from scripts.analysis import mois_k1_nine_palette_cluster_audit as m


class OriginalNinePaletteAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/"genuine-original.parquet"
        self.receipt=self.root/"receipt.json"
        self.output=self.root/"report"
        self.rows=[]
        states=[
            ("white","white"),("yellow_orange","yellow"),
            ("red_pink","pink"),("blue_purple","purple")
        ]
        for i in range(72):
            hue,major=states[i%4]
            palette={x:0.0 for x in m.PALETTE}
            palette[major]=0.68+(i%9)*0.019
            palette["bronze"]=0.04
            if major!="white":palette["white"]=0.18
            palette["blue" if major!="blue" else "orange"]=0.05
            self.rows.append(self.make_row(1000+i,hue,m.CLASSIFIED,palette))
        amb={x:0.0 for x in m.PALETTE}
        amb.update(white=0.42,purple=0.45,pink=0.13)
        self.rows.append(self.make_row(1090,"mixed_uncertain",m.AMBIGUOUS,amb))
        self.rows.append(self.make_row(1091,"mixed_uncertain",m.AMBIGUOUS,{}))
        roi={x:0.0 for x in m.PALETTE}
        roi["red"]=0.91
        self.rows.append(self.make_row(1092,"mixed_uncertain",
                                       "not_evaluable_roi_or_flip_gate",roi))
        self.rows.append(self.make_row(1093,"mixed_uncertain",
                                       "not_evaluable_roi_or_flip_gate",{}))
        self.rows.append(self.make_row(1094,"mixed_uncertain",m.ACQ_FAIL,{}))
        self.persist()

    @staticmethod
    def make_row(num,hue,status,palette):
        record=dict(photo_id=str(num),observation_id=f"uuid-{num}",
                    inat_taxon_id=str(num+9),morph=hue,
                    measurement_status=status,
                    image_sha256="" if status==m.ACQ_FAIL else "a"*64)
        for color in m.PALETTE:
            record[f"flower_fraction_{color}"]=palette.get(color,0.0)
        return record

    def persist(self):
        table=pd.DataFrame(self.rows)
        table.to_parquet(self.source,index=False)
        classified=sum(r["measurement_status"]==m.CLASSIFIED for r in self.rows)
        failed=sum(r["measurement_status"]==m.ACQ_FAIL for r in self.rows)
        self.original_four_count=classified
        self.receipt.write_text(json.dumps({
            "schema":"mois_global93m_observation_baseline_real_RGB_photo_results_v1",
            "measured_photo_parquet_sha256":m.sha(self.source),
            "attempted":len(self.rows),
            "classifiable":classified,
            "acquisition_failed":failed,
            "unclassifiable_acquired":len(self.rows)-classified-failed,
            "technical_unknowns_not_assumed_nonwhite":True,
            "raw_licensed_photo_pixels_redistributed":False,
        }))

    def run_audit(self,output=None):
        return m.run(self.source,self.receipt,output or self.output,
                     expected_source_sha=m.sha(self.source),
                     expected_receipt_sha=m.sha(self.receipt),
                     expected_attempts=len(self.rows),
                     expected_classified=self.original_four_count,
                     expected_failed=1,
                     expected_other=len(self.rows)-self.original_four_count-1,
                     k_values=(2,3,4),train_limit=68,score_limit=24)

    def test_real_original_four_hues_stay_unchanged_but_nine_composition_clusters(self):
        result=self.run_audit()
        self.assertEqual(result["source_four_hue_original_classifications_UNCHANGED"],72)
        self.assertEqual(result["source_image_acquired_but_unclassifiable_UNCHANGED"],4)
        self.assertEqual(result["source_image_acquisition_failures_UNCHANGED"],1)
        self.assertEqual(result["original_four_hue_classified_WITH_nine_palette_features"],72)
        self.assertEqual(result["source_ambiguous_palette_composition_with_usable_features"],1)
        self.assertEqual(result["other_nonclassifiable_with_9palette_NOT_treated_as_flowers"],1)
        self.assertEqual(result["source_technical_unknown_photos_promoted_to_four_hue"],0)
        self.assertIn(result["number_clusters_selected_exploratory"],(2,3,4))
        path=self.output/"original_K1_photo_9palette_and_exploratory_clusters.parquet"
        d=pd.read_parquet(path)
        self.assertEqual(set(d.morph),{"white","yellow_orange","red_pink",
                                       "blue_purple","mixed_uncertain"})
        self.assertEqual(d.source_four_hue_classified_unchanged.sum(),72)
        assigned=d.nine_palette_exploratory_cluster.astype(str).ne("")
        self.assertEqual(int(assigned.sum()),73)
        for photo in ("1091","1092","1093","1094"):
            self.assertEqual(d.loc[d.photo_id==photo,
                                   "nine_palette_exploratory_cluster"].iloc[0],"")
        row=d[d.photo_id=="1090"].iloc[0]
        self.assertTrue(row.nine_palette_exploratory_cluster.startswith("C"))
        self.assertEqual(row.cluster_assignment_provenance,
                         "exploratory_palette_ambiguous_unclassified")
        self.assertEqual(len(pd.read_csv(self.output/"cluster_k_stability_and_silhouette.csv")),3)
        self.assertEqual(len(pd.read_csv(self.output/"nine_palette_cluster_archetypes.csv")),
                         result["number_clusters_selected_exploratory"])

    def test_existing_known_four_hue_with_no_saved_palette_is_not_deleted(self):
        for color in m.PALETTE:
            self.rows[0][f"flower_fraction_{color}"]=None
        self.persist()
        result=self.run_audit()
        self.assertEqual(result["source_four_hue_original_classifications_UNCHANGED"],72)
        self.assertEqual(result["source_classified_missing_nine_palette_features"],1)
        self.assertEqual(result["original_four_hue_classified_WITH_nine_palette_features"],71)
        d=pd.read_parquet(self.output/"original_K1_photo_9palette_and_exploratory_clusters.parquet")
        r=d[d.photo_id=="1000"].iloc[0]
        self.assertEqual(r.morph,"white")
        self.assertEqual(r.nine_palette_exploratory_cluster,"")

    def test_source_photo_or_receipt_byte_tamper_rejected(self):
        original=m.sha(self.source)
        self.source.write_bytes(self.source.read_bytes()+b"invalid")
        with self.assertRaisesRegex(ValueError,"byte SHA changed"):
            m.run(self.source,self.receipt,self.output,
                  expected_source_sha=original,expected_receipt_sha=m.sha(self.receipt),
                  expected_attempts=77,expected_classified=72,expected_failed=1,
                  expected_other=4)
        self.assertFalse(self.output.exists())

    def test_no_fake_palette_on_failed_photo(self):
        self.rows[-1]["flower_fraction_white"]=1.0
        self.persist()
        with self.assertRaisesRegex(ValueError,"Acquisition failure cannot carry"):
            self.run_audit()
        self.assertFalse(self.output.exists())

    def test_mixed_roi_with_palette_cannot_be_promoted_to_colour_label(self):
        result=self.run_audit()
        self.assertEqual(result["source_technical_unknown_photos_promoted_to_four_hue"],0)
        status=pd.read_csv(self.output/"genuine_original_photo_missingness_by_status.csv")
        ambiguous=status[status.measurement_status==m.AMBIGUOUS].iloc[0]
        self.assertEqual(ambiguous.n_photo_attempts,2)
        self.assertEqual(ambiguous.exploratory_unclassified_palette_candidate,1)
        gate=status[status.measurement_status=="not_evaluable_roi_or_flip_gate"].iloc[0]
        self.assertEqual(gate.positive_palette_but_failed_other_quality_gate,1)

    def test_frozen_hue_and_status_must_be_mutually_consistent(self):
        self.rows[0]["morph"]="mixed_uncertain"
        self.persist()
        with self.assertRaisesRegex(ValueError,"status does not agree"):
            self.run_audit()

    def test_cluster_validation_uses_only_frozen_nine_palette_not_geo(self):
        frame=pd.DataFrame(self.rows)
        _,valid,palette=m.nine_palette_from_frame(frame)
        known=frame.measurement_status.eq(m.CLASSIFIED).to_numpy()
        p=palette[known & valid]
        ids=frame.photo_id.to_numpy(str)[known & valid]
        a=m.cluster_sqrt_palette(p,ids,k_values=(2,3),
                                 train_limit=64,score_limit=24)
        b=m.cluster_sqrt_palette(p,ids,k_values=(2,3),
                                 train_limit=64,score_limit=24)
        self.assertEqual(a["k"],b["k"])
        self.assertEqual(a["assign"](p),b["assign"](p))
        self.assertEqual(a["training_photo_hash"],b["training_photo_hash"])
        self.assertTrue(a["model_fit_only_source_classified"])

    def test_immutable_output(self):
        self.run_audit()
        with self.assertRaisesRegex(ValueError,"overwrite"):
            self.run_audit()


if __name__=="__main__":
    unittest.main()
