"""Synthetic SHA and URL/taxon collision regression for the GBIF CC-media rescue."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from scripts.acquisition.mois_zero9_gbif_other_cc_media_union import (
    union, CANDIDATE, STATUS
)
from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import (
    BASELINE, SCHEMA, source_sha
)


class GBIFCC304TaxaMergeTests(unittest.TestCase):
    def create_sample(self,root):
        roots=[]
        for i,(offset,n) in enumerate(zip((0,48,112,176,240),(48,64,64,64,64))):
            folder=root/f"sample{i}";folder.mkdir()
            t=pd.DataFrame({
                "inat_taxon_id":[str(offset+x+1) for x in range(n)],
                "name_status":["accepted_exact"]*n,
                "source_binomial":[f"Plantus epithet{offset+x+1}" for x in range(n)]
            })
            rows=[]
            for j in range(n):
                url="https://www.example.org/images/photo"+str(offset+j)+".jpg"
                if i==1 and j==0:
                    url="https://www.example.org/images/photo0.jpg"
                rows.append({
                    "inat_taxon_id":str(offset+j+1),
                    "gbif_occurrence_key":str(1000+offset+j),
                    "image_url":url,
                    "image_license":"https://creativecommons.org/licenses/by/4.0/"
                })
            photo=folder/CANDIDATE;status=folder/STATUS
            pd.DataFrame(rows).to_csv(photo,index=False,compression="gzip")
            t.to_csv(status,index=False)
            report={
                "schema":SCHEMA,
                "source_signed_original_atlas_sha256":BASELINE,
                "original_photo_pixels_downloaded":0,
                "nine_fraction_flower_colour_measurements":0,
                "selected_zero9_binomial_taxa_sample":n,
                "sample_source_fixed_hash_offset":offset,
                "noninat_CC_candidate_media_metadata":len(rows),
                "outputs_sha256":{photo.name:source_sha(photo),status.name:source_sha(status)}
            }
            (folder/"receipt.json").write_text(json.dumps(report)+"\n")
            roots.append(folder)
        return roots

    def test_url_dedup_quarantines_both_source_taxa(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);src=self.create_sample(root)
            out=root/"union"
            info=union(src,out)
            self.assertEqual(info["stratified_hash_draw_taxa"],304)
            self.assertEqual(info["raw_candidate_url_rows"],304)
            self.assertEqual(info["eligible_unique_external_URL_metadata_count"],302)
            self.assertEqual(info["duplicate_URL_rows_held_for_review"],2)
            self.assertEqual(info["image_hosts_seen"],1)
            self.assertEqual(info["photo_pixels_downloaded"],0)
            self.assertEqual(info["new_strict9_numerical_photo_measurements"],0)
            self.assertEqual(len(pd.read_csv(out/"mois_CC_noninat_GBIF_public_image_host_inventory.csv")),1)

    def test_wrong_signed_taxa_or_media_sha_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);roots=self.create_sample(root)
            photo=roots[0]/CANDIDATE
            with open(photo,"ab") as f:f.write(b"tampered")
            with self.assertRaisesRegex(ValueError,"byte-SHA"):
                union(roots,root/"untrusted")


if __name__=="__main__":
    unittest.main()
