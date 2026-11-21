"""No-network regressions for direct GBIF COL XR human-observation original-photo source."""
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from scripts.acquisition import mois_CCGBIF_colxr_field_photo_direct_search as m
from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import BASELINE,source_sha
from scripts.acquisition.mois_zero9_gbif_other_cc_media_union import SCHEMA as UNION


class HumanObservationCCSourceTests(unittest.TestCase):
    def row(self):
        return dict(inat_taxon_id="100",source_binomial="Plantus alba",
                    accepted_colxr_key="GBIF_100")

    def hit(self,basis="HUMAN_OBSERVATION",species="Plantus alba",license="https://creativecommons.org/licenses/by/4.0/",url="https://cdn.botanica.org/flower.jpg"):
        return {"key":201,"datasetKey":"ds1","basisOfRecord":basis,"species":species,
                "media":[{"type":"StillImage","identifier":url,"license":license}]}

    def test_field_still_image_only_with_source_species_and_cc_license(self):
        good,rejected=m.field_media(self.hit(),**{
            "inat_taxon_id":"100","binomial":"Plantus alba","accepted_key":"GBIF_100"})
        self.assertEqual(len(good),1)
        self.assertEqual(rejected,0)
        self.assertTrue(good[0]["GBIF_media_attached_to_correct_species"])
        self.assertEqual(good[0]["gbif_basisOfRecord"],"HUMAN_OBSERVATION")
        self.assertFalse(good[0]["floral_organ_visible_verified"])
        self.assertFalse(good[0]["numeric9_flower_ROI_measured"])
        self.assertTrue(good[0]["source_photo_URL_sha256"])

    def test_herbarium_conflict_and_inaturalist_urls_never_pass(self):
        for specimen in (
            self.hit(basis="PRESERVED_SPECIMEN"),
            self.hit(species="Plantus rubra"),
            self.hit(license=""),
            self.hit(url="https://inaturalist-open-data.s3.amazonaws.com/photos/123/large.jpg"),
        ):
            results,_=m.field_media(specimen,inat_taxon_id="100",
                                   binomial="Plantus alba",accepted_key="GBIF_100")
            self.assertEqual(results,[])
        nonjpeg=m.hit if False else self.hit()
        nonjpeg["media"][0]["type"]="MovingImage"
        results,_=m.field_media(nonjpeg,inat_taxon_id="100",
                               binomial="Plantus alba",accepted_key="GBIF_100")
        self.assertEqual(results,[])

    def test_signed_289_name_ledger_and_disjoint_strata(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            taxa=[]
            for i in range(304):
                taxa.append({"inat_taxon_id":str(i+1),
                             "source_binomial":"Plantus alba",
                             "name_status":"accepted_exact" if i<289 else "not_exact",
                             "accepted_colxr_key":str(i+1000) if i<289 else ""})
            df=pd.DataFrame(taxa)
            file=root/m.PALETTE_SOURCE
            df.to_csv(file,index=False)
            report={"schema":UNION,
                "source_original_42111_taxon_photo_atlas_sha256":BASELINE,
                "stratified_hash_draw_taxa":304,
                "official_COL_XR_exact_accepted_taxa_count":289,
                "photo_pixels_downloaded":0,
                "outputs_sha256":{file.name:source_sha(file)}}
            (root/"receipt.json").write_text(json.dumps(report)+"\n")
            groups=[m.signed_accepted(root,i) for i in range(4)]
            joined=pd.concat(groups)
            self.assertEqual(len(joined),289)
            self.assertTrue(joined.inat_taxon_id.is_unique)
            self.assertTrue(joined.accepted_colxr_key.is_unique)
            file.write_text("tampered")
            with self.assertRaisesRegex(ValueError,"SHA drift"):
                m.signed_accepted(root,0)


if __name__=="__main__":
    unittest.main()
