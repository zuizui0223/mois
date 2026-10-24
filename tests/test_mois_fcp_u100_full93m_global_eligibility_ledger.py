"""Unit-test full 16-way FCP-U100 taxonomic source eligibility without pixels."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_u100_full93m_global_eligibility_ledger as m


class Exact42111SourceAggregationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.groups=[]
        self.output=self.root/"all42111"
        self.n_taxa=16
        self.n_photo=0
        self.eligible=0
        self.obs=0
        self.unknown=0
        for group in range(16):
            d=self.root/f"group{group:02d}"
            d.mkdir()
            taxon=str(128+group*8)
            n=105 if group==0 else 99 if group==1 else 0 if group==2 else 100
            photos=n+(2 if group==0 else 0)
            up=2 if group==0 else 0
            e=n>=100
            self.n_photo+=photos
            self.obs+=n
            self.eligible+=e
            self.unknown+=up
            rows=[{
                "inat_taxon_id":taxon,
                "species":f"Nominal plant {taxon}",
                "full_93m_original_observation_count":str(n),
                "full_93m_original_photo_count":str(photos),
                "FCP_U100_source_eligible":str(e).lower(),
                "FCP_U100_missing_observations":str(max(0,100-n)),
                "selected_observed_months":str(1 if e else 0),
                "selected_latitude_bands":str(1 if e else 0),
            }]
            file=d/"fcp_original100_complete_nominal_taxa_source_eligibility.csv"
            with file.open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=m.HEADER,lineterminator="\n")
                w.writeheader();w.writerows(rows)
            sig={
                "group":group,"original_photos":photos,
                "original_tier_counts":{str(i):photos if i==0 else 0 for i in range(14)},
                "original_group_receipt_SHA256":"a"*64,
                "original_full_93m_source_registry_SHA256":"b"*64,
            }
            r={
                "schema":m.GROUP_SCHEMA,
                "full_original_photo_source_IDs":0, # set after all group counts known
                "full_original_named_plant_taxa":16,
                "original_taxon_modulo_source_group":group,
                "source_nominal_taxa":1,
                "source_distinct_observations":n,
                "source_original_photo_IDs":photos,
                "source_tier_archives_read_exactly_once":14,
                "no_original_image_downloads":True,
                "no_original_ROI_or_hue_inferences":True,
                "original_source_eligibility_and_all_taxa_preserved":True,
                "FCP_original_photo_observation_target_per_taxon":100,
                "all_old_FCP_MOIS_original_photo_IDs_must_be_reused_before_any_new_ROI":True,
                "FCP_original_iNaturalist_H9_phenology_annotation_verified_for_all_new_MOIS_photos":False,
                "full_WCVP_taxonomy_and_native_wild_flowering_photo_status_verified":False,
                "source_group_original_photo_ID_SHA":sig,
                "taxon_source_eligibility_SHA256":m.sha256(file),
                "U100_complete_original_photo_candidate_SHA256":"c"*64,
                "unsupported_original_photo_unknown_SHA256":"d"*64,
                "precolour_ROI_worker_shard_SHA256":{},
                "FCP_U100_source_eligible_nominal_taxa":int(e),
                "FCP_U100_source_original_photo_opportunities":100*int(e),
                "FCP_U100_original_photo_ROI_targets_loader_supported":100*int(e)-up,
                "FCP_U100_photo_unsupported_format_technical_unknown":up,
            }
            (d/"fcp_U100_full93m_group_source_receipt.json").write_text(json.dumps(r))
            self.groups.append(d)
        for d in self.groups:
            f=d/"fcp_U100_full93m_group_source_receipt.json"
            r=json.loads(f.read_text())
            r["full_original_photo_source_IDs"]=self.n_photo
            f.write_text(json.dumps(r))

    def call(self):
        return m.aggregate(self.groups,self.output,
                           expected_groups=16,expected_taxa=self.n_taxa,
                           expected_photos=self.n_photo)

    def test_exact_global_fcp_100_photo_taxon_roster_and_probability_source_counts(self):
        d=self.call()
        self.assertEqual(d["source_original_photo_IDs_reconciled"],self.n_photo)
        self.assertEqual(d["source_original_distinct_observations"],self.obs)
        self.assertEqual(d["source_FCP_style_exact100_photo_observation_taxa"],14)
        self.assertEqual(d["source_FCP_style_exact100_observation_photo_targets"],1400)
        self.assertEqual(d["source_FCP_style_supported_original_photo_extension_targets"],1398)
        self.assertEqual(d["source_FCP_style_unsupported_photo_unknowns"],2)
        self.assertEqual(d["source_taxa_with_less_than100_different_original_observations"],1)
        self.assertEqual(d["source_taxa_without_original_observations"],1)
        self.assertEqual(d["photo_ROI_pixels_opened_or_reinferred"],0)
        self.assertTrue(d["source_100taxon_scheme_is_not_worldwide_unbiased_flora_colour_rate"])
        p=self.output/"full42111_nominal_taxa_original_U100_source_eligibility.csv"
        self.assertEqual(d["complete_42111_taxa_eligibility_csv_sha256"],m.sha256(p))
        with p.open(newline="") as stream:
            r=list(csv.DictReader(stream))
        self.assertEqual(len(r),16)
        self.assertEqual([int(x["inat_taxon_id"]) for x in r],
                         sorted(int(x["inat_taxon_id"]) for x in r))

    def test_one_group_missing_cannot_represent_global_source(self):
        with self.assertRaisesRegex(ValueError,"group is absent"):
            m.aggregate(self.groups[:-1],self.output,
                        expected_groups=16,expected_taxa=16,expected_photos=self.n_photo)
        self.assertFalse(self.output.exists())

    def test_tampered_original_taxon_roster_rejected(self):
        p=self.groups[4]/"fcp_original100_complete_nominal_taxa_source_eligibility.csv"
        p.write_bytes(p.read_bytes()+b"tamper")
        with self.assertRaisesRegex(ValueError,"unsigned"):
            self.call()
        self.assertFalse(self.output.exists())

    def test_wrong_group_taxon_cross_partition_rejected(self):
        p=self.groups[5]/"fcp_original100_complete_nominal_taxa_source_eligibility.csv"
        with p.open(newline="") as f:
            r=list(csv.DictReader(f))
        r[0]["inat_taxon_id"]="128"
        with p.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=m.HEADER,lineterminator="\n")
            writer.writeheader();writer.writerows(r)
        receipt=self.groups[5]/"fcp_U100_full93m_group_source_receipt.json"
        d=json.loads(receipt.read_text())
        d["taxon_source_eligibility_SHA256"]=m.sha256(p)
        receipt.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"wrong-modulo"):
            self.call()
        self.assertFalse(self.output.exists())

    def test_unsupported_source_not_allowed_to_masquerade_as_classified(self):
        p=self.groups[0]/"fcp_U100_full93m_group_source_receipt.json"
        d=json.loads(p.read_text())
        d["FCP_U100_original_photo_ROI_targets_loader_supported"]+=1
        p.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"count mismatch"):
            self.call()

    def test_republishing_same_photo_source_is_forbidden(self):
        self.call()
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.call()


if __name__=="__main__":
    unittest.main()
