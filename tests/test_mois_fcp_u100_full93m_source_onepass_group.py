"""FCP original U100 full 93m metadata design: 14 signed source tiers/8 parts.

Synthetic source fixtures only; they produce no new flower-colour results.
Tests exact 100 observation-distinct photos, full-taxonomy source-coverage,
conditional photo inclusion probability, unknown image formats, SHA/tar and
existing-source photo/observation signatures across all disjoint tiers.
"""
from __future__ import annotations

from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.acquisition import mois_fcp_u100_full93m_source_onepass_group as m
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, obs_draw, tier_of_score, TIER_ENDS,
)


class FullOriginal93mFCPHundredSourceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.source=self.root/"registry.json"
        self.receipt=self.root/"group-00-receipt.json"
        self.tiers=self.root/"tiers"
        self.tiers.mkdir()
        self.taxafile=self.root/"taxa.csv"
        self.taxafile.write_text("test")
        self.out=self.root/"u100"
        self.alltaxa={str(i):f"Plant-{i}" for i in range(1,42112)}
        # Taxon128: 105 original biological observations, one with four
        # different original photo IDs. Taxon256: 99, taxon384: one.
        self.A=[f"source-observation-A-{i:03d}" for i in range(105)]
        self.B=[f"source-observation-B-{i:03d}" for i in range(99)]
        self.C=["source-observation-C-only"]
        top=sorted(self.A,key=m.original_observation_priority)
        self.not_eligible_photo=top[-1]
        self.technical_unknown_photo=top[0]
        self.multi_photo=top[1]
        self.rows=[]
        pid=1000
        for taxon,observations in (("128",self.A),("256",self.B),("384",self.C)):
            for observation in observations:
                count=4 if observation==self.multi_photo else 1
                for pos in range(count):
                    row={name:"" for name in COLUMNS}
                    score=obs_draw(observation)
                    tier=tier_of_score(score)
                    ext="gif" if observation==self.technical_unknown_photo else "jpg"
                    row.update(
                        photo_id=str(pid),observation_uuid=observation,
                        inat_taxon_id=taxon,observer_id=f"observer-{observation}",
                        observed_on="2025-06-15" if taxon=="128" else "2025-12-15",
                        latitude="35.1",longitude="139.2",
                        positional_accuracy="10",photo_license="cc-by",
                        image_extension=ext,
                        observation_inclusion_score=str(score),
                        tier=str(tier))
                    self.rows.append((tier,int(taxon)%128,row))
                    pid+=1
        self.nphotos=len(self.rows)
        self.obs_n=105+99+1
        part_rows={(tier,part):[] for tier in range(m.TIERS) for part in range(8)}
        for tier,part,row in self.rows:
            part_rows[(tier,part)].append(row)
        fullarchive={}
        part_receipts=[]
        groupcounts={str(tier):0 for tier in range(m.TIERS)}
        for part in range(8):
            counts=Counter()
            hashes={}
            for tier in range(m.TIERS):
                records=part_rows[(tier,part)]
                raw=io.BytesIO()
                with gzip.GzipFile(fileobj=raw,mode="wb",mtime=0) as gz:
                    with io.TextIOWrapper(gz,encoding="utf-8",newline="") as f:
                        w=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
                        w.writeheader();w.writerows(records)
                body=raw.getvalue()
                # Archive constructed by tier after collecting part bytes.
                fullarchive[(tier,part)]=body
                hashes[str(tier)]=hashlib.sha256(body).hexdigest()
                counts[str(tier)]=len(records)
                groupcounts[str(tier)]+=len(records)
            part_receipts.append({
                "part":part,
                "original_photo_IDs_in_part":sum(counts.values()),
                "counts_by_disjoint_tier":dict(counts),
                "ranked_target_csv_sha256":hashes,
            })
        archivessha={}
        for tier in range(m.TIERS):
            tar=self.tiers/f"tier-{tier:02d}-group-00.tar"
            with tarfile.open(tar,"w") as stream:
                for part in range(8):
                    data=fullarchive[(tier,part)]
                    item=tarfile.TarInfo(f"tier-{tier:02d}/part-{part:03d}.csv.gz")
                    item.size=len(data)
                    stream.addfile(item,io.BytesIO(data))
            archivessha[tar.name]=m.sha256(tar)
        self.receipt.write_text(json.dumps({
            "schema":m.RESULT_SCHEMA,
            "source_snapshot_SHA256":"a"*64,
            "source_original_photo_IDs_uniqueness_proved_by_upstream_primary_key":True,
            "all_source_photo_IDs_retained_in_one_and_only_one_tier":True,
            "same_observation_uses_same_selection_score":True,
            "images_downloaded":0,"colour_pixels_measured":0,
            "parts":part_receipts,
            "counts_by_disjoint_tier":groupcounts,
            "original_photo_IDs_in_this_group":self.nphotos,
            "archive_sha256_by_tier":archivessha,
        }))
        self.source.write_text(json.dumps({
            "schema":"mois_complete_full93m_original_photo_target_registry_v1",
            "exact_original_photo_IDs_in_all_disjoint_tiers":self.nphotos,
            "original_named_taxa_from_full_source":42111,
            "full_source_manifest_sha256":"a"*64,
            "source_only_NO_photo_pixels_or_hue_measured":True,
            "native_vs_introduced_not_hard_filtered":True,
            "published_tier_archive_SHA256":archivessha,
        }))

    def run_source(self):
        with patch.object(m,"read_taxa",return_value=self.alltaxa):
            return m.execute_group(self.source,self.receipt,self.tiers,
                                   self.taxafile,self.out,0,
                                   expected_total=self.nphotos,check_taxa_sha=False)

    def test_true_full_source_FCP_U100_exact100_observations_and_one_photo_each(self):
        r=self.run_source()
        self.assertEqual(r["source_original_photo_IDs"],self.nphotos)
        self.assertEqual(r["source_distinct_observations"],self.obs_n)
        self.assertEqual(r["source_tier_archives_read_exactly_once"],14)
        self.assertEqual(r["FCP_U100_source_eligible_nominal_taxa"],1)
        self.assertEqual(r["FCP_U100_source_original_photo_opportunities"],100)
        self.assertEqual(r["FCP_U100_photo_unsupported_format_technical_unknown"],1)
        self.assertEqual(r["FCP_U100_original_photo_ROI_targets_loader_supported"],99)
        self.assertTrue(r["all_old_FCP_MOIS_original_photo_IDs_must_be_reused_before_any_new_ROI"])
        with (self.out/"exact_100_full93m_original_photo_observation_candidates.csv").open(newline="") as f:
            selected=list(csv.DictReader(f))
        self.assertEqual(len(selected),100)
        self.assertEqual({x["inat_taxon_id"] for x in selected},{"128"})
        self.assertEqual({x["observation_id"] for x in selected},
                         set(sorted(self.A,key=m.original_observation_priority)[:100]))
        self.assertEqual(len({x["photo_id"] for x in selected}),100)
        self.assertEqual(len({x["observation_id"] for x in selected}),100)
        self.assertNotIn(self.not_eligible_photo,{x["observation_id"] for x in selected})
        self.assertIn(self.technical_unknown_photo,{x["observation_id"] for x in selected})
        self.assertIn(self.multi_photo,{x["observation_id"] for x in selected})
        for x in selected:
            count=int(x["photos_in_original_observation"])
            self.assertEqual(int(x["full_93m_distinct_original_observations_of_taxon"]),105)
            self.assertAlmostEqual(
                float(x["conditional_full_93m_photo_selection_probability"]),100/105/count)
            if x["observation_id"]==self.multi_photo:
                self.assertEqual(count,4)
                self.assertAlmostEqual(float(x["conditional_one_photo_probability"]),0.25)
            if x["observation_id"]==self.technical_unknown_photo:
                self.assertEqual(x["photo_source_disposition"],"unsupported_original_extension")
        self.assertEqual(r["source_taxa_eligibility_breakdown"]["at_least_100_observations"],1)
        self.assertEqual(r["original_four_hue_image_ROIs_inferred"],0)
        self.assertEqual(r["no_original_image_downloads"],True)

    def test_source_roster_remembers_under100_and_zero_taxa_without_fabricating_100photos(self):
        self.run_source()
        with (self.out/"fcp_original100_complete_nominal_taxa_source_eligibility.csv").open(newline="") as f:
            rows={x["inat_taxon_id"]:x for x in csv.DictReader(f)}
        self.assertEqual(len(rows),len([x for x in self.alltaxa if int(x)%128//8==0]))
        self.assertEqual(int(rows["128"]["full_93m_original_observation_count"]),105)
        self.assertEqual(rows["128"]["FCP_U100_source_eligible"],"true")
        self.assertEqual(int(rows["256"]["full_93m_original_observation_count"]),99)
        self.assertEqual(rows["256"]["FCP_U100_source_eligible"],"false")
        self.assertEqual(int(rows["256"]["FCP_U100_missing_observations"]),1)
        self.assertEqual(int(rows["384"]["full_93m_original_observation_count"]),1)
        self.assertEqual(int(rows["384"]["FCP_U100_missing_observations"]),99)
        self.assertEqual(int(rows["512"]["full_93m_original_observation_count"]),0)
        self.assertEqual(int(rows["512"]["FCP_U100_missing_observations"]),100)

    def test_immutable_full_photo_source_cannot_be_rewritten(self):
        self.run_source()
        with self.assertRaisesRegex(ValueError,"cannot be overwritten"):
            self.run_source()

    def test_signed_full_source_byte_mutation_rejected_before_final_photo_hues(self):
        path=self.tiers/"tier-00-group-00.tar"
        with path.open("ab") as f:f.write(b"changed")
        with self.assertRaisesRegex(ValueError,"archive checksum changed"):
            self.run_source()
        self.assertFalse(self.out.exists())

    def test_group_without_consistent_taxon_part_ownership_rejected(self):
        bad=dict(self.rows[0][2])
        bad["inat_taxon_id"]="129"
        import sqlite3
        with sqlite3.connect(":memory:") as db:
            m.source.create_database(db)
            with self.assertRaisesRegex(ValueError,"taxon or tier"):
                m.attach_source_photo(db,bad,int(bad["tier"]),0)


if __name__=="__main__":
    unittest.main()
