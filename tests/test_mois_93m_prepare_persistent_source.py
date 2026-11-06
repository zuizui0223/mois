"""Zero-network original 93m source archive to persistent photo chunk preparation tests."""
import csv
import gzip
import io
import json
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_93m_prepare_persistent_source as p
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS

class SignedSourceGroupPrepTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.registry=self.root/"global.json"
        self.group=self.root/"group.json"
        self.tar=self.root/"tier-12-group-00.tar"
        self.out=self.root/"source_roots"
        self.records={}
        for part,n in ((0,3),(1,4)):
            rows=[]
            for i in range(n):
                row={col:"" for col in COLUMNS}
                row.update(photo_id=str(9000+part*100+i),
                           observation_uuid=f"o-{part}-{i}",
                           inat_taxon_id="11",image_extension="jpg",
                           observed_on="2025-02-10",photo_license="cc-by",
                           tier="12")
                rows.append(row)
            f=self.root/f"part-{part:03d}.csv.gz"
            with gzip.open(f,"wt",newline="",encoding="utf-8") as stream:
                w=csv.DictWriter(stream,fieldnames=COLUMNS,lineterminator="\n")
                w.writeheader();w.writerows(rows)
            self.records[part]=(f,n)
        with tarfile.open(self.tar,"w") as tar:
            for part,(f,n) in self.records.items():
                tar.add(f,arcname=f"tier-12/part-{part:03d}.csv.gz")
        self.registry.write_text(json.dumps({
            "schema":p.REGISTRY_SCHEMA,
            "exact_original_photo_IDs_in_all_disjoint_tiers":7,
            "native_vs_introduced_not_hard_filtered":True,
            "source_only_NO_photo_pixels_or_hue_measured":True,
            "full_source_manifest_sha256":"a"*64,
            "published_tier_archive_SHA256":{self.tar.name:p.sha(self.tar)},
        }))
        self.group.write_text(json.dumps({
            "schema":"mois_onepass_full93m_disjoint_photo_target_tiers_v1",
            "source_snapshot_SHA256":"a"*64,
            "counts_by_disjoint_tier":{"12":7},
            "parts":[
                {"part":part,"counts_by_disjoint_tier":{"12":n},
                 "ranked_target_csv_sha256":{"12":p.sha(f)}}
                for part,(f,n) in self.records.items()],
        }))

    def prepare(self):
        return p.prepare_group(self.registry,self.group,self.tar,self.out,
                               tier=12,group=0,chunk_size=2,
                               expected_world_photos=7,parts_per_group=2)

    def test_two_signed_original_parts_materialize_into_four_chunks_without_image_requests(self):
        result=self.prepare()
        self.assertEqual(result["actual_signed_source_photo_IDs_accounted"],7)
        self.assertEqual(result["original_source_image_downloads"],0)
        self.assertEqual(result["nonempty_photo_parts"],2)
        self.assertEqual(result["native_taxon_season_photo_exclusions"],0)
        dirs=self.out/"source_roots.txt"
        entries=[Path(line) for line in dirs.read_text().splitlines()]
        self.assertEqual(len(entries),2)
        self.assertTrue(all(x.is_dir() for x in entries))
        counts=[]
        for d in entries:
            receipt=json.loads((d/"chunk_receipt.json").read_text())
            counts.append(receipt["actual_original_photo_IDs"])
            self.assertLessEqual(receipt["chunk_size_maximum"],2)
        self.assertEqual(counts,[3,4])
        self.assertEqual(sum(counts),7)
        collected=p.collect_roots([self.out],self.root/"all_roots.txt",
                                  world_total=7,require_world_complete=False)
        self.assertEqual(collected["signed_photo_IDs_in_root_roster"],7)
        self.assertFalse(collected["all_world_original_photo_metadata_available"])

    def test_wrong_source_archive_sha_refuses_extraction(self):
        d=json.loads(self.registry.read_text())
        d["published_tier_archive_SHA256"][self.tar.name]="f"*64
        self.registry.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"archive SHA"):
            self.prepare()
        self.assertFalse(self.out.exists())

    def test_wrong_individual_gzip_sha_fails_without_relabelling_or_outputs(self):
        d=json.loads(self.group.read_text())
        d["parts"][0]["ranked_target_csv_sha256"]["12"]="0"*64
        self.group.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"metadata gzip SHA"):
            self.prepare()
        self.assertFalse(self.out.exists())

    def test_complete_source_roots_demand_all_224_original_archives(self):
        self.prepare()
        with self.assertRaisesRegex(ValueError,"all 224"):
            p.collect_roots([self.out],self.root/"bad-global.txt",
                            world_total=7,require_world_complete=True)
        self.assertFalse((self.root/"bad-global.txt").exists())

if __name__=="__main__":
    unittest.main()
