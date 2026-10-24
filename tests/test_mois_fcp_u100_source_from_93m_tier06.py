"""FCP 100-distinct-observation scientific frame, not a fake 93m colour result."""
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.acquisition import mois_fcp_u100_source_from_93m_tier06 as fcp
from scripts.acquisition import mois_93m_taxon_month_lat_balanced_photo_sampler as source
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS


def original_row(photo, uuid, taxon, *, observed="2025-06-15", extension="jpg"):
    row={key:"" for key in COLUMNS}
    row.update(
        photo_id=str(photo), observation_uuid=str(uuid), inat_taxon_id=str(taxon),
        observer_id="independent-" + str(uuid),
        photo_license="cc-by", image_extension=extension,
        observed_on=observed,latitude="35.123",longitude="140.008",
        positional_accuracy="20",tier="6",observation_inclusion_score="0.01",
    )
    return row


class FCPHundredOriginalObservationTests(unittest.TestCase):
    def setUp(self):
        self.db=sqlite3.connect(":memory:")
        self.addCleanup(self.db.close)
        source.create_database(self.db)
        # Eligible A: 107 different observation UUIDs; one observation with
        # three extra ORIGINAL PHOTOS. 107 observations, 110 photo IDs.
        self.photos=[]
        with patch.object(source, "obs_draw", return_value=0.01):
            for i in range(107):
                uuid=f"original-observation-A-{i:03d}"
                row=original_row(10000+i,uuid,1)
                source.import_photo(self.db,row)
                self.photos.append(row)
                if i==0:
                    for j in range(3):
                        clone=original_row(60000+j,uuid,1)
                        source.import_photo(self.db,clone)
                        self.photos.append(clone)
            # B has 99 independent observations: no fake FCP U100 support.
            for i in range(99):
                row=original_row(20000+i,f"original-observation-B-{i:03d}",2,
                                 observed="2026-03-16")
                source.import_photo(self.db,row)
                self.photos.append(row)
            # C has one observed/tropical taxon, record remains missing support.
            source.import_photo(self.db,original_row(30001,"C-original",3,extension="gif"))
        self.names={"1":"flower A","2":"flower B","3":"flower C","4":"no source observations"}
        self.db.commit()

    def test_fcp_u100_100_photos_from_100_different_observations_and_no_outcome_filter(self):
        rows,inventory,eligible,months,regions=fcp.sample_fcp_u100(self.db,self.names)
        self.assertEqual(eligible,{"1":107})
        self.assertEqual(inventory["1"],(107,110))
        self.assertEqual(inventory["2"],(99,99))
        self.assertEqual(inventory["3"],(1,1))
        self.assertEqual(len(rows),100)
        self.assertEqual({r["inat_taxon_id"] for r in rows},{"1"})
        self.assertEqual(len({r["observation_id"] for r in rows}),100)
        self.assertEqual(len({r["photo_id"] for r in rows}),100)
        self.assertEqual(sum(months["1"].values()),100)
        self.assertEqual(sum(regions["1"].values()),100)
        for row in rows:
            obs_n=row["original_observations_in_taxon_tier06"]
            multiplicity=row["photos_attached_to_original_observation"]
            self.assertEqual(obs_n,107)
            self.assertAlmostEqual(
                float(row["conditional_photo_inclusion_probability_in_tier06"]),
                100/obs_n/multiplicity)
            self.assertAlmostEqual(
                float(row["nominal_source_tier_plus_photo_probability"]),
                0.0064*100/obs_n/multiplicity)
            self.assertTrue(row["photo_url_large"].endswith(".jpg"))
        self.assertFalse(any(x["inat_taxon_id"] in {"2","3","4"} for x in rows))

    def test_all_42111_source_taxa_can_keep_zero_and_under100_status(self):
        sampled,inventory,eligible,_,_=fcp.sample_fcp_u100(self.db,self.names)
        for taxon in self.names:
            n=inventory.get(taxon,(0,0))[0]
            self.assertEqual(taxon in eligible,n>=100)
        self.assertEqual({t:inventory.get(t,(0,0))[0] for t in self.names},
                         {"1":107,"2":99,"3":1,"4":0})
        self.assertEqual(fcp.FCP_MIN_CLASSIFIABLE,40)
        self.assertEqual(fcp.FCP_SOURCE_TARGET,100)
        self.assertEqual(fcp.SCHEMA,"mois_fcp_style_exact100_original_observations_per_taxon_tier06_source_v1")

    def test_legal_one_original_photo_for_each_observation_no_photographic_replication(self):
        rows,_,_,_,_=fcp.sample_fcp_u100(self.db,self.names)
        source_ids={r["photo_id"] for r in self.photos}
        self.assertTrue(set(r["photo_id"] for r in rows)<=source_ids)
        self.assertFalse(any(r["photo_source_disposition"]!="supported" for r in rows))
        obs_with_three_extra="original-observation-A-000"
        if obs_with_three_extra in {r["observation_id"] for r in rows}:
            choice=[r for r in rows if r["observation_id"]==obs_with_three_extra]
            self.assertEqual(len(choice),1)
            self.assertEqual(choice[0]["photos_attached_to_original_observation"],4)
            self.assertAlmostEqual(
                float(choice[0]["conditional_one_photo_per_observation_probability"]),.25)

    def test_inconsistent_source_taxon_frame_or_change_to_100_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"fixed at exactly 100"):
            fcp.sample_fcp_u100(self.db,self.names,target=40)
        with self.assertRaisesRegex(ValueError,"absent from 42,111-taxon source"):
            fcp.sample_fcp_u100(self.db,{"1":"A"})
        self.db.execute("UPDATE obs SET taxon='999' WHERE uuid='C-original'")
        with self.assertRaisesRegex(ValueError,"absent from 42,111-taxon source"):
            fcp.sample_fcp_u100(self.db,self.names)


if __name__ == "__main__":
    unittest.main()
