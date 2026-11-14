"""Original 93m zero-colour taxon breadth-first source selector synthetic tests."""
import unittest
from scripts.acquisition.mois_species_first_zero9_source import (
    select_rows, output_rows
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import obs_draw

def candidate(pid,taxon,obs,*,license="cc-by",ext="jpg",supported="true"):
    return {"photo_id":str(pid),"inat_taxon_id":str(taxon),"observation_uuid":str(obs),
            "photo_license":license,"image_extension":ext,"frozen_loader_supported":supported,
            "observation_inclusion_score":format(obs_draw(obs),".17g"),
            "observer_id":"observer-1","observed_on":"2026-07-01",
            "latitude":"35.0","longitude":"139.0","tier":"13"}

class SpeciesZero9SourceTests(unittest.TestCase):
    def test_only_zero_taxa_independent_observations_and_no_previously_measured_photos(self):
        rows=[candidate(1,10,"a"), candidate(2,10,"a"),candidate(3,10,"b"),
              candidate(4,10,"c"),candidate(5,11,"d"),candidate(6,12,"f")]
        result,stats=select_rows(rows,{"10","12"},{"5"},{"f"},limit=4)
        out=output_rows(result)
        self.assertEqual(set(out.inat_taxon_id),{"10"})
        self.assertEqual(out.observation_id.nunique(),3)
        self.assertEqual(len(out),3)
        self.assertEqual(stats["candidate_photo_rows"],4)
        self.assertTrue(out.photo_url_large.str.contains("inaturalist-open-data.s3.amazonaws.com").all())

    def test_limit_keeps_low_observation_hash_instead_of_first_files(self):
        obs=["z","x","y","w","v","b"]
        rows=[candidate(100+i,19,o) for i,o in enumerate(obs)]
        result,_=select_rows(rows,{"19"},set(),set(),limit=2)
        out=output_rows(result)
        expected=sorted(obs,key=lambda x:(obs_draw(x),x))[:2]
        self.assertEqual(list(out.observation_id),expected)
        self.assertEqual(list(out.source_observation_rank),[1,2])

    def test_invalid_nine_palette_source_licence_skipped(self):
        rows=[candidate(1,18,"a",license="all-rights-reserved"),
              candidate(2,18,"b",ext="gif"),
              candidate(3,18,"c",supported="false"),
              candidate(4,18,"d")]
        result,stats=select_rows(rows,{"18"},set(),set())
        self.assertEqual(len(output_rows(result)),1)
        self.assertEqual(stats["candidate_photo_rows"],1)
        self.assertEqual(stats["missing_nonCC_photo_licence"],1)
        self.assertEqual(stats["unsupported_image_extension"],2)

if __name__=="__main__":
    unittest.main()
