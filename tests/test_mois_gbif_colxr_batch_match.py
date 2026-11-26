"""Offline tests for official GBIF v2 COL XR multi-name POST matching."""
import unittest
from unittest.mock import patch
import json

from scripts.analysis.mois_gbif_colxr_batch_match import (
    match_batch, parse_row, COL_XR
)


class Source:
    inat_taxon_id="123"
    species="Campanula punctata"
    dominant_hue9="purple"
    hue_status="provisional_single"


class FakeResponse:
    def __init__(self,value):
        self.value=json.dumps(value).encode()
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def read(self,size):return self.value[:size]


class GbifBatchTests(unittest.TestCase):
    def test_post_request_order_verified(self):
        names=["Campanula punctata","Boronia safrolifera"]
        data=[{"usage":{"key":str(i),"rank":"SPECIES",
                          "status":"ACCEPTED","canonicalName":name},
               "diagnostics":{"matchType":"EXACT"},
               "classification":[{"rank":"KINGDOM","name":"Plantae"}]}
              for i,name in enumerate(names)]
        with patch("scripts.analysis.mois_gbif_colxr_batch_match.urlopen",
                   return_value=FakeResponse(data)) as mock:
            result=match_batch(names)
            self.assertEqual(len(result),2)
            self.assertEqual(result[0]["usage"]["canonicalName"],names[0])
            req=mock.call_args.args[0]
            self.assertEqual(req.get_method(),"POST")
            self.assertIn(b"scientificName",req.data)

    def test_missing_plant_kingdom_not_silently_accepted(self):
        source=Source()
        payload={"usage":{"key":"abc","canonicalName":"Campanula punctata",
                          "rank":"SPECIES","status":"ACCEPTED"},
                 "diagnostics":{"matchType":"EXACT"},
                 "classification":[{"rank":"KINGDOM","name":"Animalia"}]}
        d=parse_row(source,payload)
        self.assertEqual(d["status"],"kingdom_context_review")
        payload["classification"][0]["name"]="Plantae"
        good=parse_row(source,payload)
        self.assertEqual(good["status"],"colxr_exact_accepted_species")
        self.assertEqual(good["COL_XR_checklistKey"],COL_XR)

    def test_hard_batch_limit(self):
        with self.assertRaisesRegex(ValueError,"exceeds"):
            match_batch(["Species name"]*1001)


if __name__=="__main__":
    unittest.main()
