"""Small compressed provider fixture: test full-source join without touching AWS."""
import gzip
import io
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.acquisition import mois_full_historical_source_dates_official_inat_s3 as m


class FakeResp:
    status_code=200
    def __init__(self,compressed,etag):
        self.raw=io.BytesIO(compressed)
        self.headers={"Content-Length":str(len(compressed)),
                      "ETag":etag}
    def raise_for_status(self):pass
    def close(self):pass


def test_real_provider_stream_schema_and_compressed_content_sha(monkeypatch):
    binary=gzip.compress(
        b"photo_uuid\tphoto_id\tobservation_uuid\tobserver_id\n"
        b"uuidA\t101\tobsA\tu1\n"
        b"uuidB\t102\tobsB\tu2\n")
    entry={"url":"https://test.fake/photos.csv.gz","content_length_bytes":len(binary),"etag":'"x"'}
    reader=m.DigestReader(io.BytesIO(binary))
    monkeypatch.setattr(m,"provider_snapshot_stream",
                        lambda client,e:(FakeResp(binary,'"x"'),reader))
    rows=[]
    d=m.scan_table(None,entry,("photo_id","observation_uuid"),
                   lambda line,index:rows.append(
                       (line[index["photo_id"]],line[index["observation_uuid"]])),
                   "photos.csv.gz")
    assert rows==[(b"101",b"obsA"),(b"102",b"obsB")]
    assert d["compressed_bytes_read"]==len(binary)
    assert d["compressed_original_SHA256"]==hashlib.sha256(binary).hexdigest()
    assert d["source_verified_by_full_GZIP_CRC"]


def test_correct_uuid_photo_id_association_and_not_numeric_observation_id(monkeypatch):
    zphotos=gzip.compress(
        b"photo_uuid\tphoto_id\tobservation_uuid\tobserver_id\n"
        b"uu1\t101\tobsa-uuid\t17\n"
        b"uu2\t999999\tobswrong-uuid\t18\n"
        b"uu3\t103\tobsc-uuid\t20\n")
    zobs=gzip.compress(
        b"observation_uuid\tobserver_id\tlatitude\tlongitude\tpositional_accuracy\ttaxon_id\tquality_grade\tobserved_on\tanomaly_score\n"
        b"obsa-uuid\t17\t40\t140\t4\t300\tResearch Grade\t2025-08-14\t1\n"
        b"obswrong-uuid\t18\t22\t140\t0\t300\tResearch Grade\t2024-05-02\t1\n"
        b"obsc-uuid\t20\t22\t140\t0\t301\tResearch Grade\t2026-02-20\t1\n")
    source=pd.DataFrame({"photo_id":["101","103","1000"]})
    documents={"photos.csv.gz":zphotos,"observations.csv.gz":zobs}
    def load(client,entry):
        compressed=documents[entry["url"].rsplit("/",1)[-1]]
        return FakeResp(compressed,'"x"'),m.DigestReader(io.BytesIO(compressed))
    monkeypatch.setattr(m,"provider_snapshot_stream",load)
    probe={"source":{k:{"url":"https://example.com/"+k+".csv.gz",
                        "content_length_bytes":len(documents[k+".csv.gz"]),"etag":'"x"'}
                     for k in ("photos","observations")}}
    probe["source"]={"photo_table":probe["source"]["photos"],
                     "observation_table":probe["source"]["observations"]}
    got,ledger=m.scan_metadata(source,probe,None,m.date(2026,10,17))
    assert set(got)=={"101","103"}
    assert got["101"]["current_public_snapshot_observed_on"]=="2025-08-14"
    assert got["103"]["current_public_snapshot_taxon_id"]=="301"
    assert ledger["photo_IDs_found_in_open_snapshot"]==2
    assert ledger["matched_observation_uuid_found_in_metadata"]==2


def test_explicit_provider_date_validation():
    today=m.date(2026,10,17)
    assert m.valid_public_date(b"2026-01-18",today)=="2026-01-18"
    assert m.valid_public_date(b"2027-02-01",today) is None
    assert m.valid_public_date(b"2026-02-30",today) is None
    assert m.valid_public_date(b"not_a_date",today) is None
    assert m.valid_public_date(b"",today) is None


def test_full_source_date_rejoin_keeps_original_automatic_hue_and_unknown(monkeypatch,tmp_path):
    monkeypatch.setattr(m,"HISTORIC",200)
    monkeypatch.setattr(m,"MOIS_ROWS",201)
    monkeypatch.setattr(m,"ORIGINAL_DATED_HUES",1)
    source=pd.DataFrame([
       {"photo_id":str(i),"observation_id":str(3000+i),"inat_taxon_id":str(200+i),
        "morph":"white" if i==0 else ("mixed_uncertain" if i%2 else "blue_purple"),
        "latitude":40.1,"longitude":138.4,"observed_on":None,
        "old_mois_or_new_fcp":"historical"}
       for i in range(200)])
    other=pd.DataFrame([{"photo_id":"999","observation_id":"9999","inat_taxon_id":"999",
        "morph":"white","latitude":40.1,"longitude":138.4,
        "observed_on":"2025-05-17","old_mois_or_new_fcp":"new_fcp"}])
    combined=pd.concat([source,other],ignore_index=True)
    s3={
       str(i):{
          "current_public_snapshot_observation_uuid":f"uuid-{i}",
          "current_public_snapshot_observed_on":"2025-07-17" if i%2==0 else None,
          "current_public_snapshot_taxon_id":str(200+i),
          "current_public_quality_grade":"Research Grade"}
       for i in range(200)}
    old200=pd.DataFrame({"original_photo_id":[str(i) for i in range(200)],
          "retrieved_observed_on":["2025-07-17"]*200})
    output=tmp_path
    result=m.assemble(source,combined,s3,old200,output,{}, {})
    assert result["full_historical_attempted_source_photo_count"]==200
    assert result["photo_ids_with_verified_snapshot_observation_date"]==100
    assert result["original_before_classifiable_dated"]==1
    # half the old photos are dated, of which 99 are chromatic and 1 white;
    # missing/potentially ambiguous hue remains in the full denominator.
    assert result["new_after_classifiable_dated"]==101
    assert result["original_photo_current_open_metadata_identity_status"][
       "matched_original_photo_but_no_valid_observed_on"]==100
    assert result["original_200_date_disagreements_between_current_API_and_20260927_snapshot"]==0
    assert Path(output/m.FINAL_TABLE).exists()
    assert Path(output/m.SEASON_TABLE).exists()


def test_provider_unchanged_etag_is_hard_gate():
    class R:
        headers={"Content-Length":"1000","ETag":'"different"',
                 "Last-Modified":"Sun, 27 Sep 2026 16:45:05 GMT"}
        def raise_for_status(self):pass
    class Client:
        def head(self,*args,**kwargs):return R()
    with pytest.raises(ValueError,match="changed"):
        m.provider_snapshot_stream(Client(),{"url":"https://example.org/test.csv.gz",
            "content_length_bytes":1000,"etag":'"pinned"',
            "last_modified":"Sun, 27 Sep 2026 16:45:05 GMT"})
