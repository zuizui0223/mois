"""Original photo colour-blind 7200 historical date sampling and fake API.

No HTTP requests in pytest. Private photo/observer/site IDs never enter
aggregate public output; no missing date silently becomes non-flowering.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import zipfile
from pathlib import Path
from urllib.parse import urlparse
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_fcp_historical_date_overlay",
    ROOT/"scripts/analysis/mois_historical_photo_observation_date_reconciliation_20261016.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def sample_zip(tmp_path):
    zip_path=tmp_path/"fake_fcp_archive.zip"
    source=[];wide=[];regional=[]
    for j in range(12):
        photo=100+j;obs=200+j;tax=300+j//2
        source.append({
            "inat_taxon_id":tax,"observation_id":obs,"photo_id":photo,
            "latitude":40 if j<8 else -30,"longitude":-90 if j<8 else 120,
            "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
        })
        hue="white" if j%4==0 else "yellow_orange" if j%4<3 else "mixed_uncertain"
        label={
            "inat_taxon_id":tax,"observation_id":obs,"photo_id":photo,
            "measurement_status":"classified_four_state_morph" if hue!="mixed_uncertain"
                                 else "not_evaluable_roi_or_flip_gate",
            "morph":hue
        }
        (wide if j<7 else regional).append(label)
    def csv_gz(items):
        buf=io.StringIO()
        writer=csv.DictWriter(buf,fieldnames=list(items[0]))
        writer.writeheader();writer.writerows(items)
        return gzip.compress(buf.getvalue().encode())
    with zipfile.ZipFile(zip_path,"w") as z:
        z.writestr(m.ZIP_MEMBER,csv_gz(source))
        z.writestr(m.BREADTH,csv_gz(wide))
        z.writestr(m.REGIONAL,csv_gz(regional))
    return zip_path,source


def test_strict_outcome_blind_hash_sample_sample_size_and_source_label_provenance(tmp_path):
    zip_path,_=sample_zip(tmp_path)
    rows,n=m._historical_rows(zip_path,synthetic=True,max_photos=9)
    assert len(rows)==9 and n==12
    assert len({r["photo_id"] for r in rows})==9
    assert len({r["observation_id"] for r in rows})==9
    expected=sorted(range(100,112),
          key=lambda x:(hashlib.sha256((m.SEED+"/"+str(x)).encode()).hexdigest(),x))[:9]
    assert [r["photo_id"] for r in rows]==expected
    assert all(r["morph"] in m.PHOTO_FOUR or r["morph"]=="mixed_uncertain" for r in rows)
    with pytest.raises(ValueError,match="SHA256 mismatch"):
        m._historical_rows(zip_path,synthetic=False,max_photos=9)


def test_bounded_batched_date_current_photo_validation_and_image_nonresponse(tmp_path):
    z,source=sample_zip(tmp_path)
    source_by_observation={r["observation_id"]:r for r in source}
    calls=[];pauses=[]
    def getter(ids):
        assert 1<=len(ids)<=m.BATCH
        calls.append(ids)
        records={}
        for id_ in ids:
            donor=source_by_observation[id_]
            if id_==200:
                continue  # not returned; never invent observed_on
            row={
                "id":id_,"observed_on":"2020-04-12" if id_%2 else "2020-07-12",
                "photos":[{"id":donor["photo_id"]}],
                "taxon":{"id":donor["inat_taxon_id"]},
                "user":{"id":id_+1000},
                "annotations":[{"controlled_attribute_id":12,
                    "controlled_value_id":13,"vote_score":2}]
            }
            if id_==201:row["photos"]=[]
            if id_==202:row["observed_on"]=None
            records[id_]=row
        return records,hashlib.sha256(",".join(map(str,ids)).encode()).hexdigest()
    priv=tmp_path/"private_ids.csv";pub=tmp_path/"aggregate.json"
    r=m.run(z,priv,pub,synthetic=True,max_photos=9,
            batch_size=3,getter=getter,sleep=pauses.append)
    assert r["status"]==m.STATUS
    assert r["n_random_hash_sampled_original_observations"]==9
    assert r["sample_original_photographic_inclusion_probability"]==pytest.approx(9/12)
    assert r["n_serial_public_API_observation_batch_queries"]==3
    assert pauses==[m.GAP,m.GAP]
    assert len(r["source_original_API_response_sha256_by_batch"])==3
    assert sum(r["current_api_original_photo_and_date_status"].values())==9
    assert r["n_exact_original_photo_present_with_valid_current_API_date"]<=9
    assert r["original_four_hue_source_pixels_already_measured_not_newly_reprocessed"] is True
    assert r["original_API_observed_on_not_proven_frozen_at_archive_date"] is True
    assert r["source_not_a_probability_sample_of_all_world_flowering_plants"] is True
    assert all("public_suppressed_below" in x for x in
               r["source_observed_month_by_hemisphere_photo_white"])
    assert "observation_id" not in pub.read_text()
    assert "photo_id" not in pub.read_text()
    assert "longitude" not in pub.read_text()
    assert len(list(csv.DictReader(priv.open())))==9
    with pytest.raises(ValueError,match="Refuse overwrite"):
        m.run(z,priv,pub,synthetic=True,max_photos=9,batch_size=3,
              getter=getter,sleep=lambda _:None)


def test_api_batched_multiple_observations_extracts_only_requested_and_hashed_source():
    obs=[10,11,12]
    class FakeReply:
        status=200
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def read(self,n):
            return json.dumps({"total_results":2,"results":[
                {"id":10,"observed_on":"2020-05-13","photos":[{"id":30}]},
                {"id":11,"observed_on":"2020-05-15","photos":[{"id":31}]}
            ]}).encode()
    def opener(req,timeout):
        assert req.full_url.endswith("/10,11,12")
        return FakeReply()
    got,sha=m._response(obs,opener=opener)
    assert set(got)=={10,11} and len(sha)==64
    assert m._valid_date("2020-05-15")=="2020-05-15"
    assert m._valid_date("2027-01-01") is None
    assert m._valid_date("N/A") is None


def test_reject_unselected_returned_records_and_oversized_API_work():
    class BadReply:
        status=200
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def read(self,n):
            return json.dumps({"total_results":1,"results":[{"id":555}]}).encode()
    with pytest.raises(ValueError,match="Extraneous"):
        m._response([123],opener=lambda *a,**k:BadReply())
    with pytest.raises(ValueError,match="<=120"):
        m._response(list(range(1,130)),opener=lambda *a,**k:BadReply())


def test_multiple_species_and_local_season_mean_photo_white_and_phenology_negative():
    rows=[
        {"api_validation_status":"MATCHED_ORIGINAL_PHOTO_WITH_VALID_DATE",
         "geographic_domain":"NORTH_EXTRATROPICAL","observed_month":4,
         "measurement_status":"classified_four_state_morph",
         "morph":hue,"inat_taxon_id":tax,"api_observer_id":idx,
         "current_flowers_annotation":"NO_ANNOTATION"}
        for idx,(tax,hue) in enumerate([(1,"white"),(1,"white"),(1,"white"),
                 (2,"blue_purple"),(3,"red_pink"),(3,"white")])
    ]
    q=m._summarize(rows,min_public=1)
    assert len(q)==1
    assert q[0]["photo_white_fraction_given_classifiable"]==4/6
    assert q[0]["species_equal_white_fraction_given_classifiable"]==pytest.approx(
        (1+0+.5)/3)
    assert q[0]["positive_Flowers_annotation_source_obs"]==0
    assert m._local_season(-30,4)=="AUTUMN"
    assert m._local_season(30,4)=="SPRING"
    assert m._local_season(0,4)=="CALENDAR_Q2"
    assert m._annotation({"annotations":[]})=="NO_ANNOTATION"
