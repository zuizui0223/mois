"""Synthetic and source-SHA controls for MOIS 400 blind botanical consensus.

Synthetic flowers are never real observations or white-colour measurements.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts.analysis.mois_prepare_400_botanical_review_handoff import (
    FIELDS, token, STAGE, TISSUE,
)
from scripts.analysis.mois_reconcile_400_independent_botanical_reviews import (
    ADJ_FIELDS, SUMMARY_FIELDS, RESULT_FIELDS, reconcile, sha,
)


def write_rows(path,fieldnames,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fieldnames,lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def rows(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fixture(tmp_path):
    source=tmp_path/"source400.csv"
    doc=tmp_path/"source_receipt.json"
    handoff=tmp_path/"review_handoff.json"
    a=tmp_path/"review_A.csv"
    b=tmp_path/"review_B.csv"
    originals=[]
    for i,season in enumerate(("summer","autumn","spring","winter"),1):
        tax=101 if i<=2 else 202
        pid=str(900000+i)
        originals.append({
            "inat_taxon_id":str(tax),
            "species":("Taxon alba" if i<=2 else "Taxon blue"),
            "photo_id":pid,
            "observation_id":"o"+str(i),
            "observer_id":"u"+str(i),
            "observed_on":("2025-07-15" if season=="summer" else
                           "2025-10-01" if season=="autumn" else
                           "2025-03-15" if season=="spring" else
                           "2025-01-11"),
            "photo_license":"cc-by",
            "attribution":"photo-source-attributed-original-"+str(i),
            "photo_url_large":
               "https://inaturalist-open-data.s3.amazonaws.com/photos/"+pid+"/large.jpg",
            "site_latcell":"510","site_loncell":"300",
            "season_year":"2025","selected_season":season,
        })
    write_rows(source,list(originals[0]),originals)
    doc.write_text(json.dumps({
        "status":"ORIGINAL_PHOTO_SOURCE_ONLY_NOT_A_BLOOM_OR_COLOUR",
        "selected_photo_manifest_sha256":sha(source),
        "selected_original_photo_ids":4,
        "selected_distinct_source_observers":4,
    }))
    handoff.write_text(json.dumps({
        "status":"400_ORIGINAL_PHOTO_REVIEW_FORMS_BOTH_UNFILLED_NOT_COLOUR",
        "source_photo_manifest_SHA256":sha(source),
        "reviewer_forms_season_date_climate_colour_blind":True,
        "n_actual_independent_botanical_reviews_completed":0,
        "new_photo_pixels_measured":0,
        "n_original_photos":4,
        "n_original_named_taxa":2,
    }))
    unfilled=[]
    for r in originals:
        unfilled.append({
            "review_token":token(r["photo_id"]),"species":r["species"],
            "original_photo_url":r["photo_url_large"],
            "photo_license":r["photo_license"],"attribution":r["attribution"],
            "reviewer_id":"","open_flower_stage":"","roi_organ":"","notes":"",
        })
    write_rows(a,FIELDS,list(reversed(unfilled)))
    write_rows(b,FIELDS,[unfilled[2],unfilled[0],unfilled[3],unfilled[1]])
    return source,doc,a,b,handoff


def compute(tmp_path,A=None,B=None,third=None,name="out"):
    if A is None:
        A=fixture(tmp_path)
    src,receipt,a,b,handoff=A
    return reconcile(src,receipt,a,b,handoff,tmp_path/name,
                     adjudication=third,_synthetic_test=True)


def answer(path,pid,who,stage,organ):
    r=rows(path)
    t=token(str(pid))
    for row in r:
        if row["review_token"]==t:
            row.update(reviewer_id=who,open_flower_stage=stage,roi_organ=organ)
    write_rows(path,FIELDS,r)


def test_all_unfilled_original_forms_do_not_manufacture_white_or_nonbloom(tmp_path):
    args=fixture(tmp_path)
    r=reconcile(*args,tmp_path/"out",_synthetic_test=True)
    assert r["selected_original_photo_opportunities"]==4
    assert r["n_botanist_A_reviewed_photos"]==0
    assert r["n_botanist_B_reviewed_photos"]==0
    assert r["n_confirmed_open_petaloid_photos_FOR_LATER_ROI"]==0
    assert r["botanical_dispositions"]=={"awaiting_two_independent_reviews":4}
    assert r["photo_colour_white_frequency_estimated"] is False
    assert r["independently_confirmed_white_flower_COLOUR_COUNT"]==0
    assert r["new_colour_ROI_measurements_from_these_400"]==0
    q=rows(tmp_path/"out"/"mois_400_taxon_siteyear_season_reviewer_sample_flow.csv")
    assert len(q)==4
    assert set(q[0])==set(SUMMARY_FIELDS)
    assert all(t["source_original_photo_opportunities"]=="1" for t in q)
    assert all(t["unreviewed_or_one_reviewer_only"]=="1" for t in q)
    assert sha(tmp_path/"out"/"mois_400_source_photo_botanical_qualification_COLOUR_UNMEASURED.csv")==r["photo_qualification_csv_sha256"]
    with pytest.raises(ValueError,match="overwrite"):
        reconcile(*args,tmp_path/"out",_synthetic_test=True)


def test_two_independent_reviewers_agree_open_petal_only_eligible(tmp_path):
    src,receipt,a,b,handoff=fixture(tmp_path)
    for path,who in ((a,"botanist_A"),(b,"botanist_B")):
        answer(path,900001,who,"open_flower","petal_tepal")
        answer(path,900002,who,"fruit_only","fruit_seed")
        answer(path,900003,who,"open_flower","other_flower_organ")
        answer(path,900004,who,"uncertain","uncertain")
    r=reconcile(src,receipt,a,b,handoff,tmp_path/"out",_synthetic_test=True)
    assert r["status"]=="BOTANICAL_QUALIFICATION_COMPLETE_BUT_COLOUR_NOT_MEASURED"
    assert r["n_reviewed_by_both_independently"]==4
    assert r["n_confirmed_open_petaloid_photos_FOR_LATER_ROI"]==1
    assert r["new_colour_ROI_measurements_from_these_400"]==0
    assert r["photo_colour_white_frequency_estimated"] is False
    score={v["photo_id"]:v for v in rows(tmp_path/"out"/"mois_400_source_photo_botanical_qualification_COLOUR_UNMEASURED.csv")}
    assert score["900001"]["botanical_photo_stage_status"]=="botanical_open_petaloid_agreed"
    assert score["900001"]["eligible_independently_verified_open_petaloid_photo_for_later_roi"]=="True"
    assert score["900001"]["botanically_verified_white_flower_hue"]=="False"
    assert score["900002"]["botanical_photo_stage_status"]=="nonopen_bud_fruit_or_no_visible_flower"
    assert score["900003"]["botanical_photo_stage_status"]=="other_open_flower_organ_not_primary"
    assert score["900004"]["botanical_photo_stage_status"]=="uncertain_flower_or_tissue"


def test_one_reviewer_only_does_not_get_confirmed_flower(tmp_path):
    src,receipt,a,b,handoff=fixture(tmp_path)
    answer(a,900001,"botanist_A","open_flower","petal_tepal")
    r=reconcile(src,receipt,a,b,handoff,tmp_path/"out",_synthetic_test=True)
    assert r["n_botanist_A_reviewed_photos"]==1
    assert r["n_botanist_B_reviewed_photos"]==0
    assert r["n_confirmed_open_petaloid_photos_FOR_LATER_ROI"]==0
    assert r["botanical_dispositions"]["awaiting_second_independent_review"]==1


def test_third_reviewer_only_for_genuine_disagreement(tmp_path):
    src,receipt,a,b,handoff=fixture(tmp_path)
    answer(a,900001,"botanist_A","open_flower","petal_tepal")
    answer(b,900001,"botanist_B","fruit_only","fruit_seed")
    base=reconcile(src,receipt,a,b,handoff,tmp_path/"base",_synthetic_test=True)
    assert base["botanical_dispositions"]["reviewer_disagreement_pending"]==1
    assert base["n_confirmed_open_petaloid_photos_FOR_LATER_ROI"]==0
    third=tmp_path/"third.csv"
    write_rows(third,ADJ_FIELDS,[{
        "review_token":token("900001"),"adjudicator_id":"botanist_C",
        "adjudicated_open_flower_stage":"open_flower",
        "adjudicated_roi_organ":"petal_tepal","adjudication_notes":"synthetic agreement",
    }])
    r=reconcile(src,receipt,a,b,handoff,tmp_path/"final",
                adjudication=third,_synthetic_test=True)
    assert r["n_confirmed_open_petaloid_photos_FOR_LATER_ROI"]==1
    assert r["botanical_dispositions"]["botanical_open_petaloid_adjudicated"]==1
    r3=rows(third);r3[0]["adjudicator_id"]="botanist_A"
    write_rows(third,ADJ_FIELDS,r3)
    with pytest.raises(ValueError,match="Third-party"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"bad",
                  adjudication=third,_synthetic_test=True)


def test_wrong_photo_url_or_colour_location_column_detected(tmp_path):
    src,receipt,a,b,handoff=fixture(tmp_path)
    rows_b=rows(b)
    rows_b[0]["original_photo_url"]="https://example.invalid/fake"
    write_rows(b,FIELDS,rows_b)
    with pytest.raises(ValueError,match="metadata changed"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"bad",_synthetic_test=True)
    write_rows(b,tuple(FIELDS)+("observed_on",),[
        {**x,"observed_on":"2025-10-01"} for x in rows_b])
    with pytest.raises(ValueError,match="leakage"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"another",_synthetic_test=True)


def test_reviewer_must_have_stable_and_independent_identity(tmp_path):
    src,receipt,a,b,handoff=fixture(tmp_path)
    answer(a,900001,"same","open_flower","petal_tepal")
    answer(b,900001,"same","open_flower","petal_tepal")
    with pytest.raises(ValueError,match="independent humans"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"bad",_synthetic_test=True)
    answer(b,900001,"other","open_flower","petal_tepal")
    answer(a,900002,"third-person","open_flower","petal_tepal")
    with pytest.raises(ValueError,match="multiple identities"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"bad2",_synthetic_test=True)


def test_malformed_partial_and_unknown_review_codes_fail(tmp_path):
    src,receipt,a,b,handoff=fixture(tmp_path)
    arr=rows(a)
    arr[0]["open_flower_stage"]="open_flower"
    write_rows(a,FIELDS,arr)
    with pytest.raises(ValueError,match="Partial reviewer reply"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"bad",_synthetic_test=True)
    arr[0]["reviewer_id"]="botanist_A"
    arr[0]["roi_organ"]="fictional_white_flower"
    write_rows(a,FIELDS,arr)
    with pytest.raises(ValueError,match="Unregistered"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"bad2",_synthetic_test=True)


def test_source_rewriting_fails_receipt_and_real_guard(tmp_path):
    src,receipt,a,b,handoff=fixture(tmp_path)
    with pytest.raises(ValueError,match="Real source"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"fake_prod")
    src.write_text(src.read_text().replace("2025-07-15","2025-08-15",1))
    with pytest.raises(ValueError,match="not attested"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"bad",_synthetic_test=True)


def test_no_unsigned_adjudication_of_agreement(tmp_path):
    src,receipt,a,b,handoff=fixture(tmp_path)
    answer(a,900001,"A_bot","open_flower","petal_tepal")
    answer(b,900001,"B_bot","open_flower","petal_tepal")
    third=tmp_path/"third.csv"
    write_rows(third,ADJ_FIELDS,[{"review_token":token("900001"),"adjudicator_id":"C_bot",
        "adjudicated_open_flower_stage":"fruit_only",
        "adjudicated_roi_organ":"fruit_seed","adjudication_notes":""}])
    with pytest.raises(ValueError,match="cannot alter agreeing"):
        reconcile(src,receipt,a,b,handoff,tmp_path/"bad",
                  adjudication=third,_synthetic_test=True)
