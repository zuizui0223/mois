"""Mois photo anthesis QC must never transform automated ROI colour into flowers."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path

import pytest

from scripts.analysis.mois_validate_open_flower_photo_review import (
    REVIEW_COLUMNS, digest, decide, make_blank, audit,
)


def csv_bytes(fields, values):
    buf=io.StringIO()
    w=csv.DictWriter(buf,fieldnames=fields,lineterminator="\n")
    w.writeheader();w.writerows(values)
    return buf.getvalue().encode("utf-8")


def source(tmp_path):
    package=tmp_path/"source.zip"
    roi=[
        dict(measurement_id="FCPR-100",photo_id="100",species="Test alba",
             observation_id="obs-1",measurement_status="classified_four_state_morph",
             morph="white",flower_fraction_white="1",flower_fraction_yellow="0"),
        dict(measurement_id="FCPR-101",photo_id="101",species="Test alba",
             observation_id="obs-2",measurement_status="classified_four_state_morph",
             morph="white",flower_fraction_white=".45",flower_fraction_yellow=".35"),
        dict(measurement_id="FCPR-102",photo_id="102",species="Test alba",
             observation_id="obs-3",measurement_status="not_evaluable_roi_or_flip_gate",
             morph="mixed_uncertain",flower_fraction_white="",flower_fraction_yellow=""),
    ]
    urls=[
        dict(measurement_id=x["measurement_id"],
             image_filename=x["measurement_id"]+".jpg",
             photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{x['photo_id']}/large.jpg",
             photo_license="cc-by") for x in roi
    ]
    raw=csv_bytes(list(roi[0]),roi)
    doc={
        "status":"FIRST_GENUINE_MONTHLY_ORIGINAL_GEO_MATCHED_PIXELS_MEASURED_NOT_1000_SPECIES",
        "old_measured_colours_imported":False,
        "new_original_photo_pixels_processed":len(roi),
        "original_geographic_source_join_terminal_sha256":digest(raw),
    }
    with zipfile.ZipFile(package,"w") as z:
        z.writestr("measured/first_original_geographically_matched_photo_terminal.csv",raw)
        z.writestr("prepared/first_matched_pair_original_photo_source_acquisition.csv",
                   csv_bytes(list(urls[0]),urls))
        z.writestr("measured/first_geo_pair_real_pixel_receipt.json",
                   json.dumps(doc))
    return package


def records(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def save(path,records):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=REVIEW_COLUMNS,lineterminator="\n")
        w.writeheader();w.writerows(records)


def set_votes(record,stage="open_flower",tissue="petal_tepal",
              reviewer1="botanist_A",reviewer2="botanist_B"):
    record.update(r1_reviewer_id=reviewer1,r1_open_flower_stage=stage,
                  r1_roi_tissue=tissue,r2_reviewer_id=reviewer2,
                  r2_open_flower_stage=stage,r2_roi_tissue=tissue)


def test_original_photo_blank_review_is_no_biological_evidence(tmp_path):
    z=source(tmp_path);review=tmp_path/"blank.csv"
    frozen=make_blank(z,review)
    assert frozen["original_photos"]==3
    assert frozen["not_true_nonflower_or_white_rate"] is True
    assert frozen["original_unreviewed_flower_stages"]==3
    assert list(records(review)[0])==list(REVIEW_COLUMNS)
    r=audit(review,z,tmp_path/"results")
    assert r["status"]=="BOTANICAL_REVIEW_INCOMPLETE_NO_FLOWER_WHITE_RATE"
    assert r["n_reviewed_primary_flower_organs"]==0
    assert r["botanical_statuses"]=={"unreviewed":3}
    assert r["original_image_pixels_not_reopened"] is True
    assert r["no_white_rate_without_complete_bloom_and_organ_validation"] is True


def test_two_independent_reviewers_protect_open_petal_white_denominator(tmp_path):
    z=source(tmp_path);review=tmp_path/"review.csv"
    make_blank(z,review)
    rows=records(review)
    set_votes(rows[0])
    set_votes(rows[1])
    set_votes(rows[2],stage="fruit_only",tissue="fruit_seed")
    save(review,rows)
    q=audit(review,z,tmp_path/"out")
    assert q["status"]=="BOTANICAL_STAGE_REVIEW_COMPLETED_NOT_CAUSAL"
    assert q["n_reviewed_primary_flower_organs"]==1
    assert q["botanical_statuses"]["petal_tepal_open_flower_reviewed"]==1
    assert q["botanical_statuses"]["anthesis_reviewed_but_original_photo_colour_not_usable"]==1
    assert q["botanical_statuses"]["non_open_flower_reviewed_not_colour"]==1
    found=records(tmp_path/"out"/"audited_first_original64_botanical_flower_stage.csv")
    assert found[0]["eligible_original_open_petaloid_flower_photo_colour"]=="True"
    assert found[1]["eligible_original_open_petaloid_flower_photo_colour"]=="False"
    assert found[2]["eligible_original_open_petaloid_flower_photo_colour"]=="False"


def test_same_reviewer_cannot_supply_both_votes(tmp_path):
    z=source(tmp_path);review=tmp_path/"review.csv"
    make_blank(z,review)
    rows=records(review)
    set_votes(rows[0],reviewer1="same",reviewer2="same")
    save(review,rows)
    with pytest.raises(ValueError,match="independent"):
        audit(review,z,tmp_path/"out")


def test_disagreement_requires_independent_adjudicator(tmp_path):
    z=source(tmp_path);review=tmp_path/"review.csv"
    make_blank(z,review)
    rows=records(review)
    set_votes(rows[0])
    rows[0]["r2_open_flower_stage"]="fruit_only"
    rows[0]["r2_roi_tissue"]="fruit_seed"
    save(review,rows)
    q=audit(review,z,tmp_path/"out")
    assert q["status"]=="BOTANICAL_REVIEW_INCOMPLETE_NO_FLOWER_WHITE_RATE"
    assert q["botanical_statuses"]["needs_independent_adjudication"]==1
    rows[0]["adjudicator_id"]="third_botanist"
    rows[0]["adjudicated_stage"]="open_flower"
    rows[0]["adjudicated_roi_tissue"]="petal_tepal"
    save(review,rows)
    q=audit(review,z,tmp_path/"resolved")
    assert q["n_reviewed_primary_flower_organs"]==1
    rows[0]["adjudicator_id"]="botanist_A"
    save(review,rows)
    with pytest.raises(ValueError,match="third independent"):
        audit(review,z,tmp_path/"bad")


def test_source_identity_and_receipt_attestation_required(tmp_path):
    z=source(tmp_path);review=tmp_path/"review.csv"
    make_blank(z,review)
    rows=records(review)
    rows[0]["photo_id"]="FAKE"
    save(review,rows)
    with pytest.raises(ValueError,match="identity changed"):
        audit(review,z,tmp_path/"bad1")
    with zipfile.ZipFile(z) as orig:
        content={n:orig.read(n) for n in orig.namelist()}
    content["measured/first_original_geographically_matched_photo_terminal.csv"]+=b"\n"
    corrupt=tmp_path/"bad.zip"
    with zipfile.ZipFile(corrupt,"w") as out:
        for n,b in content.items():out.writestr(n,b)
    with pytest.raises(ValueError,match="Changed original photo"):
        make_blank(corrupt,tmp_path/"bad.csv")


def test_decision_matrix_never_calls_uncertain_buds_or_bracts_usable():
    assert decide("open_flower","petal_tepal","open_flower","petal_tepal")=="petal_tepal_open_flower_reviewed"
    assert decide("open_flower","spikelet_glume_bract","open_flower","spikelet_glume_bract")=="noncomparable_bract_fruit_foliage_or_background"
    assert decide("buds_only","petal_tepal","buds_only","petal_tepal")=="non_open_flower_reviewed_not_colour"
    assert decide("uncertain","uncertain","uncertain","uncertain")=="unresolved_flower_stage_or_roi"
    assert decide("open_flower","petal_tepal","fruit_only","fruit_seed")=="needs_independent_adjudication"
    assert decide("","","","")=="unreviewed"
