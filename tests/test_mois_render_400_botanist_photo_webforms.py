"""MOIS independently blinded, licensed real original-photo reviewer HTML contracts."""
from __future__ import annotations

import csv
import json
import re

import pytest

from scripts.analysis.mois_prepare_400_botanical_review_handoff import FIELDS
from scripts.analysis.mois_render_400_botanist_photo_webforms import build,sha


def fixture(tmp_path):
    a=tmp_path/"A.csv"
    b=tmp_path/"B.csv"
    receipt=tmp_path/"receipt.json"
    records=[]
    for i in (1,2):
        records.append({
            "review_token":f"{i:020x}","species":"Synthetic fleur",
            "original_photo_url":
              f"https://inaturalist-open-data.s3.amazonaws.com/photos/{i}/large.jpg",
            "photo_license":"cc-by",
            "attribution":f"orig photo {i} CC-BY </script> sample",
            "reviewer_id":"","open_flower_stage":"","roi_organ":"","notes":"",
        })
    def save(path,items):
        with path.open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=FIELDS,lineterminator="\n")
            w.writeheader();w.writerows(items)
    save(a,records)
    save(b,list(reversed(records)))
    receipt.write_text(json.dumps({
        "status":"400_ORIGINAL_PHOTO_REVIEW_FORMS_BOTH_UNFILLED_NOT_COLOUR",
        "reviewer_forms_season_date_climate_colour_blind":True,
        "all_review_response_fields_blank":True,
        "n_actual_independent_botanical_reviews_completed":0,
        "new_photo_pixels_measured":0,
        "n_original_photos":2,
        "n_original_named_taxa":1,
        "reviewer_form_A_SHA256":sha(a),
        "reviewer_form_B_SHA256":sha(b),
    }))
    return a,b,receipt


def payload(html):
    m=re.search(r'<script type="application/json" id="mois-original-blind-source">(.*?)</script>',
                html,re.S)
    assert m
    return json.loads(m.group(1))


def test_two_fully_independent_unfilled_photo_reviewers_html(tmp_path):
    a,b,receipt=fixture(tmp_path)
    result=build(a,b,receipt,tmp_path/"out",_synthetic_test=True)
    assert result["status"]=="TWO_STANDALONE_BROWSER_REVIEW_PORTALS_CREATED_WITH_ZERO_BOTANICAL_DECISIONS"
    assert result["original_botanical_review_response_count"]==0
    assert result["no_seasonal_white_fraction_or_climate_effect_inferred"] is True
    assert result["no_original_image_bytes_downloaded_into_artifact"] is True
    assert result["html_contains_client_only_csv_export_import_not_server_upload"] is True
    for role in ("A","B"):
        f=tmp_path/"out"/result["portals"][role]["filename"]
        h=f.read_text()
        assert sha(f)==result["portals"][role]["html_sha256"]
        assert "CSVを書き出す" in h and "CSVを読み込む" in h
        assert "referrerpolicy=\"no-referrer\"" in h
        assert '<script type="application/json"' in h
        assert 'URL.createObjectURL' in h
        assert 'new Blob([' in h
        assert 'fetch(' not in h and 'XMLHttpRequest' not in h
        d=payload(h)
        assert d["role"]==role
        assert d["fields"]==list(FIELDS)
        assert len(d["records"])==2
        assert all(not row["reviewer_id"] and not row["open_flower_stage"]
                   and not row["roi_organ"] for row in d["records"])
        assert all(set(row)==set(FIELDS) for row in d["records"])
        assert all("observed_on" not in row and "selected_season" not in row
                   and "latitude" not in row and "morph" not in row
                   and "climate" not in row for row in d["records"])
        assert "</script>" not in h[h.index('id="mois-original-blind-source">')+len('id="mois-original-blind-source">'):h.index('</script>',h.index('id="mois-original-blind-source">'))]
        assert "\\u003c/script\\u003e" in h


def test_file_hash_or_nonblank_old_annotations_are_rejected(tmp_path):
    a,b,receipt=fixture(tmp_path)
    raw=csv.DictReader(a.open())
    items=list(raw)
    items[0]["reviewer_id"]="accidental reviewer"
    with a.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(items)
    with pytest.raises(ValueError,match="BLANK"):
        build(a,b,receipt,tmp_path/"out",_synthetic_test=True)


def test_forged_old_source_or_real_400_boundary_rejected(tmp_path):
    a,b,receipt=fixture(tmp_path)
    with pytest.raises(ValueError,match="Original 400"):
        build(a,b,receipt,tmp_path/"out",_synthetic_test=False)
    d=json.loads(receipt.read_text())
    d["reviewer_forms_season_date_climate_colour_blind"]=False
    receipt.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="zero-response"):
        build(a,b,receipt,tmp_path/"bad",_synthetic_test=True)


def test_unlicensed_or_nonoriginal_host_cannot_be_embedded(tmp_path):
    a,b,receipt=fixture(tmp_path)
    text=a.read_text()
    a.write_text(text.replace("inaturalist-open-data.s3.amazonaws.com","example.invalid"))
    with pytest.raises(ValueError,match="Unverified original image"):
        build(a,b,receipt,tmp_path/"out",_synthetic_test=True)


def test_reviewer_source_order_can_differ_but_photo_set_must_match(tmp_path):
    a,b,receipt=fixture(tmp_path)
    records=list(csv.DictReader(b.open()))
    records[0]["review_token"]="000000000000000000ff"
    with b.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(records)
    with pytest.raises(ValueError,match="same original photo IDs"):
        build(a,b,receipt,tmp_path/"out",_synthetic_test=True)


def test_repeated_review_form_identity_rejected(tmp_path):
    a,b,receipt=fixture(tmp_path)
    records=list(csv.DictReader(a.open()))
    records[1]["review_token"]=records[0]["review_token"]
    with a.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(records)
    with pytest.raises(ValueError,match="not unique"):
        build(a,b,receipt,tmp_path/"out",_synthetic_test=True)


def test_no_overwrite_of_existing_biological_qualifications(tmp_path):
    a,b,receipt=fixture(tmp_path)
    build(a,b,receipt,tmp_path/"out",_synthetic_test=True)
    with pytest.raises(ValueError,match="overwrite"):
        build(a,b,receipt,tmp_path/"out",_synthetic_test=True)
