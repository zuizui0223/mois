"""V2 reviewer HTML cannot contain pre-review source metadata.

These are synthetic photos and fake answers, NOT botanical true labels.
"""
import csv
import json
import re
from pathlib import Path

import pytest
from scripts.analysis import mois_petal_background_contrast_review_v1 as base
from scripts.analysis import mois_petal_background_blind_gate_v2 as new


def roster():
    return [
        {**dict.fromkeys(base.ROSTER_COLUMNS, ""),
         "photo_id": str(120 + i), "observation_id": "obs" + str(i),
         "inat_taxon_id": "taxon" + str(i),
         "source_design": "bulk128_131072",
         "latitude_zone": ("N_arctic" if i == 0 else "S_temperate"),
         "local_season": ("winter" if i == 0 else "summer"),
         "postoutcome_source_gap": "P2_source_has_photos_no_colour_10plus",
         "screen_status": "candidate",
         "technical_diagnostic_class": "detectable_large_flower_box_ROI9_failed",
         "review_priority": "LARGE_VISIBLE_BOX_BUT_NO_STRICT9"}
        for i in range(3)
    ]


def test_blinded_html_contains_photo_id_but_no_region_season_taxon_or_model_status(tmp_path):
    html=Path("tools/review/mois_petal_background_manual_template.html")
    r=roster()
    result=new.build_blind(r,html,tmp_path/"out",source_sha=base.ROSTER_SHA)
    assert result["original_strict9_unknown_unchanged"]
    assert result["botanical_verified_photos"] == 0
    assert result["no_region_or_model_metadata_embedded_in_reviewer_html"]
    for slot in "AB":
        text=(tmp_path/"out"/f"mois_blind_contrast_reviewer_{slot}.html").read_text()
        json_text=re.search(
            r'<script type="application/json" id="MOIS_MANUAL_ROSTER_JSON">(.*?)</script>',
            text,re.S).group(1)
        values=json.loads(json_text)
        assert len(values)==3
        assert set(values[0])=={"photo_id"}
        assert set(z["photo_id"] for z in values)==set(x["photo_id"] for x in r)
        assert 'const fields=["photo_id"]' in text
        assert "MOIS_128_PETAL_BG_MANUAL_V2|" in text
        for forbidden in ("latitude_zone","local_season","inat_taxon_id",
                          "source_design","postoutcome_source_gap",
                          "technical_diagnostic_class","review_priority",
                          "N_arctic","S_temperate"):
            assert forbidden not in text


def test_blind_csv_expansion_matches_immutable_original_roster():
    original=roster()
    form=new.blank_form(original,"A","reviewer-A")
    full=new.expand_blind_review(original,form,"A")
    assert len(full)==3
    assert set(full[0])==set(base.REVIEW_COLUMNS)
    assert full[0]["latitude_zone"]=="N_arctic"
    assert full[0]["source_roster_sha256"]==base.ROSTER_SHA
    assert full[0]["reviewer_slot"]=="A"


def test_reject_unexpected_export_column_and_changed_id():
    original=roster()
    form=new.blank_form(original,"A","reviewer-A")
    form[0]["latitude_zone"]="S_temperate"
    with pytest.raises(ValueError,match="Unexpected blind"):
        new.expand_blind_review(original,form,"A")
    form=new.blank_form(original,"A","reviewer-A")
    form[0]["photo_id"]="999999"
    with pytest.raises(ValueError,match="identity"):
        new.expand_blind_review(original,form,"A")


def test_reject_duplicate_photo_or_reviewer_slot_mismatch():
    original=roster()
    form=new.blank_form(original,"A","reviewer-A")
    form[1]["photo_id"]=form[0]["photo_id"]
    with pytest.raises(ValueError,match="identity"):
        new.expand_blind_review(original,form,"A")
    form=new.blank_form(original,"A","reviewer-A")
    form[0]["reviewer_slot"]="B"
    with pytest.raises(ValueError,match="slot"):
        new.expand_blind_review(original,form,"A")


def test_empty_real_reviews_do_not_create_flower_or_contrast():
    original=roster()
    a=new.expand_blind_review(original,new.blank_form(original,"A","botanist-A"),"A")
    b=new.expand_blind_review(original,new.blank_form(original,"B","botanist-B"),"B")
    result=base.reconcile(original,a,b,expected_roster_sha=base.ROSTER_SHA)
    assert result["confirmed_optical_photo_ids"]==0
    assert result["technical_ROI9_unknown_still_unknown"]==3
    assert result["unresolved_or_disagreement"]==3


def test_private_positive_pair_passes_through_original_numeric_integrity_gate():
    original=roster()
    form_a=new.blank_form(original,"A","botanist-A")
    form_b=new.blank_form(original,"B","botanist-B")
    contrast=base.colour_contrast((255,255,255),(0,0,0))
    for form in (form_a,form_b):
        form[0].update(flower_visible="yes",petal_visible="yes",
                       background_visible="yes",photo_file_sha256="a"*64,
                       image_width="600",image_height="400",image_file_name="test.jpeg",
                       fg_points_json=json.dumps([[20,20,255,255,255],[30,30,255,255,255],[40,40,255,255,255]]),
                       bg_points_json=json.dumps([[100,100,0,0,0],[120,120,0,0,0],[140,140,0,0,0]]),
                       fg_r="255",fg_g="255",fg_b="255",bg_r="0",bg_g="0",bg_b="0",
                       **{key:str(value) for key,value in contrast.items()})
    result=new.reconcile_blind(original,form_a,form_b,source_sha=base.ROSTER_SHA)
    assert result["confirmed_optical_photo_ids"]==1
    assert result["confirmed_records"][0]["deltaE76_mean"]==pytest.approx(100,abs=.3)
    assert result["technical_ROI9_unknown_still_unknown"]==3


def test_same_person_cannot_claim_independent_review():
    original=roster()
    form_a=new.blank_form(original,"A","same")
    form_b=new.blank_form(original,"B","same")
    with pytest.raises(ValueError,match="independent reviewer"):
        new.reconcile_blind(original,form_a,form_b,source_sha=base.ROSTER_SHA)


def test_preflight_is_counts_only_and_never_infers_ecological_claims():
    result=new.preflight(roster())
    assert result["photo_ids"]==3
    assert result["regional_counts"]=={"N_arctic":1,"S_temperate":2}
    assert result["source_design_counts"]=={"bulk128_131072":3}
    assert result["botanically_verified_count"]==0
    assert result["original_strict9_increment"]==0
    assert result["island_mainland_status_in_roster"] is False
    assert result["source_purposeful_ROI9_failure_sample"]


def test_public_aggregate_never_contains_per_photo_or_raw_review():
    original=roster()
    a=new.blank_form(original,"A","a")
    b=new.blank_form(original,"B","b")
    result=new.reconcile_blind(original,a,b,source_sha=base.ROSTER_SHA)
    public=new.aggregate_only(original,result)
    assert "confirmed_records" not in public
    assert "photo_id" not in json.dumps(public)
    assert public["verified_photo_n"]==0
    assert public["regional_descriptive_estimates_eligible"] is False
