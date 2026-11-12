"""Source-frozen genuine 128 photo review: no inferred botany or colour."""
import csv
import json
from pathlib import Path

import pytest

from scripts.analysis import mois_petal_background_contrast_review_v1 as m


def roster(n=3):
    return [dict(photo_id=str(i+100), observation_id="obs"+str(i),
                 inat_taxon_id="tax"+str(i), source_design="historical_K1_FCP_source_heterogeneous",
                 latitude_zone="N_temperate", local_season="summer",
                 postoutcome_source_gap="P2_source_has_photos_no_colour_10plus",
                 screen_status="candidate",
                 technical_diagnostic_class="detectable_large_flower_box_ROI9_failed",
                 review_priority="TEST") for i in range(n)]


def completed(rows, slot, who, positive=("100",)):
    out=[]
    for r in rows:
        z={**r, "source_roster_sha256":m.ROSTER_SHA,
           "reviewer_slot":slot,"reviewer_id":who,
           "flower_visible":"no","petal_visible":"no","background_visible":"no",
           "photo_file_sha256":"","image_width":"","image_height":"",
           "image_file_name":"","fg_points_json":"","bg_points_json":"",
           **{col:"" for col in ("fg_r","fg_g","fg_b","bg_r","bg_g","bg_b",
                                 "deltaE76","deltaLstar","sRGB_luminance_ratio")},
           "notes":""}
        if r["photo_id"] in positive:
            metric=m.colour_contrast((255,255,255),(0,0,0))
            z.update(flower_visible="yes",petal_visible="yes",
                     background_visible="yes",photo_file_sha256="a"*64,
                     image_width="600",image_height="400",image_file_name="100.jpeg",
                     fg_points_json=json.dumps([[20,20,255,255,255],[30,30,255,255,255],[40,40,255,255,255]]),
                     bg_points_json=json.dumps([[100,100,0,0,0],[120,120,0,0,0],[140,140,0,0,0]]),
                     fg_r="255",fg_g="255",fg_b="255",bg_r="0",bg_g="0",bg_b="0",
                     **{k:str(v) for k,v in metric.items()})
        out.append(z)
    return out


def test_white_black_lab_d65_and_luminance():
    x=m.colour_contrast((255,255,255),(0,0,0))
    assert x["deltaE76"]==pytest.approx(100,abs=.3)
    assert x["deltaLstar"]==pytest.approx(100,abs=.3)
    assert x["sRGB_luminance_ratio"]==pytest.approx(21,abs=.01)
    y=m.colour_contrast((60,100,140),(60,100,140))
    assert y["deltaE76"]<1e-9
    assert y["sRGB_luminance_ratio"]==pytest.approx(1)


def test_all_blank_but_distinct_source_review_packets(tmp_path):
    p=tmp_path/"review_roster.csv"
    with p.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=m.ROSTER_COLUMNS)
        writer.writeheader();writer.writerows(roster())
    original_sha=m.sha(p)
    template=Path("tools/review/mois_petal_background_manual_template.html")
    result=m.build(p,template,tmp_path/"out",expected_sha=original_sha,expected_n=3)
    assert result["n_original_unknown_ROI9_photo_ids"]==3
    assert result["botanical_confirmed_count"]==0
    assert result["optical_contrast_valid_count"]==0
    assert result["original_strict9_unchanged"]
    for slot in ("A","B"):
        html=(tmp_path/"out"/("mois_manual_contrast_reviewer_"+slot+".html")).read_text()
        assert "reviewer_"+slot in html and "MOIS_MANUAL_ROSTER_JSON" in html
        assert "fetch(" not in html and "__SOURCE_ROWS__" not in html


def test_roster_source_sha_fails_closed(tmp_path):
    p=tmp_path/"r.csv"
    with p.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=m.ROSTER_COLUMNS)
        writer.writeheader();writer.writerows(roster())
    with pytest.raises(ValueError,match="SHA"):
        m.build(p,Path("tools/review/mois_petal_background_manual_template.html"),
                tmp_path/"out",expected_sha="0"*64,expected_n=3)


def test_identical_reviewers_not_independent():
    r=roster()
    with pytest.raises(ValueError,match="independent reviewer"):
        m.reconcile(r,completed(r,"A","same"),completed(r,"B","same"))


def test_two_independent_real_votes_only_one_paired_photo():
    r=roster()
    p=m.reconcile(r,completed(r,"A","botanist_A"),completed(r,"B","botanist_B"),
                  expected_roster_sha=m.ROSTER_SHA)
    assert p["confirmed_optical_photo_ids"]==1
    assert p["technical_ROI9_unknown_still_unknown"]==3
    assert p["confirmed_records"][0]["deltaE76_mean"]==pytest.approx(100,abs=.3)
    assert not p["claim_botanically_verified_colour_distribution"]
    assert not p["claim_world_native_flora_conspicuousness"]


def test_disagreement_or_different_image_sha_not_verified():
    r=roster()
    x=completed(r,"A","A")
    y=completed(r,"B","B",positive=())
    assert m.reconcile(r,x,y)["confirmed_optical_photo_ids"]==0
    y=completed(r,"B","B")
    y[0]["photo_file_sha256"]="b"*64
    result=m.reconcile(r,x,y)
    assert result["different_image_bytes"]==1
    assert result["confirmed_optical_photo_ids"]==0


def test_tampered_rgb_or_calculated_values_rejected():
    r=roster();x=completed(r,"A","A");y=completed(r,"B","B")
    x[0]["deltaE76"]="1"
    with pytest.raises(ValueError,match="derived colour metric"):
        m.reconcile(r,x,y)
    x=completed(r,"A","A")
    x[0]["fg_points_json"]=json.dumps([[20,20,255,255,255],[30,30,255,255,255],[800,800,255,255,255]])
    with pytest.raises(ValueError,match="point outside"):
        m.reconcile(r,x,y)
    x=completed(r,"A","A")
    x[0]["fg_points_json"]=x[0]["bg_points_json"]
    with pytest.raises(ValueError,match="shared foreground/background"):
        m.reconcile(r,x,y)


def test_original_source_identity_drift_rejected():
    r=roster();x=completed(r,"A","A");y=completed(r,"B","B")
    x[0]["inat_taxon_id"]="wrong"
    with pytest.raises(ValueError,match="source identity"):
        m.reconcile(r,x,y)
    x=completed(r,"A","A");x[0]["source_roster_sha256"]="b"*64
    with pytest.raises(ValueError,match="roster SHA"):
        m.reconcile(r,x,y,expected_roster_sha=m.ROSTER_SHA)


def test_unknown_and_missing_points_never_become_colour():
    r=roster()
    x=completed(r,"A","A",positive=())
    y=completed(r,"B","B",positive=())
    assert m.reconcile(r,x,y)["confirmed_optical_photo_ids"]==0
    x=completed(r,"A","A")
    x[0]["flower_visible"]="uncertain"
    with pytest.raises(ValueError,match="Unverified plant organ"):
        m.reconcile(r,x,completed(r,"B","B"))
