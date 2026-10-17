"""Source contract for independently completed ORIGINAL photo RGB pixels."""
from __future__ import annotations
from pathlib import Path
import pytest
from scripts.analysis import mois_checkpoint_real_fourth_photo_pixel_shards as m


def test_identity_boundaries_never_call_preimage_source_a_real_pixel():
    assert m.PRIOR_PHOTOS==309439
    assert m.PRIOR_CLASSIFIED==149662
    assert m.PRIOR_DATED_CLASSIFIED==149653
    assert m.SELECTION_SHA=="65fc139909f94448e042b6539c6749e70d1ab21bb8b13f50383bb7417180337a"
    assert m.METHOD["location_blind"]
    assert m.METHOD["no_old_measurement_csv_inputs"]
    assert m.MEASURED_SCHEMA=="mois_256_source_original_roi_v4_four_hue_terminal_v1"


def test_two_real_worker_shards_required_not_4096_unmeasured_source_ids(tmp_path):
    plan={0:{},1:{}}
    with pytest.raises(ValueError,match="Need at least 2"):
        m.gather_workers(plan,tmp_path,{"per_shard_sha256":{"0":"x","1":"y"}})


def test_worker_shard_source_not_present_is_not_a_colour_classification(tmp_path):
    with pytest.raises(ValueError,match="genuine"):
        m.complete_worker({},tmp_path,0,{"per_shard_sha256":{"0":"unknown"}})


def test_compressed_source_hue_count_is_149662_not_3430502():
    assert len(m.HUES)==4
    assert "mixed_uncertain" not in m.HUES
    assert m.PRIOR_PHOTOS>m.PRIOR_CLASSIFIED>m.PRIOR_DATED_CLASSIFIED
