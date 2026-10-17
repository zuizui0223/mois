"""Source pool quotas/canonical IDs and source-only full-archive expansion."""
from scripts.analysis.mois_full3430502_source_second4096 import (
    FULL_SHA,FIRST_SELECTED_SHA,SOURCE_N,SCREEN_SELECTED_N,
    UNMEASURED_FULL_N,UNMEASURED_OUTSIDE_FIRST_SAMPLE_N,
    WAVE_N,N_SHARDS,PHOTOS_PER_SHARD,measure_id
)

def test_real_first_full_source_population_is_not_a_million_flower_measurements():
    assert SOURCE_N==3430502
    assert UNMEASURED_FULL_N==3420814
    assert UNMEASURED_OUTSIDE_FIRST_SAMPLE_N==3315910
    assert SCREEN_SELECTED_N==106376
    assert WAVE_N==N_SHARDS*PHOTOS_PER_SHARD==4096
    assert len(FULL_SHA)==len(FIRST_SELECTED_SHA)==64

def test_source_photo_identifiers_do_not_infer_floral_stage_or_hue():
    assert measure_id("1").startswith("MOIS-")
    assert measure_id("123")==measure_id("123")
    assert measure_id("123")!=measure_id("124")
