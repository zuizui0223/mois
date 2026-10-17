"""Guard original source-photo third-wave photo-ID disjointness and measurement scope."""
from scripts.analysis.mois_full3430502_source_third4096 import (
    FULL_SHA,FIRST_SELECTED_SHA,SOURCE_N,UNMEASURED_FULL_N,
    UNMEASURED_OUTSIDE_FIRST_SAMPLE_N,WAVE_N,N_SHARDS,
    PHOTOS_PER_SHARD,SALT,measure_id
)


def test_source_census_not_floral_colour_labels():
    assert SOURCE_N==3430502
    assert UNMEASURED_FULL_N==3420814
    assert UNMEASURED_OUTSIDE_FIRST_SAMPLE_N==3315910
    assert 3315910-4096==3311814
    assert WAVE_N==N_SHARDS*PHOTOS_PER_SHARD==4096
    assert len(FULL_SHA)==len(FIRST_SELECTED_SHA)==64


def test_unseen_original_third_wave_never_reuses_second_wave_measurement_id_salt():
    assert "WAVE3" in SALT
    assert measure_id("123") == measure_id("123")
    assert measure_id("123") != measure_id("124")
    assert measure_id("123").startswith("MOIS-")
