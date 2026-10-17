"""Hard gate: completed 4096 image pixels replace partial 512/1792 cohorts."""
import numpy as np
import pandas as pd
import pytest
from scripts.analysis import mois_full313535_immutable_photo_hue_annual_environment as x


def test_real_photo_cohort_and_six_environment_axis_climate_recipe():
    assert x.TOTAL==313535 and x.NEW==4096
    assert x.N_CLASSES_NEW==1577
    assert x.CLASSES==151239 and x.DATES==151230
    assert x.OUT_PHOTO != x.old.P_OUT
    assert x.ORIGINAL_RELEASE.startswith("mois-primary-309439-")


def test_photo_classes_do_not_include_unknown_and_not_double_counted():
    assert "mixed_uncertain" not in x.old.HUES
    assert x.CLASSES == x.old.CLASS + x.N_CLASSES_NEW
    assert x.TOTAL == x.old.N + x.NEW


def test_solar_uses_correct_physical_unit():
    assert x.old.SHA_PHOTO != x.old.SHA_CLIMATE
    assert x.old.N==309439
    assert x.old.ENV[-1]=="wc2_annual_srad_normal_10m_Wm2"
