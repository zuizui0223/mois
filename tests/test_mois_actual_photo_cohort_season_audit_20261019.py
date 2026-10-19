"""Synthetic contracts for actual-photo arm union and photographic season support."""
from pathlib import Path

import pandas as pd
import pytest

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    assemble, load_pixels, run, season_support, summary,
)


def row(pid, hue, date, taxon="42", obs=None, lat=40.10, lon=140.10):
    return dict(photo_id=str(pid), inat_taxon_id=taxon,
                species="Campanula punctata" if taxon == "42" else "Cirsium brevicaule",
                morph=hue, observation_id=obs or "obs-" + str(pid),
                latitude=lat, longitude=lon, qualified_date=date)


def sources(tmp_path):
    core = [
        row(1, "white", "2025-06-15"),
        row(2, "blue_purple", "2025-01-18"),
        row(3, "white", "2024-12-15"),
        row(4, "mixed_uncertain", ""),
        row(5, "yellow_orange", "2025-07-15", taxon="13", lat=-25.3, lon=151.2),
        row(6, "white", "", taxon="13"),
    ]
    broad = [*core[:4], row(7, "white", "2025-09-03", taxon="13"),
             row(8, "red_pink", "2025-01-15", taxon="13")]
    wave2 = [row(9, "white", "2025-08-05", taxon="13"),
             row(10, "mixed_uncertain", "", taxon="13")]
    cached = [row(11, "white", "2025-09-12"),
              row(12, "mixed_uncertain", "2025-10-06"),
              row(13, "red_pink", "2025-06-11", taxon="13")]
    tail = [row(14, "white", "2025-01-15", taxon="13"),
            row(15, "yellow_orange", "2025-06-05", taxon="13")]
    paths = []
    for name, rows in zip(("core", "broad", "wave2", "cached", "tail"),
                          (core, broad, wave2, cached, tail)):
        p = tmp_path / (name + ".parquet")
        pd.DataFrame(rows).to_parquet(p, index=False)
        paths.append(p)
    return paths


def test_original_overlap_is_reused_not_double_counted(tmp_path):
    c, b, w, e, t = sources(tmp_path)
    d, arms = assemble(c, b, w, e, t, strict=False)
    assert len(d) == 15
    assert len(arms[1]) == 2
    assert len(arms[-1]) == 2
    assert d.photo_id.is_unique
    assert summary(d)["four_hue_classified"] == 12
    assert summary(d)["unclassifiable_or_failed"] == 3
    assert 0 < summary(d)["white_fraction_among_classifiable"] < 1


def test_duplicate_photos_are_rejected(tmp_path):
    c, b, w, e, t = sources(tmp_path)
    v = pd.read_parquet(w)
    v.loc[0, "photo_id"] = "1"
    v.to_parquet(w, index=False)
    with pytest.raises(ValueError, match="same original photo ID"):
        assemble(c, b, w, e, strict=False)


def test_known_hue_cannot_be_an_acquisition_failure(tmp_path):
    c, *_ = sources(tmp_path)
    v = pd.read_parquet(c)
    v["measurement_status"] = "classified_four_state_morph"
    v.to_parquet(c, index=False)
    with pytest.raises(ValueError, match="technical status contradicts"):
        load_pixels(c, "core")


def test_djf_and_jja_same_taxon_local_year_support(tmp_path):
    c, *_ = sources(tmp_path)
    d = load_pixels(c, "core")
    a = season_support(d)
    assert a["north"]["seasons"]["summer"]["photos"] == 1
    assert a["north"]["seasons"]["winter"]["photos"] == 2
    assert a["north"]["taxa_in_both_seasons"] == 1
    assert a["north"]["same_taxon_0p25deg_same_season_year_support"]["1"]["taxon_cell_year_pairs"] == 1
    assert a["north"]["same_taxon_0p25deg_same_season_year_support"]["2"]["taxon_cell_year_pairs"] == 0


def test_real_pixel_receipt_not_source_only_and_no_overwrite(tmp_path):
    c, b, w, e, _ = sources(tmp_path)
    out = tmp_path / "output"
    r = run(c, b, w, e, out, strict=False)
    assert r["actual_original_pixel_photo_union"]["attempted_original_photos"] == 13
    assert r["tail_included"] is False
    assert (out / "actual_photo_colour_and_season_support.json").is_file()
    assert (out / "cohort_white_balance.csv").is_file()
    with pytest.raises(ValueError, match="Never overwrite"):
        run(c, b, w, e, out, strict=False)


def test_legacy_unknown_morph_kept_unclassified_without_relabelling(tmp_path):
    c, *_ = sources(tmp_path)
    frame = pd.read_parquet(c)
    frame.loc[3, "morph"] = "unresolved_legacy_ROI"
    frame["classified"] = frame.morph.isin(
        ("white", "yellow_orange", "red_pink", "blue_purple"))
    frame.to_parquet(c, index=False)
    d = load_pixels(c, "core")
    z = summary(d)
    assert z["original_unclassified_morph_label_counts"]["unresolved_legacy_ROI"] == 1
    assert z["four_hue_classified"] == 5
    assert z["hue_counts"]["white"] == 3
