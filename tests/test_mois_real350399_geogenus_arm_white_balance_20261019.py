"""Synthetic (never ecological) source-arm matched geographic/genus arithmetic."""
import numpy as np
import pandas as pd

from scripts.analysis.mois_real350399_geogenus_arm_white_balance_20261019 import (
    HUES, group_hue, matched_contrasts, source_composition,
)


def row(pid, cohort, genus, taxon, hue, lat=35.1, lon=136.2):
    return {
        "photo_id": str(pid),
        "arm": "core" if cohort == "core" else "tail16384",
        "cohort": cohort,
        "inat_taxon_id": str(taxon),
        "species": genus + " species",
        "genus": genus,
        "latitude": lat, "longitude": lon,
        "valid_genus": bool(genus),
        "valid_taxon": bool(str(taxon)),
        "valid_geography": -90 <= lat <= 90 and -180 <= lon <= 180,
        "morph": hue,
        "classified": hue in HUES,
        "white": hue == "white",
        "observation_id": f"obs-{pid}"
    }


def sample():
    # The first two genera share geography, but have entirely different species
    # in the precolour taxonomic recruitment arms. An extra novel genus with
    # WHITE makes the raw arm comparison positive but does not alter any
    # matched genus+geocell mean.
    rows = [
        row(1, "core", "Campanula", "c1", "white"),
        row(2, "core", "Campanula", "c2", "white"),
        row(3, "core", "Cirsium", "i1", "blue_purple"),
        row(4, "core", "Cirsium", "i2", "red_pink"),
        row(5, "new_species_first", "Campanula", "c3", "red_pink"),
        row(6, "new_species_first", "Campanula", "c4", "red_pink"),
        row(7, "new_species_first", "Cirsium", "i3", "white"),
        row(8, "new_species_first", "Cirsium", "i4", "white"),
        row(9, "new_species_first", "Rosa", "r1", "white"),
        row(10, "new_species_first", "Rosa", "r2", "white"),
        row(11, "new_species_first", "Rosa", "r3", "white"),
        row(12, "new_species_first", "Rosa", "r4", "white"),
        # Different site no longer matched by genus even if taxon name repeats.
        row(13, "core", "Campanula", "c5", "white", lat=65.1),
        row(14, "new_species_first", "Campanula", "c6", "white", lat=-35.2),
    ]
    return pd.DataFrame(rows)


def test_matched_geogenus_prevents_unmatched_rose_inflation():
    d = sample()
    x = d[d.classified]
    z = matched_contrasts(x, tile=10, genus=True, min_taxa=2, draws=99)
    assert z["matched_strata"] == 2
    assert z["matched_nominal_genera"] == 2
    assert z["matched_stratum_equal_white_difference_pp"] == 0
    assert z["genus_equal_white_difference_pp"] == 0
    assert z["genus_equal_bootstrap_white_95_pp"][0] < 0
    assert z["genus_equal_bootstrap_white_95_pp"][1] > 0
    assert np.isclose(sum(z["matched_hue_differences_pp"].values()), 0, atol=1e-10)
    assert z["retained_taxon_cell_share_by_arm"]["new_species_first"] < 1
    raw = source_composition(x, d)
    assert raw["raw_new_minus_core_white_percentage_points"] > 0
    assert raw["overlapping_classifiable_named_taxon_ids"] == 0
    assert raw["overlapping_classifiable_nominal_genera"] == 2


def test_zero_support_not_interpreted_as_zero_ecological_effect():
    x = sample()
    s = matched_contrasts(x, tile=5, genus=True, min_taxa=2, draws=19)
    assert s["matched_strata"] == 2
    t = x[x.genus.eq("Rosa") & x.classified]
    result = matched_contrasts(t, tile=10, genus=True, min_taxa=1, draws=19)
    assert result["matched_strata"] == 0
    assert result["matched_hue_differences_pp"] is None
    assert result["genus_equal_bootstrap_white_95_pp"] is None


def test_repeat_original_photos_do_not_create_new_named_taxa():
    x = sample()
    duplicated = x.iloc[[0]].copy()
    duplicated.loc[:, "photo_id"] = "fresh-original-photo-id"
    duplicated.loc[:, "morph"] = "blue_purple"
    q = pd.concat((x, duplicated), ignore_index=True)
    matched = q[q.classified]
    cells, n_taxon_cells = group_hue(
        matched, ["cohort", "genus"])
    before = group_hue(x[x.classified], ["cohort", "genus"])[0]
    assert int(cells.taxa.sum()) == int(before.taxa.sum())
    assert n_taxon_cells == int(before.taxa.sum())
    c = cells[(cells.cohort == "core") & (cells.genus == "Campanula")].iloc[0]
    assert int(c.photos) > int(
        before[(before.cohort == "core") &
               (before.genus == "Campanula")].photos.iloc[0])
    assert c.p_white < 1


def test_unknown_hue_not_false_chromatic_four_hue():
    d = sample()
    d.loc[d.photo_id.eq("9"), "morph"] = "mixed_uncertain"
    d.loc[d.photo_id.eq("9"), "classified"] = False
    d.loc[d.photo_id.eq("9"), "white"] = False
    classified = d[d.classified]
    assert "9" not in set(classified.photo_id)
    result = matched_contrasts(classified, tile=10, genus=True, min_taxa=1, draws=9)
    assert np.isclose(sum(result["matched_hue_differences_pp"].values()), 0, atol=1e-9)
