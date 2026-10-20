"""Only real p0004 original photo pixels are appended; metadata never becomes colour."""
import pandas as pd
import pytest

from scripts.acquisition.mois_append_p0004_actual_original_ROIs import append


def old_source():
    rows = [
        dict(photo_id="1", inat_taxon_id="10", species="Plant alpha",
             observation_id="obs1", latitude=40., longitude=140.,
             observed_on="2024-07-01", morph="white", classified=True,
             white=True, arm="MOIS", measurement_status="classified_four_state_morph",
             technical_failure=False),
        dict(photo_id="2", inat_taxon_id="11", species="Plant beta",
             observation_id="obs2", latitude=41., longitude=141.,
             observed_on="2024-08-01", morph="mixed_uncertain", classified=False,
             white=False, arm="FCP", measurement_status="image_acquisition_failed",
             technical_failure=True),
    ]
    return pd.DataFrame(rows)


def p0004_source():
    return pd.DataFrame([
        dict(photo_id="3", inat_taxon_id="12", species="Plant gamma",
             observation_id="obs3", latitude="42.1", longitude="142.0",
             observed_on="2024-06-01", morph="yellow_orange",
             measurement_status="classified_four_state_morph",
             classified_four_hue="true",
             source_arm="global93m_probability_sample_observations"),
        dict(photo_id="4", inat_taxon_id="12", species="Plant gamma",
             observation_id="obs4", latitude="42.3", longitude="142.1",
             observed_on="2024-06-02", morph="mixed_uncertain",
             measurement_status="image_acquisition_failed",
             classified_four_hue="false",
             source_arm="global93m_probability_sample_observations"),
    ])


def counts():
    return {
        "attempted": 2, "classifiable": 1,
        "acquisition_failed": 1, "unclassifiable_acquired": 0,
        "downloaded": 1,
    }, {"new_photo_attempts_planned": 2}


def test_append_only_real_image_derived_hues_and_preserve_original_photo_ids():
    old = old_source()
    before = old.copy(deep=True)
    new = p0004_source()
    receipt, source = counts()
    merged, report = append(old, new, receipt, source, old_n=2, new_n=2)
    assert merged.photo_id.tolist() == ["1", "2", "3", "4"]
    assert merged.photo_id.is_unique
    assert merged.morph.tolist() == ["white", "mixed_uncertain",
                                    "yellow_orange", "mixed_uncertain"]
    assert merged.classified.tolist() == [True, False, True, False]
    assert merged.white.tolist() == [True, False, False, False]
    assert merged.technical_failure.tolist() == [False, True, False, True]
    assert old.equals(before)
    assert report["combined_unique_actual_original_photo_IDs"] == 4
    assert report["combined_actual_four_hue_classified"] == 2
    assert report["newly_classified_nominal_taxa"] == 1


def test_reprocessed_original_photo_is_rejected():
    old = old_source()
    new = p0004_source()
    new.loc[0, "photo_id"] = "1"
    receipt, source = counts()
    with pytest.raises(ValueError, match="reused|duplicated"):
        append(old, new, receipt, source, old_n=2, new_n=2)


def test_technical_unknown_cannot_be_labeled_chromatic_or_classified():
    old = old_source()
    new = p0004_source()
    new.loc[1, "classified_four_hue"] = "true"
    receipt, source = counts()
    with pytest.raises(ValueError, match="relabelled"):
        append(old, new, receipt, source, old_n=2, new_n=2)


def test_inconsistent_four_hue_receipt_cannot_be_counted():
    old = old_source()
    new = p0004_source()
    receipt, source = counts()
    receipt["classifiable"] = 2
    with pytest.raises(ValueError, match="disagree"):
        append(old, new, receipt, source, old_n=2, new_n=2)
