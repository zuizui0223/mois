"""Synthetic original-photo ID and old-source date preservation tests.

No synthetic data are promoted to the actual MOIS ecological result. The
real pipeline replays signed GitHub original/photo provider source artifacts.
"""
import pandas as pd
import pytest
from scripts.analysis import mois_append_finished_fcp_wave4_to_real295017 as m


COLS=["photo_id","inat_taxon_id","species","morph","source_original_measured",
      "old_mois_or_new_fcp","measurement_status","latitude","longitude",
      "photo_license","observation_id","observed_on","time_observed_at",
      "wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay",
      "soil_available_water_proxy","wc_bio1","wc_bio5","wc_bio12","wc_bio15"]


def record(photo, hue, *, taxon="400", observed_on="2025-06-14",
           original_source="fcp_wave4_partial_original_ROI_v4",
           observation_id=None):
    d={x:"" for x in COLS}
    d.update(photo_id=photo,inat_taxon_id=taxon,species="Campanula sample",
             morph=hue,source_original_measured=original_source,
             old_mois_or_new_fcp="fcp_additional_ROI",
             measurement_status=("classified_four_state_morph" if hue in m.HUES
                                 else "not_evaluable_roi_or_flip_gate"),
             latitude="39.01",longitude="140.01",
             observation_id=observation_id or "obs-"+photo,observed_on=observed_on)
    return d


def fake_sources(monkeypatch):
    monkeypatch.setattr(m,"OLD_ORIGINAL_PHOTOS",3)
    a=pd.DataFrame([
        record("11","white",original_source="other_old_mois"),
        record("12","red_pink",original_source="fcp_wave4_partial_original_ROI_v4"),
        record("13","mixed_uncertain",observed_on="",
               original_source="other_old_mois")])
    a.loc[a.photo_id.eq("13"),"measurement_status"]="not_evaluable_roi_or_flip_gate"
    b=a.copy()
    b["qualified_date"]=["2025-06-14","2025-06-14","2024-01-01"]
    b["original_source_date_provenance"]=["old","old","historic_s3"]
    fcp=pd.DataFrame([a.iloc[0].to_dict(),a.iloc[1].to_dict(),a.iloc[2].to_dict(),
                      record("14","yellow_orange",observed_on="2026-07-18"),
                      record("15","mixed_uncertain",observed_on="2025-11-19")])
    return a,b,fcp


def test_genuine_new_original_source_photo_and_day_preserved(monkeypatch):
    a,b,fcp=fake_sources(monkeypatch)
    original_a=a.copy(deep=True)
    original_b=b.copy(deep=True)
    x,y,stats=m.expand(a,b,fcp,min_partitions=181)
    assert len(x)==len(y)==5
    assert x.photo_id.is_unique
    assert y.photo_id.is_unique
    assert x.iloc[:3].reset_index(drop=True).equals(original_a)
    assert y.iloc[:3].reset_index(drop=True).equals(original_b)
    assert stats["additional_original_photo_pixels"]==2
    assert stats["additional_original_four_hue_ROI_labels"]==1
    assert stats["additional_new_photo_hues_with_original_observed_date"]==1
    assert stats["expanded_dated_original_four_hue_ROIs"]==3
    assert y.loc[y.photo_id.eq("14"),"qualified_date"].iloc[0]=="2026-07-18"
    assert y.loc[y.photo_id.eq("15"),"morph"].iloc[0]=="mixed_uncertain"


def test_exact_photo_not_double_counted_even_when_observation_ids_differ(monkeypatch):
    a,b,fcp=fake_sources(monkeypatch)
    fcp.loc[fcp.photo_id.eq("12"),"observation_id"]="iNat-current-uuid-not-original-numeric"
    _,_,counts=m.expand(a,b,fcp,min_partitions=181)
    assert counts["FCP_original_photo_source_numeric_observation_ID_discrepancies"]==1
    assert counts["additional_original_photo_pixels"]==2


def test_original_taxon_or_hue_disagreement_refused(monkeypatch):
    a,b,fcp=fake_sources(monkeypatch)
    fcp.loc[fcp.photo_id.eq("12"),"morph"]="white"
    with pytest.raises(ValueError,match="morph"):
        m.expand(a,b,fcp,min_partitions=181)
    fcp.loc[fcp.photo_id.eq("12"),"morph"]="red_pink"
    fcp.loc[fcp.photo_id.eq("12"),"inat_taxon_id"]="999"
    with pytest.raises(ValueError,match="inat_taxon_id"):
        m.expand(a,b,fcp,min_partitions=181)


def test_original_observed_on_conflict_refused(monkeypatch):
    a,b,fcp=fake_sources(monkeypatch)
    fcp.loc[fcp.photo_id.eq("12"),"observed_on"]="2025-07-14"
    with pytest.raises(ValueError,match="observed_on"):
        m.expand(a,b,fcp,min_partitions=181)


def test_new_source_is_fcp_real_completed_roi_only(monkeypatch):
    a,b,fcp=fake_sources(monkeypatch)
    fcp.loc[fcp.photo_id.eq("14"),"source_original_measured"]="not_measured"
    with pytest.raises(ValueError,match="FCP ROI-v4"):
        m.expand(a,b,fcp,min_partitions=181)


def test_new_invalid_date_and_unclassified_colours_not_turned_white(monkeypatch):
    a,b,fcp=fake_sources(monkeypatch)
    fcp.loc[fcp.photo_id.eq("14"),"observed_on"]="2026-02-31"
    with pytest.raises(ValueError,match="invalid date"):
        m.expand(a,b,fcp,min_partitions=181)
    fcp.loc[fcp.photo_id.eq("14"),"observed_on"]="2026-07-18"
    fcp.loc[fcp.photo_id.eq("15"),"morph"]="not_in_original_palette"
    with pytest.raises(ValueError,match="outcome states"):
        m.expand(a,b,fcp,min_partitions=181)


def test_no_automatic_promotion_without_new_pixel_photo(monkeypatch):
    a,b,fcp=fake_sources(monkeypatch)
    with pytest.raises(ValueError,match="No genuine new"):
        m.expand(a,b,fcp.iloc[:3],min_partitions=181)
