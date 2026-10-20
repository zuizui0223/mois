"""Synthetic tests: the values below are NOT empirical MOIS flower colours."""
import pandas as pd
import pytest
from scripts.acquisition.mois_append_first_second93m_probability_original_photos import (
    combine,
)


def row(pid,taxon,obs,morph,arm,status="classified_four_state_morph"):
    return {
        "photo_id":str(pid),"inat_taxon_id":str(taxon),
        "species":"Campanula punctata","observation_id":str(obs),
        "latitude":37.0,"longitude":137.0,"observed_on":"2025-07-02",
        "morph":morph,
        "classified":morph in ("white","yellow_orange","red_pink","blue_purple"),
        "white":morph=="white","arm":arm,"measurement_status":status,
        "technical_failure":status=="image_acquisition_failed"
    }


def source():
    old=pd.DataFrame([
        row(1,10,"oldA","white","core"),
        row(2,11,"oldB","red_pink","FCP_final_447500_actual_original_ROI")
    ])
    new=pd.DataFrame([
        row(3,12,"newA","blue_purple","global93m_probability_sample_observations"),
        row(4,13,"newB","mixed_uncertain",
            "global93m_probability_sample_observations","image_acquisition_failed")
    ])
    return old,new


def test_second_photo_wave_appends_new_real_original_ids_without_relabelling():
    old,new=source()
    result,stats=combine(old,new,expected_old=2,expected_new=2)
    assert len(result)==4 and result.photo_id.tolist()==["1","2","3","4"]
    assert stats["all_unique_original_photo_attempts"]==4
    assert stats["new_four_hue_classified"]==1
    assert stats["total_four_hue_classified"]==3
    assert stats["new_unknown_or_technical_failed"]==1
    assert stats["new_acquisition_failed"]==1
    assert stats["new_four_hue_nominal_taxa_not_in_old"]==1
    assert result.loc[result.photo_id.eq("1"),"morph"].iloc[0]=="white"
    assert result.loc[result.photo_id.eq("4"),"classified"].iloc[0]==False
    assert result.loc[result.photo_id.eq("4"),"white"].iloc[0]==False
    assert result.loc[result.photo_id.eq("3"),"arm"].iloc[0]=="global93m_second_p0004"
    assert stats["native_climate_pollination_models_not_run"] is True


def test_duplicate_second_source_photo_ids_are_rejected():
    old,new=source()
    new.loc[1,"photo_id"]="3"
    with pytest.raises(ValueError,match="incomplete or duplicated"):
        combine(old,new,expected_old=2,expected_new=2)


def test_second_source_must_exclude_every_earlier_real_image_photo_id():
    old,new=source()
    new.loc[0,"photo_id"]="1"
    with pytest.raises(ValueError,match="previously genuinely measured original ID"):
        combine(old,new,expected_old=2,expected_new=2)


def test_unknown_image_cannot_be_silently_relabelled_as_white():
    old,new=source()
    new.loc[1,"white"]=True
    with pytest.raises(ValueError,match="technically failed photographic ROIs"):
        combine(old,new,expected_old=2,expected_new=2)


def test_unknown_four_hue_cannot_be_falsely_set_classified():
    old,new=source()
    new.loc[1,"classified"]=True
    with pytest.raises(ValueError,match="Classifiability disagrees"):
        combine(old,new,expected_old=2,expected_new=2)

def test_second_source_string_coordinates_and_climate_keep_old_numeric_and_parquet_roundtrip(tmp_path):
    old, new = source()
    old["latitude"] = pd.Series([37.0, 35.5], dtype="float64")
    old["longitude"] = pd.Series([137.0, 129.5], dtype="float64")
    old["wc_bio1"] = pd.Series([12.0, 11.5], dtype="float64")
    old["soil_pH"] = pd.Series([6.0, 6.5], dtype="float64")
    new["latitude"] = ["29.7127033333", ""]
    new["longitude"] = ["128.123456789", "129.75"]
    new["wc_bio1"] = ["18.7", ""]
    new["soil_pH"] = ["", "5.25"]
    frozen = old.copy(deep=True)
    result, stats = combine(old, new, expected_old=2, expected_new=2)
    assert old.equals(frozen), "Source-signed historical data must not be mutated by a union"
    for col in ("latitude", "longitude", "wc_bio1", "soil_pH"):
        assert pd.api.types.is_numeric_dtype(result[col]), col
    assert result.latitude.iloc[2] == pytest.approx(29.7127033333)
    assert pd.isna(result.latitude.iloc[3])
    assert result.wc_bio1.iloc[2] == pytest.approx(18.7)
    assert pd.isna(result.wc_bio1.iloc[3])
    assert result.soil_pH.iloc[3] == pytest.approx(5.25)
    assert result.photo_id.astype(str).tolist() == ["1", "2", "3", "4"]
    assert stats["new_four_hue_classified"] == 1
    output = tmp_path / "preserved_original_photo_sources.parquet"
    result.to_parquet(output, index=False, compression="zstd")
    loaded = pd.read_parquet(output)
    assert loaded.photo_id.astype(str).tolist() == ["1", "2", "3", "4"]
    assert loaded.latitude.iloc[2] == pytest.approx(29.7127033333)
    assert pd.isna(loaded.latitude.iloc[3])


def test_second_photo_non_numeric_metadata_cannot_be_silently_dropped():
    old, new = source()
    new.loc[0, "latitude"] = "invalid-coordinate"
    with pytest.raises(ValueError, match="Non-numeric"):
        combine(old, new, expected_old=2, expected_new=2)
