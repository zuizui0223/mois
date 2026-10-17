"""Fixtures below are synthetic identity QA, not ecological source measurements."""
import pandas as pd
import pytest

from scripts.analysis import mois_merge_real_FCP180_and_MOIS_three4096_original_photos as m


def fake_row(pid,oid,source,morph="white",date="2025-07-11",taxon="123"):
    row={k:"" for k in m.OUTCOLS}
    row.update({"photo_id":str(pid),"observation_id":str(oid),
                "inat_taxon_id":str(taxon),"species":"Example taxon",
                "morph":morph,"observed_on":date,
                "old_mois_or_new_fcp":source,
                "measurement_status":"classified_four_state_morph",
                "latitude":42.1,"longitude":134.1})
    return row


def sample(monkeypatch):
    monkeypatch.setattr(m,"FCP_BASE_N",2)
    monkeypatch.setattr(m,"FCP_BASE_HUES",2)
    monkeypatch.setattr(m,"MOIS_THREE_N",2)
    monkeypatch.setattr(m,"MOIS_THREE_HUES",2)
    monkeypatch.setattr(m,"MOIS_THREE_NEW_N",1)
    monkeypatch.setattr(m,"MOIS_THREE_NEW_HUES",1)
    monkeypatch.setattr(m,"EXPECTED_SHARED",1)
    monkeypatch.setattr(m,"UNION_N",3)
    monkeypatch.setattr(m,"UNION_HUES",3)
    # One measured photo ID reused by two source pipelines with a numeric
    # iNaturalist observation ID versus source observation UUID respectively.
    a=pd.DataFrame([
        fake_row("232092182","136022799","fcp_additional_ROI",
                 morph="yellow_orange",date="2022-09-19"),
        fake_row("9999999","2999999","new_fcp")])[list(m.OUTCOLS)]
    b=pd.DataFrame([
        fake_row("232092182","2307d7f1-2ef4-4e25-b38f-210c6fd965c5",
                 "mois_new_original_ROI",morph="yellow_orange",
                 date="2022-09-19"),
        fake_row("8888888","3999999","mois_new_original_ROI",morph="red_pink")
    ])[list(m.OUTCOLS)]
    return a,b


def test_true_pixel_photo_id_dedup_and_source_observation_namespaces(monkeypatch):
    a,b=sample(monkeypatch)
    source,ledger=m.photo_union(a,b)
    assert len(source)==3
    assert source.photo_id.is_unique
    assert source.morph.isin(m.HUES).sum()==3
    assert ledger["new_original_photo_pixels"]==1
    assert ledger["new_original_four_hue_photo_labels"]==1
    assert ledger["source_shared_conflicts"][
        "observation_ID_source_namespace_mismatches"]==1
    assert ledger["source_shared_conflicts"]["duplicate_photo_IDs_placed_in_result_twice"]==0
    assert source.loc[source.photo_id.eq("232092182"),
           "observation_id"].iloc[0]=="136022799"


def test_conflicting_original_photo_hue_cannot_be_merged_as_independent(monkeypatch):
    a,b=sample(monkeypatch)
    b.loc[b.photo_id.eq("232092182"),"morph"]="blue_purple"
    with pytest.raises(ValueError,match="morph"):
        m.photo_union(a,b)


def test_unexpected_obs_identity_namespace_difference_rejected(monkeypatch):
    a,b=sample(monkeypatch)
    b.loc[b.photo_id.eq("232092182"),"observation_id"]="111112222"
    with pytest.raises(ValueError,match="namespace discrepancy"):
        m.photo_union(a,b)


def test_historical_public_dates_are_sidecar_and_missing_hues_remain_unknown(monkeypatch):
    monkeypatch.setattr(m,"DATE_SOURCE_N",1)
    monkeypatch.setattr(m,"NEW_DATE_ATTEMPTS",3)
    monkeypatch.setattr(m,"NEW_DATE_HUES",2)
    source=pd.DataFrame([
        fake_row("1","100","historical",morph="white",date=None),
        fake_row("2","200","mois_new_original_ROI",morph="blue_purple"),
        fake_row("3","300","new_fcp",morph="mixed_uncertain")])[list(m.OUTCOLS)]
    hist=pd.DataFrame({
        "photo_id":["1"],"original_numeric_observation_id":["100"],
        "original_taxon_id":["123"],"original_frozen_ROI_morph":["white"],
        "qualified_observed_on":["2021-09-09"]})
    out=m.restore_historic_dates(source,hist)
    assert len(out)==3
    assert out.loc[out.photo_id.eq("1"),"observed_on"].isna().all()
    assert out.loc[out.photo_id.eq("1"),"qualified_date"].iloc[0]=="2021-09-09"
    assert out.loc[out.photo_id.eq("1"),"original_source_date_provenance"].iloc[0]==(
        "historic_snapshot_verified_S3_or_original_API")
    stats=m.season_table(out)
    assert stats.original_four_hue_classifiable.sum()==2
    assert stats.original_white_photo_labels.sum()==1
    assert stats.dated_original_photo_attempts.sum()==3


def test_historic_source_photo_taxonomy_mismatch_fails(monkeypatch):
    monkeypatch.setattr(m,"DATE_SOURCE_N",1)
    source=pd.DataFrame([fake_row("1","100","historical",date=None)])[list(m.OUTCOLS)]
    hist=pd.DataFrame({
        "photo_id":["1"],"original_numeric_observation_id":["100"],
        "original_taxon_id":["99999"],"original_frozen_ROI_morph":["white"],
        "qualified_observed_on":["2021-09-09"]})
    with pytest.raises(ValueError,match="taxon identity conflict"):
        m.restore_historic_dates(source,hist)
