"""Synthetic complete photo terminal reassembly; no real source photos or ecology."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

executor_test=module("mois_synthetic_run",ROOT/"tests/test_mois_execute_sample_roi_v4_20261016.py")
m=module("mois_terminal_reassembly",ROOT/"scripts/analysis/mois_collect_sample_roi_white_20261016.py")


def data(root:Path):
    private,packets=executor_test.create(root)
    units=root/"sealed_units"
    executor_test.run_all(private,units,packets)
    return private,units,packets


def test_complete_blind_measurements_join_only_after_source_hard_gate(tmp_path):
    private,units,packet=data(tmp_path)
    out=tmp_path/"sample_white"
    rec=m.collect(private,units,out,synthetic=True)
    assert rec["status"]==m.SEAL_STATUS
    assert rec["original_photo_terminals_verified"]==3
    assert rec["original_colour_four_state_classifiable"]==2
    assert rec["original_image_acquisition_permanent_failure"]==1
    assert rec["species_with_at_least_one_classifiable_photo"]==2
    assert rec["all13_photo_site_predictors_joined"] is False
    assert rec["full_FCP_93m_photos_measured"] is False
    rows=list(csv.DictReader((out/"private_species_conditional_photo_white.csv").open()))
    assert len(rows)==2
    assert {r["genus"] for r in rows}=={"Alpha","Beta"}
    alpha=next(x for x in rows if x["genus"]=="Alpha")
    assert alpha["selected_observations"]=="2"
    assert alpha["classified_four_colour"]=="1"
    assert alpha["classified_white"]=="1"
    assert abs(float(alpha["photo_weighted_classifiability"])-.5)<1e-12
    assert abs(float(alpha["white_given_classified_photo_weighted"])-1)<1e-12
    seen=list(csv.DictReader((out/"private_joined_sample_flower_colours.csv").open()))
    assert len(seen)==3
    assert set(seen[0])==set(m.JOIN+m.PRIVATE_JOIN_EXTRA)
    assert all(len(json.loads(v["nine_colour_fractions_json"]))==9 for v in seen)
    with pytest.raises(ValueError,match="previous final"):
        m.collect(private,units,out,synthetic=True)


def test_missing_batch_does_not_become_a_biological_data_result(tmp_path):
    private,units,packet=data(tmp_path)
    sub=next(units.glob("terminal_p*_b*"))
    rename=sub.with_name(sub.name+".TEMPORARY")
    sub.rename(rename)
    with pytest.raises(ValueError,match="Incomplete sample"):
        m.collect(private,units,tmp_path/"result",synthetic=True)
    assert not (tmp_path/"result").exists()


def test_relabelled_colour_and_mixed_measurement_source_rejected(tmp_path):
    private,units,packet=data(tmp_path)
    folder=next(units.glob("terminal_p*_b*"))
    csv_path=folder/"original_roi_sample_terminal.csv"
    csv_path.write_text(csv_path.read_text().replace("classified_four_state_morph","modified_colour_outcome"))
    with pytest.raises(ValueError,match="hash"):
        m.collect(private,units,tmp_path/"result",synthetic=True)


def test_extra_batch_rejected(tmp_path):
    private,units,packet=data(tmp_path)
    (units/"terminal_p127_b9999999").mkdir()
    with pytest.raises(ValueError,match="unexpected"):
        m.collect(private,units,tmp_path/"result",synthetic=True)


def test_synthetic_terminal_cannot_pass_production_claim(tmp_path):
    private,units,packet=data(tmp_path)
    with pytest.raises(ValueError,match="method gate"):
        m.collect(private,units,tmp_path/"result",synthetic=False)
