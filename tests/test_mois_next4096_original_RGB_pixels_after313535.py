"""Fifth 4096 original RGB photo wave is not metadata-only colour evidence."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import pandas as pd
import pytest

from scripts.analysis import mois_next4096_original_photo_ids_after313535 as s
from scripts.analysis import mois_merge_fifth4096_original_ROI_after313535 as m
from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    HUES,METHOD,MEASURED_SCHEMA,
)

def source_row(n,tax,obs):
    import datetime
    date="2025-06-19"
    return {"photo_id":str(n),"inat_taxon_id":str(tax),
            "species":f"Campanula tax{tax}","observer_id":str(obs),
            "observation_id":str(obs),"observed_on":date,
            "latitude":"37.5","longitude":"137.5",
            "photo_license":"cc-by",
            "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{n}/large.jpg",
            "source_snapshot_date":"2026-09-27"}


def test_blind_new_source_gates(monkeypatch):
    monkeypatch.setattr(s,"SOURCE_N",20)
    monkeypatch.setattr(s,"N_NEW",4)
    monkeypatch.setattr(s,"TAXA",2)
    rows=[source_row(100+k,k%2,100+k) for k in range(20)]
    selected,count=s.choose(rows,set(),set())
    assert len(selected)==4
    assert len({x["observer_id"] for x in selected})==4
    assert len({x["inat_taxon_id"] for x in selected})==2
    assert not any("morph" in x for x in selected)
    assert count["provider_original_photo_ID_rows"]==20
    assert all(x["image_filename"].startswith("MOIS-") for x in selected)


def test_prior_pixel_and_targeted_reservation_must_not_be_readmitted(monkeypatch):
    monkeypatch.setattr(s,"SOURCE_N",20)
    monkeypatch.setattr(s,"N_NEW",4)
    monkeypatch.setattr(s,"TAXA",2)
    rows=[source_row(100+k,k%2,100+k) for k in range(20)]
    selected,count=s.choose(rows,{"100","101"},{"102","103"})
    assert count["photo_already_real_RGB_measured"]==2
    assert count["photo_reserved_for_season_target"]==2
    assert not {"100","101","102","103"}.intersection(
        {x["photo_id"] for x in selected})


def test_true_pixel_receipt_rejects_model_drift(tmp_path):
    selected=[]
    measured=[]
    for k in range(256):
        entry=source_row(10000+k,0,10000+k)
        entry["measurement_id"]=f"MOIS-{k:024d}"
        selected.append(entry)
        row={**entry,"morph":"white" if k%2 else "mixed_uncertain",
             "measurement_status":("classified_four_state_morph" if k%2
                                   else "unknown"),
             "classified_four_hue":"true" if k%2 else "false"}
        measured.append(row)
    terminal=tmp_path/"terminal.csv"
    pd.DataFrame(measured).to_csv(terminal,index=False)
    receipt=tmp_path/"receipt.json"
    freeze={"first_wave_original_4096_csv_sha256":"abc",
            "per_shard_sha256":{"0":"def"}}
    ok={"schema":MEASURED_SCHEMA,
        "shard":0,"terminal_rows":256,"original_source_ids":256,
        "source_sha256":"def","source_frozen_4096_selection_sha256":"abc",
        "source_model_attestation":METHOD,
        "pixels_ephemeral_no_original_image_bytes_persisted":True,
        "original_new_photo_hue_terminal_SHA256":m.file_sha(terminal),
        "four_hue_classified":128,"unclassified":128}
    receipt.write_text(json.dumps(ok))
    out,r=m.complete_worker(selected,receipt,terminal,freeze,0)
    assert len(out)==256
    bad={**ok,"source_model_attestation":{**METHOD,"location_blind":False}}
    receipt.write_text(json.dumps(bad))
    with pytest.raises(ValueError,match="not source-verified"):
        m.complete_worker(selected,receipt,terminal,freeze,0)


def test_real_images_only_not_reserved_metadata_ids():
    assert (s.N_PRIOR,s.N_CLASSIFIED,s.N_DATED,s.N_NEW)==(
        313535,151239,151230,4096)
    assert m.EXPECTED_TOTAL==317631
    assert len(HUES)==4 and "mixed_uncertain" not in HUES
