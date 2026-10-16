"""Sample-only FCP ROI-v4 blind packet handoff; entirely synthetic and offline."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parents[1]
PATH = BASE / "scripts/acquisition/mois_build_sample_roi_handoff_20261016.py"
spec=importlib.util.spec_from_file_location("mois_roi_packets",PATH)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def setup_sample(root:Path):
    sample=root/"selected.csv"
    items=[]
    for tax,genus,photo,obs,pilot in ((1,"Alpha",11,"o1",1),
                                       (1,"Alpha",12,"o2",1),
                                       (2,"Beta",21,"o9",0)):
        row=dict.fromkeys(m.PRIVATE_FIELDS,"")
        row.update({
            "inat_taxon_id":str(tax), "genus":genus,"photo_id":str(photo),
            "observation_id":obs,"photo_license":"cc-by",
            "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo}/large.jpg",
            "latitude":"10","longitude":"40", "bio12_mm":"500",
            "available_observations_in_stratum":"2" if tax==1 else "1",
            "sampled_observations_in_stratum":"2" if tax==1 else "1",
            "source_photos_in_observation":"1",
            "observation_inclusion_probability":"1",
            "photo_inclusion_probability":"1",
            "photo_design_weight":"1",
            "technical_pilot_flag":str(pilot),
            "region_stratum":"G_T1","precip_stratum":"P_MEDIUM",
        })
        items.append(row)
    with sample.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=m.PRIVATE_FIELDS);w.writeheader();w.writerows(items)
    receipt=root/"receipt.json"
    receipt.write_text(json.dumps({
        "schema":"mois_direct_fcp_source_stratified_photo_selection_v1",
        "status":"SOURCE_ADMITTED_SELECTED_PRIVATE_IDS_ONLY_NO_PHOTO_PIXELS_MEASURED",
        "selected_private_photo_csv_sha256":m.sha256_file(sample),
        "historical_colour_outcomes_seen":False,
        "pixels_opened_or_colours_measured":False,
        "original_13var_environment_join_complete":False,
        "selection_photo_cap_per_taxon":40,
        "source_species":3,
        "source_photo_ids_authenticated":9,
        "selected_original_photo_ids":3,
        "technical_pilot_selected_photo_ids":2,
    }))
    return sample,receipt,items


def test_sample_only_roi_worker_packets_are_blind_and_sha_sealed(tmp_path):
    sample,receipt,_=setup_sample(tmp_path)
    private=tmp_path/"roi_packets"
    out=m.assemble(sample,receipt,private,synthetic=True)
    assert out["status"]=="SAMPLE_ONLY_BLIND_METADATA_SEALED_PIXELS_NOT_MEASURED"
    assert out["selected_photo_ids"]==3
    assert out["technical_pilot_photo_ids"]==2
    assert out["original_pixels_measured"] is False
    assert sum(out["partition_photo_counts"])==3
    assert out["do_not_pass_to_fcp_full_census_handoff"] is True
    for i in range(128):
        worker=list(csv.DictReader((private/f"worker_partition_{i:03}.csv").open()))
        acq=list(csv.DictReader((private/f"acquisition_partition_{i:03}.csv").open()))
        assert len(worker)==len(acq)
        for w in worker:
            assert set(w)==set(m.WORKER)
            assert w["measurement_id"].startswith("FCPR-")
        for a in acq:
            assert set(a)==set(m.ACQUISITION)
    metadata=list(csv.DictReader((private/"private_metadata_join.csv").open()))
    assert len(metadata)==3
    assert set(metadata[0])==set(m.JOIN)
    for z in metadata:
        assert z["measurement_id"]==m.blind_id(z["photo_id"])
    with pytest.raises(ValueError,match="Refuse overwriting"):
        m.assemble(sample,receipt,private,synthetic=True)


def test_no_impersonating_an_unsampled_93m_source(tmp_path):
    sample,receipt,_=setup_sample(tmp_path)
    with pytest.raises(ValueError,match="full original"):
        m.assemble(sample,receipt,tmp_path/"packets")


def test_changed_photo_probability_fails_closed(tmp_path):
    sample,receipt,items=setup_sample(tmp_path)
    items[0]["photo_inclusion_probability"]="0.75"
    with sample.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=m.PRIVATE_FIELDS);writer.writeheader();writer.writerows(items)
    doc=json.loads(receipt.read_text());doc["selected_private_photo_csv_sha256"]=m.sha256_file(sample)
    receipt.write_text(json.dumps(doc))
    destination=tmp_path/"private_roi"
    with pytest.raises(ValueError,match="probabilities"):
        m.assemble(sample,receipt,destination,synthetic=True)
    assert not destination.exists()
    assert (tmp_path/"private_roi.incomplete").is_dir()


def test_changed_source_image_url_rejected(tmp_path):
    sample,receipt,items=setup_sample(tmp_path)
    items[1]["photo_url_large"]="https://evil.invalid/12.jpg"
    with sample.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=m.PRIVATE_FIELDS);writer.writeheader();writer.writerows(items)
    doc=json.loads(receipt.read_text());doc["selected_private_photo_csv_sha256"]=m.sha256_file(sample)
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="photo URL"):
        m.assemble(sample,receipt,tmp_path/"private_roi",synthetic=True)


def test_count_only_forged_receipt_rejected(tmp_path):
    sample,receipt,_=setup_sample(tmp_path)
    a=json.loads(receipt.read_text());a["status"]="BUDGET_ONLY_NOT_PHOTO_ID_SELECTION_OR_COLOUR"
    receipt.write_text(json.dumps(a))
    with pytest.raises(ValueError,match="No genuine"):
        m.assemble(sample,receipt,tmp_path/"private_roi",synthetic=True)
