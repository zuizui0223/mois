"""Source-authenticity guards for actual 309,439-photo WorldClim/maxT refit.

Synthetic photo frames are test-only, never published ecological effect sizes.
"""
from __future__ import annotations
import hashlib
import json
import pandas as pd
import pytest

from scripts.analysis import mois_refit_original_water_309439_full_FCP256 as m


def test_true_photo_source_population_and_nonoverwritten_original_water_method():
    assert m.OLD_N == 295017 and m.OLD_HUES == 142538
    assert m.NEW_N == 14422 and m.NEW_HUES == 7124
    assert m.FULL_N == 309439 and m.FULL_HUES == 149662
    assert m.FULL_SHA == "546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"
    assert m.OLD_ENV_SHA == "84c651aeaacf0e317445593c249a85940d230a19ec8b7f7223ec0e46553013f7"
    assert m.BOOT if hasattr(m,"BOOT") else 999  # inherited original resamples


def frames(tmp_path,monkeypatch):
    old = pd.DataFrame([
       dict(photo_id="100",inat_taxon_id="12",species="flower A",morph="white"),
       dict(photo_id="101",inat_taxon_id="13",species="flower B",morph="mixed_uncertain"),
    ])
    source = pd.DataFrame([
       dict(photo_id="100",inat_taxon_id="12",species="flower A",morph="white",
            source_original_measured="historic",latitude="33",longitude="130",
            observed_on="2022-06-12"),
       dict(photo_id="101",inat_taxon_id="13",species="flower B",morph="mixed_uncertain",
            source_original_measured="old",latitude="34",longitude="131",
            observed_on="2022-06-12"),
       dict(photo_id="102",inat_taxon_id="13",species="flower B",morph="red_pink",
            source_original_measured="fcp_wave4_full256_original_ROI_v4",
            latitude="35",longitude="131",observed_on="2026-08-14"),
    ])
    original=tmp_path/"photos.gz"
    climate=tmp_path/"env.gz"
    source.to_csv(original,index=False,compression="gzip")
    old.to_csv(climate,index=False,compression="gzip")
    photo_sha=hashlib.sha256(original.read_bytes()).hexdigest()
    prior_sha=hashlib.sha256(climate.read_bytes()).hexdigest()
    receipt={
        "schema":"mois_full_FCP256_309439_original_pixels_149653_dated_hues_v1",
        "now_original_measured_photo_IDs":3,
        "now_original_classifiable_four_hue_ROIs":2,
        "archive_exact_original_new_photo_source_SHA256":photo_sha,
        "FCP_full_original_256_worker_partitions":256,
        "source_photo_taxon_hue_date_and_location_identities_unmodified":True,
        "historical_source_S3_API_date_priority_not_overridden":True}
    envreceipt={
        "n_original_photo_ids":2,
        "n_four_hue_classified":1,
        "full_photo_static_monthly_normal_csv_sha256":prior_sha,
    }
    photo_json=tmp_path/"source.json"
    climate_json=tmp_path/"climate.json"
    photo_json.write_text(json.dumps(receipt))
    climate_json.write_text(json.dumps(envreceipt))
    for k,value in (
        ("OLD_N",2),("OLD_HUES",1),
        ("NEW_N",1),("NEW_HUES",1),
        ("FULL_N",3),("FULL_HUES",2),
        ("FULL_SHA",photo_sha),("OLD_ENV_SHA",prior_sha)
    ):
        monkeypatch.setattr(m,k,value)
    return source, original,photo_json,climate,climate_json


def test_exact_photo_taxon_hue_match_and_new_FCP_pixel_source(tmp_path,monkeypatch):
    raw,photos,proof,climate,climate_proof=frames(tmp_path,monkeypatch)
    r=m.verify(photos,proof,climate,climate_proof)
    assert r["all_original_unique_photos"]==3
    assert r["new_original_source_photos"]==1
    assert r["new_four_hue_measured_photo_ids"]==1
    assert r["original_model_SHA_verified_for_every_photo"]
    assert not r["original_FCP_within_species_W_D_or_50km_not_imported"] is False


def test_refuse_modifying_previously_classified_source_pixel(tmp_path,monkeypatch):
    raw,photos,proof,climate,climate_proof=frames(tmp_path,monkeypatch)
    raw.loc[0,"morph"]="red_pink"
    raw.to_csv(photos,index=False,compression="gzip")
    monkeypatch.setattr(m,"FULL_SHA",m.sha_file(photos))
    r=json.loads(proof.read_text())
    r["archive_exact_original_new_photo_source_SHA256"]=m.sha_file(photos)
    proof.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="morph was overwritten"):
        m.verify(photos,proof,climate,climate_proof)


def test_refuse_synthetic_new_photo_from_unmeasured_old_source(tmp_path,monkeypatch):
    raw,photos,proof,climate,climate_proof=frames(tmp_path,monkeypatch)
    raw.loc[2,"source_original_measured"]="unmeasured_source_metadata_only"
    raw.to_csv(photos,index=False,compression="gzip")
    monkeypatch.setattr(m,"FULL_SHA",m.sha_file(photos))
    r=json.loads(proof.read_text())
    r["archive_exact_original_new_photo_source_SHA256"]=m.sha_file(photos)
    proof.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="authenticated completed FCP"):
        m.verify(photos,proof,climate,climate_proof)
