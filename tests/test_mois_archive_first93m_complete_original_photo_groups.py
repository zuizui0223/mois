"""Synthetic contract checks for permanent PARTIAL source-original 93m photo checkpoints."""
import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.acquisition import mois_archive_first93m_complete_original_photo_groups as mod


def fixture(tmp_path):
    queue=tmp_path/"queue"
    groups=tmp_path/"groups"
    (queue/"shards").mkdir(parents=True)
    (groups/"group00").mkdir(parents=True)
    checksum={str(i):"unused_after_initial_groups" for i in range(73)}
    hue=("white","red_pink","yellow_orange","blue_purple")
    photoids=[]
    for k in range(4):
        n=256
        source=[]
        for i in range(n):
            pid=str(500000+k*256+i)
            photoids.append(pid)
            source.append({
                "measurement_id":f"M{k}-{i}",
                "photo_id":pid,"inat_taxon_id":str(i%3+1),
                "species":"Campanula punctata",
                "observation_id":f"O{k}-{i//2}",
                "observer_id":"77","observed_on":"2025-05-12",
                "latitude":"35","longitude":"137",
                "positional_accuracy":"100","photo_license":"cc-by",
                "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                "image_filename":f"M{k}-{i}.jpg",
                "source_arm":"global93m_probability_sample_observations",
                "selection_probability":"0.0002",
                "source_snapshot_date":"2026-09-27",
                "shard":str(k)
            })
        src=queue/"shards"/f"shard-{k:03d}.csv"
        pd.DataFrame(source).to_csv(src,index=False)
        checksum[str(k)]=mod.digest(src)
        raw=groups/"group00"/f"shard_{k:03d}"
        (raw/"sealed").mkdir(parents=True)
        photo=pd.DataFrame(source)
        photo["morph"]=[hue[(i+k)%4] for i in range(n)]
        photo["measurement_status"]="classified_four_state_morph"
        photo["classified_four_hue"]="true"
        photo["image_sha256"]="a"*64
        if k==1:
            photo.loc[0,"morph"]="mixed_uncertain"
            photo.loc[0,"classified_four_hue"]="false"
            photo.loc[0,"measurement_status"]="image_acquisition_failed"
            photo.loc[0,"image_sha256"]=""
        if k==2:
            photo.loc[0,"morph"]="mixed_uncertain"
            photo.loc[0,"classified_four_hue"]="false"
            photo.loc[0,"measurement_status"]="not_evaluable_roi_or_flip_gate"
        output=raw/"photos.csv"
        photo.to_csv(output,index=False)
        terminal=raw/"sealed"/f"terminal_p{k:03d}_b0000000.csv"
        photo[["measurement_id","morph","measurement_status"]].to_csv(
            terminal,index=False)
        rec={
            "schema":"mois_full_original93m_real_RGB_photo_unit_v1",
            "shard":k,"attempted":n,"downloaded":n-(k==1),
            "classifiable":n-(k in (1,2)),
            "unclassifiable_acquired":int(k==2),
            "acquisition_failed":int(k==1),
            "source_sha256":checksum[str(k)],
            "photos_sha256":mod.digest(output),
            "terminal_sha256":mod.digest(terminal),
            "method":mod.EXPECTED_FROZEN_FCP,
            "photo_colours_measured":True,
            "raw_licensed_image_bytes_persisted":False,
        }
        (raw/"receipt.json").write_text(json.dumps(rec))
    plan={
        "schema":"mois_full_source_observation_sampling_cohort_v1",
        "full_frame_photo_ids":93209351,
        "full_source_manifest_sha256":mod.FIRST_INDEX_SHA,
        "selection_probability":0.0002,
        "selection_unit":"observation_uuid",
        "sampled_before_any_hue_measurement":True,
        "source_original_images_measured":False,
        "colour_measured":False,
        "sampled_photo_ids_before_prior_dedup":18533,
        "new_photo_attempts_planned":18462,
        "previously_measured_photo_ids_skipped":71,
        "selected_new_photos_unsupported_by_frozen_loader":0,
        "work_units":73,
        "units_per_group":4,
        "photo_shard_sha256":checksum,
        "no_per_taxon_or_observation_photo_count_cap":True,
    }
    (queue/"manifest.json").write_text(json.dumps(plan))
    return queue,groups


def test_full_completed_group_photo_terminal_shas_and_zero_duplicate_ids(tmp_path):
    src,worker=fixture(tmp_path)
    photos,receipt=mod.validate(src,worker,first=0,stop=4)
    assert len(photos)==1024
    assert len(receipt["original_unit_receipts"])==4
    assert receipt["genuine_original_photo_attempts_in_partial"]==1024
    assert receipt["classified_four_hues"]==1022
    assert receipt["downloaded_original_images"]==1023
    assert receipt["technical_unknown_after_download"]==1
    assert receipt["verified_original_photo_acquisition_failures"]==1
    assert photos.photo_id.is_unique
    assert sum(receipt["hue_counts_classified"].values())==1022
    assert receipt["completion_of_remaining_14366_new_photos_NOT_claimed"] is True
    assert receipt["H3_native_island_or_any_abiotic_ecological_model_performed"] is False


def test_missing_original_outcome_csv_refuses_partial_source_release(tmp_path):
    source,groups=fixture(tmp_path)
    (groups/"group00"/"shard_002"/"photos.csv").unlink()
    with pytest.raises(ValueError,match="terminal missing"):
        mod.validate(source,groups,first=0,stop=4)


def test_original_photo_identity_cannot_be_swapped_without_detection(tmp_path):
    source,groups=fixture(tmp_path)
    p=groups/"group00"/"shard_002"/"photos.csv"
    frame=pd.read_csv(p,dtype=str)
    frame.loc[0,"photo_id"]="777777777"
    frame.to_csv(p,index=False)
    r=groups/"group00"/"shard_002"/"receipt.json"
    info=json.loads(r.read_text());info["photos_sha256"]=mod.digest(p)
    r.write_text(json.dumps(info))
    with pytest.raises(ValueError,match="photo ID is missing/changed"):
        mod.validate(source,groups,first=0,stop=4)


def test_unsiged_frozen_FCP_fingerprint_drift_is_rejected(tmp_path):
    source,groups=fixture(tmp_path)
    receipt=groups/"group00"/"shard_000"/"receipt.json"
    r=json.loads(receipt.read_text())
    r["method"]["detector_sha256"]="0"*64
    receipt.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="fixed model source digest"):
        mod.validate(source,groups,first=0,stop=4)


def test_write_checkpoint_parquet_and_immutable_shas(tmp_path):
    source,groups=fixture(tmp_path)
    # Writer's production default requires 16 complete units. Four synthetic
    # units may be tested by passing the frozen group boundary explicitly.
    output=tmp_path/"sealed"
    result=mod.run(source,groups,output,first=0,stop=4)
    assert pd.read_parquet(output/mod.OUT_PHOTOS).photo_id.is_unique
    assert sum(result["hue_counts_classified"].values())==1022
    assert mod.digest(output/mod.OUT_PHOTOS)==result["output_SHA256"][mod.OUT_PHOTOS]
    with pytest.raises(ValueError,match="overwrite"):
        mod.run(source,groups,output,first=0,stop=4)
