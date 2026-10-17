"""Original 400 MOIS plant-phenology annotation API coverage without photo-stage inference."""
from __future__ import annotations

import csv
import json
import uuid
from pathlib import Path

import pytest

from scripts.analysis.mois_observer_disjoint_photo_wave_20261017 import PHOTO_COLUMNS
from scripts.analysis.mois_audit_400_live_observation_phenology import (
    collect, original_400, sha, VALUE,
)


def fixtures(tmp_path):
    manifest=tmp_path/"source400.csv"
    receipt=tmp_path/"source_receipt.json"
    records=[]
    stages=["summer","autumn","summer","autumn"]
    names=["Alpha floral","Alpha floral","Beta floral","Beta floral"]
    for n in range(1,5):
        photo=str(900000000+n)
        orig_uuid=str(uuid.UUID(int=n))
        records.append({
            "inat_taxon_id":str(101 if n<=2 else 202),
            "species":names[n-1],
            "photo_id":photo,
            "observation_id":orig_uuid,
            "observer_id":"observer"+str(n),
            "observed_on":("2025-07-15" if stages[n-1]=="summer" else "2025-10-15"),
            "photo_license":"cc-by",
            "photo_creator_id":"creator"+str(n),
            "attribution":"synthetic-test-CCBY-"+str(n),
            "photo_url_large":
                f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo}/large.jpg",
            "source_snapshot_date":"2026-09-27",
            "site_latcell":"500",
            "site_loncell":"345",
            "season_year":"2025",
            "selected_season":stages[n-1],
            "source_rank":"hex"+str(n),
            "wave_photo_stage":"source-only pre-colour photo opportunity",
        })
    with manifest.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=PHOTO_COLUMNS,lineterminator="\n")
        w.writeheader();w.writerows(records)
    receipt.write_text(json.dumps({
        "status":"ORIGINAL_PHOTO_SOURCE_ONLY_NOT_A_BLOOM_OR_COLOUR",
        "selected_photo_manifest_sha256":sha(manifest),
        "did_not_infer_flower_opening_or_white_fraction":True,
        "did_not_open_or_reclassify_photo_pixels":True,
        "selected_original_photo_ids":4,
        "selected_distinct_source_observers":4,
        "selected_taxa":2,
    }))
    return manifest,receipt,records


def fake_provider(original, status_map=None, *,drop_v2=None,drop_v1=None,
                  drift_uuid=None,drift_photo=None,drift_date=None,
                  bad_uuid=False):
    by_uuid={r["observation_id"]:r for r in original}
    by_numeric={40000+i:r for i,r in enumerate(original,1)}
    u_to_id={r["observation_id"]:k for k,r in by_numeric.items()}
    stage=status_map or {1:13,2:14,3:21,4:None}
    def response(url,params=None):
        if "/v2/observations/" in url:
            uid=url.split("/v2/observations/")[-1].split(",")
            return {"results":[{"uuid":("bad-uuid" if bad_uuid else u),
                "id":u_to_id[u]} for u in uid if u!=drop_v2]}
        assert "/v1/observations/" in url
        ids=[int(x) for x in url.split("/v1/observations/")[-1].split(",")]
        result=[]
        for k in ids:
            if k==drop_v1:continue
            r=by_numeric[k]
            n=original.index(r)+1
            res={"uuid":drift_uuid if k==drift_uuid else r["observation_id"],
                 "id":k,
                 "observed_on":("2025-08-14" if k==drift_date else r["observed_on"]),
                 "photos":[{"id":("other" if k==drift_photo else r["photo_id"])}],
            }
            if stage[n] is not None:
                res["annotations"]=[{"controlled_attribute_id":12,
                                     "controlled_value_id":stage[n]}]
            result.append(res)
        return {"results":result}
    return response


def test_value_21_is_explicitly_no_evidence_at_observation_not_photo():
    assert VALUE[21]=="no_evidence_of_flowering"


def test_source_only_selection_uses_true_original_photo_and_independent_observers(tmp_path):
    source,receipt,records=fixtures(tmp_path)
    rows=original_400(source,receipt,_synthetic_test=True)
    assert len(rows)==4
    assert len({r["observation_id"] for r in rows})==4
    assert len({r["observer_id"] for r in rows})==4
    with pytest.raises(ValueError,match="Not exact five-taxon"):
        original_400(source,receipt,_synthetic_test=False)


def test_live_original_phenology_codes_and_missing_tags_do_not_make_flower_photos(tmp_path):
    source,receipt,original=fixtures(tmp_path)
    d=collect(source,receipt,tmp_path/"out",_synthetic_test=True,
              _request=fake_provider(original),pause=0)
    assert d["status"]=="CURRENT_PUBLIC_OBSERVATION_PHENOLOGY_EFFORT_ONLY_NOT_BOTANICAL_OPEN_FLOWERS"
    assert d["n_original_outcome_blind_selected_photo_opportunities"]==4
    assert d["n_live_current_v2_uuid_ids_resolved"]==4
    assert d["n_live_current_v1_full_annotation_records"]==4
    assert d["n_live_original_observed_on_and_photo_id_both_match"]==4
    assert d["n_source_taxon_siteyear_seasons"]==4
    assert d["n_independently_botanically_verified_open_photo_petals"]==0
    assert d["no_colour_ROI_photo_measurements_from_these_400"] is True
    assert d["no_new_white_flower_frequency_or_weather_effect_estimated"] is True
    with (tmp_path/"out"/"current_live_iNaturalist_400_taxon_siteyear_season_phenology_availability.csv").open(newline="") as f:
        groups=list(csv.DictReader(f))
    assert len(groups)==4
    assert sum(int(x["current_observation_flowering_NOT_photo"]) for x in groups)==1
    assert sum(int(x["current_observation_fruiting_NOT_photo"]) for x in groups)==1
    assert sum(int(x["current_observation_no_evidence_of_flowering_NOT_photo"]) for x in groups)==1
    assert sum(int(x["current_observation_phenology_missing_or_unknown"]) for x in groups)==1
    assert sum(int(x["original_picture_OPEN_petal_tepal_independently_reviewed"]) for x in groups)==0
    with pytest.raises(ValueError,match="overwrite"):
        collect(source,receipt,tmp_path/"out",_synthetic_test=True,
                _request=fake_provider(original),pause=0)


def test_missing_live_provider_observation_is_unknown_not_fruit(tmp_path):
    source,receipt,original=fixtures(tmp_path)
    d=collect(source,receipt,tmp_path/"out",_synthetic_test=True,
              _request=fake_provider(original,drop_v2=original[0]["observation_id"]),pause=0)
    assert d["n_original_outcome_blind_selected_photo_opportunities"]==4
    assert d["n_live_current_v2_uuid_ids_resolved"]==3
    assert d["n_live_original_observed_on_and_photo_id_both_match"]==3
    assert sum(x.get("no_current_public_observation_unknown",0) for x in
               d["current_observation_status_by_source_local_season"].values())==1


def test_observation_photo_disappears_from_live_API_but_original_source_photo_kept(tmp_path):
    source,receipt,original=fixtures(tmp_path)
    d=collect(source,receipt,tmp_path/"out",_synthetic_test=True,
              _request=fake_provider(original,drift_photo=40001),pause=0)
    assert d["n_live_original_observed_on_and_photo_id_both_match"]==3
    assert d["source_original_photo_ids_never_replaced_or_filtered"] is True


def test_drifted_live_observed_on_never_changes_original_month_sampling(tmp_path):
    source,receipt,original=fixtures(tmp_path)
    d=collect(source,receipt,tmp_path/"out",_synthetic_test=True,
              _request=fake_provider(original,drift_date=40002),pause=0)
    assert d["n_live_original_observed_on_and_photo_id_both_match"]==3
    assert sum(s.get("live_observation_day_changed_since_source_annotation_NOT_USABLE",0)
               for s in d["current_observation_status_by_source_local_season"].values())==1


def test_unexpected_provider_uuid_never_joins_foreign_plant(tmp_path):
    source,receipt,original=fixtures(tmp_path)
    with pytest.raises(ValueError,match="alien"):
        collect(source,receipt,tmp_path/"out",_synthetic_test=True,
                _request=fake_provider(original,bad_uuid=True),pause=0)


def test_original_source_edit_changes_SHA_and_blocks_any_API_requests(tmp_path):
    source,receipt,original=fixtures(tmp_path)
    source.write_text(source.read_text().replace("2025-07-15","2025-06-15",1))
    with pytest.raises(ValueError,match="SHA-certified"):
        collect(source,receipt,tmp_path/"out",_synthetic_test=True,
                _request=fake_provider(original),pause=0)


def test_exact_live_annotation_request_never_reads_image_colour_or_climate(tmp_path):
    source,receipt,original=fixtures(tmp_path)
    requests=[]
    mock=fake_provider(original)
    def track(url,params=None):
        requests.append((url,params))
        return mock(url,params)
    d=collect(source,receipt,tmp_path/"out",_synthetic_test=True,
              _request=track,pause=0,batch_size=2)
    assert len(requests)==4
    assert all("photos" not in str(req[1] or "") and "morph" not in str(req[1] or "")
               and "climate" not in str(req[1] or "") for req in requests)
    assert d["observation_level_annotation_does_not_confirm_open_petal_in_original_photo"] is True
