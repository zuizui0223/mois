"""Synthetic fake FCP per-species immutable source and site rainfall, no photos/pixels."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/acquisition/mois_sample_direct_from_fcp_source_20261016.py"
spec = importlib.util.spec_from_file_location("mois_direct", SOURCE)
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


class FakeBio12:
    provenance = {"synthetic": True}
    def value(self, lat, lon):
        if not lat and not lon:
            return ""
        if not lat or not lon:
            raise ValueError("Half-missing")
        la = float(lat)
        return "250" if la < 0 else ("1800" if la > 23 else "800")


def fixture(tmp):
    src = tmp / "fcp" / "complete_42111_metadata_source"
    src.mkdir(parents=True)
    counts = tmp / "complete_counts.csv"
    audit = tmp / "admission.json"
    rows = {
        1: [
            (11, "o1", -12, 10), (12, "o1", -12, 10),
            (13, "o2", 50, 125), (14, "o3", "", ""),
            (15, "o4", 25, -135),
        ],
        2: [(21, "o9", -15, 80), (22, "o10", -15, 80),
            (23, "o11", -15, 80), (24, "o12", -15, 80)],
        3: []
    }
    names = {1:"Alpha alba", 2:"Beta rosea", 3:"Gamma picta"}
    sha = {}
    for tax, items in rows.items():
        folder = src / "shard_000_of_128" / "species" / str(tax)
        folder.mkdir(parents=True)
        path = folder / "photo_manifest.ndjson"
        with path.open("w", encoding="utf-8") as f:
            for pid, oid, lat, lon in items:
                record = {
                    "inat_taxon_id":tax, "species":names[tax],
                    "observation_id":oid, "photo_id":str(pid),
                    "latitude":lat, "longitude":lon,
                    "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpeg",
                    "photo_license":"cc-by", "attribution":"synthetic test rights",
                    "observer_id":"u1", "observed_on":"2020-07-07",
                    "positional_accuracy_m":12, "snapshot_date":s.SNAPSHOT
                }
                f.write(json.dumps(record,sort_keys=True)+"\n")
        sha[tax] = s.sha256_file(path)
        (folder / "checkpoint.json").write_text(json.dumps({
            "status":"complete", "taxon_id":tax, "species":names[tax],
            "source_photo_rows":len(items), "metadata_sha256":sha[tax],
            "photo_cap":None,"colour_pixels_opened":False
        }))
    (src / "bulk_source_receipt.json").write_text(json.dumps({
        "schema":"fcp_open_data_bulk_complete_source_v1",
        "status":"COMPLETE_SNAPSHOT_METADATA_ALL_NOMINATED_TAXA",
        "source_snapshot_date":s.SNAPSHOT, "taxa_in_initial_frame":3,
        "n_target_taxon_licensed_photos":9, "number_independently_sealed_source_shards":128,
        "photo_cap":None, "old_photo_measurements_read":False,"image_pixels_opened":False,
    }))
    with counts.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","n_unique_target_photo_ids","source_status","source_manifest_sha256"])
        for tax in rows:
            w.writerow([tax,names[tax],len(rows[tax]),"complete",sha[tax]])
    audit.write_text(json.dumps({
        "schema":"fcp_complete_photo_ID_source_vs_independent_provider_audit_v1",
        "status":"REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS",
        "snapshot_date":s.SNAPSHOT, "fully_frozen_count_ledger_sha256":s.sha256_file(counts),
        "original_photo_IDs_covered_without_count_loss":True,
        "full_study_original_pixels_classified":False,
        "flower_colour_H2_W_or_environment_results":False,
        "species":3,"distinct_original_photo_ids":9
    }))
    return src,counts,audit


def test_exact_source_byte_and_observation_sampling(tmp_path):
    src, counts, audit = fixture(tmp_path)
    csv_out=tmp_path/"private.csv"
    receipt=tmp_path/"receipt.json"
    out=s.run(src,counts,audit,FakeBio12(),csv_out,receipt,
              calibration={"synthetic":True},synthetic=True)
    assert out["source_photo_ids_authenticated"]==9
    assert out["source_observations"]==8
    assert out["source_species"]==3
    assert out["source_species_with_photos"]==2
    assert out["selected_original_photo_ids"]==8
    assert out["pixels_opened_or_colours_measured"] is False
    records=list(csv.DictReader(csv_out.open()))
    assert len(records)==8
    assert all(x["photo_url_large"].startswith(s.PHOTO_URL_PREFIX) for x in records)
    assert {x["genus"] for x in records}=={"Alpha","Beta"}
    assert sum(x["precip_stratum"]=="P_UNKNOWN" for x in records)==1
    assert sum(float(x["photo_inclusion_probability"])==.5 for x in records)==1
    assert sum(int(x["source_photos_in_observation"])==2 for x in records)==1
    assert out["selected_private_photo_csv_sha256"]==s.sha256_file(csv_out)
    with pytest.raises(ValueError, match="Existing"):
        s.run(src,counts,audit,FakeBio12(),csv_out,receipt,
              calibration={},synthetic=True)


def test_modified_manifest_is_not_accepted(tmp_path):
    src,counts,audit=fixture(tmp_path)
    path=src/"shard_000_of_128/species/2/photo_manifest.ndjson"
    content=path.read_text()
    path.write_text(content.replace('"cc-by"','"cc0"'))
    target=tmp_path/"sample.csv"
    with pytest.raises(ValueError,match="SHA mismatch"):
        s.run(src,counts,audit,FakeBio12(),target,tmp_path/"receipt.json",
              calibration={},synthetic=True)
    assert not target.exists()
    assert not (tmp_path/"receipt.json").exists()
    assert not (tmp_path/"sample.csv.partial").exists()


def test_missing_even_zero_photo_manifest_aborts(tmp_path):
    src,counts,audit=fixture(tmp_path)
    (src/"shard_000_of_128/species/3/photo_manifest.ndjson").unlink()
    with pytest.raises(ValueError,match="Missing frozen"):
        s.run(src,counts,audit,FakeBio12(),tmp_path/"private.csv",
              tmp_path/"receipt.json",calibration={},synthetic=True)


def test_nonadmitted_counts_abort(tmp_path):
    src,counts,audit=fixture(tmp_path)
    a=json.loads(audit.read_text())
    a["status"]="VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT"
    audit.write_text(json.dumps(a))
    with pytest.raises(ValueError,match="count-only"):
        s.run(src,counts,audit,FakeBio12(),tmp_path/"private.csv",
              tmp_path/"receipt.json",calibration={},synthetic=True)


def test_invalid_site_calibration_fails_closed(tmp_path):
    cal=tmp_path/"sites.csv"
    with cal.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["latitude","longitude","bio12_mm"])
        for lat,lon in [(-30,-100),(-15,10),(10,10),(45,45)]*3:
            w.writerow([lat,lon,FakeBio12().value(str(lat),str(lon))])
    sha=s.sha256_file(cal)
    result=s.calibrate(FakeBio12(),cal,sha)
    assert result["calibration_sites_matching"]==12
    assert result["calibration_geographic_groups"]>=3
    cal.write_text(cal.read_text().replace(",250\n",",999\n",1))
    with pytest.raises(ValueError,match="SHA256 mismatch"):
        s.calibrate(FakeBio12(),cal,sha)
    bad_sha=s.sha256_file(cal)
    with pytest.raises(ValueError,match="differs"):
        s.calibrate(FakeBio12(),cal,bad_sha)


def test_huge_raster_not_needed_in_synthetic_ci():
    assert s.CAP==40
    assert s.PHOTO_URL_PREFIX.startswith("https://")
    assert s.SOURCE_SPECIES==42111
