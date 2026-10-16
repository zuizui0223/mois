#!/usr/bin/env python3
"""Fail-closed reassembly of all Mois SAMPLE ROI-v4 terminals into private traits.

This reconstructs the entire *selected* probabilistic photo sample, NOT the
full FCP 93M-photo corpus. It refuses any incomplete, duplicate, extra or
mixed-method batches, then joins the location-blind results to the protected
sampling key and calculates observed/classifiable photo-white proportions.
It does not extract the other 12 environments or run the 13-variable models.
No taxon IDs, source photo URLs or sensitive coordinates enter public receipts.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from scripts.acquisition.mois_stratified_photo_selection_20261016 import sha256_file
from scripts.acquisition.mois_build_sample_roi_handoff_20261016 import (
    N_PARTITIONS, WORKER, ACQUISITION, JOIN,
)
from scripts.acquisition.mois_execute_sample_roi_v4_20261016 import (
    STATUS, BATCH_CAP, FCP_BACKEND_GIT_BLOB, read_packet,
)

MORPHS=("white","yellow_orange","red_pink","blue_purple")
BIO=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
TERMINAL_COLUMNS=("measurement_id","morph","measurement_status")
PRIVATE_JOIN_EXTRA=("morph","measurement_status","nine_colour_fractions_json")
SUMMARY_COLUMNS=(
    "inat_taxon_id","genus","selected_observations",
    "classified_four_colour","classified_white","unclassifiable",
    "photo_weighted_classifiability","observation_weighted_classifiability",
    "white_given_classified_photo_weighted",
    "white_given_classified_observation_weighted",
)
SEAL_STATUS="COMPLETE_MOIS_SAMPLED_PHOTO_TRAITS_CONDITIONAL_CLASSIFIABILITY_NO_13ENV_FIT"
DIR_RE=re.compile(r"^terminal_p\d{3}_b\d{7}$")


def certified_packet(root:Path)->dict:
    receipt=json.loads((root/"mois_roi_packet_receipt.json").read_text(encoding="utf-8"))
    if (receipt.get("schema")!="mois_sample_scoped_original_roi_v4_blind_packets_v1"
        or receipt.get("status")!="SAMPLE_ONLY_BLIND_METADATA_SEALED_PIXELS_NOT_MEASURED"
        or receipt.get("only_mois_sample_not_all_fcp_source_photos") is not True
        or receipt.get("historical_colour_outcomes_reused") is not False
        or receipt.get("original_pixels_measured") is not False
        or receipt.get("worker_partitions")!=N_PARTITIONS):
        raise ValueError("Cannot assemble unsealed, historical or full-census FCP worker")
    files=receipt.get("file_sha256",{})
    for name, expected_sha in files.items():
        if sha256_file(root/name)!=expected_sha:
            raise ValueError("Original selected-photo sample packet changed after sealing")
    if ("private_metadata_join.csv" not in files or
        sum(receipt["partition_photo_counts"])!=receipt["selected_photo_ids"]):
        raise ValueError("Missing private identity join or invalid sample count")
    return receipt


def validate_batch(root:Path, folder:Path, part:int, batch:int,
                   expected:list[str], packet_hash:str,
                   packet:dict, batch_size:int, synthetic:bool)->tuple[list[dict],dict]:
    result_file=folder/"original_roi_sample_terminal.csv"
    receipt_file=folder/"work_unit_receipt.json"
    unit=json.loads(receipt_file.read_text(encoding="utf-8"))
    if (unit.get("schema")!="mois_sample_scoped_original_roi_v4_work_unit_v1"
        or unit.get("status")!=STATUS
        or unit.get("partition")!=part or unit.get("batch")!=batch
        or unit.get("compute_batch_size_not_biological_inclusion_cap")!=batch_size
        or unit.get("batch_photo_count")!=len(expected)
        or unit.get("full_sample_packets_sha256")!=packet_hash
        or unit.get("source_packet_worker_sha256")!=packet["file_sha256"][f"worker_partition_{part:03}.csv"]
        or unit.get("source_packet_acquisition_sha256")!=packet["file_sha256"][f"acquisition_partition_{part:03}.csv"]
        or unit.get("exact_selected_measurement_ids_sha256")!=hashlib.sha256(
            "".join(v+"\n" for v in expected).encode()).hexdigest()
        or unit.get("sealed_terminal_csv_sha256")!=sha256_file(result_file)
        or unit.get("fcp_full_census_completed_by_this_batch") is not False
        or unit.get("mois_sample_complete") is not False
        or unit.get("original_pixels_persisted") is not False
        or unit.get("historical_colour_outcomes_reused") is not False
        or unit.get("test_only_synthetic") is not synthetic
        or unit.get("original_fcp_roi_wrapper_git_blob")!=(
            "SYNTHETIC_TEST_ONLY" if synthetic else FCP_BACKEND_GIT_BLOB)):
        raise ValueError("ROI batch provenance, hash, sample scope or method gate failed")
    rows=[]
    with result_file.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        cols=set(reader.fieldnames or ())
        if not set(TERMINAL_COLUMNS).issubset(cols):
            raise ValueError("Missing original ROI response fields")
        if any(k in cols for k in ("species","inat_taxon_id","photo_id","latitude",
                                    "longitude","bio12_mm","photo_url_large")):
            raise ValueError("Location/colour outcomes and original metadata leaked into blind image evaluator")
        if not all("flower_fraction_"+c in cols for c in BIO):
            raise ValueError("Incomplete fixed nine-colour palette")
        for row in reader:
            if row["morph"] not in MORPHS+("mixed_uncertain",):
                raise ValueError("Unfrozen floral classification category")
            classified=row["measurement_status"]=="classified_four_state_morph"
            if classified != (row["morph"] in MORPHS):
                raise ValueError("Photographic four-state gate inconsistent with original ROI")
            fractions=[float(row["flower_fraction_"+c]) for c in BIO]
            if classified and (any(not math.isfinite(v) or v < -1e-9 for v in fractions)
                               or abs(sum(fractions)-1)>1e-7):
                raise ValueError("Nine-palette classified photo fractions invalid")
            row["nine_colour_fractions_json"]=json.dumps(fractions,separators=(",",":"))
            rows.append(row)
    if (len(rows)!=len(expected) or
        set(r["measurement_id"] for r in rows)!=set(expected) or
        len({r["measurement_id"] for r in rows})!=len(expected) or
        sum(r["measurement_status"]=="classified_four_state_morph" for r in rows)!=unit["photo_ids_classified"] or
        sum(r["measurement_status"]=="image_acquisition_failed" for r in rows)!=unit["photo_ids_acquisition_failed"]):
        raise ValueError("Incomplete/duplicated original ROI sample batch or counts")
    if unit["photo_ids_classified"]+unit["photo_ids_acquisition_failed"]+unit["unclassifiable_other"]!=len(expected):
        raise ValueError("Photographic terminal failure categories do not sum to actual source photos")
    return rows,unit


def collect(packet_root:Path, sealed_root:Path, destination:Path, *,
            batch_size:int=BATCH_CAP,synthetic:bool=False)->dict:
    if not 1<=batch_size<=BATCH_CAP:
        raise ValueError("Only bounded original ROI-v4 batch sizes 1-64")
    if destination.exists() or destination.with_name(destination.name+".INCOMPLETE").exists():
        raise ValueError("A previous final or incomplete private sample assembly exists")
    packet=certified_packet(packet_root)
    packet_sha=sha256_file(packet_root/"mois_roi_packet_receipt.json")
    expected_folders=set()
    stage=destination.with_name(destination.name+".INCOMPLETE")
    stage.mkdir(parents=True)
    database=stage/"_private_measured_sample_join.sqlite"
    connected=sqlite3.connect(database)
    connected.execute("CREATE TABLE colour (measurement_id TEXT PRIMARY KEY, morph TEXT, status TEXT, palette_json TEXT)")
    n_read=classified=acquisition_failed=other=0
    single_method=None
    try:
        for part in range(N_PARTITIONS):
            n=packet["partition_photo_counts"][part]
            if n==0:
                continue
            # Verify the WHOLE immutable source worker and its matching
            # acquisition file once through the original sample receipt.
            with (packet_root/f"worker_partition_{part:03}.csv").open(newline="",encoding="utf-8") as f:
                worker=list(csv.DictReader(f))
            if len(worker)!=n or list(worker[0])!=WORKER:
                raise ValueError("Source partition worker lost before terminal collection")
            selected_ids=[v["measurement_id"] for v in worker]
            for batch in range((n+batch_size-1)//batch_size):
                expected=selected_ids[batch*batch_size:(batch+1)*batch_size]
                stem=f"terminal_p{part:03}_b{batch:07}"
                expected_folders.add(stem)
                folder=sealed_root/stem
                if not folder.is_dir():
                    raise ValueError(f"Incomplete sample ROI batch: {stem}")
                rows,unit=validate_batch(packet_root,folder,part,batch,expected,
                                         packet_sha,packet,batch_size,synthetic)
                method=json.dumps(unit["method"],sort_keys=True)
                if single_method is None:
                    single_method=method
                elif method!=single_method:
                    raise ValueError("Mixed original ROI-v4 detector or palette models across sample work units")
                connected.executemany("INSERT INTO colour VALUES (?,?,?,?)",[
                    (r["measurement_id"],r["morph"],r["measurement_status"],
                     r["nine_colour_fractions_json"]) for r in rows])
                n_read+=len(rows)
                classified+=unit["photo_ids_classified"]
                acquisition_failed+=unit["photo_ids_acquisition_failed"]
                other+=unit["unclassifiable_other"]
        actual={p.name for p in sealed_root.iterdir()} if sealed_root.exists() else set()
        if actual!=expected_folders:
            raise ValueError("Unsealed, unexpected or duplicate sample ROI terminal directory")
        if n_read!=packet["selected_photo_ids"]:
            raise ValueError("Not all selected photographs have exact ROI measurements")
        connected.commit()
        stats=defaultdict(lambda:{"n":0,"class":0,"white":0,
                                 "allpw":0.,"classpw":0.,"whitepw":0.,
                                 "allow":0.,"classow":0.,"whiteow":0.})
        joined_file=stage/"private_joined_sample_flower_colours.csv"
        with (packet_root/"private_metadata_join.csv").open(newline="",encoding="utf-8") as inp, \
             joined_file.open("w",newline="",encoding="utf-8") as outf:
            reader=csv.DictReader(inp)
            if reader.fieldnames!=JOIN:
                raise ValueError("Source private metadata join drifted")
            writer=csv.DictWriter(outf,fieldnames=JOIN+list(PRIVATE_JOIN_EXTRA))
            writer.writeheader()
            seen=0
            for row in reader:
                ident=row["measurement_id"]
                found=connected.execute("SELECT morph,status,palette_json FROM colour WHERE measurement_id=?",(ident,)).fetchone()
                if found is None:
                    raise ValueError("Location-blind ROI result has no exact selected source ID")
                morph,status,palette=found
                writer.writerow({**row,"morph":morph,"measurement_status":status,
                                 "nine_colour_fractions_json":palette})
                tax=row["inat_taxon_id"]
                d=stats[(tax,row["genus"])]
                photo_w=float(row["photo_design_weight"])
                observation_pi=float(row["observation_inclusion_probability"])
                if (not math.isfinite(photo_w) or photo_w < 1 or not math.isfinite(observation_pi)
                    or not 0 < observation_pi <= 1):
                    raise ValueError("Invalid original-source photo sampling design weights")
                obs_w=1/observation_pi
                d["n"]+=1
                d["allpw"]+=photo_w;d["allow"]+=obs_w
                valid=status=="classified_four_state_morph"
                white=morph=="white"
                if valid:
                    d["class"]+=1;d["classpw"]+=photo_w;d["classow"]+=obs_w
                    if white:
                        d["white"]+=1;d["whitepw"]+=photo_w;d["whiteow"]+=obs_w
                seen+=1
            if seen!=packet["selected_photo_ids"]:
                raise ValueError("Source private photo join count does not equal frozen selection")
        summary_file=stage/"private_species_conditional_photo_white.csv"
        with summary_file.open("w",newline="",encoding="utf-8") as outf:
            writer=csv.DictWriter(outf,fieldnames=SUMMARY_COLUMNS)
            writer.writeheader()
            for (tax,genus),d in sorted(stats.items(),key=lambda z:int(z[0][0])):
                writer.writerow({
                    "inat_taxon_id":tax,"genus":genus,
                    "selected_observations":d["n"],
                    "classified_four_colour":d["class"],
                    "classified_white":d["white"],
                    "unclassifiable":d["n"]-d["class"],
                    "photo_weighted_classifiability":d["classpw"]/d["allpw"],
                    "observation_weighted_classifiability":d["classow"]/d["allow"],
                    "white_given_classified_photo_weighted":(
                        d["whitepw"]/d["classpw"] if d["classpw"] else ""),
                    "white_given_classified_observation_weighted":(
                        d["whiteow"]/d["classow"] if d["classow"] else ""),
                })
        # The joined individual-photo and species file remains RESTRICTED;
        # only global source denominators and hashes may leave its storage.
        report={
            "schema":"mois_full_selected_photo_roi_terminal_collection_v1",
            "status":SEAL_STATUS,
            "source_sample_packet_sha256":packet_sha,
            "source_photo_ids_selected":packet["selected_photo_ids"],
            "source_selected_taxa_with_photos":packet["selected_original_taxa_with_photo"],
            "original_photo_terminals_verified":n_read,
            "original_colour_four_state_classifiable":classified,
            "original_image_acquisition_permanent_failure":acquisition_failed,
            "original_other_unclassifiable":other,
            "species_with_at_least_one_classifiable_photo":sum(d["class"]>0 for d in stats.values()),
            "original_method_sha256_or_attestation":hashlib.sha256((single_method or "").encode()).hexdigest(),
            "private_joined_photo_csv_sha256":sha256_file(joined_file),
            "private_species_csv_sha256":sha256_file(summary_file),
            "full_FCP_93m_photos_measured":False,
            "all13_photo_site_predictors_joined":False,
            "global_colour_climate_model_fit":False,
            "genetic_white_frequency_or_fitness_estimated":False,
            "synthetic_only_test":bool(synthetic),
            "note":"Observed four-colour classifiability-conditional photographic summaries, not all-flora white prevalence or causal selection."
        }
        (stage/"sample_collection_receipt.json").write_text(
            json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        connected.close()
        database.unlink()
        stage.replace(destination)
        return report
    finally:
        try:connected.close()
        except Exception:pass


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-locked-packets",type=Path,required=True)
    p.add_argument("--sealed-roi-units",type=Path,required=True)
    p.add_argument("--out-private-traits",type=Path,required=True)
    p.add_argument("--batch-size",type=int,default=BATCH_CAP)
    a=p.parse_args()
    result=collect(a.sample_locked_packets,a.sealed_roi_units,a.out_private_traits,batch_size=a.batch_size)
    print(json.dumps({k:result[k] for k in (
        "status","source_photo_ids_selected","original_photo_terminals_verified",
        "original_colour_four_state_classifiable",
        "species_with_at_least_one_classifiable_photo"
    )},sort_keys=True))


if __name__=="__main__":
    main()
