#!/usr/bin/env python3
"""Resumable, source-SHA-verifying original dated-photo source SHARD executor.

Runs only a small explicit consecutive subset of original 128 disjoint FCP
SOURCE shards (not a 40-photo sample). For already completed shards, validate
their protected queue and receipt source hashes; do not duplicate source work.
Unsealed/partial outputs are quarantined by refusal, NOT silently overwritten.

Does not open image pixels, estimate Flowers annotations or white frequencies.
A successful range result is explicitly PARTIAL, never a completed global
93,209,351-photo day/phenology summary; separate full all-128 aggregation
is required to claim complete original photo observed_on availability.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
    audit_shard,sha,SCHEMA as SRC_SCHEMA,STATUS as SRC_STATUS,
    FRAME_SHA,COUNTS_SHA,FULL_WORK_PARTITIONS
)

SCHEMA="mois_original_fcp_full_photo_bounded_date_shard_progress_v1"
STATUS="SOURCE_SEALED_PARTIAL_DATE_SHARDS_NOT_WORLDWIDE_FLOWER_RESULT"


def run(shards:Path,counts:Path,private_root:Path,progress:Path,*,
        start:int,stop:int,partitions:int=FULL_WORK_PARTITIONS,
        synthetic:bool=False)->dict:
    if (start<0 or stop<=start or stop>partitions or stop-start>4
        or (not synthetic and (partitions!=FULL_WORK_PARTITIONS
             or sha(counts)!=COUNTS_SHA))):
        raise ValueError("Expected only up to four original source photo shards, not a biological photo cap")
    if progress.exists() or progress.with_name(progress.name+".INCOMPLETE").exists():
        raise ValueError("Do not overwrite previous independent real source range progress receipt")
    checked=[]
    for part in range(start,stop):
        source=shards/f"shard_{part:03}_of_{partitions:03}"
        original_receipt=source/"shard_source_receipt.json"
        private=private_root/f"dated_candidates_shard_{part:03}.ndjson"
        public=private_root/f"dated_candidates_shard_{part:03}_receipt.json"
        if not original_receipt.is_file():
            raise ValueError(f"Original FCP provider photo-ID shard {part} not yet source sealed")
        sealed=private.is_file() and public.is_file()
        if sealed:
            doc=json.loads(public.read_text(encoding="utf8"))
            if (doc.get("schema")!=SRC_SCHEMA or doc.get("status")!=SRC_STATUS
                or doc.get("shard_number")!=part
                or doc.get("source_shard_receipt_sha256")!=sha(original_receipt)
                or doc.get("source_original_counts_ledger_sha256")!=sha(counts)
                or doc.get("private_full_dated_source_photo_queue_sha256")!=sha(private)
                or doc.get("n_fcp_shards_expected_for_global_completion")!=partitions
                or doc.get("source_photo_cap_per_taxon","MISSING") is not None
                or doc.get("metadata_date_known_is_not_blooming_photo_confirmation") is not True
                or doc.get("photo_hue_classification_attempted") is not False
                or doc.get("synthetic_test_only") is not synthetic
                or not synthetic and doc.get("original_source_frame_sha256")!=FRAME_SHA):
                raise ValueError("Pre-existing date photo shard no longer authenticated")
            state="PREVIOUSLY_COMPLETED_REVERIFIED"
        else:
            if (private.exists() or public.exists()
                or private.with_name(private.name+".INCOMPLETE").exists()
                or public.with_name(public.name+".INCOMPLETE").exists()):
                raise ValueError("Partially written original dated source must be quarantined; no overwrite")
            doc=audit_shard(source,counts,private,public,synthetic=synthetic)
            state="NEW_SOURCE_DATES_VALIDATED"
        checked.append({
            "partition":part,
            "state":state,
            "original_source_shard_receipt_sha256":sha(original_receipt),
            "source_dated_photo_shard_receipt_sha256":sha(public),
            "n_original_eligible_source_photographs_examined":doc[
                "n_original_eligible_source_photographs_examined"],
            "n_original_snapshot_date_usable_photo_candidates":doc[
                "n_original_photo_observation_dates_usable"],
            "n_real_photo_bloom_image_ROI_confirmed":0,
        })
    result={
        "schema":SCHEMA,"status":STATUS,
        "source_shard_start":start,"source_shard_stop_exclusive":stop,
        "source_universe_expected_shards":partitions,
        "n_shards_verified_in_this_execution":len(checked),
        "photo_opportunities_in_this_partial_source_block":sum(
            r["n_original_eligible_source_photographs_examined"] for r in checked),
        "dated_candidates_in_this_partial_source_block":sum(
            r["n_original_snapshot_date_usable_photo_candidates"] for r in checked),
        "source_photo_sampling_cap_per_taxon":None,
        "photo_image_pixels_opened":False,
        "flowering_photos_counted_or_inferred":False,
        "worldwide_complete_42111_taxon_93million_date_census":False,
        "per_original_shard":[*checked],
        "synthetic_test_only":bool(synthetic),
        "note":"Only the requested complete-source shard range has been checked. These are dated plant photo opportunities, not actual flower images; a complete 128-shard audit and subsequent actual photo ROI results are essential."
    }
    progress.parent.mkdir(parents=True,exist_ok=True)
    stage=progress.with_name(progress.name+".INCOMPLETE")
    stage.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    stage.replace(progress)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fcp-source-shards-root",required=True,type=Path)
    p.add_argument("--original-42111-count-ledger",required=True,type=Path)
    p.add_argument("--private-dated-photos-root",required=True,type=Path)
    p.add_argument("--start-shard",required=True,type=int)
    p.add_argument("--stop-shard-exclusive",required=True,type=int)
    p.add_argument("--out-progress-receipt",required=True,type=Path)
    a=p.parse_args()
    d=run(a.fcp_source_shards_root,a.original_42111_count_ledger,
          a.private_dated_photos_root,a.out_progress_receipt,
          start=a.start_shard,stop=a.stop_shard_exclusive)
    print(json.dumps({k:d[k] for k in (
        "status","source_shard_start","source_shard_stop_exclusive",
        "photo_opportunities_in_this_partial_source_block",
        "dated_candidates_in_this_partial_source_block"
    )},sort_keys=True))


if __name__=="__main__":main()
