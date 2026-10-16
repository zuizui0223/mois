#!/usr/bin/env python3
"""Build a PRIVATE selected-observation flowering annotation overlay.

Requires a pre-authorized, independently acquired *iNaturalist observation*
API/export JSONL. Does NOT issue web requests or claim that the four-table
Open Dataset contains annotations. Missing observation UUID in JSONL means
NOT_LOOKED_UP (unknown), not NO_FLOWERS. An empty annotations array for a
verified observation means LOOKED_UP with NO_ANNOTATION.

The overlay is an independent, current-time observation label. The historical
2026-09-27 photo source does not imply annotations existed on that date.
This must never be a back-door colour-based photo inclusion decision.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from datetime import date
from pathlib import Path

OUT=["observation_id","annotation_term_id","annotation_value_ids",
     "lookup_status","annotation_source","lookup_snapshot_date"]
VALID={13,14,15,21}


def sha(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def build(sample_csv:Path, snapshot_jsonl:Path, snapshot_sha256:str,
          snapshot_date:str,out_csv:Path,out_receipt:Path,
          *,synthetic:bool=False)->dict:
    if out_csv.exists() or out_receipt.exists() or out_csv.with_suffix(out_csv.suffix+".incomplete").exists():
        raise ValueError("No overwrite of private annotation snapshots")
    if len(snapshot_sha256)!=64 or sha(snapshot_jsonl)!=snapshot_sha256:
        raise ValueError("Source iNaturalist observation JSONL has wrong SHA256")
    date.fromisoformat(snapshot_date)
    selected=set()
    with sample_csv.open(newline="",encoding="utf8") as f:
        reader=csv.DictReader(f)
        if not {"observation_id","inat_taxon_id","photo_id"}.issubset(reader.fieldnames or []):
            raise ValueError("No source-frozen selected observations in private Mois sample")
        for r in reader:
            ident=r["observation_id"].strip()
            if not ident or ident in selected:
                raise ValueError("Sampling should retain exactly one photo per source observation")
            selected.add(ident)
    checked={}
    with snapshot_jsonl.open(encoding="utf8") as f:
        for line in f:
            if not line.strip():
                raise ValueError("Empty record in claimed complete annotation source JSONL")
            item=json.loads(line)
            uuid=item.get("uuid")
            if not uuid or uuid in checked or uuid not in selected:
                raise ValueError("Selected source observation missing, duplicated or unexpected in annotation input")
            annotations=item.get("annotations")
            if not isinstance(annotations,list):
                raise ValueError("No complete annotations array for verified observation")
            ids=set()
            for a in annotations:
                if not isinstance(a,dict):
                    raise ValueError("Malformed platform annotation")
                if int(a.get("controlled_attribute_id",-1))!=12:
                    continue
                x=int(a.get("controlled_value_id",-1))
                if x not in VALID:
                    raise ValueError("Unrecognized Flowers/Fruits phenology value")
                # iNat annotations have voting. Do not use disagreements
                # with nonpositive votes as confirmed flowering evidence.
                if "vote_score" not in a or not isinstance(a["vote_score"],int):
                    raise ValueError("iNaturalist annotation evidence must contain vote_score")
                if a["vote_score"]>0:
                    ids.add(x)
            if 21 in ids and len(ids)>1:
                raise ValueError("Conflicting positive annotation vs explicit no Flowers/Fruits")
            checked[uuid]="|".join(map(str,sorted(ids)))
    out_csv.parent.mkdir(parents=True,exist_ok=True)
    out_receipt.parent.mkdir(parents=True,exist_ok=True)
    tmp=out_csv.with_suffix(out_csv.suffix+".incomplete")
    with tmp.open("w",newline="",encoding="utf8") as f:
        writer=csv.writer(f);writer.writerow(OUT)
        for ident in sorted(selected):
            found=ident in checked
            writer.writerow([ident,"12",checked.get(ident,""),
                             "LOOKED_UP" if found else "NOT_LOOKED_UP",
                             "inat_api_v1",snapshot_date])
    result={
        "schema":"mois_observation_phenology_annotation_overlay_preparation_v1",
        "status":"SELECTED_OBSERVATION_PHENOLOGY_OVERLAY_NOT_IMAGE_CLASSIFICATION",
        "source_sample_sha256":sha(sample_csv),
        "source_snapshot_jsonl_sha256":snapshot_sha256,
        "source_query_date":snapshot_date,
        "n_selected_original_observations":len(selected),
        "n_observations_annotation_LOOKED_UP":len(checked),
        "n_observations_NOT_LOOKED_UP":len(selected)-len(checked),
        "n_checked_with_accepted_Flowers_value_13":sum(13 in
            {int(x) for x in token.split("|") if x} for token in checked.values()),
        "private_output_sha256":sha(tmp),
        "annotation_source_not_original_monthly_20260927_OpenData":True,
        "observation_not_photo_level":True,
        "not_separately_verified_original_flower_pixels":True,
        "can_designate_unchecked_observations_as_not_flowering":False,
        "synthetic_only":bool(synthetic),
    }
    tmp.replace(out_csv)
    out_receipt.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--private-selected-sample",required=True,type=Path)
    p.add_argument("--authorized-observations-jsonl",required=True,type=Path)
    p.add_argument("--authorized-observations-sha256",required=True)
    p.add_argument("--annotation-snapshot-date",required=True)
    p.add_argument("--out-private-overlay",required=True,type=Path)
    p.add_argument("--out-private-receipt",required=True,type=Path)
    a=p.parse_args()
    report=build(a.private_selected_sample,a.authorized_observations_jsonl,
                 a.authorized_observations_sha256,a.annotation_snapshot_date,
                 a.out_private_overlay,a.out_private_receipt)
    print(json.dumps({k:report[k] for k in (
        "status","n_selected_original_observations",
        "n_observations_annotation_LOOKED_UP",
        "n_observations_NOT_LOOKED_UP",
        "n_checked_with_accepted_Flowers_value_13")},sort_keys=True))

if __name__=="__main__":
    main()
