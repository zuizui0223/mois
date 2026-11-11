#!/usr/bin/env python3
"""Audit all 16 tenth original-photo RGB9 ROI-v4 worker units, only real images."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_species_first_palette9 import PALETTE,class_labels,sha256
from scripts.analysis.mois_species_first_incremental_first256 import require

HUES={"white","yellow_orange","red_pink","blue_purple"}
PIX=["flower_fraction_"+p for p in PALETTE]
SCHEMA="mois_species_first_tenth4096_authentic_originalROI9_audit_v1"

def run(source,workers,output):
    source,workers,output=map(Path,(source,workers,output))
    require(not output.exists(),"Do not overwrite existing tenth real-image audit")
    manifest=json.loads((source/"manifest.json").read_text())
    require(manifest.get("schema")=="mois_full_source_observation_sampling_cohort_v1" and
            manifest.get("source_arm")=="MOIS_SPECIES_FIRST_TENTH_PRIORITIZED4096" and
            manifest.get("colour_measured") is False and
            manifest.get("work_units")==16 and
            manifest.get("worker_ids")==list(range(16)) and
            manifest.get("new_photo_attempts_planned")==4096 and
            manifest.get("source_previous_attempts")==1201882,
            "Tenth original-photo manifest is not signed fresh prioritized 4096")
    allframes=[];receipts={}
    for i in range(16):
        source_file=source/"shards"/f"shard-{i:03d}.csv"
        p=workers/f"mois-tenth4096-real-original-worker-{i//4}"/f"shard_{i:03d}"
        terminal=p/"photos.csv";receipt=p/"receipt.json"
        require(source_file.is_file() and terminal.is_file() and receipt.is_file(),
                "Incomplete tenth real original worker photo shard "+str(i))
        require(sha256(source_file)==manifest["photo_shard_sha256"][str(i)],
                "Tenth frozen original-photo source SHA drift")
        r=json.loads(receipt.read_text())
        require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1" and
                r.get("attempted")==256 and r.get("shard")==i and
                r.get("source_sha256")==sha256(source_file) and
                r.get("photos_sha256")==sha256(terminal) and
                r.get("photo_colours_measured") is True,
                "Tenth original photo inference worker receipt SHA mismatch")
        a=pd.read_csv(source_file,dtype=str,keep_default_na=False)
        b=pd.read_csv(terminal,dtype=str,keep_default_na=False)
        require(len(a)==len(b)==256 and a.photo_id.is_unique and b.photo_id.is_unique and
                set(a.photo_id)==set(b.photo_id),
                "Tenth worker photo-ID mismatch")
        for k in ("inat_taxon_id","observation_id"):
            require(a.set_index("photo_id").sort_index()[k].equals(
                    b.set_index("photo_id").sort_index()[k]),
                    "Tenth signed original observation/taxon IDs were modified")
        allframes.append(b)
        receipts[str(i)]=sha256(receipt)
    photos=pd.concat(allframes,ignore_index=True)
    require(len(photos)==4096 and photos.photo_id.is_unique and
            photos.observation_id.is_unique and photos.inat_taxon_id.is_unique,
            "Tenth photo/observation/taxon duplicates found")
    numeric=photos[PIX].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    strict=(photos.morph.isin(HUES).to_numpy(bool) & np.isfinite(numeric).all(axis=1) &
            (numeric>=0).all(axis=1) & (numeric<=1.00001).all(axis=1) &
            (np.abs(np.nan_to_num(numeric).sum(axis=1)-1)<=0.001))
    labels=np.full(len(photos),"",dtype=object)
    labels[strict]=class_labels(numeric[strict])
    photos["strict9_original_photo_numeric"]=strict
    photos["dominant_hue9_label"]=labels
    photos.loc[~strict,PIX]=""
    output.mkdir(parents=True)
    path=output/"mois_tenth4096_actual_original_nine_roi.csv.gz"
    photos.to_csv(path,index=False,compression="gzip")
    report={"schema":SCHEMA,
            "actual_independent_original_photo_attempts":4096,
            "actual_genuine_strict9_numeric_success":int(strict.sum()),
            "actual_clear_nine_palette_hues":int((labels!="").sum()),
            "actual_technical_unknowns":int((~strict).sum()),
            "photo_observation_source_taxon_ID_collisions":0,
            "sixteen_original_photo_worker_receipt_sha256":receipts,
            "frozen_original_manifest_sha256":sha256(source/"manifest.json"),
            "audited_derivative_sha256":sha256(path),
            "GBIF_occurrences_downloaded":0,
            "species_wide_flower_colour_fixed_proven":False}
    (output/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_TENTH4096_AUTHENTIC_REAL_ORIGINAL_ROI9",json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source","workers","output"):p.add_argument("--"+k,required=True)
    args=p.parse_args()
    run(args.source,args.workers,args.output)
