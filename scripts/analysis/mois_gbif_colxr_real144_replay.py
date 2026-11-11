#!/usr/bin/env python3
"""Offline reclassification of genuine GBIF v2 name match responses (no new API).

Corrects v2 diagnostics.matchType schema handling on the exact signed 144
response payload, distinguishing 138 exact accepted species, synonyms
and higher-rank matches. Preserves raw payload SHA for audit.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json
from pathlib import Path
import pandas as pd
from scripts.analysis.mois_gbif_colxr_nine_palette_name_pilot import (
    COL_XR, inspect_response, require
)

SCHEMA="mois_gbif_v2_colxr_144_response_replay_v1"
ORIGINAL_CSV_SHA="5a0c234993e4273bc137d282c653f64477680799bd921d009364100c3cd41174"
ORIGINAL_RAW_SHA="9d4dea5f609c01c7332630cf17414ff2875137acdf30c6f23fbbbf2ecd0300fe"


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def replay(source,raw_path,output,*,verify_shas=True):
    source,raw_path,output=map(Path,(source,raw_path,output))
    require(not output.exists(),"Existing frozen source replay cannot be overwritten")
    if verify_shas:
        require(sha(source)==ORIGINAL_CSV_SHA and sha(raw_path)==ORIGINAL_RAW_SHA,
                "Original GBIF API v2 sample raw SHA changed")
    old=pd.read_csv(source,dtype=str,keep_default_na=False)
    with gzip.open(raw_path,"rt",encoding="utf-8") as f:
        raw=[json.loads(line) for line in f]
    require(len(old)==len(raw)==144 and old.inat_taxon_id.is_unique,
            "Original 144 sampled source taxa or matching responses lost")
    request_index={str(d["inat_taxon_id"]):d for d in raw}
    require(len(request_index)==144,"Duplicate original GBIF response taxon IDs")
    rows=[]
    for r in old.itertuples(index=False):
        original=request_index[r.inat_taxon_id]
        require(original["original_source_species"]==r.source_species,
                "Original GBIF input name changed after real API request")
        detail=inspect_response(original["response"],r.source_species)
        rows.append(dict(inat_taxon_id=r.inat_taxon_id,
                         source_species=r.source_species,
                         dominant_hue9=r.dominant_hue9,
                         hue_status=r.hue_status,
                         requested_checklistKey=COL_XR,
                         api_version="v2",
                         **detail))
    df=pd.DataFrame(rows)
    accepted=df.loc[df.status.eq("colxr_exact_accepted_species")]
    duplicate_keys=set(accepted.loc[accepted.key.duplicated(keep=False),"key"])
    if duplicate_keys:
        df.loc[df.key.isin(duplicate_keys),"status"]="duplicate_accepted_species_key_review"
    report={"schema":SCHEMA,"original_GBIF_API_run":42930778507,
            "original_response_count":len(raw),
            "original_csv_sha256":sha(source),
            "original_raw_response_sha256":sha(raw_path),
            "gbif_checklistKey":COL_XR,
            "status_counts":{str(k):int(v) for k,v in df.status.value_counts().items()},
            "exact_accepted_species":int(df.status.eq("colxr_exact_accepted_species").sum()),
            "no_GBIF_requeries_for_correction":True,
            "occurrence_records_not_yet_fetched":True}
    output.mkdir(parents=True)
    df.to_csv(output/"mois_real144_gbif_colxr_reinterpreted_names.csv",index=False)
    df.loc[df.status.eq("colxr_exact_accepted_species")].to_csv(
        output/"mois_real144_138_exact_accepted_species_keys.csv",index=False)
    (output/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_ACTUAL_GBIF_V2_REPLAY_144_NAMES",json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-csv",required=True);p.add_argument("--raw",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args();replay(a.original_csv,a.raw,a.output)


if __name__=="__main__":
    main()
