#!/usr/bin/env python3
"""Conservative actual GBIF COL XR species-name match for MOIS nine-colour taxa.

A pilot of sampled identified nominal taxa, NOT a whole-world GBIF occurrence
join. Use official GBIF v2 match against COL XR, storing exact requested
names, raw responses, accepted taxon concepts and unresolved ambiguity.
Never silently coerce synonym/higher-rank or fuzzy matches into a species.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

SCHEMA="mois_species_first_gbif_colxr_v2_name_pilot_v1"
COL_XR="7ddf754f-d193-4cc9-b351-99906754a03b"
API="https://api.gbif.org/v2/species/match"
NAME_RULE=re.compile(r"^[A-Z][A-Za-z-]+\s+[a-z][a-z-]+(?:\s|$)")
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")


def require(ok,reason):
    if not bool(ok):raise ValueError("MOIS_GBIF_MATCH: "+reason)


def make_pilot(catalogue,per_hue=16):
    require(1<=per_hue<=100,"Unsafe per-hue pilot size")
    df=catalogue.copy()
    require({"inat_taxon_id","species","dominant_hue9","hue_status"}<=set(df.columns),
            "Source nine-colour species/name roster incomplete")
    for c in ("inat_taxon_id","species","dominant_hue9","hue_status"):
        df[c]=df[c].fillna("").astype(str).str.strip()
    require(df.inat_taxon_id.is_unique,"Nominal source species IDs repeated")
    df=df.loc[df.dominant_hue9.isin(PALETTE) &
              df.species.str.match(NAME_RULE)].copy()
    df["fixed_hash"]=df.inat_taxon_id.map(lambda x:hashlib.sha256(
        ("MOIS-GBIF-v2-NINEHUE-NAME-PILOT-20261111|"+x).encode()).hexdigest())
    return df.sort_values("fixed_hash",kind="stable").groupby(
        "dominant_hue9",sort=True,observed=True).head(per_hue).sort_values(
        ["dominant_hue9","fixed_hash"],kind="stable").drop(columns=["fixed_hash"]).reset_index(drop=True)


def inspect_response(payload,input_name):
    """Accept only exact species-level COL XR name usages with matching canon."""
    if not isinstance(payload,dict):
        return dict(status="api_schema_unrecognized",key="",matchType="",candidate_canonical="")
    usage=payload.get("usage") or {}
    diagnostics=payload.get("diagnostics") or {}
    if not isinstance(diagnostics,dict):diagnostics={}
    match_type=payload.get("matchType") or diagnostics.get("matchType","")
    if not isinstance(usage,dict):
        usage={}
    accepted=payload.get("acceptedUsage") or {}
    if not isinstance(accepted,dict):
        accepted={}
    key=accepted.get("key") or (usage.get("key") if usage.get("status")=="ACCEPTED" else None)
    can=usage.get("canonicalName") or usage.get("name") or ""
    rank=usage.get("rank") or ""
    taxon_status=usage.get("status") or ""
    accepted_can=accepted.get("canonicalName") or accepted.get("name") or ""
    expected=" ".join(input_name.strip().split()[:2]).lower()
    observed=" ".join(str(can).strip().split()[:2]).lower()
    conservative=(bool(key) and rank=="SPECIES" and
                  taxon_status=="ACCEPTED" and
                  str(match_type).upper()=="EXACT" and expected==observed and
                  len(input_name.split())==2)
    if conservative: status="colxr_exact_accepted_species"
    elif not usage: status="no_match_or_unrecognized"
    elif rank!="SPECIES": status="not_species_rank_review"
    elif taxon_status!="ACCEPTED": status="synonym_or_doubtful_review"
    elif str(match_type).upper()!="EXACT":status="fuzzy_or_nonexact_review"
    else: status="canonical_disagreement_review"
    return {
        "status":status,
        "key":str(key) if key is not None else "",
        "matchType":str(match_type),
        "candidate_canonical":str(can),
        "accepted_canonical":str(accepted_can),
        "candidate_rank":str(rank),
        "candidate_status":str(taxon_status),
        "usage_key":str(usage.get("key") or ""),
    }


def gbif_request(scientific_name,*,attempts=3):
    qs=urlencode({"scientificName":scientific_name,"checklistKey":COL_XR})
    request=Request(API+"?"+qs,headers={
        "User-Agent":"MOIS-species-first-research/1.0 (ecological taxonomic pilot)",
        "Accept":"application/json"})
    last=None
    for attempt in range(attempts):
        try:
            with urlopen(request,timeout=25) as r:
                data=r.read(500000)
            return json.loads(data)
        except HTTPError as e:
            last=e
            if e.code not in (429,500,502,503,504):break
            if attempt+1<attempts:time.sleep(2**attempt)
        except (URLError,TimeoutError,json.JSONDecodeError) as e:
            last=e
            if attempt+1<attempts:time.sleep(2**attempt)
    raise RuntimeError("GBIF v2 COL XR API request failed: "+str(last))


def run(source,output,*,per_hue=16,min_delay=0.3):
    source=Path(source);output=Path(output)
    require(not output.exists(),"Refuse overwrite of existing source-specific GBIF pilot")
    require(0.2<=min_delay<=5.0,"Do not disable polite API pacing")
    pilot=make_pilot(pd.read_csv(source,compression="infer",dtype=str),per_hue=per_hue)
    require(len(pilot)>=30,"Insufficient independently named hue taxon samples")
    output.mkdir(parents=True)
    statuses=[];responses=[]
    for i,row in enumerate(pilot.itertuples(index=False)):
        try:
            raw=gbif_request(row.species)
            interpret=inspect_response(raw,row.species)
            code="ok"
        except Exception as e:
            raw={"request_error":str(e)}
            interpret=dict(status="api_error",key="",matchType="",candidate_canonical="")
            code="error"
        record={"inat_taxon_id":row.inat_taxon_id,"source_species":row.species,
                "dominant_hue9":row.dominant_hue9,
                "hue_status":row.hue_status,
                "requested_checklistKey":COL_XR,"api_version":"v2",
                **interpret}
        statuses.append(record)
        responses.append({"inat_taxon_id":row.inat_taxon_id,
                          "original_source_species":row.species,"response":raw})
        if (i+1)%20==0:print("MOIS_GBIF_COLXR_PILOT",i+1,len(pilot),flush=True)
        time.sleep(min_delay)
    out=pd.DataFrame(statuses)
    accepted=out.loc[out.status.eq("colxr_exact_accepted_species")].copy()
    # Different iNaturalist source labels can resolve to the same GBIF accepted
    # species; do not double count species in the distribution join.
    duplicates=set(accepted.loc[accepted.key.duplicated(keep=False),"key"])
    if duplicates:
        out.loc[out.key.isin(duplicates), "status"]="duplicate_accepted_species_key_review"
    (output/"gbif_colxr_v2_species_name_pilot.csv").write_text(
        out.to_csv(index=False),encoding="utf-8")
    with gzip.open(output/"gbif_colxr_v2_raw_responses.jsonl.gz","wt",encoding="utf-8") as f:
        for row in responses:f.write(json.dumps(row,ensure_ascii=False)+"\n")
    report={"schema":SCHEMA,"input_nine_colour_roster_sha256":hashlib.sha256(
        source.read_bytes()).hexdigest(),
        "gbif_taxonomy":"COL_XR","gbif_checklist_key":COL_XR,
        "gbif_api":API,
        "nominal_taxa_attempted":len(pilot),
        "colour_classes_sampled":int(pilot.dominant_hue9.nunique()),
        "match_status_counts":{str(k):int(v) for k,v in out.status.value_counts().items()},
        "conservative_accepted_taxa":int(out.status.eq("colxr_exact_accepted_species").sum()),
        "ambiguous_review_or_errors":int(out.status.ne("colxr_exact_accepted_species").sum()),
        "published_GBIF_occurrence_download_DOI":False,
        "GBIF_occurrence_records_downloaded":0,
        "non_exhaustive_pilot":True,
        "photo_colour_is_not_geographically_fixed":True}
    (output/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_GBIF_COLXR_NAME_PILOT_REAL",json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--per-hue",type=int,default=16)
    p.add_argument("--min-delay",type=float,default=0.3)
    x=p.parse_args()
    run(x.source,x.output,per_hue=x.per_hue,min_delay=x.min_delay)


if __name__=="__main__":
    main()
