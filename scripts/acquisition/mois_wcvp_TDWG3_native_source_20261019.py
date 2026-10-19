#!/usr/bin/env python3
"""Audit Kew WCVP accepted taxa/native TDWG3 metadata against 350399 REAL MOIS photographs.

This is a SOURCE-ONLY taxon-range lookup, not a photographed-petal,
island-native, independently replicated or phenological colour effect.
The entire downloaded official WCVP archive is read once/twice via streams.
Only unambiguously name-matched, accepted species IDs are retained.

WCVP TDWG level 3 native status NEVER establishes nativeness of one
specific island, original photo point or plant population. No polygon overlay
or point-native classifications are emitted by this script.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    assemble, require, sha256,
)

FRAME_N=350399
KNOWN_HUES=164297
RANGE_FILE="mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz"
MATCH_FILE="mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"
RECEIPT="mois350399_wcvp_TDWG3_source_only_taxon_native_receipt.json"
BINOMIAL=re.compile(r"^[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+$")
WCVP_SOURCE="https://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip"


def clean_name(name):
    """No genus-only, unidentified, qualifier or trinomial coercion."""
    if not isinstance(name,str):
        return ""
    name=" ".join(name.strip().split())
    return name if BINOMIAL.fullmatch(name) else ""


def _member_name(z, suffixes):
    matches=[n for n in z.namelist()
             if not n.endswith("/") and Path(n).name.lower() in suffixes]
    require(len(matches)==1,
            "Expected exactly one official WCVP CSV member from "+str(sorted(suffixes)))
    return matches[0]


def wcvp_reader(z, member):
    with z.open(member) as raw:
        reader=csv.DictReader(io.TextIOWrapper(
            raw,encoding="utf-8-sig",errors="strict",newline=""))
        require(reader.fieldnames is not None,"WCVP header is absent")
        for row in reader:
            yield row


def source_taxon_table(real_photos, *, expected_photos=FRAME_N, expected_hues=KNOWN_HUES):
    require(len(real_photos)==expected_photos and real_photos.photo_id.is_unique
            and int(real_photos.classified.sum())==expected_hues,
            "MOIS source must be genuinely attempted 350399 pixel ROIs")
    d=real_photos[["inat_taxon_id","species","photo_id","classified","arm"]].copy()
    group=d.groupby("inat_taxon_id",sort=True)
    taxa=group.agg(
        source_species_n=("species","nunique"),
        original_photos=("photo_id","size"),
        classifiable_photos=("classified","sum"),
        source_species=("species","first")
    ).reset_index()
    taxa["match_binomial"]=taxa.source_species.map(clean_name)
    taxa["source_taxon_inconsistent"]=taxa.source_species_n.ne(1)
    return taxa


def collect_wcvp_names(z, wanted, *, min_source_rows=10000):
    """Exact species-rank names, resolving synonyms only via accepted WCVP ID."""
    member=_member_name(z,{"wcvp_names.csv"})
    header=None
    candidates=collections.defaultdict(set)
    candidate_statuses=collections.defaultdict(set)
    n_total=0
    for row in wcvp_reader(z,member):
        if header is None:
            header=set(row)
            require({"plant_name_id","taxon_name","taxon_rank","taxon_status",
                     "accepted_plant_name_id"}<=header,
                    "WCVP accepted-name authority or rank schema changed")
        n_total+=1
        name=clean_name(row["taxon_name"])
        if not name or name not in wanted:
            continue
        rank=row["taxon_rank"].strip().casefold()
        state=row["taxon_status"].strip().casefold()
        if rank!="species" or state not in ("accepted","synonym"):
            candidate_statuses[name].add("unsupported_rank_or_status")
            continue
        accepted=(row["accepted_plant_name_id"].strip()
                  or row["plant_name_id"].strip()
                  if state=="accepted" else row["accepted_plant_name_id"].strip())
        if not accepted:
            candidate_statuses[name].add("no_accepted_id")
            continue
        candidates[name].add(accepted)
        candidate_statuses[name].add(state)
    require(n_total>=min_source_rows,
            "Official archive has too few plant names; refuse incomplete WCVP source")
    return candidates,candidate_statuses,n_total,member


def resolve_photo_taxa(taxa, candidates, states, *, expected_photos=FRAME_N):
    resolved=[]
    selected=set()
    for row in taxa.to_dict("records"):
        original=str(row["source_species"])
        name=row["match_binomial"]
        values=candidates.get(name,set())
        if not str(row["inat_taxon_id"]).strip():
            match="missing_original_taxon_id"
        elif row["source_taxon_inconsistent"]:
            match="inconsistent_original_taxon_id"
        elif not name:
            match="not_exact_species_binomial"
        elif len(values)==1:
            match="unique_exact_name_to_accepted_id"
        elif len(values)>1:
            match="multiple_accepted_ids_ambiguous"
        else:
            match="no_accepted_species_match"
        accepted=next(iter(values)) if match=="unique_exact_name_to_accepted_id" else ""
        if accepted:selected.add(accepted)
        resolved.append({
            "inat_taxon_id":str(row["inat_taxon_id"]),
            "source_species":original,
            "source_species_exact_binomial":name,
            "source_original_photo_attempts":int(row["original_photos"]),
            "source_four_hue_classified_photos":int(row["classifiable_photos"]),
            "source_taxon_name_matches":int(len(values)),
            "wcvp_match_state":match,
            "wcvp_accepted_plant_name_id":accepted,
            "wcvp_name_status_seen":";".join(sorted(states.get(name,set()))),
            "wcvp_native_on_specific_island": "not_assessed",
            "source_photo_point_native_status":"not_assessed"
        })
    out=pd.DataFrame(resolved)
    require(out.inat_taxon_id.is_unique and
            int(out.source_original_photo_attempts.sum())==expected_photos,
            "Original WCVP taxa summary lost photo IDs")
    return out,selected


def accepted_name_metadata(z, selected):
    """Get official accepted names for exact accepted IDs, not guessed synonyms."""
    member=_member_name(z,{"wcvp_names.csv"})
    accepted={}
    for row in wcvp_reader(z,member):
        pid=row["plant_name_id"].strip()
        if pid not in selected:
            continue
        if row["taxon_status"].strip().casefold()!="accepted":
            continue
        name=clean_name(row["taxon_name"])
        if row["taxon_rank"].strip().casefold()!="species" or not name:
            continue
        if pid in accepted and accepted[pid]!=name:
            raise ValueError("WCVP accepted plant ID has inconsistent accepted species names")
        accepted[pid]=name
    return accepted


def _flag(x):
    s=str(x).strip().casefold()
    return {"0":False,"1":True,"false":False,"true":True,"f":False,"t":True}.get(s)


def _tdwg_state(row):
    """Distinguish uncertain, extinct and conflicting locality records."""
    doubtful=_flag(row["location_doubtful"])
    extinct=_flag(row["extinct"])
    introduced=_flag(row["introduced"])
    if doubtful is True:return "doubtful"
    if extinct is True:return "extinct"
    if doubtful is None or extinct is None or introduced is None:return "unknown_flag"
    return "introduced_at_TDWG3" if introduced else "native_at_TDWG3"


def native_distribution_TDWG3(z, selected, accepted, *, min_source_rows=10000):
    member=_member_name(z,{"wcvp_distribution.csv","wcvp_distributions.csv"})
    seen=collections.defaultdict(set)
    loc_names=collections.defaultdict(set)
    n_source=0
    matched=0
    for row in wcvp_reader(z,member):
        if n_source==0:
            require({"plant_name_id","area_code_l3","introduced","extinct",
                     "location_doubtful"}<=set(row),
                    "WCVP TDWG level 3 native/introduced schema changed")
        n_source+=1
        pid=row["plant_name_id"].strip()
        if pid not in selected or pid not in accepted:
            continue
        area=row["area_code_l3"].strip()
        if not area:continue
        seen[(pid,area)].add(_tdwg_state(row))
        loc_names[(pid,area)].add(row.get("area","").strip())
        matched+=1
    require(n_source>=min_source_rows,
            "Official WCVP distribution table is implausibly small")
    out=[]
    for (pid,area), states in sorted(seen.items()):
        state=next(iter(states)) if len(states)==1 else "contradictory_TDWG3_records"
        out.append({
            "wcvp_accepted_plant_name_id":pid,
            "wcvp_accepted_name":accepted[pid],
            "tdwg_level3_code":area,
            "tdwg_level3_name":";".join(sorted(loc_names[(pid,area)])),
            "native_introduced_state":state,
            "official_status_records_combined":len(states),
            "wcvp_status_conflict":len(states)>1,
            "native_on_specific_island":"not_established_by_TDWG3",
            "native_at_original_photo_lat_lon":"not_evaluated_no_TDWG3_polygon"
        })
    return pd.DataFrame(out), n_source, matched, member


def run(core,broad,wave2,cached8192,tail16384,wcvp,out):
    out=Path(out);require(not out.exists() or not any(out.iterdir()),
                          "Never overwrite source-locked WCVP/native audit")
    photographs,arms=assemble(core,broad,wave2,cached8192,tail16384,strict=True)
    taxa=source_taxon_table(photographs)
    require(Path(wcvp).is_file(),"Missing official original Kew WCVP archive")
    with zipfile.ZipFile(wcvp) as archive:
        require(archive.testzip() is None,"Official WCVP zip archive failed CRC integrity")
        candidates,states,name_n,name_member=collect_wcvp_names(
            archive,set(taxa.match_binomial)-{""})
        names,accepted_ids=resolve_photo_taxa(taxa,candidates,states)
        accepted=accepted_name_metadata(archive,accepted_ids)
        # An ID linked to an unavailable accepted name cannot be treated as native.
        avail=set(accepted)
        unmatched=names.wcvp_accepted_plant_name_id.ne("") & ~names.wcvp_accepted_plant_name_id.isin(avail)
        names.loc[unmatched,"wcvp_match_state"]="accepted_id_not_verified_as_species"
        names.loc[unmatched,"wcvp_accepted_plant_name_id"]=""
        native,n_distributions,n_native,distribution_member=native_distribution_TDWG3(
            archive,set(accepted),accepted)
    out.mkdir(parents=True,exist_ok=True)
    names.to_csv(out/MATCH_FILE,index=False,compression="gzip",lineterminator="\n")
    native.to_csv(out/RANGE_FILE,index=False,compression="gzip",lineterminator="\n")
    report={
        "schema":"mois350399_kew_wcvp_official_accepted_species_TDWG3_native_source_v1",
        "status":"SOURCE_MATCH_ONLY_NOT_PHOTO_POINT_NATIVE",
        "source_photo_attempts":len(photographs),
        "source_four_hue_classified":int(photographs.classified.sum()),
        "taxa_with_any_source_photos":len(taxa),
        "taxa_with_classified_photos":int(names.source_four_hue_classified_photos.gt(0).sum()),
        "wcvp_original_source":WCVP_SOURCE,
        "wcvp_archive_sha256":sha256(wcvp),
        "wcvp_archive_bytes":Path(wcvp).stat().st_size,
        "wcvp_name_source_rows_scanned":name_n,
        "wcvp_distribution_source_rows_scanned":n_distributions,
        "wcvp_distribution_matching_rows_scanned":n_native,
        "wcvp_name_archive_member":name_member,
        "wcvp_distribution_archive_member":distribution_member,
        "photo_taxon_exact_match_counts":{
            str(k):int(v) for k,v in names.wcvp_match_state.value_counts().items()
        },
        "taxa_with_verified_WCVP_accepted_species":len(accepted),
        "matched_TDWG3_region_accepted_species":int(
            native.wcvp_accepted_plant_name_id.nunique()) if len(native) else 0,
        "matched_TDWG3_native_and_introduced_states":{
            str(k):int(v) for k,v in
            native.native_introduced_state.value_counts().items()
        } if len(native) else {},
        "source_parquet_sha256":{
            arm:sha256(path) for arm,path in (
                ("core",core),("broad",broad),("wave2",wave2),
                ("cached8192",cached8192),("tail16384",tail16384)
            )
        },
        "wcvp_distributions_resolution":"TDWG_level3_only",
        "species_taxon_ID_not_equal_WCVP_plant_name_id":True,
        "ambiguous_or_missing_taxonomy_not_forced_to_native":True,
        "island_native_status_inferred":False,
        "photo_coordinate_overlaid_TDWG3_polygons":False,
        "independent_island_flora_endemism_verified":False,
        "H3_native_island_hypothesis_tested":False,
        "hue_measurements_changed":False,
        "outcome_confirmatory_analysis_run":False,
        "data_files_sha256":{name:sha256(out/name) for name in (MATCH_FILE,RANGE_FILE)}
    }
    (out/RECEIPT).write_text(json.dumps(report,indent=2,sort_keys=True,
                                       ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_Kew_WCVP350399_REAL_PHOTO_NAME_TDWG3_NATIVE_SOURCE_ONLY",
          json.dumps({k:report[k] for k in (
              "source_photo_attempts","taxa_with_any_source_photos",
              "taxa_with_classified_photos",
              "taxa_with_verified_WCVP_accepted_species",
              "matched_TDWG3_region_accepted_species",
              "photo_taxon_exact_match_counts",
              "matched_TDWG3_native_and_introduced_states")},
              sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("core","broad","wave2","cached8192","tail16384","wcvp","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.core,a.broad,a.wave2,a.cached8192,a.tail16384,a.wcvp,a.output)

if __name__=="__main__":
    main()
