#!/usr/bin/env python3
"""MOIS measured-photo expansion: append FCP's already measured 149,900 photos.

No botanical review, no image downloads, no new measurement, no changed original
MOIS regression. The two studies retain distinct estimands despite source reuse.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

PHOTO_SHA="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
FCP_SHA={
 "discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
 "validation":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
 "prospective":"57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
MEMBERS=(
 "fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz",
 "fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz",
)
COLORS=("white","yellow_orange","red_pink","blue_purple")
FIELDS=("photo_id","inat_taxon_id","species","morph","source_original_measured",
        "old_mois_or_new_fcp","measurement_status","latitude","longitude",
        "photo_license","observation_id")
FCP_N={"discovery":50000,"validation":50000,"prospective":49900}


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def require(cond:bool,msg:str):
    if not cond:raise ValueError(msg)


def normalized(row:dict,source:str,kind:str)->dict:
    r={k:str(row.get(k,"") or "").strip() for k in FIELDS}
    require(r["photo_id"].isdigit() and r["inat_taxon_id"].isdigit()
            and r["species"] and r["morph"],
            "Source original photograph/taxon/phenotype identity invalid")
    require(r["morph"] in (*COLORS,"mixed_uncertain"),
            "Unrecognized original fixed four-hue classification state")
    r["source_original_measured"]=source
    r["old_mois_or_new_fcp"]=kind
    return r


def read_historical(path:Path):
    require(sha(path)==PHOTO_SHA,"MOIS original 100543 photograph source SHA drift")
    first={}
    raw=0; conflict=0
    with zipfile.ZipFile(path) as z:
        for member in MEMBERS:
            with gzip.open(io.BytesIO(z.read(member)),"rt",encoding="utf-8",newline="") as f:
                for row in csv.DictReader(f):
                    raw+=1
                    pid=row["photo_id"]
                    if pid in first:
                        old=first[pid]
                        if (old["morph"],old["inat_taxon_id"],old["observation_id"])!=(
                            row["morph"],row["inat_taxon_id"],row["observation_id"]):
                            conflict+=1
                        continue
                    first[pid]=normalized(row,"historical_mois_20261014","historical")
    require((raw,len(first),conflict)==(127448,100543,1),
            "Original MOIS source denominator or conflict count changed")
    require(sum(v["morph"] in COLORS for v in first.values())==45953,
            "Original 45953 classifier-admitted hue labels changed")
    return first,raw,conflict


def run(photo_zip:Path,fcp_paths:dict[str,Path],outdir:Path)->dict:
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Refuse to overwrite original merged-photograph execution")
    old,old_rows,old_conflicts=read_historical(photo_zip)
    outdir.mkdir(parents=True,exist_ok=True)
    old_stats=Counter()
    all_stats=Counter()
    cohorts={}
    taxonomy=defaultdict(set)
    by_taxon=defaultdict(lambda:Counter())
    repeated_fcp=set()
    old_ids=set(old)
    seen_fcp=set()
    overlap_ids=set()
    class_shared=0;class_agreed=0;taxon_conflicts=0
    species_source=defaultdict(set)
    union_csv=outdir/"mois_fcp_combined_original_247032_photo_hues.csv.gz"
    # Byte-deterministic compressed source table.
    with union_csv.open("wb") as raw:
        with gzip.GzipFile(filename="",fileobj=raw,mode="wb",mtime=0) as gz:
            with io.TextIOWrapper(gz,encoding="utf-8",newline="") as out:
                writer=csv.DictWriter(out,fieldnames=FIELDS,lineterminator="\n")
                writer.writeheader()
                for pid,r in old.items():
                    writer.writerow(r)
                    old_stats["photos"]+=1
                    old_stats["classified"]+=r["morph"] in COLORS
                    all_stats["photos"]+=1
                    all_stats["classified"]+=r["morph"] in COLORS
                    all_stats["hue_"+r["morph"]]+=1
                    species_source[r["inat_taxon_id"]].add("historical")
                    taxonomy[r["inat_taxon_id"]].add(r["species"])
                    by_taxon[r["inat_taxon_id"]]["historical_photos"]+=1
                    by_taxon[r["inat_taxon_id"]]["historical_classified"]+=r["morph"] in COLORS
                    by_taxon[r["inat_taxon_id"]]["historical_"+r["morph"]]+=1
                for name,path in fcp_paths.items():
                    require(sha(path)==FCP_SHA[name],f"Original FCP {name} SHA mismatch")
                    count=overlap=classified_new=classified_fcp=0
                    with path.open(newline="",encoding="utf-8-sig") as f:
                        reader=csv.DictReader(f)
                        required={"photo_id","inat_taxon_id","species","morph","global_classifiable",
                                  "latitude","longitude"}
                        require(required.issubset(reader.fieldnames or ()),
                                f"Original FCP {name} no standard hue/source fields")
                        for source in reader:
                            count+=1
                            pid=source["photo_id"].strip()
                            require(pid not in seen_fcp,"FCP 3 source cohorts share a photo ID")
                            seen_fcp.add(pid)
                            r=normalized(source,"fcp_"+name,"new_fcp")
                            yes=r["morph"] in COLORS
                            require(source["global_classifiable"].strip().lower() == ("true" if yes else "false"),
                                    f"Original FCP {name} classifiability changed")
                            classified_fcp+=yes
                            if pid in old_ids:
                                overlap+=1
                                overlap_ids.add(pid)
                                oldr=old[pid]
                                if yes and oldr["morph"] in COLORS:
                                    class_shared+=1
                                    class_agreed+=r["morph"]==oldr["morph"]
                                if oldr["inat_taxon_id"]!=r["inat_taxon_id"]:
                                    taxon_conflicts+=1
                                continue
                            # Append ONLY new source photo IDs, never pooled duplicates.
                            classified_new+=yes
                            writer.writerow(r)
                            all_stats["photos"]+=1
                            all_stats["classified"]+=yes
                            all_stats["hue_"+r["morph"]]+=1
                            by_taxon[r["inat_taxon_id"]]["new_fcp_photos"]+=1
                            by_taxon[r["inat_taxon_id"]]["new_fcp_classified"]+=yes
                            by_taxon[r["inat_taxon_id"]]["new_fcp_"+r["morph"]]+=1
                            species_source[r["inat_taxon_id"]].add("new_fcp")
                            taxonomy[r["inat_taxon_id"]].add(r["species"])
                    require(count==FCP_N[name],f"FCP original {name} source length changed")
                    cohorts[name]={"original_photo_rows":count,
                                   "shared_with_original_mois":overlap,
                                   "appended_new_photo_ids":count-overlap,
                                   "fcp_original_colour_classifiable":classified_fcp,
                                   "appended_new_four_hue_classifiable":classified_new}
    require(len(seen_fcp)==149900,"Original FCP 149900 photo IDs are not distinct")
    require([cohorts[x]["shared_with_original_mois"] for x in FCP_N]==[1248,1326,837],
            "Frozen original FCP/MOIS photo-ID overlap has changed")
    require((len(overlap_ids),class_shared,class_agreed,taxon_conflicts)==(3411,1693,1693,1),
            "Measured shared-photo hue agreement/taxon anomaly differs from original audit")
    require(all_stats["photos"]==247032,
            "Original 100543+149900-3411 unique photo-ID census changed")
    require(old_stats["photos"]==100543 and old_stats["classified"]==45953,
            "Original MOIS historical evidence changed")

    species_csv=outdir/"mois_fcp_combined_original_species_hue_counts.csv"
    with species_csv.open("x",encoding="utf-8",newline="") as f:
        cols=("inat_taxon_id","species_from_sources","name_conflict",
              "historical_photos","historical_classified",
              "new_fcp_photos","new_fcp_classified",
              "combined_photos","combined_classified",
              "white","yellow_orange","red_pink","blue_purple",
              "mixed_uncertain","photographic_white_fraction_if_classified")
        w=csv.DictWriter(f,fieldnames=cols,lineterminator="\n")
        w.writeheader()
        for tax,counts in sorted(by_taxon.items(),key=lambda p:int(p[0])):
            n=counts["historical_photos"]+counts["new_fcp_photos"]
            k=counts["historical_classified"]+counts["new_fcp_classified"]
            shades={c:counts["historical_"+c]+counts["new_fcp_"+c]
                    for c in (*COLORS,"mixed_uncertain")}
            require(n==sum(shades.values()) and k==sum(shades[c] for c in COLORS),
                    "Taxon original complete hue measurements not additive")
            w.writerow({
                "inat_taxon_id":tax,"species_from_sources":" | ".join(sorted(taxonomy[tax])),
                "name_conflict":len(taxonomy[tax])>1,
                "historical_photos":counts["historical_photos"],
                "historical_classified":counts["historical_classified"],
                "new_fcp_photos":counts["new_fcp_photos"],
                "new_fcp_classified":counts["new_fcp_classified"],
                "combined_photos":n,"combined_classified":k,
                **shades,
                "photographic_white_fraction_if_classified":f'{shades["white"]/k:.8f}' if k else ""
            })
    candidate_counts=Counter()
    for counts in by_taxon.values():
        for label,n in (("historical",counts["historical_classified"]),
                        ("union",counts["historical_classified"]+counts["new_fcp_classified"])):
            for depth in (1,2,5,10,20,40,50):
                candidate_counts[f"{label}_taxa_classified_ge_{depth}"]+=int(n>=depth)
    summary={
        "schema":"mois_historical_and_fcp_highdepth_measured_photo_union_v1",
        "status":"REAL_FROZEN_FOUR_HUE_PHOTO_UNION_NO_NEW_IMAGE_PIXELS",
        "historical_original_mois_source_archive_sha256":PHOTO_SHA,
        "historical_original_photo_rows_before_dedup":old_rows,
        "historical_internal_duplicate_conflicts":old_conflicts,
        "historical_unique_original_photo_ids":old_stats["photos"],
        "historical_four_hue_classified":old_stats["classified"],
        "fcp_original_three_cohort_sha256":FCP_SHA,
        "fcp_original_highdepth_photo_ids":len(seen_fcp),
        "fcp_original_source_1499_species":1499,
        "cross_repository_shared_photo_ids":len(overlap_ids),
        "classifiable_shared_photo_ids":class_shared,
        "shared_four_hue_labels_agree":class_agreed,
        "shared_photo_different_taxon_id_warning":taxon_conflicts,
        "cohorts":cohorts,
        "new_unique_photo_ids_added":len(seen_fcp)-len(overlap_ids),
        "new_classified_hue_photo_ids_added":all_stats["classified"]-old_stats["classified"],
        "combined_unique_photo_ids":all_stats["photos"],
        "combined_classified_hue_photo_ids":all_stats["classified"],
        "combined_technical_unclassified_photo_ids":all_stats["photos"]-all_stats["classified"],
        "combined_hue_counts":{c:all_stats["hue_"+c] for c in (*COLORS,"mixed_uncertain")},
        "historical_and_union_taxon_depth":dict(candidate_counts),
        "taxa_with_any_photo":len(by_taxon),
        "taxa_with_both_photo_sources":sum(len(x)==2 for x in species_source.values()),
        "taxon_id_with_changed_species_string_count":sum(len(x)>1 for x in taxonomy.values()),
        "photo_hue_union_csv_sha256":sha(union_csv),
        "photo_hue_species_counts_sha256":sha(species_csv),
        "historical_global_genus_moisture_slopes_reestimated":False,
        "new_fcp_photo_site_month_climate_join_complete":False,
        "cross_source_measurement_independence":False,
        "photo_classifier_algorithm_source":"Identical original ROI-v4 family; 1693 of 1693 double-classifiable overlap photo labels agree",
        "sampling_note":"Original MOIS is broad and shallow, FCP adds high-depth photos in 1499 source-selected taxa. These are NOT 147489 uniformly randomized new taxa. Species/genus weighting and source-stratified inference must accompany new macroclimate fits.",
        "causal_adaptation_or_genetic_polymorphism_claim":False,
    }
    (outdir/"mois_fcp_measured_photo_union_receipt.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("MOIS_ACTUAL_247032_MEASURED_PHOTO_DATA_UNION",json.dumps({
        k:summary[k] for k in ("historical_unique_original_photo_ids",
                              "new_unique_photo_ids_added","new_classified_hue_photo_ids_added",
                              "combined_unique_photo_ids","combined_classified_hue_photo_ids")
    },sort_keys=True))
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mois-original-photo-archive",type=Path,required=True)
    for name in FCP_N:
        p.add_argument("--fcp-"+name,required=True,type=Path)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.mois_original_photo_archive,{name:getattr(a,"fcp_"+name) for name in FCP_N},a.output)
