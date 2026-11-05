#!/usr/bin/env python3
"""GloBI moth-family visitor evidence × MOIS botanical colour coverage: no photo work.

The live GloBI API is an EXPLORATORY, CAPPED source, never the frozen primary
species-level pollination truth. An actual moth visit does not prove nocturnal
pollination; a butterfly visit does not prove absence of moths.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests

MOTH_FAMILIES=("Sphingidae","Noctuidae","Geometridae","Erebidae","Crambidae","Pyralidae","Hepialidae")
BUTTERFLY_FAMILIES=("Papilionidae","Nymphalidae","Pieridae","Lycaenidae","Hesperiidae")
RELATIONS=("visitsFlowersOf","pollinates")
ENDPOINT="https://api.globalbioticinteractions.org/interaction.csv"
CANONICAL_SHA="f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c"
WCVP_SHA="9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3"
NAMES_SHA="dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2"
ISLAND_SHA="7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577"
CELL_SIZES=(.1,.25,.5,1.0)
BINOMIAL=re.compile(r"^[A-Z][A-Za-z-]+ [a-z][a-zA-Z-]+$")


def check(ok,msg):
    if not bool(ok):
        raise ValueError("MOIS_GLOBI_GUILD_INCREMENT: "+msg)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for data in iter(lambda:stream.read(1<<20),b""):h.update(data)
    return h.hexdigest()


def clean_record(row,family,relation):
    """Strict taxonomic path and interaction orientation checks, incl. API leaks."""
    check(relation in RELATIONS,"Not a flower interaction")
    lineage=[p.strip() for p in str(row.get("source_taxon_path","")).split("|")]
    if family not in lineage or "Lepidoptera" not in lineage:
        return None
    botanical=[p.strip().casefold() for p in str(row.get("target_taxon_path","")).split("|")]
    if not any(p in botanical for p in ("plantae","viridiplantae","embryophyta","angiosperms","streptophyta","magnoliopsida")):
        return None
    name=str(row.get("target_taxon_name","")).strip()
    if not BINOMIAL.fullmatch(name):
        return None
    insect=str(row.get("source_taxon_name","")).strip()
    study=str(row.get("study_title","")).strip()
    if not study or len(insect.split())<2:
        return None
    if str(row.get("interaction_type",""))!=relation:
        return None
    return {"plant_name":name,"guild":"moth_family_visit" if family in MOTH_FAMILIES else "butterfly_family_visit",
            "source_family":family,"relation":relation,"insect_name":insect,
            "study_title":study[:500]}


def query_family(family,relation,session,*,limit=250,max_pages=8):
    check(family in MOTH_FAMILIES+BUTTERFLY_FAMILIES,"Unknown taxon family")
    records=[];receipt=[];failed=False;truncated=False
    for page in range(max_pages):
        params={"sourceTaxon":family,"interactionType":relation,
                "includeObservations":"true","limit":str(limit),"offset":str(page*limit)}
        answer=None
        for attempt in range(3):
            try:
                answer=session.get(ENDPOINT,params=params,timeout=40,
                        headers={"User-Agent":"MOIS-Floral-ALAN-GloBI-Evidence/1.0"})
                if answer.status_code==200:break
            except requests.RequestException:pass
            time.sleep(.3*(attempt+1))
        if answer is None or answer.status_code!=200:
            failed=True;receipt.append({"family":family,"relation":relation,"page":page,
                                        "status":getattr(answer,"status_code",None),"source_error":"UNAVAILABLE"})
            break
        text=answer.text
        reader=csv.DictReader(io.StringIO(text))
        cols=set(reader.fieldnames or [])
        check({"source_taxon_name","source_taxon_path","target_taxon_name","target_taxon_path",
               "interaction_type","study_title"}<=cols,"GloBI source schema drift")
        page_records=list(reader)
        receipt.append({"family":family,"relation":relation,"page":page,
                        "http":answer.status_code,"sha256":hashlib.sha256(answer.content).hexdigest(),
                        "n_provider_rows":len(page_records),
                        "no_provider_release_pin":True})
        for row in page_records:
            candidate=clean_record(row,family,relation)
            if candidate is not None:records.append(candidate)
        if len(page_records)<limit:
            break
        if page==max_pages-1:
            truncated=True
    return records,receipt,{"family":family,"relation":relation,"verified_plant_moth_links":len(records),
                             "query_failed":failed,"page_cap_possible_truncation":truncated,
                             "n_pages":len(receipt)}


def source_taxonomy(names,accepted):
    key=names.loc[names.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
        ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    check(key.inat_taxon_id.is_unique,"Original taxon ID duplicate")
    a=accepted[["wcvp_accepted_plant_name_id","wcvp_accepted_name"]].drop_duplicates()
    check(a.groupby("wcvp_accepted_plant_name_id").wcvp_accepted_name.nunique().le(1).all(),
          "WCVP accepted ID name conflict")
    a=a.drop_duplicates("wcvp_accepted_plant_name_id")
    return key.merge(a,on="wcvp_accepted_plant_name_id",how="left",validate="many_to_one")


def coverage(globi,canonical,taxonomy,direct):
    check(not canonical.photo_id.duplicated().any() and not canonical.observation_id.duplicated().any(),
          "Repeated source photo or ecological observation")
    source=canonical.loc[canonical.source_subarm.eq("historical_main") &
           canonical.strict9 &
           canonical.observed_on.astype(str).str[:4].between("2019","2024") &
           canonical.latitude.between(-60,60) & canonical.longitude.between(-180,180),
           ["photo_id","observation_id","inat_taxon_id","observed_on","latitude","longitude","original_roi_four_hue"]].copy()
    source=source.merge(taxonomy,on="inat_taxon_id",how="left",validate="many_to_one")
    check(source.photo_id.is_unique,"Source taxonomy caused repeated photo IDs")
    direct=direct.loc[(direct.trait_name=="flower_primary_color")&
                 direct.quality.isin(("high","medium"))&
                 direct.resolution_status.eq("resolved")&
                 ~direct.classification.isin(("cultivar_contamination","partial_self_compatibility")),
                 ["accepted_species","normalized_value"]].copy()
    check(direct.accepted_species.is_unique,"Direct colour species duplicate")
    source=source.merge(direct,left_on="wcvp_accepted_name",right_on="accepted_species",how="left",validate="many_to_one")
    documented=np.where(source.normalized_value.eq("white"),"white",
               np.where(source.normalized_value.isin(("yellow_orange","red_pink",
                    "blue_purple","green_brown_inconspicuous")),"nonwhite","unknown"))
    source["direct_colour"]=documented
    source["year"]=source.observed_on.astype(str).str[:4]
    candidates=globi.drop_duplicates(["plant_name","guild"])[["plant_name","guild"]]
    # Species may have both butterfly and moth visits. Keep both and never turn
    # "no GloBI record" into diurnal/moth absence.
    selected=source.merge(candidates,left_on="wcvp_accepted_name",
                          right_on="plant_name",how="inner",validate="many_to_many")
    unique=selected[["wcvp_accepted_name","guild"]].drop_duplicates()
    yearly=[]
    for guild,d in selected.groupby("guild"):
        vals=d.loc[d.direct_colour.eq("white")]
        yearly.append({"guild":guild,"distinct_accepted_species":d.wcvp_accepted_name.nunique(),
                       "original_strict9_photo_records":len(d),
                       "independent_direct_white_species":vals.wcvp_accepted_name.nunique(),
                       "independent_direct_white_photo_records":len(vals),
                       "direct_nonwhite_species":d.loc[d.direct_colour.eq("nonwhite"),"wcvp_accepted_name"].nunique()})
    support=[]
    for label,part in [("2020",selected.loc[selected.year.eq("2020")]),
                       ("2024",selected.loc[selected.year.eq("2024")]),
                       ("2019-2024 pooled",selected)]:
        for deg in CELL_SIZES:
            t=part.loc[part.direct_colour.ne("unknown")].copy()
            if t.empty:
                support.append({"period":label,"grid_degrees":deg,"cells_any_both_white_groups":0,
                                "cells_both_white_guilds_ge3_species":0,"cells_four_groups_ge1_species":0,
                                "cells_four_groups_ge3_species":0})
                continue
            t["cell"]=(np.floor((t.longitude+180)/deg).astype(int).astype(str)+"_"+
                       np.floor((t.latitude+90)/deg).astype(int).astype(str))
            ct=t.groupby(["cell","guild","direct_colour"]).wcvp_accepted_name.nunique().unstack(["guild","direct_colour"]).fillna(0)
            def count(g,w):return ct.get((g,w),pd.Series(0,index=ct.index))
            a=count("moth_family_visit","white");b=count("butterfly_family_visit","white")
            c=count("moth_family_visit","nonwhite");d=count("butterfly_family_visit","nonwhite")
            support.append({"period":label,"grid_degrees":deg,
               "cells_any_both_white_groups":int(((a>=1)&(b>=1)).sum()),
               "cells_both_white_guilds_ge3_species":int(((a>=3)&(b>=3)).sum()),
               "cells_four_groups_ge1_species":int(((a>=1)&(b>=1)&(c>=1)&(d>=1)).sum()),
               "cells_four_groups_ge3_species":int(((a>=3)&(b>=3)&(c>=3)&(d>=3)).sum())})
    return pd.DataFrame(yearly),pd.DataFrame(support),{
        "source_historical_photo_count_2019_2024":len(source),
        "source_moth_or_butterfly_linked_photo_count":len(selected),
        "source_accepted_species_with_linked_flower_interaction":int(unique.wcvp_accepted_name.nunique()),
        "species_linked_to_both_guilds_not_assumed_day_only":
           int(unique.groupby("wcvp_accepted_name").guild.nunique().ge(2).sum()),
        "no_white_imputation":True,"no_moth_absence_from_no_record":True,
        "no_nocturnal_pollination_truth_inferred":True
    }


def run(canonical,names,accepted,direct,out,*,page_size=250,max_pages=8):
    files={"canonical":(Path(canonical),CANONICAL_SHA),
           "names":(Path(names),NAMES_SHA),"accepted":(Path(accepted),WCVP_SHA),
           "direct":(Path(direct),ISLAND_SHA)}
    for key,(path,digest) in files.items():
        check(sha(path)==digest,"Input SHA drift: "+key)
    out=Path(out);check(not out.exists(),"Refuse to overwrite old external data")
    session=requests.Session()
    rows=[];receipts=[];quality=[]
    for family in MOTH_FAMILIES+BUTTERFLY_FAMILIES:
        for relation in RELATIONS:
            events,checks,q=query_family(family,relation,session,limit=page_size,max_pages=max_pages)
            rows.extend(events);receipts.extend(checks);quality.append(q)
            print("MOIS_GLOBI_FAMILY_SOURCE",json.dumps(q,sort_keys=True),flush=True)
    check(len(rows)>0,"No actual external moth/butterfly flower interactions returned")
    records=pd.DataFrame(rows).drop_duplicates()
    check(records.plant_name.str.fullmatch(BINOMIAL).all(),"Invalid botanical binomial")
    c=pd.read_csv(files["canonical"][0],compression="gzip",low_memory=False,
             dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str})
    n=pd.read_csv(files["names"][0],compression="gzip",dtype=str)
    a=pd.read_csv(files["accepted"][0],compression="gzip",dtype=str,
                   usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    t=pd.read_csv(files["direct"][0],compression="gzip",dtype=str,
                   usecols=["accepted_species","trait_name","normalized_value",
                            "quality","resolution_status","classification"])
    taxonomy=source_taxonomy(n,a)
    totals,support,meta=coverage(records,c,taxonomy,t)
    out.mkdir(parents=True)
    # Publish only ecological aggregates and source-query hashes; source original
    # study-level claim rows remain local and are NOT auto published to repo.
    totals.to_csv(out/"guild_record_photo_coverage.csv",index=False)
    support.to_csv(out/"guild_record_spatial_support.csv",index=False)
    # Preserve only public GloBI interaction source claims, never the restricted
    # Island species-colour rows. This frozen cache prevents repeat API queries.
    records.to_csv(out/"GloBI_bounded_interaction_source_cache.csv.gz",index=False,
                   compression={"method":"gzip","mtime":0})
    pd.DataFrame(quality).to_csv(out/"source_query_status.csv",index=False)
    pd.DataFrame(receipts).to_csv(out/"source_query_SHA_pages.csv",index=False)
    result={"schema":"MOIS_ALAN_GloBI_live_API_moth_family_source_feasibility_v1",
            "source":"GloBI live API on analysis date; exploratory NOT version-pinned confirmatory",
            "search_families_moth":MOTH_FAMILIES,"search_families_butterfly":BUTTERFLY_FAMILIES,
            "relation_types":RELATIONS,
            "validated_record_count":len(records),"number_of_source_queries":len(quality),
            "some_queries_failed":any(q["query_failed"] for q in quality),
            "some_queries_may_be_truncated":any(q["page_cap_possible_truncation"] for q in quality),
            "not_all_insect_family_names_in_API_results_are_correct":True,
            "record_taxon_lineage_checked_at_species_level":True,
            "plant_species_photo_colour_not_relabelled":True,
            "moth_family_visit_not_documented_effective_nocturnal_pollination":True,
            "butterfly_family_visit_not_proven_no_moth_visit":True,
            "API_data_not_original_stable_GloBI_0p9_source":True,
            "MOIS_source_SHAs":{k:v[1] for k,v in files.items()},
            **meta,
            "GloBI_source_claims_cached_in_GitHub_artifact":True,
            "source_examples_not_raw_study_records_published":False,
            "restricted_Island_species_colour_rows_published":False,
            "new_original_image_downloads":0,"new_image_inferences":0,
            "MNAR_sensitivity_scenarios_run":False,
            "result_SHA256":{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},
           }
    (out/"guild_source_integration_receipt.json").write_text(json.dumps(result,indent=2,
                ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_GLOBI_ADDITIONAL_TRUE_FLOWER_GUILD_SOURCE_RESULT",
          json.dumps({k:result[k] for k in ("validated_record_count",
          "source_accepted_species_with_linked_flower_interaction",
          "source_moth_or_butterfly_linked_photo_count",
          "some_queries_may_be_truncated","some_queries_failed")},sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("canonical","names","accepted","direct","out"):p.add_argument("--"+arg,required=True)
    p.add_argument("--page-size",type=int,default=250)
    p.add_argument("--max-pages",type=int,default=8)
    args=p.parse_args();run(args.canonical,args.names,args.accepted,args.direct,args.out,
        page_size=args.page_size,max_pages=args.max_pages)


if __name__=="__main__":main()
