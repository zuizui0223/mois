#!/usr/bin/env python3
"""Incrementally extend capped GloBI moth/butterfly pages; exact MOIS source photo IDs.

Existing signed data only for flower photographs. A new API record is documented
flower INTERACTION, not proof of nocturnal effective pollination or absence.
All source pages hashed. Source-bound species colors are never redistributed.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,time,concurrent.futures
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from scripts.analysis import mois_alan_globi_moth_plant_source_increment as old

SHA={
 "globi":"6116acab80f2277af7ec109715e5815bf6c4461767ae61580da8d842c9e230e4",
 "status":"183d1a2da5a5bf159876e77441a292218ca2f9bab202e9f5feb68ee8e24b9eda",
 "canonical":old.CANONICAL_SHA,
 "names":old.NAMES_SHA,
 "accepted":old.WCVP_SHA,
 "direct":old.ISLAND_SHA,
 "viirs2020":"3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3",
}
# Actual status SHA pinned here, no silent fallback to a different source.
SHA["status"]="183d1a2da5a5bf159876e77441a292218ca2f9bab202e9f5feb68ee8e24b9eda"
SIZES=(.05,.1,.25,.5,1.)
PERIODS=("2020","2024","2019-2024")
GUILDS=("moth_only","butterfly_only")

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_ALAN_MORE_DATA: "+msg)

def filesha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def classify_groups(links):
    require({"plant_name","guild"}<=set(links),"GloBI schema changed")
    x=links[["plant_name","guild"]].drop_duplicates()
    require(x.guild.isin(["moth_family_visit","butterfly_family_visit"]).all(),
            "Encountered unverified interaction guild")
    counts=x.assign(value=1).pivot(index="plant_name",columns="guild",values="value").fillna(0)
    m=counts.get("moth_family_visit",pd.Series(0,index=counts.index)).astype(bool)
    b=counts.get("butterfly_family_visit",pd.Series(0,index=counts.index)).astype(bool)
    require((m|b).all(),"No guild source evidence")
    return pd.Series(np.select([m&~b,b&~m,m&b],["moth_only","butterfly_only","both"],default="unknown"),
                    index=counts.index,name="source_evidence_group")

def extract_new_page(source_family,rel,offset,limit,session):
    url=old.ENDPOINT
    params={"sourceTaxon":source_family,"interactionType":rel,
            "includeObservations":"true","offset":str(offset),"limit":str(limit)}
    for attempt in range(3):
        try:
            res=session.get(url,params=params,timeout=45,headers={"User-Agent":"MOIS-ALAN-AdditionalGloBI-LimitedStableAudit/2.0"})
            if res.status_code!=200:
                if attempt<2 and res.status_code in (429,500,502,503,504):
                    time.sleep(.7*(attempt+1));continue
                return [],{"family":source_family,"relation":rel,"offset":offset,
                           "status":res.status_code,"failed":True,"n_rows":0},False
            reader=csv.DictReader(io.StringIO(res.text))
            expected={"source_taxon_name","source_taxon_path","target_taxon_name","target_taxon_path",
                      "interaction_type","study_title"}
            require(expected<=set(reader.fieldnames or ()),"Unexpected GloBI header")
            rows=list(reader)
            accepted=[]
            for r in rows:
                val=old.clean_record(r,source_family,rel)
                if val is not None:accepted.append(val)
            info={"family":source_family,"relation":rel,"offset":offset,
                  "limit":limit,"status":200,"sha256":hashlib.sha256(res.content).hexdigest(),
                  "n_rows":len(rows),"n_accepted_before_dedup":len(accepted),"failed":False}
            return accepted,info,len(rows)<limit
        except (requests.RequestException,UnicodeDecodeError,ValueError) as exc:
            if attempt>=2:
                return [],{"family":source_family,"relation":rel,"offset":offset,
                           "status":None,"failed":True,"error":str(exc)[:250],"n_rows":0},False
            time.sleep(.7*(attempt+1))
    raise AssertionError("Retry exhausted")

def extend(original,status,*,limit=1000,max_additional_pages=12,workers=6):
    """Only requests the pages beyond the previous 8*250 source windows."""
    require(limit>=250 and max_additional_pages>=1 and 1<=workers<=8,"Invalid courteous bounded API request")
    require(set(old.MOTH_FAMILIES+old.BUTTERFLY_FAMILIES)==set(status.family),"Previous family list changed")
    todo=status.loc[status.page_cap_possible_truncation.astype(bool),
                    ["family","relation"]].drop_duplicates().to_dict("records")
    require(len(todo)>0,"No genuinely capped GloBI queries to extend")
    def task(t):
        session=requests.Session()
        records=[];receipts=[];end=False
        for page in range(max_additional_pages):
            offset=2000+page*limit
            rs,meta,finished=extract_new_page(t["family"],t["relation"],offset,limit,session)
            records+=rs;receipts.append(meta)
            if meta["failed"] or finished:
                end=True;break
        return records,receipts,{"family":t["family"],"relation":t["relation"],
                "requested_extra_pages":len(receipts),"reached_end_or_failed":end,
                "truncated_after_new_page_cap":not end and
                 bool(receipts and receipts[-1].get("n_rows",0)==limit),
                "provider_failed":bool(receipts and receipts[-1]["failed"])}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        results=list(pool.map(task,todo))
    more=[r for r,_,_ in results]
    received=[p for _,rr,_ in results for p in rr]
    done=[s for _,_,s in results]
    added=pd.DataFrame([x for batch in more for x in batch],
        columns=list(original.columns))
    combined=pd.concat([original,added],ignore_index=True).drop_duplicates()
    require(combined.plant_name.notna().all() and
            combined.guild.isin(("moth_family_visit","butterfly_family_visit")).all(),
            "Moth and butterfly study evidence source corrupted")
    return combined,added,pd.DataFrame(received),pd.DataFrame(done)

def assess(s,source_names, direct,positions):
    x=positions.copy()
    x=x.loc[x.source_subarm.eq("historical_main") & x.strict9 &
            x.latitude.between(-60,60) & x.longitude.between(-180,180) &
            x.observed_on.astype(str).str[:4].between("2019","2024")].copy()
    require(not x.photo_id.duplicated().any(),"Photo ID duplicated before lineage join")
    x=x.merge(source_names,on="inat_taxon_id",how="left",validate="many_to_one")
    x["group"]=x.wcvp_accepted_name.map(classify_groups(s))
    valid=direct.loc[direct.trait_name.eq("flower_primary_color") &
         direct.quality.isin(("high","medium")) &
         direct.resolution_status.eq("resolved") &
         ~direct.classification.isin(("cultivar_contamination","partial_self_compatibility")),
         ["accepted_species","normalized_value"]]
    require(valid.accepted_species.is_unique,"Direct species colour duplicated")
    x=x.merge(valid,left_on="wcvp_accepted_name",right_on="accepted_species",
              how="left",validate="many_to_one")
    mapping={"white":1,"yellow_orange":0,"red_pink":0,
             "blue_purple":0,"green_brown_inconspicuous":0}
    x["direct_white"]=x.normalized_value.map(mapping)
    x["year"]=x.observed_on.astype(str).str[:4]
    x=x.loc[x.group.isin(GUILDS)].copy()
    summaries=[]
    for period in PERIODS:
        y=x if period=="2019-2024" else x.loc[x.year.eq(period)]
        for deg in SIZES:
            t=y.copy()
            t["cell"]=(np.floor((t.longitude+180)/deg).astype(int).astype(str)+"_"+
                       np.floor((t.latitude+90)/deg).astype(int).astype(str))
            t=t.groupby(["cell","wcvp_accepted_name","group"],observed=True,as_index=False).agg(
                 direct_white=("direct_white","first"))
            counts=t.groupby(["cell","group"],observed=True).wcvp_accepted_name.nunique().unstack(fill_value=0)
            for k in GUILDS:
                if k not in counts:counts[k]=0
            w=t.loc[t.direct_white.notna()].pivot_table(index="cell",columns=["group","direct_white"],
                  values="wcvp_accepted_name",aggfunc="nunique",fill_value=0)
            def g(label,white):
                return w.get((label,float(white)),pd.Series(0,index=w.index))
            a,b,c,d=g("moth_only",1),g("butterfly_only",1),g("moth_only",0),g("butterfly_only",0)
            summary={"period":period,"grid_degrees":deg,
                     "source_photo_rows_with_exclusive_visit_evidence":int(len(y)),
                     "cells_with_both_visit_groups_ge1":int(((counts.moth_only>=1)&(counts.butterfly_only>=1)).sum()),
                     "cells_with_both_visit_groups_ge3":int(((counts.moth_only>=3)&(counts.butterfly_only>=3)).sum()),
                     "cells_both_groups_independently_documented_white_ge1":int(((a>=1)&(b>=1)).sum()),
                     "cells_four_direct_colour_and_visit_groups_ge1":int(((a>=1)&(b>=1)&(c>=1)&(d>=1)).sum()),
                     "cells_four_direct_colour_and_visit_groups_ge2":int(((a>=2)&(b>=2)&(c>=2)&(d>=2)).sum()),
                     "cells_four_direct_colour_and_visit_groups_ge3":int(((a>=3)&(b>=3)&(c>=3)&(d>=3)).sum())}
            summaries.append(summary)
    return pd.DataFrame(summaries),{
        "source_photo_rows":len(x),
        "strict9_sampled_accepted_species_with_exclusive_visit_evidence":
          int(x.wcvp_accepted_name.nunique()),
        "white_direct_species_in_both_sides":
          int(x.loc[x.direct_white.eq(1)].wcvp_accepted_name.nunique()),
        "unknown_species_colour_not_imputed":True,
        "recorded_exclusive_visit_evidence_NOT_true_exclusive_pollination":True,
        "nocturnal_time_of_day_unverified":True,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("prior","status","canonical","names","accepted","direct","viirs2020","out"):
        p.add_argument("--"+name.replace("_","-"),dest=name,required=True,type=Path)
    p.add_argument("--extra-pages",type=int,default=12)
    p.add_argument("--page-limit",type=int,default=1000)
    p.add_argument("--workers",type=int,default=6)
    args=p.parse_args()
    out=args.out
    require(not out.exists(),"Refuse overwrite")
    paths={"globi":args.prior,"status":args.status,"canonical":args.canonical,
           "names":args.names,"accepted":args.accepted,
           "direct":args.direct,"viirs2020":args.viirs2020}
    observed={k:filesha(p) for k,p in paths.items()}
    require(all(observed[k]==SHA[k] for k in SHA),"Frozen science input digest disagrees")
    prev=pd.read_csv(args.prior,compression="gzip")
    stat=pd.read_csv(args.status)
    merged,new,requests_log,done=extend(prev,stat,limit=args.page_limit,
                max_additional_pages=args.extra_pages,workers=args.workers)
    # Acceptance criteria: real provider pages must actually return accepted,
    # source-provenance plant-moth edges.  Never invent more data.
    require((~requests_log.failed).any(),"Provider never returned an additional valid page")
    names=pd.read_csv(args.names,compression="gzip",dtype=str)
    accepted=pd.read_csv(args.accepted,compression="gzip",dtype=str,
             usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    mapping=old.source_taxonomy(names,accepted)
    photo=pd.read_csv(args.canonical,compression="gzip",low_memory=False,
             usecols=["photo_id","observation_id","inat_taxon_id","source_subarm",
                      "strict9","strict9_white","observed_on","latitude","longitude"],
             dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str})
    require(len(photo)==1032882 and int(photo.strict9.sum())==401703 and
            photo.photo_id.is_unique and photo.observation_id.is_unique,
            "Frozen original photo identity changed")
    traits=pd.read_csv(args.direct,compression="gzip",dtype=str,
             usecols=["accepted_species","trait_name","normalized_value",
                      "quality","resolution_status","classification"])
    old_support,oldmeta=assess(prev,mapping,traits,photo)
    new_support,newmeta=assess(merged,mapping,traits,photo)
    comparisons=new_support.merge(old_support,on=["period","grid_degrees"],suffixes=("_new","_old"))
    out.mkdir(parents=True)
    merged.to_csv(out/"GloBI_extended_source_cache.csv.gz",index=False,
                 compression={"method":"gzip","mtime":0})
    pd.DataFrame(requests_log).to_csv(out/"additional_GloBI_pages_sha256.csv",index=False)
    pd.DataFrame(done).to_csv(out/"additional_source_retrieval_status.csv",index=False)
    comparisons.to_csv(out/"original_photo_support_prior_vs_extended.csv",index=False)
    pgroup=classify_groups(prev).value_counts().to_dict()
    ngroup=classify_groups(merged).value_counts().to_dict()
    report={
      "schema":"mois_ALAN_more_actual_GloBI_species_source_20261105_v1",
      "input_sha256":observed,"original_GloBI_links":len(prev),
      "new_validated_source_rows_before_dedup":len(new),
      "combined_distinct_GloBI_study_links":len(merged),
      "new_distinct_source_records":len(merged)-len(prev),
      "source_name_group_prior":pgroup,
      "source_name_group_extended":ngroup,
      "source_api_pages_attempted":len(requests_log),
      "source_api_pages_failed":int(requests_log.failed.sum()),
      "source_queries_still_capped":int(done.truncated_after_new_page_cap.sum()),
      "source_WCVP_and_photo_source_unchanged":True,
      "photo_pixels_downloaded_or_inferred":0,
      "existing_strictROI9_success_unchanged":401703,
      "documented_colour_missing_never_nonwhite":True,
      "moth_visitor_NOT_effective_nocturnal_pollination":True,
      "butterfly_without_moth_record_NOT_verified_diurnal_only":True,
      "API_source_is_unversioned_2026_11_05_exploratory":True,
      "no_MNAR_sensitivity_scenario":True,
      "old_observation_support":oldmeta,
      "new_observation_support":newmeta,
      "outputs_SHA256":{p.name:filesha(p) for p in sorted(out.iterdir())}}
    (out/"extended_source_receipt.json").write_text(
      json.dumps(report,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_GLOBI_SOURCE_EXPANDED",json.dumps({
      "old":len(prev),"new":len(merged),"new_rows":len(merged)-len(prev),
      "source_queries_still_capped":report["source_queries_still_capped"],
      "white_fourway_25km_old":int(old_support.loc[
         (old_support.period=="2019-2024")&(old_support.grid_degrees==.25),
         "cells_four_direct_colour_and_visit_groups_ge1"].iloc[0]),
      "white_fourway_25km_new":int(new_support.loc[
         (new_support.period=="2019-2024")&(new_support.grid_degrees==.25),
         "cells_four_direct_colour_and_visit_groups_ge1"].iloc[0])},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
