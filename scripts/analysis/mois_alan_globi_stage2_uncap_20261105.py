#!/usr/bin/env python3
"""Resume the actually capped GloBI pages using stage-1 signed source offsets.

GloBI returns limit+1 rows in many pages. The stage-1 cap count mistakenly used
equality to limit, so it recorded zero remaining caps. This repair determines
remaining families from the actual page records and continues from last offset.
No photo downloads, no new ROI classification, no inferred nocturnality.
"""
from __future__ import annotations
import argparse,concurrent.futures,json
from pathlib import Path
import pandas as pd
import requests
from scripts.analysis import mois_alan_expand_real_globi_and_photo_support as first
from scripts.analysis import mois_alan_globi_moth_plant_source_increment as old

SOURCE="6e5ba3ca7466677615d68ac8a5657a82aca12c32c4386970b9ab8830ab283f31"
PAGE="8c1185b9d43587e2eadd566eec3920699eb48d00ea0b711dd9d5d581f5aeeed1"

def find_capped(log):
    required={"family","relation","offset","limit","n_rows","failed"}
    first.require(required<=set(log),"Provider stage1 page log drift")
    tail=log.sort_values("offset").groupby(["family","relation"],as_index=False).tail(1)
    return tail.loc[(~tail.failed.astype(bool))&
                    (tail.n_rows>=tail["limit"]),["family","relation","offset","limit"]].to_dict("records")

def add_pages(previous,log,*,maximum=24,limit=1000,workers=5):
    first.require(1<=workers<=8 and maximum>=1 and limit>=250,"Invalid API acquisition settings")
    tasks=find_capped(log)
    first.require(len(tasks)>0,"No truly capped provider pages found")
    def get(t):
        session=requests.Session()
        collected=[];seen=[];finished=False
        start=int(t["offset"])+int(t["limit"])
        for j in range(maximum):
            values,meta,complete=first.extract_new_page(t["family"],t["relation"],
                      start+j*limit,limit,session)
            collected+=values;seen.append(meta)
            if meta["failed"] or complete:
                finished=True;break
        return collected,seen,{"family":t["family"],"relation":t["relation"],
          "pages":len(seen),"last_offset":seen[-1]["offset"],
          "last_n_rows":seen[-1]["n_rows"],
          "query_failed":bool(seen[-1]["failed"]),
          "still_capped":not finished and seen[-1]["n_rows"]>=limit}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        result=list(executor.map(get,tasks))
    more=pd.DataFrame([r for batch,_,_ in result for r in batch],columns=previous.columns)
    merged=pd.concat([previous,more],ignore_index=True).drop_duplicates()
    page_log=pd.DataFrame([p for _,batch,_ in result for p in batch])
    status=pd.DataFrame([s for _,_,s in result])
    return merged,more,page_log,status

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("prior","page_log","canonical","names","accepted","direct","output"):
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--max-pages",type=int,default=24)
    a=p.parse_args();out=a.output
    first.require(not out.exists(),"Never overwrite a signed source output")
    first.require(first.filesha(a.prior)==SOURCE and first.filesha(a.page_log)==PAGE,
                  "Previous GloBI stage source mismatch")
    previous=pd.read_csv(a.prior,compression="gzip")
    logs=pd.read_csv(a.page_log)
    combined,added,pages,status=add_pages(previous,logs,maximum=a.max_pages)
    for key in ("canonical","names","accepted","direct"):
        first.require(first.filesha(getattr(a,key))==first.SHA[key],
                      "Previous source hash drift: "+key)
    x=pd.read_csv(a.canonical,compression="gzip",low_memory=False,
         dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
         usecols=["photo_id","observation_id","inat_taxon_id",
           "strict9","strict9_white","latitude","longitude","observed_on","source_subarm"])
    first.require(len(x)==1032882 and x.photo_id.is_unique and x.observation_id.is_unique,
                  "Source photo/observation identities drifted")
    names=pd.read_csv(a.names,compression="gzip",dtype=str)
    accepted=pd.read_csv(a.accepted,compression="gzip",dtype=str,
                   usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    taxonomy=old.source_taxonomy(names,accepted)
    direct=pd.read_csv(a.direct,compression="gzip",dtype=str,
          usecols=["accepted_species","trait_name","normalized_value","quality",
                   "resolution_status","classification"])
    old_stats,_=first.assess(previous,taxonomy,direct,x)
    new_stats,_=first.assess(combined,taxonomy,direct,x)
    compare=new_stats.merge(old_stats,on=["period","grid_degrees"],
                 validate="one_to_one",suffixes=("_new","_old"))
    out.mkdir(parents=True)
    combined.to_csv(out/"GloBI_stage2_signed_extension.csv.gz",index=False,
           compression={"method":"gzip","mtime":0})
    pages.to_csv(out/"GloBI_stage2_provider_page_SHA256.csv",index=False)
    status.to_csv(out/"GloBI_stage2_remaining_pagination.csv",index=False)
    compare.to_csv(out/"source_same_photo_colour_support_stage1_vs_stage2.csv",index=False)
    report={"schema":"MOIS_ALAN_GloBI_stage2_precise_pagination_20261105",
       "input_stage1_source_SHA256":SOURCE,"input_stage1_page_SHA256":PAGE,
       "original_stage1_interaction_rows":len(previous),
       "new_validated_provider_rows_before_dedup":len(added),
       "new_unique_interaction_rows":len(combined)-len(previous),
       "combined_unique_interaction_rows":len(combined),
       "initial_truncated_queries_corrected":len(find_capped(logs)),
       "additional_GloBI_page_requests":len(pages),
       "provider_errors":int(pages.failed.sum()),
       "remaining_truncated_queries":int(status.still_capped.sum()),
       "original_photo_source_unchanged":True,
       "source_nocturnality_or_pollination_effectiveness_not_verified":True,
       "no_record_not_ecological_absence":True,
       "independent_species_colour_unknown_not_imputed":True,
       "photo_image_downloads_or_inferences":0,
       "MNAR_scenarios_run":False,
       "output_SHA256":{f.name:first.filesha(f) for f in out.iterdir()}}
    (out/"source_stage2_receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_ALAN_GLOBI_PAGINATED_STAGE2",json.dumps({
        k:report[k] for k in ("original_stage1_interaction_rows","new_unique_interaction_rows",
        "combined_unique_interaction_rows","initial_truncated_queries_corrected",
        "remaining_truncated_queries")},sort_keys=True),flush=True)
    for deg in (.25,.5,1.):
        row=compare.loc[compare.period.eq("2019-2024")&compare.grid_degrees.eq(deg)].iloc[0]
        print("MOIS_ALAN_EVIDENCE_FOUR_CLASS_SUPPORT",json.dumps({
         "degrees":deg,
         "old_four_ge1":int(row.cells_four_direct_colour_and_visit_groups_ge1_old),
         "new_four_ge1":int(row.cells_four_direct_colour_and_visit_groups_ge1_new),
         "new_four_ge2":int(row.cells_four_direct_colour_and_visit_groups_ge2_new),
         "new_four_ge3":int(row.cells_four_direct_colour_and_visit_groups_ge3_new)},sort_keys=True),flush=True)

if __name__=="__main__":main()
