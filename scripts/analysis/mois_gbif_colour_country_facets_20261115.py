#!/usr/bin/env python3
"""Real public GBIF COL XR georeferenced country/continent facets for 180 source-blind species.
A sampled species' photographed hue is NOT measured at GBIF occurrence locations.
These are country records, not native ranges or global 100km occupancy.
"""
from __future__ import annotations
import argparse,hashlib,io,json,time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import URLError,HTTPError
from zipfile import ZipFile
import pandas as pd

SHA_ZIP="abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a"
SHA_CSV="5a77f3c968cf221287cd98a855559b3c3e8a8cf1322c08add28878adeb293df8"
COLXR="7ddf754f-d193-4cc9-b351-99906754a03b"
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
SEED="MOIS-COLXR-20261115-COUNTRY20-EACH9-HUE-BLIND"
API="https://api.gbif.org/v1/occurrence/search"
def require(x,msg):
    if not bool(x):raise ValueError("MOIS_GBIF_REAL_COUNTRY_FACETS: "+msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def choose(df,n=20,start=0):
    require(n>=1 and n<=100 and start>=0 and start+n<=150,"GBIF sampling quota or disjoint offset out of bounds")
    d=df.copy()
    d["blind_rank"]=[sha((SEED+"|"+key).encode()) for key in d.key]
    d=d.sort_values(["blind_rank","key"],kind="stable")
    d["within_hue_rank"]=d.groupby("source_hue9",sort=True).cumcount()
    pick=d.loc[(d.within_hue_rank>=start)&(d.within_hue_rank<start+n)].drop(columns=["blind_rank","within_hue_rank"])
    require(len(pick)==9*n and pick.groupby("source_hue9").size().eq(n).all(),
            "nine-hue source quota incomplete")
    return pick.sort_values(["source_hue9","key"],kind="stable")

def roster(path,n=20,start=0):
    raw=Path(path).read_bytes();require(sha(raw)==SHA_ZIP,"source artifact ZIP digest drift")
    with ZipFile(io.BytesIO(raw)) as z:b=z.read("current_exact_accepted_COLXR_species.csv.gz")
    require(sha(b)==SHA_CSV,"accepted source CSV SHA changed")
    d=pd.read_csv(io.BytesIO(b),compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==37104 and d.key.is_unique and d.inat_taxon_id.is_unique,
            "accepted COL XR source species denominator wrong")
    require(d.status.eq("colxr_exact_accepted_species").all()
            and d.GBIF_kingdom.eq("Plantae").all()
            and d.candidate_rank.eq("SPECIES").all()
            and d.COL_XR_checklistKey.eq(COLXR).all()
            and set(d.source_hue9)==set(HUES),"COLXR accepted/species labels drift")
    return choose(d,n,start)
def extract(response,fields=("country","continent")):
    require(isinstance(response,dict) and isinstance(response.get("count"),int)
            and response["count"]>=0 and not response.get("results"),
            "GBIF did not return search-only facet metadata")
    f={str(v.get("field","")).lower():v for v in response.get("facets",[])}
    out={}
    for key in fields:
        require(key in f,"GBIF missing requested facet "+key)
        rows=f[key].get("counts",[])
        require(isinstance(rows,list) and len(rows)<200,
                "GBIF facet potentially capped; no truncated country-range inference")
        require(all(str(x.get("name","")) and type(x.get("count"))==int and x["count"]>=0 for x in rows),
                "invalid country/continent counts")
        require(len({str(x["name"]) for x in rows})==len(rows),"GBIF facet duplicate country")
        require(sum(x["count"] for x in rows)<=response["count"],
                "GBIF facet tally exceeds total occurrences")
        out[key]=rows
    return out
def fetch(key,*,pause=1.25):
    params={"checklistKey":COLXR,"taxonKey":str(key),"hasCoordinate":"true",
            "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT","limit":"0",
            "facet":["country","continent"],"facetLimit":"200"}
    url=API+"?"+urlencode(params,doseq=True)
    req=Request(url,headers={"Accept":"application/json",
       "User-Agent":"MOIS-9hue-GBIF-country-facets/20261115 (+github.com/zuizui0223/mois)"})
    for attempt in range(4):
        time.sleep(pause)
        try:
            with urlopen(req,timeout=40) as f:raw=f.read(4_000_000)
            require(len(raw)<4_000_000,"GBIF response size cap")
            obj=json.loads(raw)
            return url,obj,raw
        except HTTPError as e:
            if e.code not in (408,429,500,502,503,504):raise
            if attempt==3:raise
        except (URLError,TimeoutError):
            if attempt==3:raise
        time.sleep(min(16,2**attempt))
    raise RuntimeError("Unreachable GBIF retry exit")
def run(source,out,n=20,start=0):
    out=Path(out);require(not out.exists(),"do not overwrite a historical source")
    selected=roster(source,n,start)
    out.mkdir(parents=True);(out/"api_raw").mkdir()
    selected.to_csv(out/"source_selected_20_per_photo_hue.csv",index=False)
    records=[];country=[];continent=[];checks=[]
    for idx,row in enumerate(selected.itertuples(index=False)):
        url,response,raw=fetch(row.key)
        grouped=extract(response)
        cs=grouped["country"];ct=grouped["continent"]
        base={"inat_taxon_id":str(row.inat_taxon_id),"gbif_colxr_species_key":str(row.key),
              "species":row.species,"photo_hue9":row.source_hue9,
              "photo_hue_status":row.hue_status}
        records.append({**base,"GBIF_georeferenced_occurrence_total":int(response["count"]),
                        "GBIF_country_coded_occurrences":sum(v["count"] for v in cs),
                        "GBIF_continent_coded_occurrences":sum(v["count"] for v in ct),
                        "countries_GBIF_1plus":len(cs),
                        "countries_GBIF_3plus":sum(v["count"]>=3 for v in cs),
                        "countries_GBIF_10plus":sum(v["count"]>=10 for v in cs),
                        "continents_GBIF_1plus":len(ct),
                        "native_range_verified":False,
                        "photo_hue_measured_at_GBIF_record":False,
                        "species_wide_fixed_hue_proven":False})
        country.extend([{**base,"country":v["name"],"GBIF_georeferenced_records":v["count"]} for v in cs])
        continent.extend([{**base,"continent":v["name"],"GBIF_georeferenced_records":v["count"]} for v in ct])
        name="api_raw/"+str(idx).zfill(3)+"_"+row.key+".json"
        (out/name).write_bytes(raw)
        checks.append({"file":name,"sha256":sha(raw),"url":url,"GBIF_count":response["count"]})
        if idx%10==9 or idx==len(selected)-1:
            print("GBIF_HUE_COUNTRY_QUERY",idx+1,"of",len(selected),"last",row.key,flush=True)
            pd.DataFrame(records).to_csv(out/"partial_species_result.csv",index=False)
    r=pd.DataFrame(records)
    require(len(r)==n*9 and r.gbif_colxr_species_key.is_unique,"incomplete GBIF responses")
    r.to_csv(out/"species_GBIF_country_continent_facets.csv",index=False)
    pd.DataFrame(country).to_csv(out/"species_x_country_real_GBIF_records.csv",index=False)
    pd.DataFrame(continent).to_csv(out/"species_x_continent_real_GBIF_records.csv",index=False)
    (out/"requests_SHA_and_URL.json").write_text(json.dumps(checks,sort_keys=True,indent=2)+"\n")
    response={"schema":"mois_COLXR_nine_photo_hue_real_GBIF_country_facets_v1",
              "source_accepted_COLXR_species":37104,
              "sampled_accepted_species":len(r),"sample_n_per_photo_hue":n,"within_hue_rank_start":start,
              "source_blind_selection_seed":SEED,
              "GBIF_checklistKey":COLXR,
              "GBIF_filter":"hasCoordinate=true;hasGeospatialIssue=false;occurrenceStatus=PRESENT",
              "sampled_species_with_GBIF_georeferenced_records":int(r.GBIF_georeferenced_occurrence_total.gt(0).sum()),
              "georeferenced_occurrences_total_for_sample":int(r.GBIF_georeferenced_occurrence_total.sum()),
              "species_country_rows":len(country),"species_continent_rows":len(continent),
              "all_country_facets_complete_under_cap":True,
              "GBIF_authenticated_bulk_download_submitted":False,
              "GBIF_bulk_DOI_created":False,"native_status_verified":False,
              "photo_hue_is_occurrence_level_truth":False,
              "global_37104_species_distributions_completed":False,
              "output_SHA256":{p.name:sha(p.read_bytes()) for p in out.iterdir() if p.is_file()},
              "raw_response_SHA256":{x["file"]:x["sha256"] for x in checks}}
    (out/"receipt.json").write_text(json.dumps(response,sort_keys=True,indent=2)+"\n")
    print("MOIS_GBIF_REAL_COUNTRY_FACETS_DONE",json.dumps({k:v for k,v in response.items()
          if k not in ("output_SHA256","raw_response_SHA256")},sort_keys=True),flush=True)
    return response
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--per-hue",type=int,default=20)
    p.add_argument("--rank-start",type=int,default=0)
    x=p.parse_args();run(x.source,x.output,x.per_hue,x.rank_start)
