#!/usr/bin/env python3
"""One-time audit of 38 previously failed photo IDs using the unchanged FCP ROI-v4 method.

The frozen algorithm receives provider LARGE derivatives, not full-resolution originals.
This file never updates the main 42,111-taxon atlas; it outputs only non-promoted candidates.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,re,zipfile
from pathlib import Path
from urllib.parse import urlsplit
import pandas as pd

PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
HUES={"white","yellow_orange","red_pink","blue_purple"}
SHA={
    "qa_zip":"e585e3052028e9a8a998c9f7eafaef56cd79e350f5dfaf1cee58497d80252848",
    "identities_zip":"60e8f57896449d614f98923fcf326b2bc850f654e9bbce92137efced6e6b79fe",
    "latest_zip":"2686b566485c8f13593e4ba463ae4b84d6d3b497c37ba6eabee4e0a9bd6a1696",
    "latest_csv":"6d08bfc203f317f1fc4a8d3d5f22015a7d0b393cba8ce9fcb18dd0ea74fe1909",
    "original_parquet":"3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835",
    "source_shard":"6c035941f6b7b77e0369ec8bbb62f6f1abf8046ad0187634a284d9aabfb9ce4b",
    "triage_git_blob":"f44e054bf6b6ebc60ba0ac69b80abd5aa67a9d34",
}

def require(cond,msg):
    if not bool(cond):
        raise ValueError("MOIS_VISUAL38_FROZEN_REMEASURE: "+msg)

def sha(blob):
    if isinstance(blob,(str,Path)):blob=Path(blob).read_bytes()
    return hashlib.sha256(blob).hexdigest()

def member(p,expected_sha,name):
    p=Path(p)
    require(sha(p)==expected_sha,"signed source ZIP digest changed: "+str(p))
    with zipfile.ZipFile(p) as z: return z.read(name)

def source_large_url(medium,pid):
    a=urlsplit(medium)
    require(a.scheme=="https" and a.netloc=="inaturalist-open-data.s3.amazonaws.com"
            and not a.query and not a.fragment,"untrusted photo source URL")
    require(re.fullmatch(r"/photos/"+re.escape(str(pid))+r"/medium\.(?:jpg|jpeg|png)",a.path),
            "medium URL disagrees with photo ID")
    return medium.replace("/medium.","/large.")

def read_inputs(input_dir):
    folder=Path(input_dir)
    a=member(folder/"latest.zip",SHA["latest_zip"],"mois_42111_post20261114_new_independent_original_photo9.csv.gz")
    require(sha(a)==SHA["latest_csv"],"latest 42,111 atlas changed")
    atlas=pd.read_csv(io.BytesIO(a),compression="gzip",dtype={"inat_taxon_id":str},low_memory=False).fillna("")
    require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique,"source taxon atlas missing")
    raw=member(folder/"qa.zip",SHA["qa_zip"],"image_source_manifest.csv")
    source=pd.read_csv(io.BytesIO(raw),dtype=str).fillna("")
    identity=member(folder/"identities.zip",SHA["identities_zip"],"all_1181146_photo_identity_records.csv.gz")
    with zipfile.ZipFile(folder/"identities.zip") as z:
        receipt=json.loads(z.read("receipt.json"))
    require(receipt["source_parquet_sha256"]==SHA["original_parquet"]
            and sha(identity)==receipt["output_sha256"]["all_1181146_photo_identity_records.csv.gz"],
            "original 1.18m identity export altered")
    photos=pd.read_csv(io.BytesIO(identity),compression="gzip",dtype={"photo_id":str,"inat_taxon_id":str},
                       usecols=["photo_id","inat_taxon_id","strict9"])
    require(len(photos)==1181146 and photos.photo_id.is_unique and int(photos.strict9.sum())==451643,
            "historical photo/strict9 identity records changed")
    return atlas,source,photos

def original_triage(path):
    raw=Path(path).read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    require(blob==SHA["triage_git_blob"],"90-photo source classification changed")
    a=pd.read_csv(io.BytesIO(raw),dtype=str).fillna("")
    a=a.loc[a.route.eq("petaloid_roi_recheck")].copy()
    require(len(a)==38 and a.photo_id.is_unique and a.inat_taxon_id.is_unique
            and a.open_flower_visible.eq("yes").all(),"38-photo visual source selection changed")
    return a

def prepare(input_dir,triage_path,output):
    output=Path(output); require(not output.exists(),"refuse overwrite")
    atlas,manifest,photos=read_inputs(input_dir)
    chosen=original_triage(triage_path)
    d=chosen.merge(manifest,on=["photo_id","inat_taxon_id","review_sequence"],validate="one_to_one")
    require(len(d)==38 and d.canonical_observation_uuid.is_unique,"photo provenance not independent")
    require(d.retrieval_status.eq("downloaded").all()
            and d.photo_license_code.isin(["cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa"]).all(),
            "unsupported source licence/photo")
    require(d.photo_id.is_unique and d.inat_taxon_id.is_unique,"repeated known photo or taxon")
    by_taxon=atlas.set_index("inat_taxon_id")
    require(d.inat_taxon_id.isin(by_taxon.index).all()
            and by_taxon.loc[d.inat_taxon_id,"after_strict9"].eq(0).all(),
            "target was newly classified or lost from most recent signed atlas")
    reference=photos.set_index("photo_id").loc[d.photo_id]
    require((reference.inat_taxon_id.to_numpy()==d.inat_taxon_id.to_numpy()).all()
            and (~reference.strict9).all(),"historical photo was not a failed ROI9")
    d=d.rename(columns={"image_sha256":"qa_medium_image_sha256",
                        "byte_count":"qa_medium_byte_count"})
    d["observation_id"]=d.canonical_observation_uuid
    d["photo_url_large"]=[source_large_url(u,pid) for u,pid in zip(d.review_url,d.photo_id)]
    d["photo_license"]=d.photo_license_code
    d["measurement_id"]="MOIS-REMEASURE38-"+d.photo_id
    d["image_filename"]=d.measurement_id+".jpg"
    d["shard"]=0
    d["measurement_kind"]="repeat_existing_photo"
    d["source_resolution"]="provider_large_derivative"
    output.mkdir(parents=True); (output/"shards").mkdir()
    shard=output/"shards/shard-000.csv"
    d.to_csv(shard,index=False)
    require(sha(shard)==SHA["source_shard"],
            "frozen original 38-photo SHA differs; refuse changed cohort")
    receipt={"schema":"mois_full_source_observation_sampling_cohort_v1",
             "colour_measured":False,"work_units":1,"worker_ids":[0],
             "photo_shard_sha256":{"0":sha(shard)},
             "measurement_kind":"repeat_existing_photo_not_new_acquisition",
             "repeat_photo_attempts_planned":38,
             "new_photo_attempts_planned":0,"new_unique_photo_ids":0,
             "original_resolution_downloads":0,
             "input_zip_sha256":{"atlas":SHA["latest_zip"],
                                 "photos":SHA["identities_zip"],"qa":SHA["qa_zip"]},
             "selection":"AI visual petaloid candidate, not random; no altered thresholds",
             "numerical_results_admitted_to_main_catalogue":False}
    (output/"manifest.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_VISUAL38_FROZEN_SOURCE",json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

def numeric_vector(row):
    supported=row["measurement_status"]=="classified_four_state_morph"
    require(supported==(row["morph"] in HUES),"morph measurement status changed")
    if not supported:return None
    v=[float(row["flower_fraction_"+c]) for c in PALETTE]
    require(all(math.isfinite(x) and 0<=x<=1.000001 for x in v)
            and abs(sum(v)-1)<=0.001,"invalid real numerical 9-component ROI")
    return v

def audit(input_dir,source_dir,worker_dir,output):
    out=Path(output);require(not out.exists(),"refuse overwrite")
    atlas,_,_=read_inputs(input_dir)
    s=Path(source_dir); w=Path(worker_dir)/"shard_000"
    m=json.loads((s/"manifest.json").read_text()); f=s/"shards/shard-000.csv"
    require(sha(f)==SHA["source_shard"]==m["photo_shard_sha256"]["0"],"source changed")
    chosen=pd.read_csv(f,dtype=str).fillna("")
    receipt=json.loads((w/"receipt.json").read_text())
    require(receipt["attempted"]==38 and receipt["source_sha256"]==sha(f)
            and receipt["photos_sha256"]==sha(w/"photos.csv"),
            "worker numerical photo source signature changed")
    terminal=w/"sealed/terminal_p000_b0000000.csv"
    require(receipt["terminal_sha256"]==sha(terminal),"sealed frozen worker changed")
    result=pd.read_csv(w/"photos.csv",dtype=str).fillna("")
    require(len(result)==38 and result.photo_id.is_unique
            and set(result.photo_id)==set(chosen.photo_id),"worker returned wrong photo set")
    result=result.set_index("photo_id").loc[chosen.photo_id].reset_index()
    for field in ["measurement_id","inat_taxon_id","observation_id","photo_url_large","photo_license"]:
        require((result[field].to_numpy()==chosen[field].to_numpy()).all(),
                "worker source identity mismatch: "+field)
    rows=[]
    for r in result.to_dict("records"):
        v=numeric_vector(r)
        label=""
        if v is not None:
            order=sorted(range(9),key=lambda i:v[i],reverse=True)
            if v[order[0]]>=0.60 and v[order[0]]-v[order[1]]>=0.20:
                label=PALETTE[order[0]]
        if r["measurement_status"]!="image_acquisition_failed":
            require(re.fullmatch(r"[0-9a-f]{64}",r["image_sha256"]) is not None,
                    "photo pixel SHA missing")
        row={k:r[k] for k in ["photo_id","inat_taxon_id","species","observation_id",
                              "photo_url_large","photo_license","attribution",
                              "image_sha256","morph","measurement_status"]}
        row.update(old_strict9=False,new_strict9=v is not None,
                   dominant_hue9=label,measurement_kind="existing_photo_recheck",
                   new_unique_photo=False,botanical_expert_validated=False)
        row.update({"strict9_"+p:("" if v is None else v[i]) for i,p in enumerate(PALETTE)})
        rows.append(row)
    positives=[r for r in rows if r["new_strict9"]]
    clear=[r for r in positives if r["dominant_hue9"]]
    current=atlas.set_index("inat_taxon_id")
    require(current.loc[[r["inat_taxon_id"] for r in rows],"after_strict9"].eq(0).all(),
            "one photo's taxon already has numeric ROI")
    patch=[]
    for row in positives:
        patch.append({"inat_taxon_id":row["inat_taxon_id"],"species":row["species"],
                      "prior_numeric_photo_count":0,
                      "candidate_numeric_photo_count":1,
                      "prior_photo_attempts_unchanged":int(current.at[row["inat_taxon_id"],"after_photo_attempts"]),
                      "candidate_hue_status":"provisional_single" if row["dominant_hue9"] else "mixed_image_unresolved",
                      "candidate_dominant_hue9":row["dominant_hue9"],
                      **{"candidate_fraction_"+p:row["strict9_"+p] for p in PALETTE}})
    out.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(out/"38_existing_photo_remeasurement_results.csv",index=False)
    pd.DataFrame(patch,columns=["inat_taxon_id","species","prior_numeric_photo_count",
        "candidate_numeric_photo_count","prior_photo_attempts_unchanged","candidate_hue_status",
        "candidate_dominant_hue9",*["candidate_fraction_"+p for p in PALETTE]]).to_csv(
        out/"candidate_taxon_patch_NOT_PROMOTED.csv",index=False)
    status_counts=pd.Series([r["measurement_status"] for r in rows]).value_counts().to_dict()
    summary={"schema":"mois_visual38_existing_v4_20261115_audit_v1",
             "existing_photo_ids_rechecked":38,
             "source_resolution":"provider_large_derivative_not_original_resolution",
             "new_unique_photo_ids":0,"new_original_image_ids":0,
             "source_latest_accepted_photo_attempts":1231341,
             "source_latest_strict9_numerical":467565,
             "source_latest_numerical_taxa":39911,
             "genuine_new_measurements_of_previously_failed_photo_ids":len(positives),
             "new_clear_sampled_hues_candidates":len(clear),
             "technical_unknown_or_failed":38-len(positives),
             "measurement_status_counts":{k:int(v) for k,v in status_counts.items()},
             "main_atlas_promoted":False,"botanical_flower_ROI_independently_verified":False,
             "frozen_source_sha256":sha(f),
             "worker_receipt_sha256":sha(w/"receipt.json"),
             "original_worker_terminal_sha256":sha(terminal),
             "method":receipt["method"],
             "output_sha256":{p.name:sha(p) for p in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print("MOIS_VISUAL38_REAL_PHOTO_REPEAT_AUDIT",json.dumps(summary,sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("stage",choices=("prepare","audit"))
    a.add_argument("--inputs",required=True)
    a.add_argument("--triage")
    a.add_argument("--source")
    a.add_argument("--worker")
    a.add_argument("--output",required=True)
    x=a.parse_args()
    if x.stage=="prepare":prepare(x.inputs,x.triage,x.output)
    else:audit(x.inputs,x.source,x.worker,x.output)
