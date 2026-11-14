#!/usr/bin/env python3
"""Signed MOIS zero-ROI9 failure review queue, without inventing flower labels.

Uses SHA-pinned 42111-taxon repaired atlas, 1181146 original photo ID ledger,
and signed initial taxonomy. Creates an auditable manual botanical review queue.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

CORRECTED_ZIP_SHA = "ef6578a8934280de654d9ff3adf9c8458a041feca463dc52b3cd87f501ba527e"
CORRECTED_CSV_SHA = "11c0bcf0cb54e3a665f764ac74906b15a275d0f57f3b07a30d3073791b985a34"
BASELINE_ZIP_SHA = "6cc7689938dcf13baf6498b954286924aa97f2b63b3418c41181b18e346c0d2b"
BASELINE_CSV_SHA = "5f51b9115719ff3c94a46851f9a3972742f7c52192ceb31b3da110248028fad8"
ORIGINAL_PHOTOS_SHA = "3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835"
SCHEMA = "mois_zero9_96_source_photo_review_v1"


def require(value, reason):
    if not bool(value):
        raise ValueError("MOIS_ZERO9_REVIEW: " + reason)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_zip_atlas(path, expected_zip_sha, inside, expected_csv_sha):
    require(sha(path)==expected_zip_sha, "source artifact ZIP digest changed")
    with ZipFile(path) as z:
        csv=z.read(inside)
    require(hashlib.sha256(csv).hexdigest()==expected_csv_sha,
            "source atlas CSV digest changed")
    return pd.read_csv(io.BytesIO(csv),compression="gzip",
                       dtype={"inat_taxon_id":str},low_memory=False)


def select_zero_taxa(atlas):
    require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique,
            "signed 42,111 nominal source taxa missing")
    zero=atlas.loc[atlas.after_strict9.eq(0)].copy()
    require(len(zero)==2200 and zero.after_photo_attempts.ge(1).all(),
            "signed no-ROI9 count or photo attempts changed")
    zero["review_tier"]=zero.after_photo_attempts.apply(
        lambda n: "P0_10plus" if n>=10 else ("P1_3to9" if n>=3 else "P2_1to2"))
    high=zero.loc[zero.after_photo_attempts.ge(10)].copy()
    require(len(high)==96 and int(high.after_photo_attempts.sum())==1278,
            "96 repeated high-failure taxa count changed")
    return zero.sort_values(["after_photo_attempts","inat_taxon_id"],
                            ascending=[False,True]),high


def shortlist_original_images(photos, high, initial, per_taxon=3):
    require(per_taxon>=1,"invalid per-taxon review sample capacity")
    p=photos.loc[photos.inat_taxon_id.astype(str).isin(set(high.inat_taxon_id))].copy()
    for k in ("photo_id","inat_taxon_id","observation_id","source_subarm","source_design",
              "observed_on","original_roi_four_hue"):
        p[k]=p[k].fillna("").astype(str).str.strip()
    require(len(p)>0 and p.photo_id.is_unique and p.photo_id.str.fullmatch(r"\d+").all(),
            "original photo IDs malformed")
    require((~p.strict9).all() and set(p.inat_taxon_id)==set(high.inat_taxon_id),
            "unexpected original true ROI9 or missing source taxon")
    counts=p.groupby("inat_taxon_id").size()
    b=initial.set_index("inat_taxon_id")
    require((counts==b.loc[counts.index,"after_photo_attempts"]).all(),
            "source photo IDs do not reconstruct the historic atlas")
    p["uuid_known"]=p.observation_id.ne("")
    p["observation_key"]=p.observation_id.where(p.uuid_known,"photo:"+p.photo_id)
    p["hash_rank"]=[hashlib.sha256(
        ("MOIS_96_PHOTO_REVIEW_20261114|"+tax+"|"+pid).encode()).hexdigest()
                    for tax,pid in zip(p.inat_taxon_id,p.photo_id)]
    p=p.sort_values(["inat_taxon_id","hash_rank"],kind="stable").drop_duplicates(
        ["inat_taxon_id","observation_key"],keep="first")
    p["review_order"]=p.groupby("inat_taxon_id").cumcount()+1
    p=p.loc[p.review_order.le(per_taxon)].copy()
    require(p.inat_taxon_id.nunique()==96 and int(p.review_order.eq(1).sum())==96,
            "must have exactly one primary photo per high-failure taxon")
    p=p.merge(high[["inat_taxon_id","species","after_photo_attempts",
                    "full_93m_original_photo_count","full_93m_original_observation_count"]],
              on="inat_taxon_id",validate="many_to_one")
    p["photo_page_url"]="https://www.inaturalist.org/photos/"+p.photo_id
    # Intentionally unfilled: technical ROI failure does not determine if open flowers exist.
    p["manual_open_flower_present"]=""
    p["manual_stage"]=""
    p["manual_failure_cause"]=""
    p["reviewer"]=""
    p["notes"]=""
    p=p.drop(columns=["hash_rank"])
    return p.sort_values(["after_photo_attempts","inat_taxon_id","review_order"],
                         ascending=[False,True,True],kind="stable"),int(counts.sum())


def html_review(primary):
    # No borrowed photo pixels or unsupported invented extension/licence.
    from html import escape
    top="""<!DOCTYPE html><html lang="ja"><meta charset="utf-8"><title>MOIS 96種 開花器官レビュー</title>
<style>body{font:16px system-ui;max-width:880px;margin:1em auto;padding:1em;background:#f6f8f7}
section{background:white;margin:12px 0;padding:13px;border:1px solid #ddd;border-radius:8px}
a{color:#14657c}h2{font-size:1.07em}select,input,button{padding:7px;margin:3px}button{background:#215945;color:white;border:none;border-radius:5px}</style>
<h1>MOIS：9色ROI失敗96分類群・最初の写真確認</h1>
<p>原写真ページを開いて、開花の有無を独立に確認してください。このフォームの全分類は未入力です。
9色ROI失敗は花不在の証拠ではありません。外部写真の再配布はしていません。</p>
<button onclick="exportCsv()">CSVを書き出す</button><span id="n"></span><div id="items"></div>
<script>const entries=__DATA__;const result={};let root=document.getElementById("items");
function s(options,onChange){let x=document.createElement("select");options.forEach(t=>{let o=document.createElement("option");o.value=t;o.textContent=t||"未入力";x.append(o)});x.onchange=()=>onChange(x.value);return x}
function csv(x){x=String(x??"");return /[",\\n\\r]/.test(x)?'"'+x.replaceAll('"','""')+'"':x}
function refresh(){document.getElementById('n').textContent=" 入力 "+Object.values(result).filter(x=>x.stage).length+" / "+entries.length}
entries.forEach(e=>{let item=document.createElement("section");let h=document.createElement("h2");h.textContent=e.species+" ("+e.attempts+"回失敗)";item.append(h);
let link=document.createElement("a");link.textContent="原写真 "+e.photo_id+" をiNaturalistで開く";link.href=e.url;link.target="_blank";link.rel="noopener noreferrer";item.append(link);
let r={...e,stage:"",cause:"",notes:""};result[e.photo_id]=r;
item.append(document.createElement('p'),s(['','open_flower','bud_only','fruit_or_vegetative','nonplant','image_unavailable','uncertain'],v=>{r.stage=v;refresh()}),
s(['','open_flower_but_ROI9_failed','not_an_open_flower','unresolved'],v=>{r.cause=v;refresh()}));
let note=document.createElement("input");note.placeholder="備考";note.oninput=()=>r.notes=note.value;item.append(note);root.append(item)});refresh();
function exportCsv(){const cols=["inat_taxon_id","species","photo_id","observation_id","attempts","url","stage","cause","notes"];
const lines=[cols.join(","),...Object.values(result).map(row=>cols.map(k=>csv(row[k])).join(","))].join("\\r\\n");
let blob=new Blob(["\\ufeff",lines],{type:"text/csv;charset=utf-8"});let a=document.createElement("a");a.href=URL.createObjectURL(blob);a.download="mois_96_open_flower_review.csv";a.click();URL.revokeObjectURL(a.href)}
</script></html>"""
    entry=[{"inat_taxon_id":str(x.inat_taxon_id),"species":str(x.species),
            "photo_id":str(x.photo_id),"observation_id":str(x.observation_id),
            "attempts":int(x.after_photo_attempts),"url":str(x.photo_page_url)}
           for x in primary.itertuples()]
    payload=json.dumps(entry,ensure_ascii=False).replace("<",r"\u003c").replace(">",r"\u003e")
    return top.replace("__DATA__",payload)


def run(repaired_zip, initial_zip, original_parquet, output, per_taxon=3):
    current=read_zip_atlas(repaired_zip,CORRECTED_ZIP_SHA,
        "mois_42111_taxa_corrected_strict9_counts.csv.gz",CORRECTED_CSV_SHA)
    initial=read_zip_atlas(initial_zip,BASELINE_ZIP_SHA,
        "mois_species_first_42111_nine_palette_catalogue.csv.gz",BASELINE_CSV_SHA)
    zero, high=select_zero_taxa(current)
    require(sha(original_parquet)==ORIGINAL_PHOTOS_SHA,
            "frozen real-photo source parquet SHA changed")
    photos=pd.read_parquet(original_parquet,columns=[
        "photo_id","observation_id","inat_taxon_id","strict9","source_subarm",
        "source_design","observed_on","original_roi_four_hue"])
    require(len(photos)==1181146 and photos.photo_id.is_unique,
            "frozen 1,181,146 original photos missing or duplicated")
    review, initial_attempts=shortlist_original_images(photos,high,initial,per_taxon)
    output=Path(output)
    require(not output.exists(),"refuse to overwrite an existing audit")
    output.mkdir(parents=True)
    cols=["inat_taxon_id","species","review_tier","after_photo_attempts",
          "full_93m_original_photo_count","full_93m_original_observation_count"]
    zero[cols].to_csv(output/"mois_2200_zero9_qa_priority.csv",index=False)
    review.to_csv(output/"mois_96_review_original_photo_IDs.csv",index=False)
    clear4=review.original_roi_four_hue.isin(("white","yellow_orange","red_pink","blue_purple"))
    followup=review.loc[clear4].copy()
    # Technical 4-hue proxy only; not ROI9 rescue or botanical flower truth.
    require(len(followup)==31 and followup.inat_taxon_id.nunique()==29,
            "signed four-hue but strict9-unknown cohort changed")
    require(followup.source_subarm.eq("historical_main").all(),
            "old four-hue diagnostic source arm changed")
    followup.to_csv(output/"mois_29_four_hue_but_no_numeric9_followup.csv",index=False)
    (output/"mois_96_blind_bloom_review.html").write_text(
        html_review(review.loc[review.review_order.eq(1)]),encoding="utf-8")
    receipt={"schema":SCHEMA,"source_taxa":42111,"no_numerical_ROI9_taxa":2200,
        "ten_plus_failed_taxa":96,"ten_plus_total_original_attempts":1278,
        "ten_plus_frozen_initial_attempts":initial_attempts,
        "review_photo_rows":len(review),"review_primary_taxa":96,
        "technical_four_hue_positive_photo_rows":len(followup),
        "technical_four_hue_positive_taxa":int(followup.inat_taxon_id.nunique()),
        "technical_four_hue_positive_source_arm":"historical_main",
        "uuid_known_review_rows":int(review.uuid_known.sum()),
        "four_hue_diagnostic_nonempty":int(review.original_roi_four_hue.ne("").sum()),
        "new_original_photo_downloads":0,"new_roi9_successes":0,
        "open_flower_diagnosis_from_photos_completed":False,
        "original_photo_parquet_SHA256":sha(original_parquet),
        "repaired_atlas_artifact_SHA256":sha(repaired_zip),
        "input_baseline_atlas_artifact_SHA256":sha(initial_zip)}
    receipt["output_SHA256"]={p.name:sha(p) for p in output.iterdir()}
    (output/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_SIGNED_ZERO9_REVIEW_QUEUE_SUCCESS",json.dumps(receipt,sort_keys=True))
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("repaired-zip","initial-zip","original-parquet","output"):
        p.add_argument("--"+x,required=True)
    p.add_argument("--max-per-taxon",type=int,default=3)
    a=p.parse_args()
    run(a.repaired_zip,a.initial_zip,a.original_parquet,a.output,a.max_per_taxon)


if __name__=="__main__":
    main()
