#!/usr/bin/env python3
"""Prepare independent, outcome-blind 1200/240 existing-photo corolla review.

Human labels are intentionally BLANK. HTML loads one public iNaturalist image
at a time; no source pixels are downloaded by this generator or redistributed.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
from collections import Counter

from scripts.analysis.mois_observation_probability_audit_20261029 import (
    SCHEMA as SOURCE_SCHEMA, allocate, hash_rank, require, sha_file
)

REVIEW_SCHEMA = "mois_existing_observation1200_blind_human_corolla_review_v1"
DOUBLE_SUCCESS = 80
DOUBLE_UNKNOWN = 160

HTML = r'''<!doctype html>
<html lang="ja"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>MOIS 既存写真・花冠の独立目視監査</title>
<style>
:root{font-family:system-ui,-apple-system,"Noto Sans JP",sans-serif;color:#182938;background:#f6f8fa}
*{box-sizing:border-box}main{max-width:1000px;margin:auto;padding:18px}
h1{font-size:1.35rem}p{line-height:1.6;color:#485d6e}article{background:white;border:1px solid #d5dce4;border-radius:12px;padding:16px}
.grid{display:grid;grid-template-columns:1.25fr 1fr;gap:20px}@media(max-width:750px){.grid{grid-template-columns:1fr}}
#photo{width:100%;max-height:66vh;object-fit:contain;background:#f0f2f5}
label{font-weight:600;display:block;margin:10px 0 4px}select,input,textarea,button{font:inherit;max-width:100%;border:1px solid #bdc8d2;border-radius:7px;padding:10px;width:100%;background:white}
textarea{height:80px}.buttons{display:grid;grid-template-columns:1fr 1fr;gap:9px;margin-top:14px}
button{cursor:pointer}button.primary{background:#153e5b;color:white}.meta{font-size:.84rem;color:#4b647a}
.progress{height:9px;background:#ddd;border-radius:8px;overflow:hidden}.progress div{height:100%;background:#296c96;width:0}
</style></head><body><main>
<h1>MOIS｜花冠色の独立監査 <span id="role"></span></h1>
<p>撮影時の環境・技術的ROI判定結果・機械の色分類は伏せています。花が実際に開いているか、花冠の可視色を目で確認してください。自生・栽培はこの画像だけで推定しないでください。</p>
<div class="progress"><div id="bar"></div></div><p id="progress" class="meta"></p>
<article class="grid"><section>
<p id="photoid" class="meta"></p>
<img id="photo" alt="iNaturalistの元写真">
<p id="imagestatus" class="meta">写真を閲覧するときのみ公開画像を取得します。</p>
</section><section>
<label>独立レビュー担当コード</label><input id="reviewer" placeholder="例: reviewer-A">
<label>開花した花冠が見えるか（必須）</label>
<select id="flower"><option value="">未確認</option><option value="yes">はい</option><option value="no">いいえ</option><option value="uncertain">不確実</option><option value="inaccessible">画像取得不可</option></select>
<label>花冠の可視色（開花を確認できた場合）</label>
<select id="colour"><option value="">未確認</option><option value="white">白優占</option><option value="yellow_orange">黄・橙優占</option><option value="red_pink">赤・桃優占</option><option value="blue_purple">青・紫優占</option><option value="mixed">多色（優占不明）</option><option value="other">その他</option><option value="uncertain">色が不確実</option><option value="not_applicable">対象花冠なし</option></select>
<label>画像の主な問題（任意）</label>
<select id="issue"><option value="">未選択</option><option value="no_open_flower">開花なし</option><option value="small_or_occluded">小型・遮蔽</option><option value="wrong_organ">花冠以外が目立つ</option><option value="illumination">露出・照明</option><option value="blur">ピント・ブレ</option><option value="none">特になし</option><option value="other">その他・不明</option></select>
<label>メモ（任意）</label><textarea id="notes"></textarea>
<div class="buttons"><button id="prev">前の写真</button><button id="next" class="primary">保存して次へ</button>
<button id="export">CSVを書き出し</button><button id="import">以前のCSVを読み込む</button></div>
<input id="importfile" type="file" accept=".csv,text/csv" hidden>
<p class="meta">回答はこのブラウザ内だけに保存されます。CSVを書き出してバックアップしてください。自動送信なし。RGB写真はHTML内に同梱されません。</p>
</section></article></main>
<script>
"use strict";
const ROWS=__ROWS__;
const ROLE=__ROLE__;
const VERSION=__VERSION__;
document.getElementById("role").textContent=ROLE;
const IDS=new Set(ROWS.map(x=>x.photo_id));
let idx=0, state={};
try{state=JSON.parse(localStorage.getItem(VERSION)||"{}")||{};}catch(_){state={};}
const fields=["flower","colour","issue","notes"];
const el=(x)=>document.getElementById(x);
function store(){try{localStorage.setItem(VERSION,JSON.stringify(state));}catch(_){}}
function complete(v){return !!v && !!v.flower && (v.flower!=="yes" || !!v.colour);}
function progress(){
 const done=ROWS.filter(x=>complete(state[x.photo_id])).length;
 el("bar").style.width=(100*done/ROWS.length)+"%";
 el("progress").textContent=ROLE+" | 判定完了 "+done+"/"+ROWS.length+" | 表示中 "+(idx+1)+"/"+ROWS.length;
}
function save(){
 const r=ROWS[idx];if(!r)return;
 let v={reviewer:el("reviewer").value};
 for(const f of fields)v[f]=el(f).value;
 state[r.photo_id]=v;store();progress();
}
function show(i){
 idx=Math.min(ROWS.length-1,Math.max(0,i));
 const r=ROWS[idx],v=state[r.photo_id]||{};
 el("reviewer").value=v.reviewer||el("reviewer").value||"";
 for(const f of fields)el(f).value=v[f]||"";
 el("photoid").textContent="元写真 ID "+r.photo_id+"（地理・ROI成否は非表示）";
 const pic=el("photo"),msg=el("imagestatus");
 pic.removeAttribute("src");
 let n=0,ext=["jpeg","jpg","png"];
 function url(){return "https://inaturalist-open-data.s3.amazonaws.com/photos/"+r.photo_id+"/medium."+ext[n];}
 pic.onerror=function(){n++;if(n<ext.length)pic.src=url();else msg.textContent="この写真URLを取得できません。inaccessibleを選択してください。";};
 pic.onload=function(){msg.textContent="公開写真表示中。花冠の確認のみ行ってください。";};
 msg.textContent="既存公開写真を読み込み中…";pic.src=url();progress();
}
el("prev").onclick=function(){save();show(idx-1);};
el("next").onclick=function(){save();show(idx+1);};
for(const f of fields)el(f).addEventListener("change",save);
el("reviewer").addEventListener("change",save);
function csvEscape(s){s=String(s==null?"":s).replace(/[\r\n]/g," ");return '"'+s.replace(/"/g,'""')+'"';}
const HEADER=["photo_id","observation_id","review_role","reviewer_id","human_review_completed","human_open_corolla_visible","human_corolla_colour","failure_cause","notes"];
el("export").onclick=function(){
 save();let lines=[HEADER.join(",")];
 for(const r of ROWS){
  const v=state[r.photo_id]||{};
  lines.push([r.photo_id,r.observation_id,ROLE,v.reviewer||"",complete(v)?"True":"False",
     v.flower||"",v.colour||"",v.issue||"",v.notes||""].map(csvEscape).join(","));
 }
 const b=new Blob(["\ufeff"+lines.join("\r\n")+"\r\n"],{type:"text/csv;charset=utf-8"});
 const u=URL.createObjectURL(b),a=document.createElement("a");
 a.href=u;a.download="mois_audit1200_"+ROLE+"_human_review.csv";a.click();URL.revokeObjectURL(u);
};
function parseCSV(txt){
 txt=txt.replace(/^\uFEFF/,"");let rows=[],row=[],v="",q=false;
 for(let i=0;i<txt.length;i++){
  const c=txt[i];
  if(q){if(c==='"'){if(txt[i+1]==='"'){v+='"';i++;}else q=false;}else v+=c;}
  else if(c==='"')q=true;else if(c===','){row.push(v);v="";}
  else if(c==='\n'){row.push(v.replace(/\r$/,""));rows.push(row);row=[];v="";}
  else v+=c;
 }
 if(v.length||row.length){row.push(v);rows.push(row);}
 return rows;
}
el("import").onclick=function(){el("importfile").click();};
el("importfile").onchange=async function(){
 const file=this.files&&this.files[0];if(!file)return;
 const rows=parseCSV(await file.text());if(!rows.length)return;
 const h=rows.shift(),ix=Object.fromEntries(h.map((v,i)=>[v,i]));
 if(!("photo_id" in ix)||!("review_role" in ix)){alert("対応するCSVではありません");return;}
 let n=0;
 for(const r of rows){
  const p=r[ix.photo_id];if(!IDS.has(p))continue;
  if(r[ix.review_role]!==ROLE){alert("別レビュー担当のCSVは読み込めません");return;}
  state[p]={reviewer:r[ix.reviewer_id]||"",flower:r[ix.human_open_corolla_visible]||"",
   colour:r[ix.human_corolla_colour]||"",issue:r[ix.failure_cause]||"",notes:r[ix.notes]||""};n++;
 }
 store();show(idx);alert(n+"件を復元しました");
};
show(0);
</script></body></html>'''

def write_review_html(path, roster, role, roster_sha):
    require(all(r["photo_id"].isdecimal() and r["observation_id"] for r in roster),
            "Bad original photo/observation identity for HTML")
    # No machine hue, source arm, latitude, season, flower truth or photo bytes.
    data = [{"photo_id": r["photo_id"], "observation_id": r["observation_id"]}
            for r in roster]
    html = (HTML.replace("__ROWS__", json.dumps(data, separators=(",", ":")))
            .replace("__ROLE__", json.dumps(role))
            .replace("__VERSION__", json.dumps("mois1200|" + role + "|" + roster_sha)))
    path.write_text(html, encoding="utf-8")

def create_review_bundle(source_dir):
    root=Path(source_dir)
    roster_path=root/"audit1200_source_observation_roster.csv"
    receipt_path=root/"audit1200_frozen_design_receipt.json"
    original=json.loads(receipt_path.read_text())
    require(original.get("schema")==SOURCE_SCHEMA and original.get("sample_total")==1200 and
            original.get("sample_success_count")==400 and
            original.get("sample_unknown_count")==800,
            "Expected authentic 1200 photo source receipt")
    require(original["output_sha256"]["audit1200_source_observation_roster.csv"]==
            sha_file(roster_path), "Primary 1200 audit roster source SHA drifted")
    with roster_path.open(newline="") as f:
        roster=list(csv.DictReader(f))
    require(len(roster)==1200 and len({r["photo_id"] for r in roster})==1200 and
            len({r["observation_id"] for r in roster})==1200,
            "Primary audit originals duplicated")
    quota={}
    by_status=Counter(r["roi_status"] for r in roster)
    require(by_status=={"success":400,"unknown":800},"Audit outcomes changed")
    for status,n in (("success",DOUBLE_SUCCESS),("unknown",DOUBLE_UNKNOWN)):
        counts=Counter(r["stratum"] for r in roster if r["roi_status"]==status)
        quota.update(allocate(counts,n))
    chosen=[]
    for stratum in sorted(quota):
        rows=sorted((r for r in roster if r["stratum"]==stratum),
                    key=lambda r: hash_rank("second_reviewer",r["observation_id"]))
        chosen.extend(rows[:quota[stratum]])
    chosen=sorted(chosen,key=lambda r: hash_rank("second_review_order",r["observation_id"]))
    require(len(chosen)==240 and
            len({r["photo_id"] for r in chosen})==240 and
            Counter(r["roi_status"] for r in chosen)=={"success":80,"unknown":160},
            "Double-blind 20 percent overlap invalid")
    second=root/"audit1200_secondary240_review_BLANK.csv"
    with second.open("w",newline="") as f:
        fields=("photo_id","observation_id","reviewer_id","human_review_completed",
                "human_open_corolla_visible","human_corolla_colour","failure_cause","notes")
        wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader()
        for r in chosen:
            wr.writerow({"photo_id":r["photo_id"],"observation_id":r["observation_id"],
                         "reviewer_id":"","human_review_completed":"False",
                         "human_open_corolla_visible":"","human_corolla_colour":"",
                         "failure_cause":"","notes":""})
    primary_html=root/"audit1200_primary_blind_review.html"
    secondary_html=root/"audit1200_secondary240_blind_review.html"
    write_review_html(primary_html,roster,"R1",sha_file(roster_path))
    write_review_html(secondary_html,chosen,"R2",sha_file(roster_path))
    m={"schema":REVIEW_SCHEMA,"source_roster_sha256":sha_file(roster_path),
       "n_primary_photo_IDs":1200,"n_second_reviewer_photo_IDs":240,
       "double_review_success":80,"double_review_unknown":160,
       "same_real_photo_IDs_not_new_images":True,
       "human_labels_added":0,"machine_hue_hidden_in_UI":True,
       "roi_status_source_strata_geo_hidden_in_UI":True,
       "photo_binary_embedded_or_distributed":False,
       "viewer_reads_public_image_only_when_requested":True,
       "output_sha256":{p.name:sha_file(p) for p in
             (primary_html,secondary_html,second)}}
    (root/"audit1200_blind_review_UI_receipt.json").write_text(
        json.dumps(m,indent=2,sort_keys=True)+"\n")
    print("MOIS_1200_BLIND_BOTANICAL_REVIEW_UI_READY",
          json.dumps({k:m[k] for k in (
              "n_primary_photo_IDs","n_second_reviewer_photo_IDs",
              "human_labels_added")},sort_keys=True),flush=True)
    return m

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit-dir",type=Path,required=True)
    a=p.parse_args()
    create_review_bundle(a.audit_dir)

if __name__=="__main__":
    main()
