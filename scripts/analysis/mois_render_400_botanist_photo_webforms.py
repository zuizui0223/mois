#!/usr/bin/env python3
"""Build standalone HTML botanical review pages from SHA-frozen blind CSVs.

Each page shows ONLY a licensed original photo, species, attribution, blind
token and reviewer response controls. No original dates, source regions,
season, automatic flower-colour outcomes or climate fields are embedded.
Reviewer-specific original input order is preserved; responses can be exported
and re-imported as the exact frozen CSV schema, without uploading to a server.

This creates no human review; both pages initially contain ZERO answers.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from scripts.analysis.mois_prepare_400_botanical_review_handoff import (
    FIELDS, STAGE, TISSUE,
)

REVIEW_STATUS="400_ORIGINAL_PHOTO_REVIEW_FORMS_BOTH_UNFILLED_NOT_COLOUR"
URL=re.compile(r"https://inaturalist-open-data[.]s3[.]amazonaws[.]com/photos/([0-9]+)/large[.](?:jpg|jpeg|png)",re.I)


def sha(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trusted(path:Path,expected:str,*,_synthetic_test:bool=False)->list[dict]:
    if not _synthetic_test and sha(path)!=expected:
        raise ValueError("The review page must use exactly the original frozen blank reviewer sheet")
    with path.open(encoding="utf-8",newline="") as f:
        r=csv.DictReader(f)
        if tuple(r.fieldnames or ())!=FIELDS:
            raise ValueError("Review sheet included forbidden season/climate/automatic hue columns")
        rows=list(r)
    if len(rows)!=(2 if _synthetic_test else 400):
        raise ValueError("Original review photo count drifted")
    if len({x["review_token"] for x in rows})!=len(rows):
        raise ValueError("Original photo reviewer tokens are not unique")
    for r in rows:
        if (not URL.fullmatch(r["original_photo_url"])
            or len(r["review_token"])!=20
            or not re.fullmatch(r"[a-f0-9]{20}",r["review_token"])
            or not r["species"] or not r["photo_license"]):
            raise ValueError("Unverified original image URL, licence or blind token")
        if any((r["reviewer_id"],r["open_flower_stage"],r["roi_organ"],r["notes"])):
            raise ValueError("Original reviewer input must be BLANK before authoring the photo review page")
    return rows


def page(role:str,rows:list[dict])->str:
    # Prevent user-provided attribution from escaping a script element.
    payload=json.dumps({"role":role,"fields":list(FIELDS),"records":rows,
         "stages":list(STAGE),"organs":list(TISSUE)},ensure_ascii=False,separators=(",",":"))
    payload=payload.replace("<","\\u003c").replace(">","\\u003e").replace("&","\\u0026")
    template=r"""<!doctype html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="referrer" content="no-referrer">
<title>MOIS 独立植物学写真レビュー __ROLE__</title>
<style>
:root { color-scheme:light; font-family:system-ui,-apple-system,"Helvetica Neue",sans-serif; }
body { max-width:960px; margin:0 auto; padding:20px; background:#f7f8f7; color:#18221d; }
h1 {font-size:1.5rem;margin:0 0 .4rem}
p { line-height:1.55;color:#405149 }
.panel { background:#fff; border:1px solid #ccd4cf; border-radius:12px; padding:16px; margin:14px 0; }
.photo {display:block;max-width:100%;width:auto;height:auto;max-height:62vh;object-fit:contain;margin:12px auto;}
nav {display:flex;gap:12px;flex-wrap:wrap;align-items:center;}
button {min-height:43px;border:1px solid #708479;border-radius:8px;background:#fff;color:#172d20;padding:9px 16px;cursor:pointer;}
button:hover {background:#e5ede6}
select,input[type=text],textarea {width:100%;font:inherit;border:1px solid #83958b;border-radius:7px;padding:10px;box-sizing:border-box;margin:6px 0 13px;}
label {display:block;font-weight:650}
.muted {font-size:.92rem;color:#4d6256}
.photo-credit {font-size:.85rem;overflow-wrap:anywhere}
.row {display:grid;grid-template-columns:1fr 1fr;gap:16px}
@media(max-width:600px){.row{grid-template-columns:1fr}}
.counter{font-variant-numeric:tabular-nums;font-weight:700}
</style>
</head>
<body>
<h1>独立植物学レビュー — 判定者 __ROLE__</h1>
<p>元写真を直接確認し、<strong>開いている花が写っているか</strong>と<strong>判定対象の組織</strong>を記録してください。季節・撮影年月・位置・自動花色・気候条件は表示されません。別の判定者と回答を相談せずに作業します。</p>
<div class="panel">
<label for="reviewer">判定者ID（A/B間で異なるID）</label>
<input id="reviewer" type="text" placeholder="自分の判定者IDを入力">
<div class="muted">このページは回答を外部サーバーへ送信しません。作業途中でも「CSVを書き出す」で保存し、次回「CSVを読み込む」で再開できます。写真はiNaturalist公開写真サーバーから表示されます。</div>
</div>
<div class="panel">
<nav><button id="prev" type="button">← 前の写真</button><span id="progress" class="counter"></span><button id="next" type="button">次の写真 →</button></nav>
<div id="species" style="font-weight:700;margin-top:18px"></div>
<img id="photo" class="photo" loading="lazy" referrerpolicy="no-referrer" alt="観察元写真">
<div id="photo-info" class="photo-credit"></div>
<div class="row">
<div><label for="stage">開花状態</label><select id="stage"></select></div>
<div><label for="organ">写真内の確認対象組織</label><select id="organ"></select></div>
</div>
<label for="notes">備考（任意）</label>
<textarea rows="2" id="notes" placeholder="判断に迷った点、花弁の場所など"></textarea>
</div>
<div class="panel">
<nav>
<button id="export" type="button">回答CSVを書き出す</button>
<label for="import" style="font-weight:400;cursor:pointer">保存済みCSVを読み込む <input id="import" type="file" accept=".csv,text/csv" style="display:block"></label>
</nav>
<p class="muted">未回答の行は必ず空欄のまま残ります。果実・つぼみ・苞葉などを白花／有色花とは扱いません。花色測定は別工程です。</p>
</div>
<script type="application/json" id="mois-original-blind-source">__PAYLOAD__</script>
<script>
"use strict";
const source=JSON.parse(document.getElementById('mois-original-blind-source').textContent);
const fields=source.fields,records=source.records,role=source.role;
const idxByToken=new Map(records.map((r,i)=>[r.review_token,i]));
const reply=new Map();
let index=0;
const el=id=>document.getElementById(id);
const stage=el('stage'),organ=el('organ'),reviewer=el('reviewer'),notes=el('notes');
function options(node,arr){
  node.replaceChildren();
  let blank=document.createElement('option');
  blank.value="";blank.textContent="未入力";node.appendChild(blank);
  for(const v of arr){
     const o=document.createElement('option');o.value=v;o.textContent=v;node.appendChild(o);
  }
}
options(stage,source.stages);options(organ,source.organs);
function completed(t){
  const r=reply.get(t);
  return Boolean(r && r.open_flower_stage && r.roi_organ && reviewer.value.trim());
}
function store(){
  const t=records[index].review_token;
  const s=stage.value,o=organ.value,n=notes.value;
  if(s||o||n)reply.set(t,{open_flower_stage:s,roi_organ:o,notes:n});
  else reply.delete(t);
  updateProgress();
}
function updateProgress(){
  const n=records.filter(r=>completed(r.review_token)).length;
  el('progress').textContent=(index+1)+" / "+records.length+" • 完了 "+n+"件";
  el('prev').disabled=index===0;el('next').disabled=index===records.length-1;
}
function render(){
  const item=records[index],saved=reply.get(item.review_token)||{};
  el('species').textContent=item.species;
  const im=el('photo');im.src=item.original_photo_url;im.alt=item.species+" の原本写真";
  const target=el('photo-info');target.replaceChildren();
  const link=document.createElement('a');link.href=item.original_photo_url;link.target='_blank';
  link.rel='noopener noreferrer';link.textContent='元写真を開く';
  target.append(link,document.createTextNode(" · "+item.photo_license+" · "+item.attribution+
          " · review token "+item.review_token));
  stage.value=saved.open_flower_stage||"";
  organ.value=saved.roi_organ||"";
  notes.value=saved.notes||"";
  updateProgress();
}
stage.addEventListener('change',store);
organ.addEventListener('change',store);
notes.addEventListener('change',store);
reviewer.addEventListener('input',updateProgress);
el('prev').addEventListener('click',()=>{store();index=Math.max(0,index-1);render();});
el('next').addEventListener('click',()=>{store();index=Math.min(records.length-1,index+1);render();});
function csvValue(v){
 const text=String(v??'');
 return '"'+text.replaceAll('"','""')+'"';
}
function exportCSV(){
 store();
 if(!reviewer.value.trim()){alert("判定者IDを入力してください");return;}
 let count=0;
 const lines=[fields.map(csvValue).join(',')];
 for(const row of records){
   const saved=reply.get(row.review_token)||{};
   const done=Boolean(saved.open_flower_stage && saved.roi_organ);
   count+=Number(done);
   const x={...row,
     reviewer_id:done?reviewer.value.trim():"",
     open_flower_stage:done?saved.open_flower_stage:"",
     roi_organ:done?saved.roi_organ:"",
     notes:done?(saved.notes||""):""};
   lines.push(fields.map(key=>csvValue(x[key])).join(','));
 }
 const blob=new Blob(["\uFEFF"+lines.join('\r\n')+'\r\n'],{type:"text/csv;charset=utf-8"});
 const link=document.createElement('a');link.download='mois_reviewer_'+role+'_answered_'+count+'_of_'+records.length+'.csv';
 const url=URL.createObjectURL(blob);link.href=url;document.body.appendChild(link);
 link.click();link.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
el('export').addEventListener('click',exportCSV);
function parseCSV(text){
 const result=[],row=[];let value="",quoted=false;
 text=text.replace(/^\uFEFF/,'');
 for(let i=0;i<text.length;i++){
    const c=text[i];
    if(quoted){
      if(c==='"' && text[i+1]==='"'){value+='"';i++;}
      else if(c==='"'){quoted=false;}
      else value+=c;
    }else if(c==='"'){quoted=true;}
    else if(c===','){row.push(value);value="";}
    else if(c==='\n' || c==='\r'){
      if(c==='\r' && text[i+1]==='\n')i++;
      row.push(value);result.push(row.splice(0));value="";
    }else value+=c;
 }
 if(quoted)throw Error("未終了のCSV引用符");
 if(value || row.length){row.push(value);result.push(row);}
 return result;
}
el('import').addEventListener('change',async e=>{
 try{
   const file=e.target.files&&e.target.files[0];if(!file)return;
   const all=parseCSV(await file.text());
   if(JSON.stringify(all.shift())!==JSON.stringify(fields))throw Error("別形式のCSVです");
   if(all.length!==records.length)throw Error("元の400枚全件を残したCSVが必要です");
   const seen=new Set(),received=new Map(),people=new Set();
   for(const cells of all){
      if(cells.length!==fields.length)throw Error("CSV列の数が違います");
      const x=Object.fromEntries(fields.map((k,i)=>[k,cells[i]]));
      const at=idxByToken.get(x.review_token);
      if(at===undefined || seen.has(x.review_token))throw Error("未知・重複トークン");
      seen.add(x.review_token);
      const orig=records[at];
      if(["species","original_photo_url","photo_license","attribution"].some(k=>orig[k]!==x[k]))
        throw Error("原写真ID・出典が変わっています");
      const filled=Boolean(x.open_flower_stage||x.roi_organ||x.reviewer_id);
      if(filled && (!x.open_flower_stage||!x.roi_organ||!x.reviewer_id))
        throw Error("部分入力の回答が含まれています");
      if(x.open_flower_stage && !source.stages.includes(x.open_flower_stage))throw Error("不明な開花状態");
      if(x.roi_organ && !source.organs.includes(x.roi_organ))throw Error("不明な組織");
      if(x.reviewer_id)people.add(x.reviewer_id);
      if(filled)received.set(x.review_token,{open_flower_stage:x.open_flower_stage,
        roi_organ:x.roi_organ,notes:x.notes});
   }
   if(people.size>1)throw Error("複数の判定者IDが混在しています");
   if(people.size===1 && reviewer.value && reviewer.value.trim()!==[...people][0])
      throw Error("現在の判定者IDとCSVの判定者が異なります");
   reply.clear();for(const [k,v] of received)reply.set(k,v);
   if(people.size===1)reviewer.value=[...people][0];
   index=0;render();alert("回答"+received.size+"件を復元しました");
 }catch(err){alert("CSVの読み込みを拒否: "+err.message);}
});
render();
</script>
</body></html>"""
    return template.replace("__ROLE__",role).replace("__PAYLOAD__",payload)


def build(form_a:Path,form_b:Path,receipt:Path,out:Path,
          *,_synthetic_test:bool=False)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Do not overwrite sealed botanist review portal handoffs")
    d=json.loads(receipt.read_text(encoding="utf-8"))
    if (d.get("status")!=REVIEW_STATUS or
        d.get("reviewer_forms_season_date_climate_colour_blind") is not True or
        d.get("all_review_response_fields_blank") is not True or
        d.get("n_actual_independent_botanical_reviews_completed")!=0 or
        d.get("new_photo_pixels_measured")!=0):
        raise ValueError("Not the original, zero-response observer-blind source handoff")
    if not _synthetic_test and (
            d.get("n_original_photos")!=400 or d.get("n_original_named_taxa")!=5 or
            sha(form_a)!=d.get("reviewer_form_A_SHA256") or
            sha(form_b)!=d.get("reviewer_form_B_SHA256")):
        raise ValueError("Original 400 review form bytes or SHA changed")
    inputs={}
    for who,path,key in (("A",form_a,"reviewer_form_A_SHA256"),
                         ("B",form_b,"reviewer_form_B_SHA256")):
        inp=trusted(path,d[key],_synthetic_test=_synthetic_test)
        inputs[who]=inp
    if {r["review_token"] for r in inputs["A"]}!={r["review_token"] for r in inputs["B"]}:
        raise ValueError("Reviewer A and B do not share exactly the same original photo IDs")
    out.mkdir(parents=True,exist_ok=True)
    results={}
    for who,entries in inputs.items():
        dest=out/f"mois_400_review_{who}_original_photo_browser_UNREVIEWED.html"
        dest.write_text(page(who,entries),encoding="utf-8")
        results[who]={"filename":dest.name,"html_sha256":sha(dest),
                       "original_source_photos":len(entries),
                       "filled_review_responses":0}
    report={
        "schema":"mois_400_original_licensed_photo_independent_plant_phenology_browser_v1",
        "status":"TWO_STANDALONE_BROWSER_REVIEW_PORTALS_CREATED_WITH_ZERO_BOTANICAL_DECISIONS",
        "original_400_reviewer_handoff_SHA256":sha(receipt),
        "reviewer_original_csv_SHA256":{"A":sha(form_a),"B":sha(form_b)},
        "portals":results,
        "live_images_loaded_from_public_licensed_original_S3_on_view":True,
        "no_original_image_bytes_downloaded_into_artifact":True,
        "no_photo_observation_date_geography_climate_or_automated_colour_embedded":True,
        "both_reviewers_have_independently_shuffled_original_photo_order":True,
        "original_botanical_review_response_count":0,
        "no_automatically_confirmed_flowering_photos":True,
        "no_seasonal_white_fraction_or_climate_effect_inferred":True,
        "html_contains_client_only_csv_export_import_not_server_upload":True,
    }
    (out/"MOIS_400_BOTANIST_IMAGE_REVIEW_PORTAL_RECEIPT.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reviewer-A-unfilled",type=Path,required=True)
    p.add_argument("--reviewer-B-unfilled",type=Path,required=True)
    p.add_argument("--original-handoff-receipt",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(build(a.reviewer_A_unfilled,a.reviewer_B_unfilled,
                           a.original_handoff_receipt,a.outdir),sort_keys=True))


if __name__=="__main__":
    main()
