#!/usr/bin/env python3
"""Freeze a 128-existing-photo botanist review UI from real 512 ROI9 failures.

The output is an annotation tool, NOT a classification result. Preselect
technical detector negatives/errors (uncertain flower state), large positive
boxes (possible segmentation/flip-gate failure), plus independent source
arms/latitudes. No new original photos are acquired and no botanical
labels or pigmented/white classifications are invented. HTML loads original
iNaturalist medium variants as the user browses, one photo at a time;
NO image bytes are embedded, packaged, stored by the generator or posted
back to a server. Human labels can be exported as separate CSV locally.
"""
from __future__ import annotations
from collections import Counter,defaultdict
import argparse,csv,hashlib,json
from pathlib import Path

SCHEMA="mois_prior9unknown128_independent_botanical_photo_roi_gate_review_v1"
SOURCE_SHA="4114dfc1642355dc7b1dc2e38527fe3a9ff6477b932fe84a6b69b9273f709218"
N_TOTAL=512
TARGET=128
SELECT_SEED="MOIS_BOTANICAL_QC_128_OF_512_ORIGINAL9_ROI_FAILURES_20261027"

def require(ok,message):
    if not bool(ok):raise ValueError("MOIS_128_BOTANICAL_REVIEW: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def rank(item):
    return hashlib.sha256((SELECT_SEED+"|"+item["photo_id"]).encode()).hexdigest()

def choose(rows,*,target=TARGET):
    require(len(rows)>=target and target>0 and target<=N_TOTAL,
            "Not enough original photo cases for botanist review")
    photo=set();taxa=set();obs=set()
    for r in rows:
        a,b,c=r["photo_id"],r["observation_id"],r["inat_taxon_id"]
        require(a.isdecimal() and a not in photo and
                b and b not in obs and c and c not in taxa and
                r["old_ROI9_missing"]=="True" and
                r["screen_status"] in ("candidate","negative","uncertain","error"),
                "Review photos must have unique original source identities and technical unknown")
        photo.add(a);taxa.add(c);obs.add(b)
    selected=[];used=set()
    def take(items,cap):
        for row in sorted(items,key=rank):
            if len(selected)>=target or cap<=0:break
            if row["photo_id"] in used:continue
            selected.append(row);used.add(row["photo_id"]);cap-=1
    # Independently check detector negatives/errors before treating no-box
    # as evidence of no flower. This is deliberately NOT a random flora sample.
    negative=[r for r in rows if r["screen_status"] in ("negative","error","uncertain")]
    require(len(negative)<=48,
            "Detector misses/technical errors far exceed planned review budget; revise design")
    take(negative,48)
    large=[r for r in rows if r["technical_diagnostic_class"]==
           "detectable_large_flower_box_ROI9_failed"]
    take(large,min(64,target-len(selected)))
    # Fill remaining from all source designs and geography, selecting
    # underrepresented strata before repeatedly sampling a region.
    buckets=defaultdict(list)
    for row in rows:
        if row["photo_id"] in used:continue
        key=(row["source_design"],row["latitude_zone"],
             row["postoutcome_source_gap"])
        buckets[key].append(row)
    for bucket in buckets.values():
        bucket.sort(key=rank)
    cursor={key:0 for key in buckets}
    while len(selected)<target:
        did=False
        for key in sorted(buckets,key=lambda k:(
                sum(x["source_design"]==k[0] for x in selected),
                sum(x["latitude_zone"]==k[1] for x in selected),
                k)):
            arr=buckets[key];pos=cursor[key]
            if pos<len(arr):
                row=arr[pos];cursor[key]+=1
                if row["photo_id"] not in used:
                    selected.append(row);used.add(row["photo_id"]);did=True
                    if len(selected)==target:break
        if not did:break
    require(len(selected)==target and len(used)==target and
            len(set(x["inat_taxon_id"] for x in selected))==target,
            "128 botanist-review IDs or original plant taxon IDs duplicated/missing")
    return selected

HTML = r'''<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>MOIS | 原画像の花器官・ROI失敗 確認128枚</title>
<style>
:root{font-family:system-ui,-apple-system,"Noto Sans JP",sans-serif;color:#182532;background:#f7f9fa}
*{box-sizing:border-box}
main{max-width:1000px;margin:0 auto;padding:1rem}
h1{font-size:1.25rem;margin:0 0 .5rem}
p,.muted{line-height:1.55;color:#506175}
article{border:1px solid #d6dee6;border-radius:14px;background:#fff;padding:1rem}
.grid{display:grid;grid-template-columns:minmax(0,1.4fr) minmax(0,1fr);gap:1rem}
@media(max-width:760px){.grid{grid-template-columns:1fr}}
img{max-width:100%;max-height:58vh;border-radius:8px;object-fit:contain;background:#edf1f3}
.preview{display:flex;align-items:center;justify-content:center;min-height:270px}
label{display:block;margin:.8rem 0 .3rem;font-weight:600}
select,textarea,button{font:inherit;width:100%;padding:.65rem;border:1px solid #bcc8d0;border-radius:8px;background:white}
textarea{min-height:4.5rem}
button{cursor:pointer;background:#eef4f8}
button.primary{background:#183c59;color:#fff}
.controls{display:grid;grid-template-columns:repeat(3,1fr);gap:.5rem;margin-top:1rem}
.details{font-variant-numeric:tabular-nums;font-size:.91rem;color:#34495c}
code{font-size:.86rem;overflow-wrap:anywhere}
.progress{background:#e1e6ea;border-radius:99px;height:8px;overflow:hidden;margin:.8rem 0}
.progress>div{height:100%;width:0;background:#188771}
.notice{font-size:.9rem;color:#6d5364}
</style></head><body><main>
<h1>MOIS — 9色ROI判定不能画像の花器官確認</h1>
<p>既存写真128枚・128分類群の<strong>人手による検証</strong>です。YOLOの花候補枠は植物学的に確認した花弁ではありません。色・ROIの判定を自動入力しません。</p>
<div class="progress"><div id="bar"></div></div>
<p id="counter" aria-live="polite"></p>
<div class="grid">
<article><div class="preview"><img id="photo" alt="未判定の原画像"/></div>
<p id="loadstatus" class="notice">画像はiNaturalist公開画像を表示します。元画像バイナリはこのHTMLに含みません。</p>
<div class="details" id="details"></div></article>
<article>
<label for="flower">開花中の花が写っている？</label>
<select id="flower"><option value="">未評価</option><option value="yes">はい</option><option value="no">いいえ</option><option value="uncertain">不確実</option></select>
<label for="corolla">花冠・花弁の色を写真から評価できる？</label>
<select id="corolla"><option value="">未評価</option><option value="yes">はい</option><option value="no">いいえ</option><option value="uncertain">不確実</option></select>
<label for="problem">最も疑われる測定失敗要因</label>
<select id="problem">
<option value="">未評価</option>
<option value="no_flower">花器官なし／開花していない</option>
<option value="occluded_or_tiny">花器官が小さい・遮蔽</option>
<option value="background_or_wrong_object">検出対象が花ではない／背景の影響</option>
<option value="flower_visible_roi_gate">花は見えるがROI・flipゲート不成功</option>
<option value="illumination_or_color">照明・露出・色の品質問題</option>
<option value="photo_inaccessible">画像取得不能</option>
<option value="other">その他</option>
<option value="uncertain">判断不可</option>
</select>
<label for="memo">メモ（任意）</label>
<textarea id="memo" placeholder="重要な観察や判定の理由"></textarea>
<p class="muted">この評価はブラウザ上にのみ一時保存されます。<strong>「CSV保存」で定期的に書き出してください。</strong> 写真ID・植物ID・技術的検出結果は変更されません。</p>
<div class="controls"><button id="prev">前へ</button><button id="next">次へ</button><button class="primary" id="export">CSV保存</button></div>
</article></div>
<p class="notice">確認対象は技術的検出陰性・エラーおよび花候補大のROI失敗画像を意図的に多めに含みます。統計的な全球花色・開花率の推定標本ではありません。</p>
</main>
<script>
"use strict";
const roster=__SOURCE_ROWS__;
const version="mois-botanical-128-v1:__SOURCE_SHA__";
let current=0;
let state={};
try{state=JSON.parse(localStorage.getItem(version)||"{}")||{};}catch(e){state={};}
function persist(){try{localStorage.setItem(version,JSON.stringify(state));}catch(e){}}
function record(){
 const r=roster[current];if(!r)return;
 state[r.photo_id]={
  flower:document.getElementById("flower").value,
  corolla:document.getElementById("corolla").value,
  failure:document.getElementById("problem").value,
  memo:document.getElementById("memo").value
 };
 persist();count();
}
function count(){
 const n=roster.filter(x=>(state[x.photo_id]||{}).flower).length;
 document.getElementById("bar").style.width=(100*n/roster.length)+"%";
 document.getElementById("counter").textContent="判定済み "+n+"/"+roster.length+" 件 — "+(current+1)+"枚目";
}
function view(index){
 current=Math.max(0,Math.min(roster.length-1,index));
 const r=roster[current],val=state[r.photo_id]||{};
 document.getElementById("flower").value=val.flower||"";
 document.getElementById("corolla").value=val.corolla||"";
 document.getElementById("problem").value=val.failure||"";
 document.getElementById("memo").value=val.memo||"";
 document.getElementById("details").textContent=
  "Photo ID: "+r.photo_id+" | Taxon ID: "+r.inat_taxon_id+
  " | "+r.latitude_zone+" / "+r.local_season+
  " | "+r.source_design+" | 事前診断: "+r.technical_diagnostic_class;
 const pic=document.getElementById("photo");
 const status=document.getElementById("loadstatus");
 pic.removeAttribute("src");
 let variant=0;
 const exts=["jpeg","jpg","png"];
 const src=function(){return "https://inaturalist-open-data.s3.amazonaws.com/photos/"+r.photo_id+"/medium."+exts[variant];};
 pic.onerror=function(){
  variant+=1;
  if(variant<exts.length){pic.src=src();}
  else {status.textContent="写真を表示できません。取得不能と記録してください。";}
 };
 pic.onload=function(){status.textContent="既存の公開元画像。花・花冠は目視評価してください。";};
 status.textContent="画像を読み込み中…";pic.src=src();count();
}
document.querySelectorAll("select,textarea").forEach(function(e){
 e.addEventListener("change",record);
 e.addEventListener("input",record);
});
document.getElementById("prev").onclick=function(){record();view(current-1);};
document.getElementById("next").onclick=function(){record();view(current+1);};
function csvEscape(value){
 return '"'+String(value==null?"":value).replace(/"/g,'""')+'"';
}
document.getElementById("export").onclick=function(){
 record();
 const hdr=["photo_id","observation_id","inat_taxon_id","source_design",
  "latitude_zone","local_season","prior_original_strict9_unknown",
  "detector_screen_status","detector_technical_class",
  "human_flower_visible","human_corolla_colour_evaluable",
  "human_failure_reason","human_notes"];
 const records=[hdr];
 for(const r of roster){
  const v=state[r.photo_id]||{};
  records.push([r.photo_id,r.observation_id,r.inat_taxon_id,r.source_design,
    r.latitude_zone,r.local_season,"True",r.screen_status,
    r.technical_diagnostic_class,v.flower||"",v.corolla||"",v.failure||"",v.memo||""]);
 }
 const content=records.map(row=>row.map(csvEscape).join(",")).join("\r\n");
 const blob=new Blob(["\ufeff"+content],{type:"text/csv;charset=utf-8"});
 const u=URL.createObjectURL(blob),a=document.createElement("a");
 a.href=u;a.download="mois_128_original_roi9_botanical_human_review.csv";
 a.click();setTimeout(function(){URL.revokeObjectURL(u);},1000);
};
view(0);
</script></body></html>'''

def create(input_csv,out,*,validate_sha=True,target=TARGET):
    input_csv,out=map(Path,(input_csv,out))
    require(not out.exists(),"Do not overwrite an existing botanical review roster")
    if validate_sha:
        require(sha(input_csv)==SOURCE_SHA,
                "512 original photo status data or release SHA changed")
    with input_csv.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    require(len(rows)==N_TOTAL,
            "Original 512 technical ROI unknown photos incomplete")
    picks=choose(rows,target=target)
    out.mkdir(parents=True)
    fields=("photo_id","observation_id","inat_taxon_id","source_design",
        "latitude_zone","local_season","postoutcome_source_gap",
        "screen_status","technical_diagnostic_class","review_priority")
    chosen=[]
    for row in picks:
        result={k:row.get(k,"") for k in fields}
        result["review_priority"]=(
            "DETECTOR_NEGATIVE_OR_ERROR" if row["screen_status"] in
                ("negative","error","uncertain") else
            "LARGE_VISIBLE_BOX_BUT_NO_STRICT9" if row["technical_diagnostic_class"]==
                "detectable_large_flower_box_ROI9_failed" else
            "OTHER_ROI9_UNKNOWN_FOR_BALANCED_REVIEW"
        )
        chosen.append(result)
    def write(name,records,keys):
        path=out/name
        with path.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n")
            w.writeheader();w.writerows(records)
        return path
    roster=write("manual128_existing_original_photoID_review_roster.csv",
                 chosen,tuple(chosen[0]))
    labelled=[
        {**r,"human_flower_visible":"","human_corolla_colour_evaluable":"",
         "human_failure_reason":"","human_notes":"","human_review_completed":False}
        for r in chosen]
    label=write("manual128_botanical_annotation_BLANK_not_yet_reviewed.csv",
                labelled,tuple(labelled[0]))
    payload=json.dumps(chosen,ensure_ascii=False,separators=(",",":")).replace(
        "<","\\u003c").replace(">","\\u003e")
    page=HTML.replace("__SOURCE_ROWS__",payload).replace("__SOURCE_SHA__",sha(roster))
    html=out/"mois_128_botanical_flower_ROI9_failure_annotation.html"
    html.write_text(page,encoding="utf-8")
    classes=Counter(r["review_priority"] for r in chosen)
    source=Counter(r["source_design"] for r in chosen)
    report={
        "schema":SCHEMA,
        "original_existing_technical_ROI9_unknown_source_photos":N_TOTAL,
        "botanical_review_original_photo_IDs_selected":len(chosen),
        "distinct_original_taxa":len(set(r["inat_taxon_id"] for r in chosen)),
        "purposefully_oversampled_detector_misses_and_large_box_failures":True,
        "case_type_counts":dict(classes),
        "source_acquisition_design_counts":dict(source),
        "all_botanical_open_flower_or_petal_labels_still_blank":True,
        "original_strict9_unknown_not_modified":True,
        "new_photo_IDs_added_to_census":0,
        "image_bytes_persisted_or_redistributed_by_tool":False,
        "HTML_only_displays_public_source_URLs_lazily_one_photo_at_a_time":True,
        "browser_form_export_is_local_not_server_submitted":True,
        "input_original_512_technical_diagnostic_SHA256":sha(input_csv),
        "outputs_SHA256":{p.name:sha(p) for p in (roster,label,html)}
    }
    (out/"botanical128_review_design_receipt.json").write_text(
        json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("MOIS_REAL_PRIOR_UNKNOWN128_BOTANICAL_REVIEW_ROSTER_READY",
          json.dumps({k:report[k] for k in (
            "botanical_review_original_photo_IDs_selected",
            "distinct_original_taxa","case_type_counts",
            "source_acquisition_design_counts")},sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--diagnosis-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    x=p.parse_args()
    create(x.diagnosis_csv,x.output)

if __name__=="__main__":
    main()
