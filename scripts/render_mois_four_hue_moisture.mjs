#!/usr/bin/env node
// Replicate Figure 4 EXACTLY from the original PR193 1,280-row CSV;
// no new photographic, bootstrap, or multiple-testing calculation.
import fs from "node:fs";
import path from "node:path";
import {fileURLToPath} from "node:url";
const ROOT=path.resolve(path.dirname(fileURLToPath(import.meta.url)),"..");
function extract(csv){
 const lines=csv.trim().split(/\r?\n/),names=lines[0].split(",");
 if(lines.length!==1281)throw Error("Expected 1280 original row results");
 const rows=lines.slice(1).map(line=>Object.fromEntries(line.split(",").map((v,i)=>[names[i],v])));
 const hue=["white","yellow_orange","red_pink","blue_purple"],axes=["wc_bio12","wc_vapr_annual_mean"];
 const panels=axes.map((axis)=>{
  const entries=hue.map(h=>{
    const hits=rows.filter(q=>q.mode==="climate"&&q.quality==="all"&&q.region==="world"&&q.adjustment==="genus_geo"&&q.weight==="genus_equal"&&q.anchor===axis&&q.hue===h);
    if(hits.length!==1)throw Error("Missing original source hue "+axis+" "+h);
    const q=hits[0];
    const d={hue:h,value:+q.direction_pp_per_env_global_sd,lo:+q.maxT_ci_lo_pp,hi:+q.maxT_ci_hi_pp,significant:q.maxT_excludes_zero==="True",species:+q.n_species,genera:+q.n_genera};
    if(d.species!==20694||d.genera!==1934||!(d.lo<=d.value&&d.value<=d.hi))throw Error("Invalid hue interval source");
    if(d.significant!==(d.lo>0||d.hi<0))throw Error("Source FWER support indicator mismatch");
    if(d.significant !== (h==="white"))throw Error("Unexpected individual hue effect: "+axis+" "+h);
    return d;
  });
  if(Math.abs(entries.reduce((a,b)=>a+b.value,0))>1e-9)throw Error("Compositional hue slopes fail to sum to zero");
  return {axis,entries};
 });
 return {schema:"mois_original_four_hue_moisture_figure4_v1",recordSource:"source FCP PR193 exact 1280 direction ledger; no new inferential fit",nOriginalRows:1280,nSpecies:20694,nGenera:1934,panels};
}
function render(q){
 const o=[],str=x=>String(x).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll(">","&gt;");
 const R=(x,y,w,h,col,rx=0)=>o.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${rx}" fill="${col}"/>`);
 const T=(x,y,s,size=16,bold=400,color="#243649",anchor="start")=>o.push(`<text x="${x}" y="${y}" font-family="Arial,Helvetica,sans-serif" font-size="${size}" font-weight="${bold}" fill="${color}" text-anchor="${anchor}">${str(s)}</text>`);
 const L=(x,y,X,Y,col="#C7D2D9",w=1)=>o.push(`<line x1="${x}" y1="${y}" x2="${X}" y2="${Y}" stroke="${col}" stroke-width="${w}"/>`);
 const C=(x,y,r,col,stroke="#FFFFFF",sw=1)=>o.push(`<circle cx="${x}" cy="${y}" r="${r}" fill="${col}" stroke="${stroke}" stroke-width="${sw}"/>`);
 const colors={white:"#77838D",yellow_orange:"#D49D32",red_pink:"#BB6775",blue_purple:"#7064A2"};
 const labels={white:"White",yellow_orange:"Yellow / orange",red_pink:"Red / pink",blue_purple:"Blue / purple"};
 o.push(`<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="f4title f4desc" width="1450" height="850" viewBox="0 0 1450 850">`);
 o.push(`<title id="f4title">All four photographed flower-colour responses to global precipitation and vapour pressure</title>`);
 o.push(`<desc id="f4desc">Original 20,694 species genus-equal climate analysis. White has a positive simultaneous confidence interval excluding zero for both water variables. The other three hues show individual intervals including zero. Source-controlled, not a universal causal pigment effect.</desc>`);
 R(0,0,1450,850,"#FFFFFF");R(0,0,1450,14,"#345670");
 T(54,66,"FOUR PHOTO-HUE RESPONSES TO MOISTURE",28,700);T(54,101,"Genus-equal world climate26 · 20,694 species / 1,934 genera · original 24-contrast maxT intervals",17,400,"#526778");
 const panelConfig=[{x:54,name:"A  Annual precipitation (BIO12)"},{x:745,name:"B  Atmospheric vapour pressure"}];
 q.panels.forEach((panel,i)=>{
  const px=panelConfig[i].x,lo=px+235,hi=px+630,axMin=-5,axMax=8,ax=v=>lo+(v-axMin)/(axMax-axMin)*(hi-lo);
  R(px,149,653,560,"#F8FAFC",13);
  T(px+21,189,panelConfig[i].name,21,700);
  T(px+21,219,"Photo-colour share: change (percentage points per source SD)",15,400,"#526778");
  [-4,-2,0,2,4,6,8].forEach(tick=>{const xx=ax(tick);L(xx,282,xx,596,tick===0?"#596B77":"#D9E1E6",tick===0?2:1);T(xx,627,String(tick),14,400,"#526778","middle");});
  panel.entries.forEach((z,j)=>{
   const y=315+j*83,color=colors[z.hue],sig=z.significant;
   T(px+22,y-5,labels[z.hue],18,sig?700:500);
   L(ax(z.lo),y,ax(z.hi),y,color,5);
   L(ax(z.lo),y-10,ax(z.lo),y+10,color,2);L(ax(z.hi),y-10,ax(z.hi),y+10,color,2);
   C(ax(z.value),y,7,color);
   T(px+22,y+26,(z.value>=0?"+":"")+z.value.toFixed(3)+" ["+z.lo.toFixed(3)+", "+z.hi.toFixed(3)+"]",14,sig?600:400,"#526778");
  });
  T(px+22,676,"Individual CI excludes zero: white only",15,600,"#335874");
 });
 L(53,744,1398,744);
 T(54,773,"Photo probabilities sum to 1: colour-class slopes sum to 0; negative point estimates for coloured states do not individually pass maxT.",15,400,"#526778");
 T(54,803,"Original photo labels and bootstrap uncertainty unchanged; study-wide/causal significance is NOT established.",15,400,"#526778");
 o.push("</svg>");return o.join("\n")+"\n";
}
const original=fs.readFileSync(path.join(ROOT,"results/fcp_macro_direction_20261013/direction_ledger.csv"),"utf8");
const image=render(extract(original));
const output=path.join(ROOT,"figures/mois_original_four_hue_moisture.svg");
if (process.argv.includes("--check")){
 if(!fs.existsSync(output)||fs.readFileSync(output,"utf8")!==image) throw new Error("Figure 4 differs from original 1,280 coefficient source");
 console.log("PASS Figure 4: source-backed four-hue coefficients, within-family uncertainty, 20,694 species");
}else if(process.argv.includes("--write")){
 fs.writeFileSync(output,image,"utf8");console.log("Wrote source-consistent Figure 4",output);
}else {console.log("Usage: node scripts/render_mois_four_hue_moisture.mjs --check | --write");}
