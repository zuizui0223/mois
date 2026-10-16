#!/usr/bin/env node
// Deterministic source-backed synthesis of photo-white parity and the ecological taxonomic-scale rainfall signal.
// Numerical inputs are frozen aggregate original-source receipts; no new source photo or regression fits.
import fs from "node:fs";
import path from "node:path";
import {fileURLToPath} from "node:url";
const ROOT=path.resolve(path.dirname(fileURLToPath(import.meta.url)),"..");
function extract(fr,sc,si){
if(fr.n_classified_photos!==45953 || fr.n_classified_white_photos!==21574 ||fr.n_species_with_classified_photo!==25064||fr.n_genera_with_classified_photo!==5369||si.source_support.n_primary_species!==20694 || si.source_support.n_primary_genera!==1934 ||sc.source_cohort.n_primary_taxa!==20694)throw Error("Source cohort invariants");
const vals=[
 {label:"Equal photograph",denominator:fr.n_classified_photos,white:fr.point_white_proportion.photo_equal},
 {label:"Equal sampled species",denominator:fr.n_species_with_classified_photo,white:fr.point_white_proportion.species_equal},
 {label:"Equal sampled genus",denominator:fr.n_genera_with_classified_photo,white:fr.point_white_proportion.genus_equal},
 {label:"Genus with ≥5 species",denominator:fr.min_classified_species_thresholds.find(q=>q.min_classified_species_per_genus===5).n_genera,white:fr.min_classified_species_thresholds.find(q=>q.min_classified_species_per_genus===5).genus_equal_white}
];
const slopes=[{axis:"BIO12 · rainfall",key:"wc_bio12"},{axis:"Annual vapour",key:"wc_vapr_annual_mean"}].flatMap(({axis,key})=>{
const e=sc.estimates.axes[key];
return [{label:axis+" · within",axis,scale:"within genus",value:e.within_nominal_genus.point_pp_per_original_photo_sd,CI:e.within_nominal_genus.whole_genus_percentile_CI95_pp},
 {label:axis+" · between",axis,scale:"genus means",value:e.between_genus_means.point_pp_per_original_photo_sd,CI:e.between_genus_means.whole_genus_percentile_CI95_pp}];});
const sig=[{label:"Rainfall · BIO12",key:"wc_bio12"},{label:"Annual vapour",key:"wc_vapr_annual_mean"}].map(({label,key})=>{
 const d=si.cases["3"].axes[key];
 return {label,total:si.cases["3"].n_genera,pos:d.n_positive_contributing_genera,neg:d.n_negative_contributing_genera,zero:d.n_no_within_genus_photo_white_variation,positiveAmongVariable:d.fraction_positive_among_variable_genera};
});
if(Math.abs(slopes[0].value-2.254545235)>1e-8||Math.abs(slopes[1].value-8.554761564)>1e-8||sig[0].pos!==830||sig[0].neg!==719||sig[0].zero!==385)throw Error("Scope-paired axes or signed-genus invariants");
return {white:vals,slopes,signs:sig,nSpecies:sc.source_cohort.n_primary_taxa,nG:sc.source_cohort.n_primary_genera,nAll:fr.n_unique_photo_opportunities,source:"original photographed plant source, 2026-10-15; no new photo/fitness measurements"};
}
function render(q){
const out=[];
const esc=s=>String(s).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll(">","&gt;").replaceAll('"',"&quot;");
const rect=(x,y,w,h,fill,rx=0)=>out.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${rx}" fill="${fill}"/>`);
const line=(x,y,X,Y,c="#D7E1E6",w=1,dash="")=>out.push(`<line x1="${x}" y1="${y}" x2="${X}" y2="${Y}" stroke="${c}" stroke-width="${w}"${dash?` stroke-dasharray="${dash}"`:""}/>`);
const txt=(x,y,s,size=18,weight=400,fill="#243746",anchor="start")=>out.push(`<text x="${x}" y="${y}" font-family="Arial,Helvetica,sans-serif" font-size="${size}" font-weight="${weight}" fill="${fill}" text-anchor="${anchor}">${esc(s)}</text>`);
const cir=(x,y,r,color)=>out.push(`<circle cx="${x}" cy="${y}" r="${r}" fill="${color}" stroke="#FFFFFF" stroke-width="2"/>`);
const blue="#3B728C",purple="#8F719A",teal="#25766D",coral="#B76A54",gray="#A7B3BD",dark="#263849",muted="#546878";
out.push(`<svg xmlns="http://www.w3.org/2000/svg" width="1560" height="870" viewBox="0 0 1560 870" role="img" aria-labelledby="f5-title f5-desc">`);
out.push(`<title id="f5-title">Global photographed flower white balance and contrasting taxonomic scale moisture associations</title>`);
out.push(`<desc id="f5-desc">Panel A shows three differently weighted near-half photographic white proportions and a stricter genus sample. Panel B shows positive within-genus and between-genus photographed white slopes for precipitation and atmospheric vapour, with ordinary whole-genus bootstrap 95% intervals. Panel C shows positive, negative and invariant nominal genus contributions underlying worldwide positive slopes. Neither photograph frequencies nor coefficients measure genetic colour frequencies, natural selection or plant fitness.</desc>`);
rect(0,0,1560,870,"#FFFFFF");rect(0,0,1560,14,"#2F5870");
txt(49,69,"GLOBAL NEAR-PARITY, DIFFERENT ECOLOGICAL SCALES",29,700,dark);
txt(49,110,"Source-classified flower photographs · distinct statistical targets · no inferred pigment fitness",18,400,muted);
const panels=[{x:48,w:455},{x:530,w:590},{x:1146,w:365}];
panels.forEach(p=>rect(p.x,151,p.w,601,"#F8FAFC",12));
// PANEL A
const A=panels[0],ax=A.x+27,bw=401;
txt(A.x+23,196,"A   White : coloured observations",23,700,dark);
txt(A.x+23,229,"Different weighting, not independent datasets",15,400,muted);
rect(A.x+25,248,14,14,blue,3);txt(A.x+47,261,"White",15,500,dark);
rect(A.x+141,248,14,14,purple,3);txt(A.x+163,261,"Other visible hues",15,500,dark);
q.white.forEach((d,i)=>{
const top=305+i*100,fw=bw*d.white;
txt(A.x+26,top-12,d.label,19,600,dark);
txt(A.x+428,top-12,(100*d.white).toFixed(2)+"%",19,700,blue,"end");
rect(ax,top,fw,32,blue,5);rect(ax+fw,top,bw-fw,32,purple,5);
line(ax+bw*.5,top-5,ax+bw*.5,top+38,"#FFFFFF",2);
txt(A.x+28,top+54,`n = ${d.denominator.toLocaleString("en-US")}`,15,400,muted);
});
txt(A.x+25,729,"The genus ≥5 cohort excludes under-sampled genera.",14,500,muted);
// PANEL B
const B=panels[1],pmin=B.x+252,pmax=B.x+555,pw=pmax-pmin,coord=v=>pmin+(v/14)*pw;
txt(B.x+23,196,"B   Moisture: within vs. between genera",23,700,dark);
txt(B.x+23,229,`${q.nSpecies.toLocaleString("en-US")} species in ${q.nG.toLocaleString("en-US")} genera · identical source taxa`,16,400,muted);
for(const val of [0,2,4,6,8,10,12,14]){const xx=coord(val);line(xx,280,xx,651,val===0?"#597182":"#DEE6EB",val===0?2:1);txt(xx,688,String(val),15,400,muted,"middle");}
q.slopes.forEach((d,i)=>{
const y=315+98*i,color=i%2===0?blue:teal;
txt(B.x+24,y-10,(d.axis.startsWith("BIO12")?"Rain":"Vapour")+" · "+d.scale,18,600,dark);
const lo=coord(d.CI[0]),hi=coord(d.CI[1]),x=coord(d.value);
line(lo,y,hi,y,color,5);line(lo,y-11,lo,y+11,color,2);line(hi,y-11,hi,y+11,color,2);cir(x,y,7,color);
txt(B.x+24,y+35,`${d.value>=0?"+":""}${d.value.toFixed(3)} [${d.CI[0].toFixed(3)}, ${d.CI[1].toFixed(3)}]`,15,400,muted);
});
txt(B.x+25,719,"White-photo change (percentage points / source SD)",16,500,muted);
// PANEL C
const C=panels[2],cx=C.x+22,cw=320;
txt(C.x+22,196,"C   Genus contributions vary",22,700,dark);
txt(C.x+22,229,"Same 1,934 qualified nominal genera",15,400,muted);
q.signs.forEach((d,i)=>{
const y=326+i*194;
txt(C.x+23,y-35,d.label,21,600,dark);
let posW=cw*d.pos/d.total,negW=cw*d.neg/d.total,zeroW=cw*d.zero/d.total;
rect(cx,y,posW,55,teal,4);rect(cx+posW,y,negW,55,coral,0);rect(cx+posW+negW,y,zeroW,55,gray,4);
txt(cx+posW/2,y+35,String(d.pos),18,700,"#FFFFFF","middle");
txt(cx+posW+negW/2,y+35,String(d.neg),18,700,"#FFFFFF","middle");
txt(cx+posW+negW+zeroW/2,y+35,String(d.zero),18,700,dark,"middle");
txt(C.x+23,y+88,`${(100*d.positiveAmongVariable).toFixed(2)}% positive among variable genera`,15,500,muted);
});
rect(C.x+20,664,12,12,teal,2);txt(C.x+38,677,"Positive",15,500,dark);
rect(C.x+131,664,12,12,coral,2);txt(C.x+149,677,"Negative",15,500,dark);
rect(C.x+249,664,12,12,gray,2);txt(C.x+266,677,"No variation",15,500,dark);
txt(C.x+22,721,"Signs are genus contributions, not fitness.",14,500,muted);
// figure-wide footer
line(48,782,1510,782,"#C9D5DE",1.5);
txt(50,815,"A: all classifiable photographs, taxa or genera;  B: ordinary 95% genus-bootstrap intervals on identical climate26 taxa (not original maxT).",16,400,muted);
txt(50,843,"C: signs of genus cross-product contributions; observations are nonrandom and do not identify stable genetic morph ratios or causal benefits.",16,400,muted);
out.push("</svg>");return out.join("\n")+"\n";
}
const paths=["results/mois_original_photo_species_genus_balance_20261014.json","results/mois_genus_between_within_moisture_20261015.json","results/mois_congeneric_moisture_genus_sign_mosaic_20261015.json"];
const data=paths.map(p=>JSON.parse(fs.readFileSync(path.join(ROOT,p),"utf8")));
const svg=render(extract(...data));
const output=path.join(ROOT,"figures/mois_white_moisture_ecological_scales.svg");
if(process.argv.includes("--check")){
 if(!fs.existsSync(output)||fs.readFileSync(output,"utf8")!==svg)throw Error("Figure 5 disagrees with source receipts, regenerated SVG, or frozen cohort invariants");
 console.log("PASS Figure 5 photographed global near-parity, identical-taxon slopes, heterogeneous genus contributions");
}else if(process.argv.includes("--write")){
 fs.writeFileSync(output,svg,"utf8");console.log("WROTE source-verified Figure 5",output);
}else console.log("Usage: node scripts/render_mois_figure5_ecological_scales.mjs --check | --write");
