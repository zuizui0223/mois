#!/usr/bin/env node
// Reporting-only Figure 2 generator. Source JSON identities and numbers are fixed.
import fs from "node:fs";
import path from "node:path";
import {fileURLToPath} from "node:url";
const ROOT=path.resolve(path.dirname(fileURLToPath(import.meta.url)),"..");
function extract(a,b) {
 const m=b.matched_source_photo_and_species,bin=a.rain_bin_all_photo_source_counts,z=a.original_photo_wet_dry_unweighted;
 const pp=x=>100*x;
 const item=(i)=>{const n=bin.find(q=>q.rain_decile===i);return {total:n.n,classified:n.n_class,white:z[String(i)].white_fraction};};
 const o={total:a.all_original_unique_photo_ids,classified:a.source_all_classifiable,unknown:a.source_all_unclassifiable,dry:item(0),wet:item(9),n_photos:m.n_common_original_photos,n_species:m.n_common_species,n_genera:m.n_source_genera,slopes:[
  {name:"Between congeneric species",point:pp(m.between_congener_BIO12_photo_white_slope),range:m.paired_genus_bootstrap95_between.map(pp)},
  {name:"Within named species",point:pp(m.within_same_species_BIO12_photo_white_slope),range:m.paired_genus_bootstrap95_within.map(pp)},
  {name:"Between minus within",point:pp(m.between_minus_within_slope),range:m.paired_genus_bootstrap95_difference.map(pp)}
 ]};
 if(o.total!==100543||o.classified!==45953||o.unknown!==54590||o.dry.total!==10084||o.wet.total!==9990||o.n_photos!==21448||o.n_species!==6143)throw Error("Source audit mismatch");
 return o;
}
function render(o) {
 const out=['<svg xmlns="http://www.w3.org/2000/svg" width="1440" height="800" viewBox="0 0 1440 800" role="img" aria-labelledby="fig2-title fig2-desc">','<title id="fig2-title">Photographic selection and matched rainfall contrasts</title>','<desc id="fig2-desc">Original data show 100543 photos, 45953 classifiable. Wet dry classes differ in classifiability. The original matched photo contrast shows a positive between species slope but not within species.</desc>'];
 const R=(x,y,w,h,fill,rx=0)=>out.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${rx}" fill="${fill}"/>`);
 const T=(x,y,str,size=18,weight=400,anchor="start",color="#1E3447")=>out.push(`<text x="${x}" y="${y}" font-family="Arial,sans-serif" font-size="${size}" font-weight="${weight}" fill="${color}" text-anchor="${anchor}">${str}</text>`);
 const L=(x,y,X,Y,c="#CAD6DD",w=1)=>out.push(`<line x1="${x}" y1="${y}" x2="${X}" y2="${Y}" stroke="${c}" stroke-width="${w}"/>`);
 const round=(x,n=1)=>x.toFixed(n);
 const perc=(x)=>round(x*100,1)+"%";
 const a=o.classified/o.total;
 R(0,0,1440,800,"#FFFFFF");R(0,0,1440,14,"#315773");
 T(50,58,"VISIBLE FLOWER-COLOUR OPPORTUNITY AND MOISTURE",27,700);
 T(50,92,"Source-locked reporting only | original photographs, no re-estimation",17,400,"start","#5A6A78");
 [[45,130,425,565],[506,130,425,565],[968,130,425,565]].forEach(v=>R(...v,"#F6F9FB",14));
 T(69,174,"A  Photographic classification",20,700);
 T(69,210,o.total.toLocaleString("en-US")+" unique original photos",17);
 R(69,242,370,44,"#DCE5EA",8);R(69,242,370*a,44,"#347E9B",8);
 T(69,322,perc(a)+" classifiable",20,700);T(69,355,o.classified.toLocaleString("en-US")+" four-hue photo labels",16);
 T(69,417,perc(o.unknown/o.total)+" unclassifiable",20,700);
 T(69,451,o.unknown.toLocaleString("en-US")+" labels remain missing",16);
 L(69,512,443,512);
 T(69,551,"Unknown flower colour is not",17);T(69,576,"coded white or nonwhite.",17);
 T(530,174,"B  Wet–dry selection contrast",20,700);
 T(530,210,"All-photo rainfall-opportunity deciles",16,400,"start","#5A6A78");
 [{title:"Driest",...o.dry},{title:"Wettest",...o.wet}].forEach((q,i)=>{
  const y=260+i*196;
  T(530,y,q.title+" (n="+q.total.toLocaleString("en-US")+")",19,700);
  R(530,y+17,355,18,"#DCE5EA",6);R(530,y+17,355*q.classified/q.total,18,"#347E9B",6);
  T(530,y+59,perc(q.classified/q.total)+" photos classifiable",16);
  R(530,y+80,355,18,"#DCE5EA",6);R(530,y+80,355*q.white,18,"#9784B7",6);
  T(530,y+123,perc(q.white)+" white among classifiable",16);
 });
 T(990,174,"C  Same original photographs",20,700);
 T(990,210,o.n_photos.toLocaleString("en-US")+" photos | "+o.n_species.toLocaleString("en-US")+" species",16);
 T(990,242,"BIO12 white slope (pp / source SD)",15);
 const x0=1024,x1=1368,min=-2,max=6.5,px=v=>x0+(v-min)*(x1-x0)/(max-min);
 [-2,0,2,4,6].forEach(v=>{const xx=px(v);L(xx,282,xx,630,v===0?"#6A7884":"#D7DFE6",v===0?2:1);T(xx,654,String(v),13,400,"middle","#637281");});
 o.slopes.forEach((q,i)=>{
  const y=342+i*105,c=["#307DA2","#AD6B6E","#315773"][i],lo=q.range[0],hi=q.range[1];
  R(985,y-43,365,27,"#F6F9FB",4);
  T(989,y-24,q.name,17,700);
  L(px(lo),y,px(hi),y,c,5);L(px(lo),y-10,px(lo),y+10,c,2);L(px(hi),y-10,px(hi),y+10,c,2);
  out.push(`<circle cx="${px(q.point)}" cy="${y}" r="7" fill="${c}"/>`);
  R(985,y+16,305,26,"#F6F9FB",4);
  T(989,y+35,(q.point>=0?"+":"")+round(q.point,3)+" ["+round(lo,3)+", "+round(hi,3)+"]",15,400,"start","#526576");
 });
 T(991,675,"Original paired genus-bootstrap 95% CIs",14,400,"start","#526576");
 L(49,721,1395,721);
 T(50,750,"Conditional photographed flower colour, not pigment chemistry or natural-population morph frequencies.",15,400,"start","#566A7B");
 T(50,778,"Retrospective between/within slopes use identical photo IDs, but different statistical estimands.",15,400,"start","#566A7B");
 out.push("</svg>");return out.join("\n")+"\n";
}
const m=JSON.parse(fs.readFileSync(path.join(ROOT,"results/source_receipts/photo_colour_missingness_mnar.json"),"utf8"));
const p=JSON.parse(fs.readFileSync(path.join(ROOT,"results/source_receipts/matched_congeneric_vs_withinspecies_rainfall.json"),"utf8"));
const expected=render(extract(m,p));
const dest=path.join(ROOT,"figures/mois_photo_selection_and_matched_rainfall.svg");
if(process.argv.includes("--check")){
 if(!fs.existsSync(dest)||fs.readFileSync(dest,"utf8")!==expected)throw Error("Figure 2 does not match archived source receipts");
 console.log("PASS mois Figure 2: original photo denominators and same-photo slopes");
}else if(process.argv.includes("--write")){
 fs.writeFileSync(dest,expected,"utf8");console.log("Wrote",dest);
}else{console.log("Use --check or --write");}
