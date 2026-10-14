#!/usr/bin/env node
// Deterministic source-only transfer Figure 3. Original photo/abiotic fits are not rerun.
import fs from "node:fs";
import path from "node:path";
import {fileURLToPath} from "node:url";
const ROOT=path.resolve(path.dirname(fileURLToPath(import.meta.url)),"..");
function extract(a,b,c) {
 const pa=a.primary36,rb=b.main36_same_20786_evaluated_species,rc=c.primary_result.genus_equal.total.contrasts;
 const p=(name,d)=>({name,point:d.genus_equal_gain,ci:d.CI95});
 const q=(name,d)=>({name,point:d.estimate,ci:d.CI95});
 const result=[
 {title:"A | Change the prediction target",subtitle:"Different biological/photographic target cohorts",data:[
  p("Absolute species hue (20,786 spp.)",pa.absolute_species_colour),
  p("Within-genus deviation (16,801 spp.)",pa.within_unseen_genus_species_differences),
  p("Within-species photos (5,939 spp.)",pa.within_species_photo_differences)
 ]},
 {title:"B | Exclude other geographic support",subtitle:"Same 20,786 original target species",data:[
  p("Unseen genus only",rb.genus_only),
  p("Unseen genus + sector",rb.genus_and_region),
  p("Plus all source photos within 500 km",rb.genus_region_and_500km)
 ]},
 {title:"C | Model-form sensitivity at 500 km",subtitle:"Same 20,786 source species, genus-equal loss",data:[
  q("Ridge: adding all36",rc.ridge_environment_gain),
  q("Nonlinear: adding all36",rc.nonlinear_environment_gain),
  q("Nonlinear full over ridge full",rc.full_model_improvement)
 ]}
 ];
 if(result[0].data[0].point<0.004||Math.abs(result[1].data[2].point-0.00017405463489672535)>1e-10||Math.abs(result[2].data[1].point-0.00488179501096242)>1e-10)throw Error("Source transport drift");
 return result;
}
function render(panels) {
 const v=['<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="820" viewBox="0 0 1500 820" role="img" aria-labelledby="title description">','<title id="title">How environmental prediction transfers across taxa and geography</title>','<desc id="description">Original source genus-equal squared-error improvement from adding environment. Grouped results differ by biological target, geography removal, and model form; no direct causal fitness or universal environmental response is inferred.</desc>'];
 const rect=(x,y,w,h,fill,r=0)=>v.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${r}" fill="${fill}"/>`);
 const t=(x,y,text,n=16,b=400,c="#223747",anchor="start")=>v.push(`<text x="${x}" y="${y}" font-family="Arial,sans-serif" font-size="${n}" font-weight="${b}" text-anchor="${anchor}" fill="${c}">${text}</text>`);
 const l=(x,y,X,Y,c="#D4DEE5",w=1)=>v.push(`<line x1="${x}" y1="${y}" x2="${X}" y2="${Y}" stroke="${c}" stroke-width="${w}"/>`);
 rect(0,0,1500,820,"#FFFFFF");rect(0,0,1500,15,"#305872");
 t(47,63,"WHERE DOES ENVIRONMENTAL HUE INFORMATION TRANSFER?",29,700);
 t(47,100,"Original source genus-equal loss reduction • positive = improved prediction",18,400,"#54687A");
 panels.forEach((q,k)=>{
  const x=45+490*k,y=140;
  rect(x,y,463,535,"#F8FAFC",12);t(x+21,183,q.title,20,700);
  t(x+21,215,q.subtitle,14,400,"#5C6E7F");
  const xmin=-0.002,max=0.008,xa=x+62,xb=x+432,m=n=>xa+(n-xmin)*(xb-xa)/(max-xmin);
  [-0.002,0,0.002,0.004,0.006,0.008].forEach(z=>{const xx=m(z);l(xx,280,xx,601,z===0?"#718190":"#D8E1E7",z===0?2:1);t(xx,626,z.toFixed(3),12,400,"#617387","middle");});
  q.data.forEach((d,j)=>{
   const yy=330+j*114,col=["#2A7599","#9C6F9A","#36566F"][j];
   t(x+21,yy-35,d.name,15,600);
   l(m(d.ci[0]),yy,m(d.ci[1]),yy,col,5);
   l(m(d.ci[0]),yy-9,m(d.ci[0]),yy+9,col,2);l(m(d.ci[1]),yy-9,m(d.ci[1]),yy+9,col,2);
   v.push(`<circle cx="${m(d.point)}" cy="${yy}" r="7" fill="${col}"/>`);
   t(x+22,yy+30,(d.point>=0?"+":"")+d.point.toFixed(6)+" ["+d.ci[0].toFixed(5)+", "+d.ci[1].toFixed(5)+"]",13,400,"#586D7F");
  });
 });
 l(46,709,1453,709);
 t(48,740,"Panel A uses different source target cohorts; Panel B fixes the 20,786 evaluation species while training support changes.",16,400,"#596B79");
 t(48,772,"Panel C compares distinct geography baselines: nonlinear environmental gain does not imply full nonlinear model superiority.",16,400,"#596B79");
 v.push("</svg>");return v.join("\n")+"\n";
}
const paths=[
 "results/source_receipts/unseen_genus_colour_environment.json",
 "results/source_receipts/unseen_genus_other_geography_500km.json",
 "results/source_receipts/nonlinear_geographic_transport.json"
];
const data=paths.map(q=>JSON.parse(fs.readFileSync(path.join(ROOT,q),"utf8")));
const fig=render(extract(...data));
const dest=path.join(ROOT,"figures/mois_environmental_transfer_and_model_form.svg");
if(process.argv.includes("--check")){
 if(!fs.existsSync(dest)||fs.readFileSync(dest,"utf8")!==fig)throw Error("Figure 3 not reproduced from source-verified results");
 console.log("PASS Figure 3 source-locked transfer comparisons");
}else if(process.argv.includes("--write")){fs.writeFileSync(dest,fig,"utf8");console.log("wrote",dest);}
else console.log("Use --check or --write");
