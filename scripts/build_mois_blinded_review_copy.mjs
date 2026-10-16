#!/usr/bin/env node
// Regenerate a source-derived GEB double-anonymous working review copy.
// Redacts publicly identifying research-repository/commit information.
// IMPORTANT: it does not provide the required anonymous reviewer DATA link.
import fs from "node:fs";
import path from "node:path";
import {fileURLToPath} from "node:url";
const ROOT=path.resolve(path.dirname(fileURLToPath(import.meta.url)),"..");
function anonymize(src){
 const title=src.split("\n")[0];
 if(!title.startsWith("# The geographic and measurement dependence"))throw Error("Source manuscript title unexpected");
 const running="Floral colour and moisture patterns";
 if(running.length>=40)throw Error("Running title too long");
 const a=src.indexOf("## Abstract\n"),b=src.indexOf("## Data and code availability\n"),ref=src.indexOf("## References\n"),tbl=src.indexOf("## Main tables\n");
 if(Math.min(a,b,ref,tbl)<0||!(a<b&&b<ref&&ref<tbl))throw Error("Manuscript editorial boundaries changed");
 const body=src.slice(a,b).trim();
 const refs=src.slice(ref,tbl).trim();
 const tables=src.slice(tbl).trim();
 const available=`## Data and code availability

**BLINDED-REVIEW SOURCE ACCESS — NOT YET PROVIDED.** Archived source photographic hue and environmental inputs, analysis code, original model results, and independent receipt verifiers will be made available in an appropriately licensed and reviewer-anonymized stable archive. The private reviewer access link and the permanent dataset record **have not been established**; the manuscript must not be submitted until these are supplied. GitHub source mirrors alone do not satisfy the journal's stable archiving requirement. Original photo/licence restrictions may require separately justified access conditions.

`;
 let out=title+"\n\n**Running title:** "+running+"\n\n"+body+"\n\n"+refs+"\n\n"+available+tables+"\n";
 // Strip text links that reveal author repositories; keep scholarly DOI addresses.
 out=out.replace(/\[([^\]]+)\]\(https:\/\/github\.com\/[^)\s]+\)/gi,"$1 (source held for blinded review)");
 out=out.replace(/https:\/\/github\.com\/[^\s)]+/gi,"[repository withheld for blinded review]");
 out=out.replace(/`(?:results|scripts|figures|docs)\/[^`]+`/g,"[source file withheld for blinded review]");
 out=out.replace(/`[a-f0-9]{40,64}`/gi,"[source fingerprint withheld for blinded review]");
 out=out.replace(/\b[a-f0-9]{64}\b/gi,"[source fingerprint withheld for blinded review]");
 out=out.replace(/\bFCP\b/g,"the historical source study");
 out=out.replace(/\bmois\b/gi,"the source study");
 out=out.replace(/\bPR #\d+\b/g,"archived source analysis");
 out=out.replace(/`?source_pr\d+[^`\s]*`?/gi,"[source file withheld for blinded review]");
 const lower=out.toLowerCase();
 for(const forbidden of ["zuizui0223","github.com","/fcp/","mois_geb","fcp_macro","run 385"]){if(lower.includes(forbidden))throw Error("Author repository identifier remains "+forbidden);}
 if(!out.includes("## Data and code availability")||!out.includes("## Main tables")||!out.includes("## Figure legends"))throw Error("Blinded structural output damaged");
 if(!out.includes("NOT YET PROVIDED"))throw Error("Missing immutable reviewer data gate");
 if((out.match(/^## Data and code availability/gm)||[]).length!==1)throw Error("Duplicate blinded data statement");
 return out;
}
const source=fs.readFileSync(path.join(ROOT,"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md"),"utf8");
const blind=anonymize(source);
const numbered=blind.split("\n").map((x,i)=>String(i+1).padStart(4,"0")+"  "+x).join("\n");
const md=path.join(ROOT,"docs/MOIS_GEB_BLINDED_REVIEW_COPY_20261014.md");
const txt=path.join(ROOT,"docs/MOIS_GEB_LINE_NUMBERED_REVIEW_COPY_20261014.txt");
if(process.argv.includes("--check")){
 if(!fs.existsSync(md)||fs.readFileSync(md,"utf8")!==blind)throw Error("Blinded manuscript differs from regenerated author-facing source");
 if(!fs.existsSync(txt)||fs.readFileSync(txt,"utf8")!==numbered)throw Error("Line-numbered review copy differs from frozen source");
 console.log("PASS double-anonymous working copy from author-facing manuscript: "+numbered.split("\n").length+" numbered lines");
}else if(process.argv.includes("--write")){
 fs.writeFileSync(md,blind,"utf8");fs.writeFileSync(txt,numbered,"utf8");console.log("Wrote blinded manuscript and line-numbered companion");
}else{console.log("Usage: node scripts/build_mois_blinded_review_copy.mjs --check | --write");}
