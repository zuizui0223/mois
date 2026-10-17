#!/usr/bin/env node
"use strict";
/* Browser-free quality checks: reviewer CSV roundtrip and scientific blinding.
 * This checks the pure browser parser and syntax; actual image access remains
 * conditional on each original iNaturalist URL being accessible to reviewers.
 */
const fs = require("node:fs");
const crypto = require("node:crypto");
const vm = require("node:vm");
const assert = require("node:assert/strict");

const args=process.argv.slice(2);
if(args.length!==3) throw Error("Usage: node test_mois_400_browser_review_ui.js browser.html A.csv B.csv");
const [htmlPath, aPath, bPath]=args;
const html=fs.readFileSync(htmlPath,"utf8");
assert(html.startsWith("<!DOCTYPE html>"));
assert(!html.includes('src="https://') && !html.includes("fetch("));
assert(!html.includes("XMLHttpRequest") && !html.includes("navigator.sendBeacon"));
assert(!html.includes("selected_season"));
assert(!html.includes("observed_on"));
assert(!html.includes("latitude") && !html.includes("longitude"));
const scriptMatch=html.match(/<script>\s*([\s\S]*?)\s*<\/script>/);
assert(scriptMatch && scriptMatch[1]);
new vm.Script(scriptMatch[1]); // compile all viewer JS without a browser DOM
const code=scriptMatch[1];
const start=code.indexOf("function parseCSV(");
const end=code.indexOf("function notify(");
assert(start>=0 && end>start);
const context={};
vm.runInNewContext(code.slice(start,end)+"\nthis.parseCSV=parseCSV;this.encodeCSV=encodeCSV;",context);
const expected=["review_token","species","original_photo_url","photo_license","attribution",
                "reviewer_id","open_flower_stage","roi_organ","notes"];
const frozen={
  A:"576a55e7e613faf8c93652b4ac2d72b79f8a38bec4bfe8267c8d67bbc3fed322",
  B:"bfbbfc92106585ca260fbd5a33573ad2598209d19907a52bbce2b941b1164203",
};
for (const [who,p] of [["A",aPath],["B",bPath]]) {
  const raw=fs.readFileSync(p,"utf8");
  assert.equal(crypto.createHash("sha256").update(raw).digest("hex"),frozen[who],
               "original fixed review form SHA mismatch");
  const parsed=context.parseCSV(raw);
  assert.equal(parsed.length,401);
  assert.equal(JSON.stringify(parsed[0]),JSON.stringify(expected));
  assert(parsed.every(row=>row.length===9));
  assert.equal(new Set(parsed.slice(1).map(row=>row[0])).size,400);
  assert(parsed.slice(1).every(row=>row[5]===""&&row[6]===""&&row[7]===""));
  assert(parsed.slice(1).every(row=>new URL(row[2]).hostname==="inaturalist-open-data.s3.amazonaws.com"));
  const serialized=context.encodeCSV(parsed);
  assert.equal(JSON.stringify(context.parseCSV(serialized)),JSON.stringify(parsed),
               "CSV must preserve original URL/licences/quote fields on download");
}
const complicated=[expected,["t","Test with comma, dash","https://inaturalist-open-data.s3.amazonaws.com/photos/1/large.jpg",
             "cc-by","© Author, CC-BY","A","open_flower","petal_tepal",
             "Reviewed two petals, line one\nLine two \"quoted\""]];
assert.equal(JSON.stringify(context.parseCSV(context.encodeCSV(complicated))),
             JSON.stringify(complicated),"RFC4180 quote/newline roundtrip is broken");
console.log("MOIS_400_BROWSER_ORIGINAL_CSV_BOTH_REVIEWERS_NO_PHENOLOGY_OR_COLOUR_REPLAY_PASS");
