#!/usr/bin/env python3
"""MOIS × DoPI v2: independent positive flower-visiting guild evidence.

Original public DoPI v2, DOI 10.5281/zenodo.22876365 (European literature),
is a 384,896,230-byte, 676,900-row CSV, SHA256 pinned below. Import only
positive documented study×plant×visitor-guild records, with real source
pollination evidence grades:
  1 confirmed pollen transfer or seed set;
  2 evidence of pollen carriage, not proven successful transfer;
  3 flower visit without confirmed pollen carriage/transfer;
  4 robber/herbivory/nonpollinating: EXCLUDE.

Five guilds fixed from independent visitor taxonomy as in verified fd
release: Bombus, other_bee, butterfly_family, Sphingidae_hawkmoth,
Syrphidae_hoverfly. An unrecorded guild is UNKNOWN, never an absence.
Photo nine hues are exact MOIS accepted species nominal ROI labels, not
measured flower-petal colour.

Analysis:
 A. source and exact-binomial support by guild, source publication and hue
    under two predeclared evidence scopes: pollen 1/2, visits 1/2/3.
 B. within ORIGINAL publication, only studies with >=2 matched MOIS
    plants with *positive* evidence of BOTH compared guilds. Compare
    9 hues' proportions between the two positive-plant subsets.
 C. within that same study, conditional on shared positively documented
    plant genera on both sides (>=2 shared genera), report genus-equal
    differences and their publication support, NOT fabricated guild absence.
 D. 5 fixed guild pairs × 9 hues × 2 grade scopes = 90 exploratory tests,
    BH FDR across the entire supported family, not per chosen colour.
    A study with <10 independently documented comparisons => no p-value.
 D is an observational species-pool concordance comparison, not floral
    preference, biological effectiveness, or pollinator-driven evolution.
 All study and plant IDs remain internal and are NOT written to outputs.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter,defaultdict
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy import stats

MOIS_ZIP_SHA="d024e4f497c67c0eb3073244dc84189e50fe1c22deec2231f27adef8ab5d7035"
MOIS_CSV="MOIS_37104_accepted_species_photo9_and_current_environment.csv.gz"
MOIS_CSV_SHA="af74536f54317f3fc923b8f5f00eb61f5cba3803a19b9be021b073a7a1ea2d46"
DOPI_SHA="b95d3084399ed2c0978af2485801108a27d14e843b9fed2611934ea582b745c8"
DOPI_BYTES=384896230
DOPI_ROWS=676900
DOPI_DOI="10.5281/zenodo.22876365"
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
GUILDS=("Bombus","other_bee","butterfly_family",
        "Sphingidae_hawkmoth","Syrphidae_hoverfly")
PAIRS=(
    ("Sphingidae_hawkmoth","Bombus"),
    ("butterfly_family","Bombus"),
    ("Sphingidae_hawkmoth","butterfly_family"),
    ("Bombus","Syrphidae_hoverfly"),
    ("other_bee","Bombus"),
)
SCOPES=("positive_pollen_evidence_grade12",
        "positive_visit_or_pollen_grade123")
BEE_FAMILIES={"apidae","andrenidae","halictidae","colletidae",
              "megachilidae","melittidae","stenotritidae"}
BUTTERFLY_FAMILIES={"nymphalidae","papilionidae","pieridae",
                    "lycaenidae","hesperiidae","riodinidae"}
SEED=20261118
N_BOOT=3000
MIN_PLANTS_PER_GUILD_STUDY=2
MIN_STUDIES=10
MIN_SHARED_GENERA=2

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_DOPI_POSITIVE_VISITORS: "+str(msg))

def file_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def byte_sha(blob):
    return hashlib.sha256(blob).hexdigest()

def get_mois_source(zippath):
    require(file_sha(zippath)==MOIS_ZIP_SHA,
            "Original 37104 accepted species source adapter ZIP SHA changed")
    with ZipFile(zippath) as z:
        source=z.read(MOIS_CSV)
        receipt=json.loads(z.read("source_receipt.json"))
    require(byte_sha(source)==MOIS_CSV_SHA==receipt["source_CSV_SHA256"]
            and receipt["global_accepted_species"]==37104,
            "MOIS accepted species source CSV digest changed")
    data=pd.read_csv(io.BytesIO(source),compression="gzip",
                     usecols=["key","species","genus","source_hue9"],
                     dtype=str,keep_default_na=False)
    require(len(data)==37104 and data.key.is_unique
            and data.species.is_unique and set(data.source_hue9)==set(HUES),
            "Accepted species 9-hue/unique exact name denominator drifted")
    return {
        species:(hue,genus)
        for species,hue,genus in zip(
            data.species,data.source_hue9,data.genus)
    },receipt

def guild(order,family,species):
    o=(order or "").strip().lower()
    f=(family or "").strip().lower()
    s=(species or "").strip().lower()
    if s.startswith("bombus ") and f in BEE_FAMILIES:return "Bombus"
    if f in BEE_FAMILIES:return "other_bee"
    if f=="sphingidae":return "Sphingidae_hawkmoth"
    if f in BUTTERFLY_FAMILIES:return "butterfly_family"
    if f=="syrphidae":return "Syrphidae_hoverfly"
    return None

def normalize_study(rec):
    doi=(rec.get("doi") or "").strip().lower()
    if doi and doi not in {"n/a","na","unknown","none"}:return doi
    pub=(rec.get("publication") or "").strip().lower()
    if pub and pub not in {"n/a","na","unknown","none"}:return pub
    return None

def scan_real_source(path,mois_species):
    require(Path(path).stat().st_size==DOPI_BYTES,
            "DoPI original source byte length mismatch")
    require(file_sha(path)==DOPI_SHA,
            "DoPI v2 original source SHA mismatch")
    positive=defaultdict(set)   # (scope, publication, observed guild)->set(accepted binomials)
    count=Counter()
    all_recorded_guild_plants=defaultdict(set)
    all_positive_studies=defaultdict(set)
    with Path(path).open("rt",encoding="utf-8-sig",errors="replace",
                         newline="") as file:
        reader=csv.DictReader(file)
        cols=set(reader.fieldnames or [])
        required={"plant_species","pollinator_species","pollinator_family",
                  "pollinator_order","doi","publication","pollination"}
        require(required<=cols,
                "DoPI v2 original source taxonomic/evidence column drift")
        for row in reader:
            count["raw_source_rows"]+=1
            plant=(row.get("plant_species") or "").strip()
            if plant not in mois_species:continue
            count["source_rows_with_exact_MOIS_accepted_plant"]+=1
            study=normalize_study(row)
            if study is None:
                count["mois_matched_rows_without_publication"]+=1
                continue
            k=guild(row.get("pollinator_order"),row.get("pollinator_family"),
                    row.get("pollinator_species"))
            if k is None:continue
            count["MOIS_matched_rows_with_one_of_five_guilds"]+=1
            grade=(row.get("pollination") or "").strip()
            if grade not in ("1","2","3"):
                count["excluded_nonpollinating_or_unknown_evidence_rows"]+=1
                continue
            count["MOIS_matched_valid_positive_rows"]+=1
            count["positive_grade_"+grade]+=1
            all_recorded_guild_plants[k].add(plant)
            all_positive_studies[k].add(study)
            positive[(SCOPES[1],study,k)].add(plant)
            if grade in ("1","2"):
                positive[(SCOPES[0],study,k)].add(plant)
    require(count["raw_source_rows"]==DOPI_ROWS,
            "Original DoPI 676900 source rows different")
    require(count["source_rows_with_exact_MOIS_accepted_plant"]>100,
            "Real DoPI plants too scarce after exact MOIS species match")
    require(count["MOIS_matched_valid_positive_rows"]>50,
            "No source-backed positive guild-pollination records available")
    return positive,count,all_recorded_guild_plants,all_positive_studies

def data_support(positive,lookup,count,all_plants,all_studies):
    aggregate=[]
    for scope in SCOPES:
        for name in GUILDS:
            units=[(s,p) for (sp,s,g),p in positive.items()
                   if sp==scope and g==name]
            spp=set().union(*(sp for study,sp in units)) if units else set()
            genres={lookup[species][1] for species in spp}
            by_hue=Counter(lookup[species][0] for species in spp)
            for hue in HUES:
                aggregate.append({
                    "evidence_scope":scope,"positive_visitor_guild":name,
                    "photo9_hue":hue,
                    "n_distinct_exact_MOIS_species_in_positive_records":len(spp),
                    "n_original_published_studies_positive":len(units),
                    "n_distinct_positive_species_genera":len(genres),
                    "n_species_photo9_hue":int(by_hue[hue]),
                    "proportion_of_positive_guild_species_with_photo_hue":(
                         by_hue[hue]/len(spp) if spp else None),
                })
    return aggregate

def paired_study_hue(positive,lookup):
    rng=np.random.default_rng(SEED)
    studies_by_scope={
        scope:sorted({study for sc,study,g in positive if sc==scope})
        for scope in SCOPES
    }
    results=[]
    for scope in SCOPES:
        studies=studies_by_scope[scope]
        for A,B in PAIRS:
            study_inputs=[]
            for study in studies:
                x=positive.get((scope,study,A),set())
                y=positive.get((scope,study,B),set())
                if min(len(x),len(y))<MIN_PLANTS_PER_GUILD_STUDY:
                    continue
                shared_genus=(
                    {lookup[s][1] for s in x}
                    &{lookup[s][1] for s in y})
                study_inputs.append((study,x,y,shared_genus))
            for hue in HUES:
                # Study is only one statistical unit, however many times
                # the same species or visitor was observed in that study.
                diff=[]
                both_genus=[]
                positive_a=positive_b=0
                for study,x,y,shared in study_inputs:
                    xa=sum(lookup[s][0]==hue for s in x)
                    ya=sum(lookup[s][0]==hue for s in y)
                    positive_a+=int(xa>0)
                    positive_b+=int(ya>0)
                    diff.append(xa/len(x)-ya/len(y))
                    if len(shared)>=MIN_SHARED_GENERA:
                        pergenus=[]
                        for genus in shared:
                            xx=[s for s in x if lookup[s][1]==genus]
                            yy=[s for s in y if lookup[s][1]==genus]
                            pergenus.append(
                                sum(lookup[s][0]==hue for s in xx)/len(xx)
                               -sum(lookup[s][0]==hue for s in yy)/len(yy))
                        both_genus.append(float(np.mean(pergenus)))
                n=len(diff)
                if n>=MIN_STUDIES:
                    v=np.asarray(diff,dtype=float)
                    orig=float(v.mean())
                    simulations=v[rng.integers(n,size=(N_BOOT,n))].mean(axis=1)
                    ci=np.quantile(simulations,[.025,.975])
                    # Exploratory study-cluster t p-value, not causal guild
                    # selection; BH across 90 planned hue-pair-grade panels.
                    if v.std(ddof=1)>0:
                        p=float(stats.ttest_1samp(v,0).pvalue)
                    else:
                        p=1.0 if abs(orig)<1e-12 else 0.0
                    supported=True
                else:
                    orig=ci=p=None
                    supported=False
                entry={
                    "evidence_scope":scope,
                    "guild_A_positive_plants":A,
                    "guild_B_positive_plants":B,
                    "photo9_hue":hue,
                    "n_same_published_studies_both_positive_guilds":n,
                    "n_studies_with_hue_positive_in_guild_A":positive_a,
                    "n_studies_with_hue_positive_in_guild_B":positive_b,
                    "n_same_study_shared_genus_panels":len(both_genus),
                    "positive_guildA_minus_guildB_photo_hue_fraction":orig,
                    "bootstrap95_by_publication_low":(
                        float(ci[0]) if supported else None),
                    "bootstrap95_by_publication_high":(
                        float(ci[1]) if supported else None),
                    "source_publication_t_p":p,
                    "eligible_for_exploratory_test":supported
                }
                if len(both_genus)>=MIN_STUDIES:
                    a=np.asarray(both_genus,float)
                    if a.std(ddof=1)>0:
                        g_p=float(stats.ttest_1samp(a,0).pvalue)
                    else:g_p=1.0 if abs(a.mean())<1e-12 else 0.0
                    entry["within_publication_shared_genus_equal_effect"]=float(a.mean())
                    entry["within_publication_shared_genus_t_p"]=g_p
                else:
                    entry["within_publication_shared_genus_equal_effect"]=None
                    entry["within_publication_shared_genus_t_p"]=None
                results.append(entry)
    require(len(results)==len(SCOPES)*len(PAIRS)*len(HUES),
            "Original predeclared photo9 guild test family count changed")
    eligible=[j for j,r in enumerate(results) if r["eligible_for_exploratory_test"]]
    # BH across ALL 90 hypotheses; unsupported panels receive no p/q.
    # Declaring unsupported hypotheses untested cannot boost a lone nominal
    # p-value to an ex-ante confirmed interaction effect.
    if eligible:
        from scipy.stats import false_discovery_control
        q=false_discovery_control(
             np.asarray([results[j]["source_publication_t_p"] for j in eligible]),
             method="bh")
        for j,v in zip(eligible,q):
            results[j]["BH_q_among_supported_90_screen"]=float(v)
    eligible_study_pairs=defaultdict(set)
    for (scope,study,guild_name),plants in positive.items():
        if len(plants)>=MIN_PLANTS_PER_GUILD_STUDY:
            eligible_study_pairs[(scope,study)].add(guild_name)
    return results,{
      "planned_tests":len(results),
      "tested_supported_panels":len(eligible),
      "unsupported_panels":len(results)-len(eligible),
      "n_publications_with_at_least_two_supported_guilds_by_evidence_scope":{
          scope:int(sum(
              len(guilds)>=2 for (sc,study),guilds in eligible_study_pairs.items()
              if sc==scope))
          for scope in SCOPES
      },
      "exploratory_FDR_scope":"90 predeclared photo-colour comparisons, tested support subset only",
    }

def run(mois_zip,dopi_csv,out):
    root=Path(out)
    require(not root.exists(),"Never overwrite original DoPI evidence results")
    lookup,source_receipt=get_mois_source(mois_zip)
    positive,count,all_plants,all_studies=scan_real_source(dopi_csv,lookup)
    descriptive=data_support(positive,lookup,count,all_plants,all_studies)
    comparison,coverage=paired_study_hue(positive,lookup)
    root.mkdir(parents=True)
    source_table=root/"positive_guild_exact_MOIS_photo9_species_hue_support.csv"
    pd.DataFrame(descriptive).to_csv(source_table,index=False)
    pair_table=root/"within_publication_positive_guild_pair_nine_hue_tests.csv"
    pd.DataFrame(comparison).to_csv(pair_table,index=False)
    guide={
      "source":"DoPI v2 European positively documented visitor records × MOIS 37104 exact accepted photo9 species",
      "original_DoPI_DOI":DOPI_DOI,
      "original_DoPI_rows":DOPI_ROWS,"original_DoPI_csv_size_bytes":DOPI_BYTES,
      "source_DoPI_SHA256":DOPI_SHA,
      "source_MOIS_signed_adapter_ZIP_SHA256":MOIS_ZIP_SHA,
      "source_MOIS_species_CSV_SHA256":MOIS_CSV_SHA,
      "original_source_accepted_photo9_species":len(lookup),
      "all_source_row_screen_counts":dict(count),
      "coverage_by_guild":{
          name:{
             "n_distinct_MOIS_plants_positive":len(all_plants[name]),
             "n_distinct_original_published_studies_positive":len(all_studies[name]),
          } for name in GUILDS
      },
      "planned_paired_test_scopes":list(SCOPES),
      "planned_guild_comparisons":[list(x) for x in PAIRS],
      "planned_nine_photo_hues":list(HUES),
      "paired_comparison_support":coverage,
      "only_positive_visitor_or_pollen_records_not_absences":True,
      "grade3_visitation_not_proven_pollen_transfer":True,
      "grade4_nonpollinating_excluded":True,
      "attribution_is_DoPI_original_publications_and_zenodo_release":True,
      "same_publication_visitor_lists_not_ecological_preference_or_fitness":True,
      "Europe_source_not_global_local_pollinator_community":True,
      "photo_colour_not_verified_actual_petal_colour":True,
      "species_names_and_original_publication_IDs_not_exported":True,
    }
    report=root/"source_and_inferential_limits.json"
    report.write_text(json.dumps(guide,indent=2,allow_nan=False)+"\n")
    receipt={
      **guide,
      "output_SHA256":{p.name:file_sha(p) for p in
                       (source_table,pair_table,report)}
    }
    (root/"receipt.json").write_text(
         json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    group_summary=[]
    for scope in SCOPES:
        for name in GUILDS:
            subset=[r for r in descriptive if
                    r["evidence_scope"]==scope and
                    r["positive_visitor_guild"]==name]
            group_summary.append({
                "scope":scope,"guild":name,
                "n_positive_species":subset[0]["n_distinct_exact_MOIS_species_in_positive_records"],
                "n_studies":subset[0]["n_original_published_studies_positive"],
                "white_fraction":next(r["proportion_of_positive_guild_species_with_photo_hue"]
                     for r in subset if r["photo9_hue"]=="white")
            })
    print("MOIS_DOPI_ORIGINAL_676900_NAMED_GUILDS_POSITIVE_EVIDENCE_SUCCESS",
          json.dumps({
             "source_rows":count["raw_source_rows"],
             "source_rows_exact_MOIS":count["source_rows_with_exact_MOIS_accepted_plant"],
             "source_matched_positive":count["MOIS_matched_valid_positive_rows"],
             "positive_guilds":group_summary,
             "comparison_coverage":coverage,
             "comparison_supported":[{k:r[k] for k in (
                  "evidence_scope","guild_A_positive_plants",
                  "guild_B_positive_plants","photo9_hue",
                  "n_same_published_studies_both_positive_guilds",
                  "positive_guildA_minus_guildB_photo_hue_fraction",
                  "source_publication_t_p","BH_q_among_supported_90_screen"
             )} for r in comparison if
                r["eligible_for_exploratory_test"] and
                r["photo9_hue"] in ("white","red","purple")],
          },allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--mois",required=True)
    p.add_argument("--dopi",required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    run(args.mois,args.dopi,args.out)
