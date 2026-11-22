#!/usr/bin/env python3
"""Seebens v3.1 first observations of aliens x MOIS photo-hue cohorts.

This does NOT measure introduction dates, within-species flower colour change,
selection strength, naturalization rates, or adaptive radiation. Curated first
records are used only for research under the source citation conditions.
"""
import argparse, hashlib, io, json, re
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import pycountry

SOURCE_SHA="7897fbd88d300f48d63039767c5d99b4ddc4799cc88604be861456b5ef6ecbd4"
GBIF_SHA="422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
ALIASES={"United States":"US","United States of America":"US",
         "United Kingdom":"GB","Great Britain":"GB","Russia":"RU",
         "South Korea":"KR","North Korea":"KP","Czech Republic":"CZ",
         "Iran":"IR","Vietnam":"VN","Tanzania":"TZ","Bolivia":"BO",
         "Venezuela":"VE","Laos":"LA","Taiwan":"TW","Kosovo":None}
PERIOD_EDGES=[1500,1800,1900,1950,2000,2024]
PERIOD_NAMES=["1500–1799","1800–1899","1900–1949","1950–1999","2000–2023"]

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def guard(ok,message):
    if not bool(ok):
        raise ValueError("FIRST_RECORD_SOURCE_GATE: "+message)

def countries(names):
    m={}
    for c in pycountry.countries:
        for name in (c.name,getattr(c,"official_name",None),getattr(c,"common_name",None)):
            if name:
                m[name]=c.alpha_2
    m.update(ALIASES)
    return names.map(m)

def sources(gbif,first_records):
    guard(digest(gbif)==GBIF_SHA,"GBIF original ZIP SHA changed")
    guard(digest(first_records)==SOURCE_SHA,"Seebens original XLSX SHA changed")
    with ZipFile(gbif) as z:
        r=json.loads(z.read("receipt.json"))
        guard(r["all_16_shards_source_verified"] and
              r["GBIF_occurrences_passed_all_filters"]==396142825,"GBIF original lineage changed")
        s=pd.read_csv(io.BytesIO(z.read("all_37104_accepted_species_grid_coverage_and_photo9.csv.gz")),
                      compression="gzip",dtype={"key":str,"species":str})
    guard(len(s)==37104 and s.key.is_unique and s.species.is_unique
          and set(s.source_hue9)==set(HUES),"MOIS accepted taxon/hue source drift")
    f=pd.read_excel(first_records,sheet_name="FirstRecords",dtype=str,engine="openpyxl")
    needed={"Region","LifeForm","TaxonName","FirstRecord","FirstRecord_orig",
            "PresentStatus","Phylum","Source","DataUsage"}
    guard(len(f)==77375 and len(f.columns)==20 and needed.issubset(f.columns),
          "Seebens 77,375 original records schema changed")
    return s,f

def link(s,f):
    f=f.copy()
    f["LifeForm"]=f.LifeForm.fillna("").str.strip()
    f["TaxonName"]=f.TaxonName.fillna("").str.strip()
    f["Region"]=f.Region.fillna("").str.strip()
    f["iso2"]=countries(f.Region)
    f["year"]=pd.to_numeric(f.FirstRecord,errors="coerce")
    f["valid_year"]=f.year.between(1500,2023)&f.year.eq(np.floor(f.year))
    orig=f.FirstRecord_orig.fillna("").str.strip().str.replace(r"\.0$","",regex=True)
    f["original_exact"]=orig.str.fullmatch(r"(1[5-9]\d{2}|20[0-2]\d)") & f.year.eq(pd.to_numeric(orig,errors="coerce"))
    plant=f[f.LifeForm.str.contains("plant",case=False,regex=False)].copy()
    s=s.copy();s["genus"]=s.species.str.split().str[0]
    x=plant.merge(s[["key","species","genus","source_hue9","hue_status","occupied_grid_cells"]],
                  left_on="TaxonName",right_on="species",validate="many_to_one")
    guard(len(x)>100,"Insufficient exact plant species matches")
    duplicates=int(x.duplicated(["TaxonName","Region"]).sum())
    x=x.sort_values(["TaxonName","Region","year","original_exact"],
                    ascending=[True,True,True,False],kind="stable")
    x=x.drop_duplicates(["TaxonName","Region"],keep="first").reset_index(drop=True)
    report={"all_source_rows":len(f),
            "plant_rows":len(plant),
            "source_lifeforms":{str(k):int(v) for k,v in f.LifeForm.value_counts().items()},
            "plant_unique_taxa":int(plant.TaxonName.nunique()),
            "MOIS_matched_unique_taxa":int(x.key.nunique()),
            "MOIS_matched_taxon_region":len(x),
            "MOIS_matched_duplicate_taxon_region_records":duplicates,
            "MOIS_matched_country":int(x.iso2.notna().sum()),
            "MOIS_matched_country_unique_species":int(x.loc[x.iso2.notna(),"key"].nunique()),
            "MOIS_matched_country_valid_year":int((x.iso2.notna()&x.valid_year).sum()),
            "MOIS_matched_country_original_exact_year":int((x.iso2.notna()&x.valid_year&x.original_exact).sum()),
            "MOIS_matched_country_imputed_or_approx_year":int((x.iso2.notna()&x.valid_year&~x.original_exact).sum()),
            "MOIS_country_firstrecord_year_range":[int(x.loc[x.iso2.notna()&x.valid_year,"year"].min()),
                                                  int(x.loc[x.iso2.notna()&x.valid_year,"year"].max())],
            "MOIS_matched_lifeforms":{str(k):int(v) for k,v in x.LifeForm.value_counts().items()},
            "MOIS_matched_PresentStatus":{str(k):int(v) for k,v in x.PresentStatus.fillna("(missing)").value_counts().items()},
            "subnational_or_unresolved_region_top20":[{"region":str(k),"n":int(v)} for k,v in
                                                   x.loc[x.iso2.isna(),"Region"].value_counts().head(20).items()]}
    return s,x,report

def dataset(x, quality="all"):
    y=x.loc[x.iso2.notna()&x.valid_year&x.original_exact].copy()
    if quality=="no_bronze":
        y=y.loc[y.source_hue9!="bronze"]
    elif quality=="consensus_only":
        y=y.loc[y.hue_status=="sampling_consensus_single"]
    y["site_genus"]=y.iso2+"|"+y.genus
    group=y.groupby("site_genus")
    eligible=group.agg(species=("key","nunique"),years=("year","nunique"))
    fixed=eligible.index[(eligible.species>=2)&(eligible.years>=2)]
    y=y[y.site_genus.isin(fixed)].copy().reset_index(drop=True)
    return y, int(len(fixed))

def permutations(y,s,draws,model,seed):
    if len(y)<100 or y.site_genus.nunique()<20:
        return {"status":"insufficient >=100 records or >=20 country-genus groups",
                "n_records":len(y),"n_country_genus_groups":int(y.site_genus.nunique())}, None
    ix=pd.Index(s.key.astype(str)).get_indexer(y.key.astype(str))
    guard((ix>=0).all(),"MOIS species identity lost in matching")
    donor=np.unique(ix)
    labels=s.source_hue9.to_numpy(str)
    gcode=pd.factorize(s.genus,sort=True)[0]
    bins=np.digitize(s.occupied_grid_cells.to_numpy(int),[2,6,21,101])
    strata={}
    for i in donor:
        code=(int(gcode[i]),int(bins[i])) if model=="genus_breadth" else (int(gcode[i]),)
        strata.setdefault(code,[]).append(int(i))
    exchange=[np.array(ids,dtype=int) for ids in strata.values()
              if len(ids)>1 and len(set(labels[ids]))>1]
    groups, _ = pd.factorize(y.site_genus,sort=True)
    gs=np.bincount(groups)
    years=y.year.to_numpy(float)
    yc=(years-np.bincount(groups,weights=years)[groups]/gs[groups])/100
    wx=yc/gs[groups]
    denom=np.sum(wx*yc)
    guard(denom>0,"No within-group temporal predictor variation")
    white=lambda hue: float(np.dot(wx,hue[ix]=="white")/denom)
    observed=white(labels)
    rng=np.random.default_rng(seed)
    simulated=np.empty(draws)
    for k in range(draws):
        hue=labels.copy()
        for members in exchange:
            hue[members]=labels[members][rng.permutation(len(members))]
        simulated[k]=white(hue)
    lo=(1+int(np.sum(simulated<=observed+1e-12)))/(draws+1)
    hi=(1+int(np.sum(simulated>=observed-1e-12)))/(draws+1)
    movable=np.unique(np.concatenate(exchange)) if exchange else np.array([],dtype=int)
    return {"status":"completed","hue":"white",
            "interpretation":"later first observed alien species have a higher white-photo-label proportion if slope>0",
            "outcome":"within-country-genus slope in white source photo-hue label / century of first OBSERVATION",
            "source_record_unit":"species x country after unique earliest curated first record",
            "within_country_genus_groups":int(y.site_genus.nunique()),
            "source_plant_country_rows":len(y),"source_unique_species":int(y.key.nunique()),
            "source_unique_genera":int(y.genus.nunique()),"source_countries":int(y.iso2.nunique()),
            "donor_species":len(donor),"variable_donor_groups":len(exchange),
            "movable_species":len(movable),"movable_share":len(movable)/len(donor),
            "dated_country_rows_with_movable_species":int(np.isin(ix,movable).sum()),
            "observed_slope_per_century":observed,"null_mean":float(simulated.mean()),
            "observed_minus_null_slope":float(observed-simulated.mean()),
            "null_sd":float(simulated.std(ddof=1)),
            "null_reference_95pct_interval":np.percentile(simulated,[2.5,97.5]).tolist(),
            "rank_p_lower":lo,"rank_p_upper":hi,"rank_p_two_sided":min(1.,2*min(lo,hi)),
            "null_shared_species_hues_across_all_countries":True,
            "no_effect_confidence_interval_estimated":True},simulated

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--seebens",required=True)
    p.add_argument("--gbif",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--permutations",type=int,default=999)
    a=p.parse_args()
    out=Path(a.out);guard(not out.exists(),"Output may not overwrite another analysis")
    s,f=sources(a.gbif,a.seebens)
    s,x,audit=link(s,f)
    out.mkdir(parents=True)
    year=x.loc[x.iso2.notna()&x.valid_year].copy()
    year["period"]=pd.cut(year.year,PERIOD_EDGES,labels=PERIOD_NAMES,right=False)
    years=year.groupby(["period","source_hue9","original_exact"],observed=True).agg(
        n_species_country=("key","size"),n_species=("key","nunique"),
        n_countries=("iso2","nunique")).reset_index()
    years.to_csv(out/"historical_hue_cohort_counts.csv",index=False)
    stats={};frames=[]
    for sample in ("all","no_bronze","consensus_only"):
        y,_=dataset(x,sample)
        print("MOIS_TEMPORAL_SUBSET",sample, len(y),y.key.nunique(),y.site_genus.nunique(),flush=True)
        for model in ("genus_breadth","genus_only"):
            r,draw=permutations(y,s,a.permutations,model,20261122+(0 if model=="genus_breadth" else 137))
            stats[sample+"_"+model]=r
            if draw is not None:
                frames.append(pd.DataFrame({"sample":sample,"model":model,
                                            "draw_id":np.arange(1,len(draw)+1),
                                            "slope":draw}))
            print("MOIS_FIRSTRECORD_MODEL",sample,model,json.dumps(r,allow_nan=False),flush=True)
    if frames:
        pd.concat(frames,ignore_index=True).to_csv(
            out/"synchronous_species_hue_permutation_draws.csv.gz",index=False,
            compression={"method":"gzip","mtime":0})
    summary={"schema":"mois_seebens_v31_firstrecord_exact_country_genus_photo_hue_20261122",
             "source_Seebens_v3_1_sha256":SOURCE_SHA,
             "source_original_MOIS_GBIF_sha256":GBIF_SHA,
             "source_record_url":"https://zenodo.org/records/10039630",
             "source_literature":"Seebens et al. 2017 Nature Communications 8 14435; Seebens et al. 2018 PNAS 115 E2264-E2273",
             "source_use":"non-commercial research with attribution",
             "source_audit":audit,"source_firstrecord_not_introduction":True,
             "source_uncertain_original_dates_removed_from_inference":True,
             "ISO_countries_exact_name_or_explicit_alias_only":True,
             "photographic_hues_not_independently_botanically_confirmed":True,
             "no_time_varying_pheno_change_tested":True,
             "no_causal_post_introduction_selection_or_adaptive_radiation_tested":True,
             "main_primary_sample":"all_genus_breadth","main_primary_outcome":"white label x first-record calendar year within country-genus",
             "sensitivities":["all_genus_only","no_bronze_genus_breadth",
                              "consensus_only_genus_breadth"],
             "temporal_results":stats,"n_permutations":a.permutations}
    (out/"temporal_linkage_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    members={p.name:digest(p) for p in sorted(out.iterdir())}
    (out/"receipt.json").write_text(json.dumps(
        {"status":"complete","source_sha256":{"seebens":SOURCE_SHA,"MOIS":GBIF_SHA},
         "script_SHA256":digest(__file__),"outputs_sha256":members},indent=2)+"\n")
    print("MOIS_FIRSTRECORD_LINK_COMPLETE",json.dumps({"source_audit":audit,
          "temporal_primary":stats["all_genus_breadth"]},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
