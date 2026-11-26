#!/usr/bin/env python3
"""FirstRecords v4.0 × MOIS: observed inter-regional first-record accumulation.

The outcome is a chronology of documentary first records in distinct source
locationIDs. It is NOT measured biological dispersal, introduction time,
post-introduction flower-colour evolution, selection, or speciation.
Original data: FirstRecords 2026 v4.0, CC BY 4.0, Zenodo 18759840.
"""
from __future__ import annotations
import argparse,hashlib,io,json,sys
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

SRC_MD5={
 "FirstRecords_dataset_public_v4.0.csv":"cbba8f2c2c518519674fd4f0540809d4",
 "FirstRecords_taxonomy_table_public_v4.0.csv":"d3e21741364b33fabb3cef3b249543bd",
 "FirstRecords_location_table_public_v4.0.csv":"cb69dd1ad613cbfe24db316b89f5884a",
}
GBIF_SHA="422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
SEED=20261122
HUES={"white","yellow","orange","red","pink","magenta","purple","blue","bronze"}
MAIN_YEAR_MIN=1800
LAST_COMPLETE_YEAR=2023
YEAR_ENDPOINTS={25:1998,50:1973,100:1923}

def require(ok,msg):
    if not bool(ok): raise ValueError("MOIS_FIRSTRECORDS_V4_GATE: "+msg)

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def source_tables(path,gbif_zip):
    root=Path(path)
    files={}
    for name,md5 in SRC_MD5.items():
        p=root/name
        require(p.is_file(),"Missing original FirstRecords v4 file "+name)
        require(hashlib.md5(p.read_bytes()).hexdigest()==md5,"Original Zenodo source checksum drift "+name)
        files[name]=p
    r=pd.read_csv(files["FirstRecords_dataset_public_v4.0.csv"],
                  sep=";",encoding="utf-8-sig",encoding_errors="replace",
                  dtype=str,keep_default_na=False,low_memory=False)
    tax=pd.read_csv(files["FirstRecords_taxonomy_table_public_v4.0.csv"],
                    encoding="utf-8-sig",encoding_errors="replace",
                    dtype=str,keep_default_na=False,low_memory=False)
    loc=pd.read_csv(files["FirstRecords_location_table_public_v4.0.csv"],
                    encoding="utf-8-sig",encoding_errors="replace",
                    dtype=str,keep_default_na=False,low_memory=False)
    require((len(r),len(tax),len(loc))==(90023,56971,557),
            "FirstRecords v4.0 public row counts differ from authenticated source")
    need={"locationID","taxonID","confidenceFirstRecordEvent","firstRecordEvent",
          "verbatimFirstRecordEvent","occurrenceStatus","establishmentMeans",
          "degreeOfEstablishment","accessRights"}
    require(need.issubset(r.columns),"Original v4 first record fields changed")
    require({"taxonID","species","genus","GBIFstatus","GBIFmatchType",
             "GBIFtaxonRank","taxaGroup"}.issubset(tax.columns),
            "Original taxonomy fields changed")
    require({"locationID","location","continent"}.issubset(loc.columns),
            "Original location fields changed")
    require(digest(gbif_zip)==GBIF_SHA,"MOIS original photo9 species ZIP checksum drift")
    with ZipFile(gbif_zip) as z:
        receipt=json.loads(z.read("receipt.json"))
        require(receipt["all_16_shards_source_verified"] is True and
                receipt["GBIF_occurrences_passed_all_filters"]==396142825,
                "Original MOIS geographic source provenance changed")
        flora=pd.read_csv(io.BytesIO(z.read(
            "all_37104_accepted_species_grid_coverage_and_photo9.csv.gz")),
            compression="gzip",dtype={"key":str,"species":str})
    require(len(flora)==37104 and flora.key.is_unique and flora.species.is_unique
            and set(flora.source_hue9)==HUES,"MOIS original 37,104 species source drift")
    flora["genus"]=flora.species.str.split().str[0]
    return r,tax,loc,flora,{
        "v4_source_zip_or_file_MD5":SRC_MD5,
        "v4_dataset_sha256":{name:digest(path) for name,path in files.items()},
        "MOIS_gbif_original_SHA256":GBIF_SHA,
        "v4_table_rows":{"records":len(r),"taxonomy":len(tax),"locations":len(loc)},
        "v4_record_unparseable_UTF8_replacement_count":sum(
            files["FirstRecords_dataset_public_v4.0.csv"].read_bytes()
            .decode("utf-8-sig",errors="replace").count("\ufffd") for _ in [0])
    }

def strict_links(r,tax,loc,flora):
    # LocationID is the only supported nonoverlap unit, never fuzzy country
    # inference. Verbatim location variants can repeat a canonical ID.
    normal=loc[["locationID","location","continent"]].copy()
    normal["locationID"]=normal.locationID.str.strip()
    loc_collisions=(normal.groupby("locationID")[["location","continent"]].nunique()>1).any(axis=1)
    badloc=set(loc_collisions[loc_collisions].index)
    norm=normal.loc[~normal.locationID.isin(badloc)].drop_duplicates("locationID")
    norm=norm.rename(columns={"location":"official_location","continent":"continent"})
    require(norm.locationID.is_unique,"LocationID official normalisation still not unique")
    original_tax_count=len(tax)
    tax=tax[(tax.taxaGroup=="Vascular plants") & (tax.GBIFstatus=="ACCEPTED") &
            (tax.GBIFmatchType=="EXACT") & (tax.GBIFtaxonRank=="SPECIES")].copy()
    tax["species"]=tax.species.str.strip()
    tax=tax[(tax.species!="") & (~tax.species.str.contains("\ufffd",regex=False))]
    tax_conflicts=tax.groupby("taxonID").species.nunique()
    conflicting_ids=set(tax_conflicts[tax_conflicts>1].index)
    tax=tax[~tax.taxonID.isin(conflicting_ids)].drop_duplicates("taxonID")
    taxonomy=tax[["taxonID","species"]].merge(
        flora[["key","species","genus","source_hue9","hue_status",
               "occupied_grid_cells"]],on="species",
        how="inner",validate="many_to_one")
    source=r.copy()
    source["first_year"]=pd.to_numeric(source.firstRecordEvent,errors="coerce")
    source["raw_year"]=source.verbatimFirstRecordEvent.str.strip()
    source["exact_year"]=source.raw_year.str.fullmatch(r"(1[5-9]\d{2}|20[0-2]\d)") & (
            source.first_year==pd.to_numeric(source.raw_year,errors="coerce"))
    source["high_confidence"]=source.confidenceFirstRecordEvent.eq("high confidence")
    source["record_is_present"]=source.occurrenceStatus.eq("present")
    source["record_is_alien"]=source.establishmentMeans.str.lower().eq("introduced")
    source["record_is_established"]=source.degreeOfEstablishment.str.lower().isin(
        ["established","invasive"])
    source["calendar_valid"]=source.first_year.between(1500,LAST_COMPLETE_YEAR) & (
        source.first_year==np.floor(source.first_year))
    source["public_allowed"]=source.accessRights.str.startswith("Free")
    source["id_clean"]=~(source.taxon.str.contains("\ufffd",regex=False) |
                         source.location.str.contains("\ufffd",regex=False))
    filters=("exact_year","high_confidence","record_is_present","record_is_alien",
             "record_is_established","calendar_valid","public_allowed","id_clean")
    masks=np.ones(len(source),dtype=bool)
    filter_stats={"raw_records":len(source)}
    for column in filters:
        masks&=source[column].fillna(False).to_numpy(bool)
        filter_stats["after_"+column]=int(masks.sum())
    source=source.loc[masks].copy()
    matched=source.merge(taxonomy,on="taxonID",how="inner",validate="many_to_one")
    matched=matched.merge(norm,on="locationID",how="inner",validate="many_to_one")
    require(matched.continent.ne("").all(),"Missing continent for accepted source location")
    matched["year"]=matched.first_year.astype(int)
    require((matched.year>=1500).all() and (matched.year<=2023).all(),
            "Invalid first-record calendar year slipped through")
    initial_dup=int(matched.duplicated(["key","locationID"]).sum())
    matched=matched.sort_values(["key","locationID","year"],kind="stable")
    matched=matched.drop_duplicates(["key","locationID"],keep="first")
    require(not matched.duplicated(["key","locationID"]).any(),"Duplicated source region after species mapping")
    report={
        "strict_filter_stages":filter_stats,
        "v4_original_taxa_count":original_tax_count,
        "taxa_matching_exact_accepted_species":int(tax.taxonID.nunique()),
        "v4_taxonID_with_conflicting_accepted_species_removed":len(conflicting_ids),
        "location_ID_with_conflicting_official_normalization_removed":len(badloc),
        "source_location_IDs_after_validated_dedup":len(norm),
        "exact_name_matched_MOIS_taxa":int(taxonomy.key.nunique()),
        "first_record_rows_strict_before_species_location_dedup":len(matched)+initial_dup,
        "species_location_duplicate_rows_removed":initial_dup,
        "strict_MOIS_first_record_species_locations":len(matched),
        "strict_MOIS_unique_species":int(matched.key.nunique()),
        "strict_MOIS_unique_source_locationID":int(matched.locationID.nunique()),
        "strict_MOIS_continent_counts":{str(k):int(v) for k,v in matched.continent.value_counts().items()},
        "first_record_year_min":int(matched.year.min()) if len(matched) else None,
        "first_record_year_max":int(matched.year.max()) if len(matched) else None,
        "no_unconfirmed_naturalisation_used":True,
        "no_imputed_year_used":True,
        "no_source_geographic_centroid_inferred":True,
    }
    return matched,flora,report

def region_curves(matched,flora):
    records=matched[["key","locationID","year","continent"]].copy()
    require(not records.duplicated(["key","locationID"]).any(),
            "Need one first record per species x source nonoverlapping location")
    first=records.groupby("key",sort=False).year.min().rename("first_year")
    records=records.join(first,on="key")
    records["lag"]=records.year-records.first_year
    require(records.lag.ge(0).all(),"All source region first records must postdate source species baseline")
    def counted(mask):
        return mask.groupby(records.key).sum().astype(int)
    results=pd.DataFrame(index=first.index)
    results["first_year"]=first
    results["n_recorded_regions"]=records.groupby("key").locationID.nunique()
    results["same_year_regions"]=counted(records.lag.eq(0))
    for horizon in (25,50,100):
        results["new_regions_"+str(horizon)+"yr"]=counted(
            records.lag.between(1,horizon,inclusive="both"))
        results["any_new_regions_"+str(horizon)+"yr"]=(
            results["new_regions_"+str(horizon)+"yr"]>0).astype(int)
    first_regions=records.loc[records.lag.eq(0)].copy()
    cont=first_regions.groupby("key").continent.agg(
        lambda c: str(c.iloc[0]) if c.nunique()==1 else "MULTIPLE_CONTINENTS")
    results["first_record_continent"]=cont
    result=results.reset_index().merge(
        flora[["key","species","genus","source_hue9","hue_status",
               "occupied_grid_cells"]],on="key",validate="one_to_one",how="left")
    require(result.source_hue9.notna().all(),"Species lost original photo hue")
    require((result.n_recorded_regions >= result.same_year_regions).all(),
            "Same-year records exceed known source regions")
    for horizon in (25,50,100):
        require((result["new_regions_"+str(horizon)+"yr"] <=
                 result.n_recorded_regions-result.same_year_regions).all(),
                "Counted downstream first-record regions exceed distinct known locations")
    return result

def model_sample(curves,horizon,remove_bronze=False,single_start=False,
                 consensus_only=False):
    require(horizon in YEAR_ENDPOINTS,"Unexpected first-record observation horizon")
    cohort=curves[curves.first_year.between(MAIN_YEAR_MIN,YEAR_ENDPOINTS[horizon])].copy()
    if remove_bronze:cohort=cohort.loc[cohort.source_hue9!="bronze"].copy()
    if single_start:cohort=cohort.loc[cohort.same_year_regions==1].copy()
    if consensus_only:
        cohort=cohort.loc[cohort.hue_status=="sampling_consensus_single"].copy()
    cohort["first_year_25bin"]=((cohort.first_year-MAIN_YEAR_MIN)//25).astype(int)
    cohort["cohort_control"]=(
        cohort.first_year_25bin.astype(str)+"|"+cohort.first_record_continent)
    return cohort.reset_index(drop=True)

def test_hue_association(cohort,flora,horizon,model,draws,seed):
    """FWL colour-to-new-regions effect after cohort×first-continent demeaning.
    One original species photo hue per global null permutation. No attempt to
    call this unconfounded dispersal or infer within-species selection.
    """
    require(model in {"genus_breadth","genus_only"},"Invalid null model")
    n=len(cohort)
    units=int(cohort.cohort_control.nunique())
    idx=pd.Index(flora.key.astype(str)).get_indexer(cohort.key.astype(str))
    require((idx>=0).all(),"Species in observed data missing from source")
    species_ids=np.unique(idx)
    full_labels=flora.source_hue9.to_numpy(str)
    gen=flora.genus.to_numpy(str)
    breadth=np.digitize(flora.occupied_grid_cells.to_numpy(int),[2,6,21,101])
    strata={}
    for i in species_ids:
        k=(gen[i],int(breadth[i])) if model=="genus_breadth" else (gen[i],)
        strata.setdefault(k,[]).append(i)
    movable=[np.array(m,dtype=int) for m in strata.values()
             if len(m)>=2 and len(set(full_labels[m]))>1]
    mobile_ids=np.concatenate(movable) if movable else np.array([],dtype=int)
    sample_audit={"source_species":len(species_ids),
                  "source_species_regions_with_valid_window":n,
                  "year_first_continent_strata":units,
                  "exchangeable_species_count":len(mobile_ids),
                  "exchangeable_species_fraction":len(mobile_ids)/max(n,1),
                  "exchangeable_label_strata":len(movable),
                  "firstrecord_year_range":[int(cohort.first_year.min()),int(cohort.first_year.max())] if n else [],
                  "n_unique_genera":int(cohort.genus.nunique()),
                  "n_source_firstrecord_continents":int(cohort.first_record_continent.nunique()),
                  "has_at_least_two_hue_categories":bool(cohort.source_hue9.nunique()>1),
                  "unconditional_fifty_year_country_spread_not_measured":True}
    if n<200 or units<10 or len(mobile_ids)<50:
        return {**sample_audit,"status":"insufficient primary exchangeability or sample"},None
    gcode=pd.factorize(cohort.cohort_control,sort=True)[0]
    gcount=np.bincount(gcode).astype(float)
    y=cohort["any_new_regions_"+str(horizon)+"yr"].to_numpy(float)
    ycenter=y-np.bincount(gcode,weights=y)[gcode]/gcount[gcode]
    def score(current_labels):
        white=(current_labels[idx]=="white").astype(float)
        xc=white-np.bincount(gcode,weights=white)[gcode]/gcount[gcode]
        denom=np.dot(xc,xc)
        return float(np.dot(xc,ycenter)/denom) if denom>1e-12 else np.nan
    observed=score(full_labels)
    require(np.isfinite(observed),"No within-cohort year/continent colour variation")
    rng=np.random.default_rng(seed)
    null=np.empty(draws)
    for b in range(draws):
        labs=full_labels.copy()
        for group in movable:
            labs[group]=full_labels[group][rng.permutation(len(group))]
        null[b]=score(labs)
    finite=np.isfinite(null)
    require(finite.sum()==draws,"Some permuted score undefined, cannot report rank-p")
    lower=(1+int(np.sum(null<=observed+1e-12)))/(draws+1)
    upper=(1+int(np.sum(null>=observed-1e-12)))/(draws+1)
    result={**sample_audit,
            "status":"completed",
            "source_phenotype":"one unverified photograph-derived colour per accepted species",
            "outcome":"binary first record in another public FirstRecords location 1-"+str(horizon)+"yr after species earliest public first record",
            "horizon_years":horizon,
            "n_permutations":draws,
            "n_white":int(cohort.source_hue9.eq("white").sum()),
            "n_nonwhite":int(cohort.source_hue9.ne("white").sum()),
            "observed_any_new_recorded_region_fraction":float(y.mean()),
            "observed_white_partial_probability_difference":observed,
            "null_mean":float(null.mean()),
            "null_adjusted_difference":float(observed-null.mean()),
            "null_reference95":np.percentile(null,[2.5,97.5]).tolist(),
            "p_rank_lower":lower,"p_rank_upper":upper,
            "p_rank_two_sided":min(1.,2*min(lower,upper)),
            "followup_right_censoring_guard":int(cohort.first_year.max())+horizon<=LAST_COMPLETE_YEAR,
            "same_species_hue_coherent_across_regions":True,
            "no_biological_dispersal_or_colonization_rate_observed":True,
            "not_a_causal_selection_test":True,
            "no_effect_confidence_interval":True}
    return result,null

def run(args):
    out=Path(args.out)
    require(not out.exists(),"Refuse to overwrite prior signed regional chronology result")
    records,tax,locations,flora,sources=source_tables(args.v4,args.gbif)
    strict,flora,link_audit=strict_links(records,tax,locations,flora)
    curves=region_curves(strict,flora)
    out.mkdir(parents=True)
    # CC BY4: provide derived species-level chronology for inspection; do not
    # mirror publisher's original 90,023-row source data.
    curves.to_csv(out/"MOIS_species_recorded_region_accumulation.csv.gz",
                  index=False,compression={"method":"gzip","mtime":0})
    samples={}
    draws=[]
    configs=[
        ("main_50",50,False,False,False,"genus_breadth"),
        ("no_bronze_50",50,True,False,False,"genus_breadth"),
        ("single_initial_region_50",50,False,True,False,"genus_breadth"),
        ("genus_only_50",50,False,False,False,"genus_only"),
        ("horizon25",25,False,False,False,"genus_breadth"),
        ("horizon100",100,False,False,False,"genus_breadth"),
        ("photo_consensus_50",50,False,False,True,"genus_breadth")
    ]
    for index,(label,horizon,bronze,single,consensus,model) in enumerate(configs):
        data=model_sample(curves,horizon,bronze,single,consensus)
        result,sim=test_hue_association(data,flora,horizon,model,args.permutations,
                                        SEED+index*137)
        result.update({"comparison_label":label,"remove_bronze":bronze,
                       "require_single_initial_location":single,
                       "require_photo_sampling_consensus":consensus,"null_model":model})
        samples[label]=result
        if sim is not None:
            draws.append(pd.DataFrame({"comparison":label,"draw_id":np.arange(1,len(sim)+1),
                                       "white_partial_probability_difference":sim}))
        print("MOIS_FIRSTRECORDS_V40_SPREAD_RESULT",json.dumps(result,ensure_ascii=False),flush=True)
    if draws:
        pd.concat(draws,ignore_index=True).to_csv(out/"species_hue_synchronous_null_draws.csv.gz",
            index=False,compression={"method":"gzip","mtime":0})
    # No raw copyright source mirror: only aggregate time-bin composition by
    # photographic category, and species-first-record summary derived under CC BY.
    curves["initial_period"]=pd.cut(curves.first_year,bins=[1500,1800,1850,1900,1950,1975,2000,2024],
        labels=["1500-1799","1800-1849","1850-1899","1900-1949","1950-1974","1975-1999","2000-2023"],
        right=False)
    grouped=curves.groupby(["initial_period","source_hue9"],observed=True).agg(
        n_species=("key","size"),mean_recorded_regions=("n_recorded_regions","mean"),
        median_recorded_regions=("n_recorded_regions","median")).reset_index()
    grouped.to_csv(out/"original_first_record_cohort_photo9_composition.csv",index=False)
    metadata={
        "source":"FirstRecords v4.0 2026 DOI 10.5281/zenodo.18759840",
        "source_cc_by_4_0_attribution_required":True,
        "not_introduction_or_clonal_spread_dates":True,
        "not_flower_color_evolution_or_phylogenetic_radiation":True,
        "comparisons_were_specified_before_testing":True,
        "primary":"main_50",
        "followup_window_years":50,
        "strict_confidence_original_year":True,
        "source_data_sha":sources,
        "lineage_audit":link_audit,
        "species_curve_summary":{
            "n_species_with_exact_accepted_mois_link":len(curves),
            "with_two_or_more_distinct_source_locations":int(curves.n_recorded_regions.ge(2).sum()),
            "with_same_year_multiple_regions":int(curves.same_year_regions.ge(2).sum()),
            "min_source_first_year":int(curves.first_year.min()) if len(curves) else None,
            "max_source_first_year":int(curves.first_year.max()) if len(curves) else None,
            "horizon_50_eligible_species":int(curves.first_year.between(1800,1973).sum())
        },
        "comparison_results":samples,
        "rank_p_not_corrected_for_all_discovery_searches":True,
        "no_fossil_refugia_biomass_pollinator_causal_inference":True,
    }
    (out/"summary.json").write_text(json.dumps(metadata,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    digest_dict={p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out/"receipt.json").write_text(json.dumps({
        "original_source_SHA":sources,
        "code_sha256":digest(__file__),
        "files_sha256":digest_dict,
        "synchronous_species_permutation":True,
        "raw_original_first_records_not_republished":True
    },ensure_ascii=False,indent=2)+"\n")
    print("MOIS_FIRSTRECORDS_V40_REGION_ACCUMULATION_COMPLETE",
          json.dumps({"source":link_audit,
                      "curve":metadata["species_curve_summary"],
                      "primary":samples["main_50"]},ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--v4",required=True)
    p.add_argument("--gbif",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--permutations",type=int,default=999)
    args=p.parse_args()
    require(args.permutations>=99,"Cannot run tiny null and report real p")
    run(args)
