#!/usr/bin/env python3
"""MOIS real 313535-source ROI-hue missingness: solar-stratified wet/dry tipping.

Use ONLY the source-released nine rainfall x solar photo-attempt cells that
include prior original four-hue successes, white counts and unknown counts.
NO images downloaded, NO genetic/petal/phenology imputation or novel model.
These are illustrative explicit counterfactual assumptions for unmeasurable
photo hue, NOT estimated missing-white outcomes. The annual precipitation
(+1.99 pp/SD) genus-999 maxT model is on a DIFFERENT 309439-photo population
and is NOT numerically corrected by these wet/dry sensitivity contrasts.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import pandas as pd

SCHEMA="mois_source_313535_photo_binary_white_ROI_MNAR_wetdry_solar3_tipping_v1"
PARENT_SOURCE_TAG="mois-original313535-four-hues-water-solar-genus-20261018-v1"
PARENT_CSV_SHA="04728fc4426b160a73cb2cff56094f373beac1b5825c861675016712425daafb"
PARENT_N=313535
PARENT_CLASSIFIED=151239
PARENT_WHITE=71427
REQUIRED=(
 "rainfall_tertile","solar_tertile","original_processed_photo_attempts",
 "real_original_classifiable_ROI_photos",
 "technically_unknown_hue_photo_attempts","white_photo_ROIs",
 "yellow_orange_photo_ROIs","red_pink_photo_ROIs","blue_purple_photo_ROIs"
)
RAIN=("low","middle","high")
SOLAR=("low","middle","high")
DELTA_PPT=(0,5,10,15,20)

def require(ok,msg):
    if not bool(ok):
        raise ValueError("MOIS_GEB_SOURCE_MNAR: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for x in iter(lambda:f.read(1<<20),b""):h.update(x)
    return h.hexdigest()

def validate(cells,*,verify_global=True):
    require(set(REQUIRED)<=set(cells.columns),"Original rain/solar attempt+unknown 3x3 columns absent")
    require(len(cells)==9 and cells[["rainfall_tertile","solar_tertile"]].drop_duplicates().shape[0]==9,
            "Exactly nine original pre-classification rain by solar cells required")
    require(set(cells.rainfall_tertile)==set(RAIN)
            and set(cells.solar_tertile)==set(SOLAR),
            "Unexpected original rain/solar tertile levels or absent cells")
    a=cells.copy()
    num=REQUIRED[2:]
    for k in num:
        a[k]=pd.to_numeric(a[k],errors="raise")
        require(a[k].notna().all() and a[k].ge(0).all(),
                "Missing or negative original photograph counts in "+k)
        require(a[k].eq(a[k].astype(int)).all(),
                "Noninteger original photograph counts in "+k)
    require((a.original_processed_photo_attempts
              ==a.real_original_classifiable_ROI_photos+
                 a.technically_unknown_hue_photo_attempts).all(),
            "True original photo attempts do not equal measured+unknown")
    require((a.real_original_classifiable_ROI_photos ==
             a[["white_photo_ROIs","yellow_orange_photo_ROIs",
                "red_pink_photo_ROIs","blue_purple_photo_ROIs"]].sum(axis=1)).all(),
            "Known hue photo counts do not exhaust original source four categories")
    require(a.real_original_classifiable_ROI_photos.gt(0).all()
            and a.technically_unknown_hue_photo_attempts.gt(0).all(),
            "Cannot define missing-label sensitivity without observed hue and unknowns")
    if verify_global:
        total=int(a.original_processed_photo_attempts.sum())
        classified=int(a.real_original_classifiable_ROI_photos.sum())
        white=int(a.white_photo_ROIs.sum())
        require((total,classified,white)==(313204,151107,71368),
                "Actual original 313535 photo source climate-overlap totals mutated")
        require(PARENT_N-total==331 and
                PARENT_CLASSIFIED-classified==132 and
                PARENT_WHITE-white==59,
                "Full source photos with missing climate strata drifted")
    return a

def row_for(cells,rain,solar):
    result=cells[(cells.rainfall_tertile==rain)&
                 (cells.solar_tertile==solar)]
    require(len(result)==1,"Missing one original water/solar sampling cell")
    return result.iloc[0]

def extract(cell):
    total=int(cell.original_processed_photo_attempts)
    good=int(cell.real_original_classifiable_ROI_photos)
    white=int(cell.white_photo_ROIs)
    unknown=total-good
    require(0<=white<=good<total,"Invalid true white/known/source photo denominator")
    return {"attempted":total,"known":good,"white":white,
            "unknown":unknown,"classifiable":good/total,
            "white_known":white/good,"unknown_fraction":unknown/total}

def mean_white_under_unknown(g,white_unknown):
    """Explicit hypothetical missing-state label, NOT a measured phenotype."""
    require(0<=white_unknown<=1,"Unmeasured white-label sensitivity fraction must be in [0,1]")
    return (g["white"]+g["unknown"]*white_unknown)/g["attempted"]

def one_solar(dry,wet,solar):
    d=extract(dry);w=extract(wet)
    observed=100*(w["white_known"]-d["white_known"])
    # Assumption: dry unknown photos would be white-labelled at the SAME
    # rate as classified dry photos (MAR only within dry rainfall x solar).
    # This assumption is deliberately tested, not claimed to be true.
    assumed_dry_unknown=d["white_known"]
    wet_unknown_tipping=(
        mean_white_under_unknown(d,assumed_dry_unknown)*w["attempted"]-w["white"]
    )/w["unknown"]
    tipping_delta=wet_unknown_tipping-w["white_known"]
    require(0<=wet_unknown_tipping<=1 and tipping_delta<=0,
            "Actual source unexpectedly already reverses wet-dry under dry MAR")
    out={
        "solar_tertile":solar,
        "original_photo_attempts_dry":d["attempted"],
        "original_photo_attempts_wet":w["attempted"],
        "original_four_hue_known_dry":d["known"],
        "original_four_hue_known_wet":w["known"],
        "original_four_hue_unknown_dry":d["unknown"],
        "original_four_hue_unknown_wet":w["unknown"],
        "original_four_hue_ROI_success_dry":d["classifiable"],
        "original_four_hue_ROI_success_wet":w["classifiable"],
        "photo_white_share_among_classifiable_dry":d["white_known"],
        "photo_white_share_among_classifiable_wet":w["white_known"],
        "photographic_white_wet_minus_dry_completecase_pp":observed,
        "wet_missing_white_probability_needed_to_tie_if_dry_missing_is_MAR":
            wet_unknown_tipping,
        "wet_missing_minus_wet_observed_probability_needed_to_tie_pp":
            100*tipping_delta,
        "wet_unknown_probability_tipping_is_UNMEASURED_NOT_ESTIMATED":True,
        # Conservative *hypothetical* unconditional-photo label limits:
        "unidentified_wet_minus_dry_low_all_wet_unknown_0_dry_unknown_1_pp":
            100*(mean_white_under_unknown(w,0)-mean_white_under_unknown(d,1)),
        "unidentified_wet_minus_dry_high_all_wet_unknown_1_dry_unknown_0_pp":
            100*(mean_white_under_unknown(w,1)-mean_white_under_unknown(d,0)),
    }
    scenarios=[]
    for step in DELTA_PPT:
        delta=step/100
        # One-sided: all observed dry missing resemble classified dry;
        # wet missing are increasingly white-depleted.
        p_wet=max(0,w["white_known"]-delta)
        one=100*(mean_white_under_unknown(w,p_wet)-
                 mean_white_under_unknown(d,d["white_known"]))
        # Two-sided sensitivity, not an empirical estimate:
        # all wet missing are delta less white; dry missing delta more white.
        p_dry=min(1,d["white_known"]+delta)
        sym=100*(mean_white_under_unknown(w,p_wet)-
                 mean_white_under_unknown(d,p_dry))
        for name,value,ud,uw in (
            ("wet_missing_white_minus_delta_dry_missing_MAR",one,
             d["white_known"],p_wet),
            ("symmetric_wet_depletion_dry_enrichment",sym,p_dry,p_wet)
        ):
            scenarios.append({
                "solar_tertile":solar,
                "counterfactual_MNAR_scenario":name,
                "assumed_unmeasured_delta_percentage_points":step,
                "assumed_white_fraction_for_dry_missing":ud,
                "assumed_white_fraction_for_wet_missing":uw,
                "hypothetical_photo_white_wet_minus_dry_pp":value,
                "direction_remains_positive_under_assumption":value>0,
                "counterfactual_not_a_measurement":True
            })
    return out,scenarios

def calculate(table,*,check_source=True):
    d=validate(table,verify_global=check_source)
    rows=[];scenarios=[]
    for solar in SOLAR:
        result,alts=one_solar(row_for(d,"low",solar),row_for(d,"high",solar),solar)
        rows.append(result);scenarios.extend(alts)
    reg=pd.DataFrame(rows)
    sensitivity=pd.DataFrame(scenarios)
    # Fixed pre-colour (source attempted) weights among wet/dry in each
    # solar band. This standardization intentionally excludes middle rain
    # and does not adjust taxon, genus or geographic support.
    support=(reg.original_photo_attempts_dry+reg.original_photo_attempts_wet)
    weights=support/support.sum()
    observed=float((weights*reg.photographic_white_wet_minus_dry_completecase_pp).sum())
    wet_unknown_mean=float((weights*(1-reg.original_four_hue_ROI_success_wet)).sum())
    one_sided_drop_to_zero=observed/wet_unknown_mean
    require(observed>0 and 0<one_sided_drop_to_zero<100,
            "Actual source unexpected wet-dry support or white-label direction")
    weighted=[]
    for (scenario,step),g in sensitivity.groupby(
         ["counterfactual_MNAR_scenario","assumed_unmeasured_delta_percentage_points"],
         observed=True):
        mapping=g.set_index("solar_tertile")
        changed=float(sum(weights.iloc[i]*
          mapping.loc[solar,"hypothetical_photo_white_wet_minus_dry_pp"]
          for i,solar in enumerate(SOLAR)))
        weighted.append({
            "counterfactual_MNAR_scenario":scenario,
            "assumed_unmeasured_delta_percentage_points":int(step),
            "photo_attempt_equal_solar_support_standardized_wet_minus_dry_pp":changed,
            "direction_remains_positive_under_assumption":changed>0,
            "solar_support_weights_formed_pre_classification":True,
            "not_adjusted_for_taxonomic_composition_or_worldwide_inclusion":True,
        })
    summary={
      "schema":SCHEMA,
      "source_release":PARENT_SOURCE_TAG,
      "true_parent_original_photo_attempts":PARENT_N,
      "true_parent_original_four_hue_classifiable":PARENT_CLASSIFIED,
      "three_by_three_original_photo_attempts_with_both_env_axes":int(
          d.original_processed_photo_attempts.sum()),
      "three_by_three_original_photo_hue_known":int(
          d.real_original_classifiable_ROI_photos.sum()),
      "three_by_three_original_photo_ROI9_UNKNOWN_not_same_definition_as_four_hue":True,
      "three_by_three_original_photo_four_hue_unknown":int(
          d.technically_unknown_hue_photo_attempts.sum()),
      "source_photos_excluded_from_3by3_because_climate_missing":331,
      "number_wet_dry_comparison_solar_strata":3,
      "standardization":"sum attempted wet+dry photos by solar tertile / all wet+dry attempted photos, prespecified before four-hue availability",
      "completecase_attempt_weighted_solar_stratified_wet_minus_dry_white_pp":observed,
      "one_sided_assumed_white_fraction_depletion_of_wet_unmeasured_to_erase_contrast_pp":
          one_sided_drop_to_zero,
      "one_sided_zero_assumes_dry_missing_has_observed_dry_white_fraction":True,
      "MNAR_scenario_results_do_not_correct_genus_adjusted_309439_worldclim_model":True,
      "missing_photo_hue_latent_unverified_flower_presence":True,
      "strict_nine_palette_photo_ROI_not_inferred_from_four_hue":True,
      "NO_unknown_original_photo_assigned_white_or_chromatic":True,
      "new_RGB_photo_downloads":0,
      "new_statistical_environment_model_fits":0,
    }
    return reg,sensitivity,pd.DataFrame(weighted),summary

def run(source,out):
    source,out=Path(source),Path(out)
    require(not out.exists(),"Refuse overwrite of existing source-conditional MNAR audit")
    require(sha(source)==PARENT_CSV_SHA,
            "Original real313535 nine-cell photo-climate table SHA changed")
    original=pd.read_csv(source,low_memory=False)
    strata,scenarios,weighted,summary=calculate(original)
    out.mkdir(parents=True)
    outputs={
        "real313535_wetdry_hue_unknown_tipping_by_solar.csv":strata,
        "real313535_explicit_MNAR_assumptions_by_solar.csv":scenarios,
        "real313535_precolour_attempt_weighted_MNAR_sensitivity.csv":weighted
    }
    summary["signed_source_CSV_SHA256"]=sha(source)
    summary["output_SHA256"]={}
    for name,df in outputs.items():
        dest=out/name
        df.to_csv(dest,index=False,float_format="%.12g")
        summary["output_SHA256"][name]=sha(dest)
    (out/"real313535_fourhue_missingness_MNAR_tipping_receipt.json").write_text(
        json.dumps(summary,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_REAL313535_SOURCE_MNAR_SENSITIVITY_NOT_GLOBAL_CAUSAL",
          json.dumps({
            "white_pp_wetdry_source_attempt_solar_adjusted":
                summary["completecase_attempt_weighted_solar_stratified_wet_minus_dry_white_pp"],
            "wet_unknown_white_probability_depletion_to_tie_pp":
                summary["one_sided_assumed_white_fraction_depletion_of_wet_unmeasured_to_erase_contrast_pp"],
            "actual_source_3by3_attempts":
                summary["three_by_three_original_photo_attempts_with_both_env_axes"]
          },sort_keys=True),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--signed-source-cells",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    x=p.parse_args()
    run(x.signed_source_cells,x.output)

if __name__=="__main__":
    main()
