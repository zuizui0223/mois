#!/usr/bin/env python3
"""MOIS: automatic, no-human ROI9 selection-bias diagnostics on frozen pixels.

All 1,181,146 genuine existing original photo attempts are reused, with one
outcome-blind canonical photo per observation. These are *photographic labels*
conditional on frozen ROI9, not verified open petals or native community colour.
The selected four-arm source has no global-flora first-stage inclusion weights.

For each arm x broad latitude x hemisphere-local-season stratum:
  - exact observed white fraction among strict9 successes,
  - assumption-conditional MAR poststratification by all observations,
  - MNAR delta-odds sensitivity including unmeasured outcome uncertainty,
  - identification bounds, nonoverlap and sparsity.
Strict9 failures are never relabelled or represented as measured flowers.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_roi9_observation_frame_20261029 import (
    ARMS, COLS, EXPECTED_PHOTOS, EXPECTED_SUCCESS, sha_file,
    observation_frame, coarse_zone, require,
)
from scripts.analysis.mois_1181146_original_roi9_missingness_gap_decision import (
    latitude_band, recorded_local_season,
)

SCHEMA = "mois_real1181146_observation_roi9_automatic_bias_bounds_v1"
HUES = ("white", "yellow_orange", "red_pink", "blue_purple")
DELTA_GRID = (-2., -1., -.5, 0., .5, 1., 2.)
N_OBS = 1032882
N_SUCCESS = 401703
N_UNKNOWN = 631179
SOURCE_ARMS = ARMS

def logit_tilt(p, delta):
    """Numerically stable odds tilt; sensitivity, NOT inferred unknown colour."""
    p = float(p)
    require(0. <= p <= 1. and math.isfinite(delta),
            "Invalid unknown-colour sensitivity parameter")
    if p == 0. or p == 1.:
        return p
    log_odds = math.log(p / (1. - p)) + float(delta)
    return 1. / (1. + math.exp(-max(-700., min(700., log_odds))))

def source_groups(f):
    """Group data at source selection strata, keeping missing colour unknown."""
    keys = ["source_design", "coarse_zone", "local_season"]
    require(f.photo_id.is_unique and f.observation_id.is_unique,
            "Canonical photos/observations are not distinct")
    require(f.strict9.notna().all() and pd.api.types.is_bool_dtype(f.strict9),
            "ROI9 technical status is not genuine boolean")
    require((~f.strict9 | f.hue.isin(HUES)).all(),
            "Strict9 success with invalid recorded original colour label")
    f = f.copy()
    f["_white"] = (f.strict9 & f.hue.eq("white"))
    f["_missing"] = ~f.strict9
    t = f.groupby(keys, observed=True, dropna=False, sort=True).agg(
        n_observations=("observation_id", "size"),
        strict9_success=("strict9", "sum"),
        strict9_missing=("_missing", "sum"),
        observed_white=("_white", "sum"),
        n_taxa=("inat_taxon_id", "nunique"),
    ).reset_index()
    for c in ("n_observations", "strict9_success", "strict9_missing",
              "observed_white", "n_taxa"):
        t[c] = t[c].astype(int)
    require(len(t)>0 and
            (t.n_observations == t.strict9_success + t.strict9_missing).all() and
            (t.observed_white<=t.strict9_success).all(), "Broken photo totals")
    t["roi_success_rate"] = t.strict9_success / t.n_observations
    t["white_on_success"] = t.observed_white / t.strict9_success.replace(0,np.nan)
    t["has_observed_colour_support"] = t.strict9_success.gt(0)
    t["sparse_strict9_lt30"] = t.strict9_success.lt(30)
    # Preserve exact groups where no valid ROI9 exists, rather than "predicting" a colour.
    return t

def group_summary(t, *, scope):
    """Assumption-labelled standardization and missing-not-at-random bounds."""
    N = int(t.n_observations.sum())
    K = int(t.strict9_success.sum())
    U = int(t.strict9_missing.sum())
    W = int(t.observed_white.sum())
    require(N == K+U and K>=W and N>0, "Invalid grouped ROI9 counts")
    available = t.loc[t.has_observed_colour_support]
    N_supported = int(available.n_observations.sum())
    N_unsupported = N-N_supported
    S = int((t.strict9_success<30).sum())
    unknown_unidentified = int(t.loc[~t.has_observed_colour_support, "strict9_missing"].sum())
    mar_supported_total = float((available.n_observations *
                                available.white_on_success).sum())
    mar = None
    if N_unsupported == 0:
        mar = mar_supported_total/N
    mar_supported_fraction = (
        mar_supported_total/N_supported if N_supported else None
    )
    # Across technical missing photos no actual Y is observed: these bounds
    # require only the restricted photographic-label definition, not MAR.
    lower = W/N
    upper = (W+U)/N
    require(0<=lower<=upper<=1, "Impossible identification bound")
    out = {
        "scope":scope,
        "observations":N, "roi9_success":K, "roi9_unknown":U,
        "strict9_classifiability":K/N,
        "white_strict9_success":W,
        "complete_case_white_rate":W/K if K else None,
        "nonparametric_white_fraction_identification_bounds":[lower,upper],
        "strata":int(len(t)),
        "sparse_strata_success_lt30":S,
        "strata_without_any_valid_roi9":int((~t.has_observed_colour_support).sum()),
        "observations_without_stratum_colour_support":N_unsupported,
        "fraction_observation_frame_with_colour_support":N_supported/N,
        "MAR_stratum_poststandardized_photo_white_rate":mar,
        "MAR_supported_only_photo_white_rate":mar_supported_fraction,
        "MAR_unidentified_strata_interval": [
            mar_supported_total/N,(mar_supported_total+N_unsupported)/N
        ],
        "MNAR_delta_log_odds_scenarios": {},
    }
    # For very sparse cells we use a Jeffreys reference success fraction only
    # for sensitivity shifts. The exact unregularized MAR point is separate.
    for delta in DELTA_GRID:
        assumptions=[]
        for r in available.itertuples(index=False):
            p_smooth = (r.observed_white+0.5)/(r.strict9_success+1.0)
            p_missing = logit_tilt(p_smooth,delta)
            assumptions.append(r.observed_white + r.strict9_missing * p_missing)
        # Scenarios include exact already-measured white success labels;
        # technically unknown photographs get only *hypothetical* colours.
        total = float(sum(assumptions))
        assert total <= N+1e-8
        key=("%+.1f" % delta)
        out["MNAR_delta_log_odds_scenarios"][key] = {
            "scenario_only_not_measured":True,
            "scenario_white_fraction_if_all_strata_supported":(
                total/N if N_unsupported==0 else None
            ),
            "interval_if_no_colour_support":[
                total/N,(total+N_unsupported)/N
            ],
            "assumed_log_odds_shift_for_unknown_vs_success":delta,
        }
    out["fraction_unknown_needed_for_half_white_label"] = (
        (0.5*N-W)/U if U else None
    )
    out["half_white_possible_with_arbitrary_missing_colours"] = (
        lower<=0.5<=upper
    )
    return out

def taxon_balance(f):
    """Taxon-equal summary only for species with measured source hues."""
    z=f.loc[f.inat_taxon_id.ne("")].copy()
    z["_white"] = z.strict9 & z.hue.eq("white")
    t=z.groupby(["source_design","inat_taxon_id"], observed=True,sort=True).agg(
        attempted=("photo_id","size"),
        success=("strict9","sum"),
        known_white=("_white","sum"),
    ).reset_index()
    summary=[]
    for arm in SOURCE_ARMS:
        arm_t=t.loc[t.source_design.eq(arm)]
        supported=arm_t.loc[arm_t.success.gt(0)]
        out={
            "source_design":arm,
            "observed_taxa":int(len(arm_t)),
            "taxa_with_strict9":int(len(supported)),
            "taxa_no_strict9":int((arm_t.success==0).sum()),
            "taxa_without_strict9_cannot_be_corrected_to_truth":True,
            "equal_taxon_white_rate_among_supported_taxa":(
                float((supported.known_white/supported.success).mean())
                if len(supported) else None
            ),
        }
        summary.append(out)
    return summary

def run(photo_path, receipt_path, output, *, verify_counts=True):
    photo_path,receipt_path,output=map(Path,(photo_path,receipt_path,output))
    require(not output.exists(), "Never replace original bias analysis output")
    frozen=json.loads(receipt_path.read_text())
    require(frozen.get("schema") ==
        "mois_signed1181146_trueROI9_bulk128_original_metadata_georestore_v2"
        and frozen.get("same_original_image_ID_and_original_9palette_unmodified") is True
        and frozen.get("source_georestored_original_parquet_SHA256")==sha_file(photo_path),
        "Input is not the exact corrected v2 source Parquet and signed receipt")
    cols = list(COLS) + ["original_roi_four_hue"]
    d=pd.read_parquet(photo_path,columns=cols)
    if verify_counts:
        require(len(d)==EXPECTED_PHOTOS and int(d.strict9.sum())==EXPECTED_SUCCESS,
                "Original full frozen numeric strict9 census changed")
    f, census = observation_frame(d)
    # The canonical choice has original row indices: preserve original colours
    # only for real strict9 successes, retaining hue text as unknown otherwise.
    f["hue"] = d.loc[f.index,"original_roi_four_hue"].fillna("").astype(str).to_numpy()
    if verify_counts:
        require(len(f)==N_OBS and int(f.strict9.sum())==N_SUCCESS and
                int((~f.strict9).sum())==N_UNKNOWN,
                "Real canonical observation/ROI9 counts changed")
    f["latitude_zone"] = latitude_band(f.latitude)
    f["coarse_zone"] = f.latitude_zone.map(coarse_zone)
    f["local_season"] = recorded_local_season(f.observed_on,f.latitude)
    strata=source_groups(f)
    records=[group_summary(strata,scope="four_arm_selected_observation_pool_NONPROBABILITY")]
    for arm in SOURCE_ARMS:
        subset=strata.loc[strata.source_design.eq(arm)]
        require(len(subset)>0,"Original source arm disappeared: "+arm)
        records.append(group_summary(subset,scope=arm))
    for zone in sorted(strata.coarse_zone.unique()):
        records.append(group_summary(strata.loc[strata.coarse_zone.eq(zone)],
                                     scope="zone:"+str(zone)))
    for season in sorted(strata.local_season.unique()):
        records.append(group_summary(strata.loc[strata.local_season.eq(season)],
                                     scope="season:"+str(season)))
    for arm in SOURCE_ARMS:
        for zone in sorted(strata.loc[strata.source_design.eq(arm),"coarse_zone"].unique()):
            records.append(group_summary(
                strata.loc[strata.source_design.eq(arm)&strata.coarse_zone.eq(zone)],
                scope="arm_zone:"+arm+"|"+str(zone)))
    for arm in SOURCE_ARMS:
        for season in sorted(strata.loc[strata.source_design.eq(arm),"local_season"].unique()):
            records.append(group_summary(
                strata.loc[strata.source_design.eq(arm)&strata.local_season.eq(season)],
                scope="arm_season:"+arm+"|"+str(season)))
    overall=records[0]
    require(overall["observations"]==len(f) and overall["roi9_success"]==
            census["canonical_success_count"] and
            overall["roi9_unknown"]==census["canonical_unknown_count"] and
            sum(strata.n_observations)==len(f), "Incorrect overall ROI9 denominator")
    output.mkdir(parents=True)
    stratum_csv=output/"mois_auto_roi9_source_geoseason_strata.csv"
    strata.to_csv(stratum_csv,index=False)
    records_file=output/"mois_auto_roi9_bias_by_scope.json"
    records_file.write_text(json.dumps(records,indent=2,ensure_ascii=False,
                                      sort_keys=True,allow_nan=False)+"\n")
    taxon_summary=taxon_balance(f)
    taxon_file=output/"mois_auto_roi9_taxon_supported_only_summary.json"
    taxon_file.write_text(json.dumps(taxon_summary,indent=2,sort_keys=True,
                                    allow_nan=False)+"\n")
    manifest={
        "schema":SCHEMA,
        "status":"REAL_AUTOMATIC_SOURCE_OBSERVATION_BIAS_DIAGNOSTIC_COMPLETED",
        "input_original_georestored_parquet_SHA256":sha_file(photo_path),
        "source_restoration_receipt_SHA256":sha_file(receipt_path),
        "source_photo_attempts":census["source_original_photo_attempts"],
        "source_distinct_observations":census["distinct_represented_observations"],
        "suppressed_same_observation_extra_photos":
            census["additional_photos_of_already_represented_observations"],
        "frozen_roi9_success_observations":census["canonical_success_count"],
        "frozen_roi9_unknown_observations":census["canonical_unknown_count"],
        "actual_human_labels":0,
        "new_image_downloads":0,"new_image_inferences":0,
        "strict9_unknowns_upgraded_to_known":0,
        "four_arm_mixture_not_native_global_probability_sample":True,
        "automatic_MAR_correction_is_assumption_dependent":True,
        "MNAR_true_missing_colour_not_identified":True,
        "botanical_flower_or_native_status_NOT_verified":True,
        "white_label_definition":"original_four_hue=white among frozen original strict9; NOT white petal truth",
        "geography_season_strata":int(len(strata)),
        "overall":overall,
        "output_SHA256":{p.name:sha_file(p) for p in
                         (stratum_csv,records_file,taxon_file)},
    }
    (output/"mois_auto_roi9_bias_receipt.json").write_text(
        json.dumps(manifest,indent=2,sort_keys=True,ensure_ascii=False,
                   allow_nan=False)+"\n"
    )
    print("MOIS_REAL_AUTO_NOHUMAN_ROI9_BIAS_COMPLETE",
          json.dumps({
            "obs":len(f),"success":int(f.strict9.sum()),
            "strata":len(strata),
            "complete_case_white":overall["complete_case_white_rate"],
            "MAR":overall["MAR_stratum_poststandardized_photo_white_rate"],
            "MNAR_worst_case":overall["nonparametric_white_fraction_identification_bounds"],
            "no_human_labels":True,
          },sort_keys=True),flush=True)
    return manifest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-parquet",type=Path,required=True)
    p.add_argument("--georestored-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.photo_parquet,args.georestored_receipt,args.output)

if __name__=="__main__":
    main()
