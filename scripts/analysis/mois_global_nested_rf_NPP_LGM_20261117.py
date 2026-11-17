#!/usr/bin/env python3
"""MOIS genuine global nested RF: contemporary climate vs productivity vs LGM.

One accepted species = one row and one vote; 3 genus-disjoint folds shared
between *all* models. Photographic modal 9-hue is an unverified image-ROI proxy.
LGM climate is climate at present GBIF range cells, NOT reconstructed LGM range.
No photo status enters model X (avoid making label quality into a predictor).

Both inputs must match the exact signed 2026-11-16/17 Actions artifacts.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import (
    roc_auc_score, log_loss, balanced_accuracy_score, average_precision_score,
    brier_score_loss,
)

MAIN_SHA = "1da364ff63c506f9e2c06d1d71c9609a2c48fec7a428a21ba2c3c4c6ff8f4800"
MAIN_FILE = "GBIF37104_photo9_species_NPP23years_LGM_AGB_NDVI.parquet"
ENV_SHA = "d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
ENV_FILE = "BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
COLORS = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
SEED = 20261117
FOLDS = 3
NPP = "NPP_mean_GBIF_2001_2023_kgC_m2_yr"
NPP_COVER = "NPP_fraction_occupied_cells_18plus_years"
LGM_T = "threeGCM_delta_temperature_C_median"
LGM_P = "threeGCM_delta_log_precip_median"
TEMP = [f"BIO{i:02d}" for i in range(1, 12)] + [
    "tavg_monthly_mean", "tavg_monthly_range", "tavg_monthly_sd",
    "tavg_warmest_month_C",
]
MOIST = [f"BIO{i:02d}" for i in range(12, 20)] + [
    "prec_monthly_mean", "prec_monthly_range", "prec_monthly_sd",
    "prec_annual_sum_of_monthly_normal_mm", "vapr_monthly_mean",
    "vapr_monthly_range", "vapr_monthly_sd", "vpd_kPa_approx",
]
SOLAR_WIND = [
    "srad_monthly_mean", "srad_monthly_range", "srad_monthly_sd",
    "wind_monthly_mean", "wind_monthly_range", "wind_monthly_sd",
]
GEOGRAPHY_EFFORT = ["mean_abs_latitude_over_occupied_cells", "log_GBIF_occurrences"]
BASE = TEMP + MOIST + SOLAR_WIND + GEOGRAPHY_EFFORT
PRODUCTIVITY = [
    NPP, "MODIS_NDVI_2001_2021_mean", "MODIS_NDVI_2001_2021_sd",
    "log_ESACCI_AGB_woody_Mg_ha_2020",
    "ESA_AGB_woody_positive_area_fraction_2020",
    "ESA_AGB_2020_minus_2010_Mg_ha",
]
HISTORY = [LGM_T, LGM_P]
STAGES = {
    "A_current_full_climate": BASE,
    "B_plus_NPP_NDVI_AGB": BASE + PRODUCTIVITY,
    "C_plus_LGM_3GCM": BASE + PRODUCTIVITY + HISTORY,
}
COMPARISONS = [
    ("A_current_full_climate", "B_plus_NPP_NDVI_AGB"),
    ("B_plus_NPP_NDVI_AGB", "C_plus_LGM_3GCM"),
]
N_ROI_COLORS = len(COLORS)


def require(expr, msg):
    if not bool(expr):
        raise ValueError("MOIS_GLOBAL_RF_NPP_LGM: " + msg)


def sha256(data):
    if isinstance(data, (str, Path)):
        data = Path(data).read_bytes()
    return hashlib.sha256(data).hexdigest()


def read_signed_parquet(path, expected_outer, member):
    import pyarrow.parquet as pq
    require(sha256(path) == expected_outer, "source ZIP hash mismatch: " + str(path))
    with ZipFile(path) as z:
        require(member in z.namelist(), "missing exact parquet " + member)
        inner = z.read(member)
        receipt = json.loads(z.read("receipt.json"))
    require(
        sha256(inner) == receipt["output_SHA256"][member],
        "Parquet SHA mismatch: " + member,
    )
    frame = pq.read_table(io.BytesIO(inner)).to_pandas()
    require(len(frame) == 37104, "accepted-species source count != 37104")
    frame["key"] = frame["key"].astype(str)
    require(frame.key.is_unique, "duplicate accepted GBIF species key")
    return frame, receipt


def prepare(main, env):
    for d in (main, env):
        require("key" in d.columns and d.key.is_unique, "accepted species IDs not unique")
        require(set(d.source_hue9) == set(COLORS), "wrong original 9-colour labels")
    require(set(main.key) == set(env.key), "main and BIO19 do not have same species keys")
    right = main[[
        "key", "source_hue9", "genus", "hue_status",
        NPP, NPP_COVER, LGM_T, LGM_P, "lgm_climate_coverage",
    ]].copy()
    d = env.merge(right, on="key", how="inner", validate="one_to_one",
                  suffixes=("", "_from_NPP_LGM"))
    require(len(d) == 37104, "join lost accepted species")
    require(
        d.source_hue9.equals(d.source_hue9_from_NPP_LGM)
        and d.genus.astype(str).equals(d.genus_from_NPP_LGM.astype(str)),
        "original photographic labels or genus differ between signed sources",
    )
    d = d.drop(columns=["source_hue9_from_NPP_LGM", "genus_from_NPP_LGM"])
    d["log_GBIF_occurrences"] = np.log1p(d.records_passed_spatial_QA.astype(float))
    d["log_ESACCI_AGB_woody_Mg_ha_2020"] = np.log1p(
        d.ESA_AGB_woody_Mg_ha_2020.astype(float).clip(lower=0)
    )
    require(
        set(STAGES["C_plus_LGM_3GCM"]).issubset(d.columns),
        "a contemporaneous, productivity, or ancient feature is unavailable",
    )
    complete = [
        "BIO01", "BIO19", NPP, NPP_COVER, LGM_T, LGM_P,
        "lgm_climate_coverage", "MODIS_NDVI_2001_2021_mean",
        "ESA_AGB_woody_Mg_ha_2020", "mean_abs_latitude_over_occupied_cells",
        "records_passed_spatial_QA", "occupied_grid_cells", "hue_status",
    ]
    eligible = d[complete].notna().all(axis=1)
    eligible &= d[NPP_COVER].ge(.75) & d.lgm_climate_coverage.ge(.75)
    eligible &= d.ESA_AGB_woody_Mg_ha_2020.ge(0)
    eligible &= d.BIO12.ge(0) & d.occupied_grid_cells.gt(0)
    d = d.loc[eligible].copy().reset_index(drop=True)
    require(len(d) > 30000 and d.genus.nunique() > 5000,
            "incomplete global matching support")
    require(d.key.is_unique and d.source_hue9.isin(COLORS).all(),
            "invalid accepted species sample")
    return d


def prepare_folds(d):
    y9 = d.source_hue9.to_numpy(str)
    genera = d.genus.to_numpy(str)
    splitter = StratifiedGroupKFold(n_splits=FOLDS, shuffle=True, random_state=SEED)
    folds = list(splitter.split(np.zeros(len(d)), y9, genera))
    all_test = []
    for train, test in folds:
        require(
            set(genera[train]).isdisjoint(set(genera[test])),
            "same genus appears in train and test",
        )
        require(set(y9[train]) == set(COLORS), "training fold missing a photo colour")
        all_test.extend(test.tolist())
    require(sorted(all_test) == list(range(len(d))), "test species are not disjoint")
    return folds


def fit_fold(kind, X, y, train, test, seed, multiclass=False, quick=False):
    imp = SimpleImputer(strategy="median", keep_empty_features=True)
    xt = imp.fit_transform(X[train])
    xv = imp.transform(X[test])
    require(np.isfinite(xt).all() and np.isfinite(xv).all(),
            "training-only median imputation produced nonfinite values")
    if kind == "RF":
        model = RandomForestClassifier(
            n_estimators=9 if quick else 125,
            min_samples_leaf=5 if quick else 30,
            class_weight="balanced_subsample",
            max_samples=.82, bootstrap=True, n_jobs=2, random_state=seed,
        )
    elif kind == "HGB" and not multiclass:
        model = HistGradientBoostingClassifier(
            max_iter=12 if quick else 150,
            max_leaf_nodes=15, min_samples_leaf=5 if quick else 90,
            learning_rate=.05, l2_regularization=3, random_state=seed,
        )
    else:
        raise ValueError("unrecognized model task")
    model.fit(xt, y[train])
    pr = model.predict_proba(xv)
    if multiclass:
        out = np.zeros((len(test), len(COLORS)), dtype=float)
        for j, hue in enumerate(COLORS):
            hit = np.flatnonzero(model.classes_ == hue)
            if hit.size:
                out[:, j] = pr[:, hit[0]]
        require((out.sum(axis=1) > .999).all(), "nine-hue probabilities not normalized")
        return out
    require(0 in model.classes_ and 1 in model.classes_, "binary training label lost")
    return pr[:, np.flatnonzero(model.classes_ == 1)[0]]


def loss_binary(y, p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(y * np.log(p) + (1 - y) * np.log1p(-p))


def loss_nine(y, probs):
    ids = np.asarray([COLORS.index(k) for k in y], dtype=int)
    return -np.log(np.clip(probs[np.arange(len(y)), ids], 1e-12, 1.0))


def perf_binary(y, p):
    return {
        "ROC_AUC": float(roc_auc_score(y, p)),
        "avg_precision": float(average_precision_score(y, p)),
        "balanced_accuracy_at_0p5": float(balanced_accuracy_score(y, p >= .5)),
        "log_loss": float(np.mean(loss_binary(y, p))),
        "Brier": float(brier_score_loss(y, p)),
    }


def perf_nine(y, p):
    guess = np.asarray(COLORS)[np.argmax(p, axis=1)]
    return {
        "balanced_accuracy": float(balanced_accuracy_score(y, guess)),
        "top1_accuracy": float(np.mean(guess == y)),
        "log_loss": float(np.mean(loss_nine(y, p))),
    }


def grouped_bootstrap(y, genera, results, kind, outcome, repeats):
    # Genus cluster bootstrap on exactly the same held-out species and predictions.
    # Difference = new stage minus old stage; positive AUC or negative log loss is better.
    ids, group = np.unique(genera, return_inverse=True)
    rng = np.random.default_rng(SEED + (0 if kind == "RF" else 100) +
                                (0 if outcome == "white" else 200))
    rows = []
    for old, new in COMPARISONS:
        a = results[(outcome, kind, old)]
        b = results[(outcome, kind, new)]
        if outcome == "white":
            la, lb = loss_binary(y, a), loss_binary(y, b)
            observed_auc_diff = float(roc_auc_score(y, b) - roc_auc_score(y, a))
        else:
            la, lb = loss_nine(y, a), loss_nine(y, b)
            observed_auc_diff = np.nan
        observed_nll_diff = float(np.mean(lb - la))
        ci_auc, ci_nll = [], []
        for _ in range(repeats):
            group_counts = np.bincount(
                rng.integers(0, len(ids), size=len(ids)), minlength=len(ids)
            )
            weights = group_counts[group].astype(float)
            ci_nll.append(float(np.average(lb - la, weights=weights)))
            if outcome == "white" and len(np.unique(y[weights > 0])) == 2:
                ci_auc.append(float(
                    roc_auc_score(y, b, sample_weight=weights)
                    - roc_auc_score(y, a, sample_weight=weights)
                ))
        rows.append({
            "kind": kind, "outcome": outcome, "baseline": old, "added": new,
            "new_minus_old_AUC": observed_auc_diff,
            "delta_AUC_genus_bootstrap_CI95": (
                list(map(float, np.percentile(ci_auc, [2.5, 97.5])))
                if ci_auc else None
            ),
            "new_minus_old_logloss": observed_nll_diff,
            "delta_logloss_genus_bootstrap_CI95":
                list(map(float, np.percentile(ci_nll, [2.5, 97.5]))),
            "n_genera": int(len(ids)), "resamples": int(repeats),
        })
    return rows


def evaluate(d, dest, quick=False):
    dest = Path(dest)
    require(not dest.exists(), "do not overwrite prior analysis output")
    dest.mkdir(parents=True)
    folds = prepare_folds(d)
    n = len(d)
    y = d.source_hue9.to_numpy(str)
    white = (y == "white").astype(int)
    kinds = [("white", "RF"), ("white", "HGB"), ("nine", "RF")]
    preds = {}
    fold_assign = np.full(n, -1, dtype=int)
    fold_scores = []
    for fold, (train, test) in enumerate(folds):
        fold_assign[test] = fold
        for outcome, kind in kinds:
            target = white if outcome == "white" else y
            for stage, features in STAGES.items():
                X = d[features].to_numpy(float)
                pp = fit_fold(kind, X, target, train, test,
                              seed=SEED + fold, multiclass=(outcome == "nine"),
                              quick=quick)
                key = (outcome, kind, stage)
                if key not in preds:
                    preds[key] = np.full((n, N_ROI_COLORS), np.nan) if outcome == "nine" else np.full(n, np.nan)
                preds[key][test] = pp
                perf = (perf_binary if outcome == "white" else perf_nine)(target[test], pp)
                fold_scores.append({
                    "fold": int(fold), "test_n": int(len(test)),
                    "train_n": int(len(train)), "outcome": outcome,
                    "kind": kind, "stage": stage,
                    **{k: float(v) for k, v in perf.items()}
                })
                print(f"MOIS_RF_GLOBAL fold={fold} model={kind}/{outcome} "
                      f"stage={stage} n_test={len(test)} "
                      f"ll={perf['log_loss']:.6f}", flush=True)
    require((fold_assign >= 0).all(), "fold coverage incomplete")
    scores = []
    for (outcome, kind, stage), pp in preds.items():
        require(np.isfinite(pp).all(), f"nonfinite out-of-fold prediction {outcome} {kind} {stage}")
        perf = (perf_binary if outcome == "white" else perf_nine)(
            white if outcome == "white" else y, pp
        )
        scores.append({
            "outcome": outcome, "kind": kind, "stage": stage, "n_species": n,
            "n_genera": int(d.genus.nunique()),
            "n_predictors": len(STAGES[stage]),
            **perf,
        })
    nboot = 25 if quick else 500
    comparisons = []
    for outcome, kind in kinds:
        comparisons.extend(grouped_bootstrap(
            white if outcome == "white" else y,
            d.genus.to_numpy(str), preds, kind, outcome, nboot
        ))
    # Diagnostic only: this stratum is colour-biased, not an unbiased gold standard.
    diagnostics = []
    for label, mask in [
        ("all", np.ones(n, dtype=bool)),
        ("lt30_absolute_lat", d.mean_abs_latitude_over_occupied_cells.to_numpy() < 30),
        ("ge30_absolute_lat", d.mean_abs_latitude_over_occupied_cells.to_numpy() >= 30),
        ("photo_sampling_consensus_single",
         d.hue_status.to_numpy(str) == "sampling_consensus_single"),
    ]:
        if mask.sum() < 90 or len(np.unique(white[mask])) != 2:
            continue
        for stage in STAGES:
            pp = preds[("white", "RF", stage)][mask]
            pf = perf_binary(white[mask], pp)
            diagnostics.append({
                "stratum": label, "stage": stage, "n_species": int(mask.sum()),
                "n_white": int(white[mask].sum()),
                **pf,
            })
    pd.DataFrame(scores).to_csv(dest / "nested_RF_HGB_pooled_OOF_performance.csv", index=False)
    pd.DataFrame(fold_scores).to_csv(dest / "nested_RF_HGB_genus_disjoint_fold_metrics.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(dest / "nested_RF_photo_quality_latitude_diagnostics.csv", index=False)
    (dest / "paired_genus_bootstrap_increments.json").write_text(
        json.dumps(comparisons, indent=2, allow_nan=False) + "\n"
    )
    # Species-level predictions are in the signed Actions artifact only;
    # the GitHub repository receives no images or location-level individual data.
    ptab = pd.DataFrame({
        "specieskey": d.key.astype(str), "genus": d.genus.astype(str),
        "photo_hue9_proxy": y, "photo_evidence_status": d.hue_status.astype(str),
        "heldout_genus_fold": fold_assign,
        "mean_abs_latitude_over_present_GBIF_range":
            d.mean_abs_latitude_over_occupied_cells,
    })
    for (outcome, kind, stage), pp in preds.items():
        if outcome == "white":
            ptab[f"p_white_{kind}_{stage}"] = pp
        else:
            for j, hue in enumerate(COLORS):
                ptab[f"p_{hue}_RF9_{stage}"] = pp[:, j]
    ptab.to_csv(dest / "all_species_out_of_fold_predictions.csv.gz",
                index=False, compression="gzip")
    source = {
        "n_global_source_accepted_species": 37104,
        "n_eligible_same_species": n,
        "n_genera": int(d.genus.nunique()),
        "n_high_confidence_photo_status":
            int((d.hue_status == "sampling_consensus_single").sum()),
        "same_genus_disjoint_3fold_and_species_all_stages": True,
        "no_photo_status_used_as_predictor": True,
        "WorldClim21_current_period": "1970-2000",
        "WorldClim14_LGM_deltas_paired_against": "1960-1990",
        "LGM_3GCM_at_current_GBIF_ranges_not_reconstructed_LGM_occurrence": True,
        "NPP_is_2001_2023_annual_carbon_flux_not_flower_investment": True,
        "AGB_is_woody_vegetation_not_focal_flower_biomass": True,
        "no_fossil_or_pollen_validation": True,
        "photographic_colour_is_not_independently_verified_petal_truth": True,
        "input_zip_SHA256": {"current_BIO19_woody_AGB": ENV_SHA,
                             "NPP_LGM_GBIF_species": MAIN_SHA},
        "GBIF_recorded_not_true_absence": True,
        "model": {"RF_trees": 9 if quick else 125, "min_samples_leaf": 5 if quick else 30,
                  "HGB_max_iter": 12 if quick else 150},
        "paired_bootstrap_cluster": "accepted genus",
        "bootstrap_repetitions": nboot,
        "n_classes": len(COLORS),
        "files_SHA256": {
            p.name: sha256(p) for p in sorted(dest.iterdir())
            if p.is_file()
        },
    }
    (dest / "receipt.json").write_text(json.dumps(source, indent=2) + "\n")
    print("MOIS_GLOBAL_NESTED_RF_DONE", json.dumps({
        "n_species": n, "n_genera": int(d.genus.nunique()),
        "scores": scores, "paired_delta": comparisons,
        "photo_diagnostics": diagnostics,
    }, allow_nan=False), flush=True)
    return scores, comparisons


def synthetic_test():
    # Synthetic arithmetic and fold integrity smoke test; never scientific data.
    rng = np.random.default_rng(17)
    d = pd.DataFrame({
        "key": [str(i) for i in range(810)],
        "genus": [f"g{i//3:04d}" for i in range(810)],
        "source_hue9": np.array(COLORS)[np.arange(810) % 9],
        "hue_status": np.where(np.arange(810) % 5 == 0, "sampling_consensus_single",
                               "provisional_single"),
        "mean_abs_latitude_over_occupied_cells": rng.uniform(0, 60, 810),
    })
    for c in set(STAGES["C_plus_LGM_3GCM"]):
        if c not in d.columns:
            d[c] = rng.normal(size=810)
    import tempfile
    with tempfile.TemporaryDirectory() as path:
        out = Path(path) / "out"
        scores, comparisons = evaluate(d, out, quick=True)
        require(len(scores) == 9 and len(comparisons) == 6,
                "expected nine pooled scores and six bootstrap comparisons")
        require((out / "receipt.json").exists(), "missing signed output receipt")
    print("MOIS_GLOBAL_NESTED_RF_SYNTHETIC_PASS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--main")
    ap.add_argument("--env")
    ap.add_argument("--output")
    ap.add_argument("--selftest", action="store_true")
    x = ap.parse_args()
    if x.selftest:
        synthetic_test()
    else:
        require(x.main and x.env and x.output,
                "--main, --env, --output are required for real analysis")
        main, mr = read_signed_parquet(x.main, MAIN_SHA, MAIN_FILE)
        env, er = read_signed_parquet(x.env, ENV_SHA, ENV_FILE)
        require(mr["source_37k_species"] == er["source_global_accepted_species"] == 37104,
                "37,104 accepted key count changed")
        d = prepare(main, env)
        evaluate(d, x.output)
