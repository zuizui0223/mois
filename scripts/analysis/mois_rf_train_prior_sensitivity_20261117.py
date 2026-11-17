#!/usr/bin/env python3
"""Train-fold colour-prior sensitivity for 2026-11-17 MOIS global RF artifact.

Strict signed source, independent genus-fold validation, no heldout labels in
prior construction. This is a diagnostic reweighting, NOT learned calibration
or a verified botanical petal-colour prediction model.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

SHA = "8fc8d980aff1e80e36c5633f71c31d9e1d7405c1cfbb06b0e051e1b6591472b0"
COLOURS = ["white","yellow","orange","red","pink","magenta","purple","blue","bronze"]
STAGES = ["A_current_full_climate", "B_plus_NPP_NDVI_AGB", "C_plus_LGM_3GCM"]


def read_source(path):
    raw = Path(path).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == SHA, "Original numeric ZIP digest changed"
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        receipt = json.loads(z.read("receipt.json"))
        for name, digest in receipt["files_SHA256"].items():
            assert hashlib.sha256(z.read(name)).hexdigest() == digest, name
        d = pd.read_csv(
            io.BytesIO(z.read("all_species_out_of_fold_predictions.csv.gz")),
            compression="gzip", dtype={"specieskey":str}
        )
    assert len(d) == 36581 and d.specieskey.is_unique
    assert d.genus.nunique() == 6562
    assert d.groupby("genus").heldout_genus_fold.nunique().max() == 1
    return d


def analyze(d, bootstrap=2000, seed=20261118):
    labels = d.photo_hue9_proxy.to_numpy(str)
    binary = (labels == "white").astype(int)
    classid = np.array([COLOURS.index(x) for x in labels], dtype=int)
    folds = d.heldout_genus_fold.to_numpy(int)
    n = len(d)
    p0 = np.empty((n, 9), dtype=float)
    b0 = np.empty(n, dtype=float)
    adjusted9 = {s:np.zeros((n,9)) for s in STAGES}
    adjustedb = {s:np.zeros(n) for s in STAGES}
    for fold in np.unique(folds):
        test = folds == fold
        train = ~test
        prior = np.array([(labels[train]==x).mean() for x in COLOURS])
        assert np.isfinite(prior).all() and prior.min() > 0
        p0[test] = prior
        b0[test] = prior[0]
        for stage in STAGES:
            source9 = d[[f"p_{x}_RF9_{stage}" for x in COLOURS]].to_numpy()[test]
            reweighted = source9 * prior
            adjusted9[stage][test] = reweighted/reweighted.sum(axis=1,keepdims=True)
            sourceb = d.loc[test, f"p_white_RF_{stage}"].to_numpy()
            z = sourceb * prior[0]
            adjustedb[stage][test] = z/(z+(1-sourceb)*(1-prior[0]))
    logloss9 = lambda prob: -np.log(np.clip(prob[np.arange(n), classid], 1e-12, 1))
    loglossb = lambda prob: -(binary*np.log(np.clip(prob,1e-12,1)) +
                                   (1-binary)*np.log(np.clip(1-prob,1e-12,1)))
    baseline = {"nine":logloss9(p0), "white":loglossb(b0)}
    losses = {}
    report = []
    for stage in STAGES:
        probs9 = adjusted9[stage]
        probsb = adjustedb[stage]
        losses[("nine",stage)] = logloss9(probs9)
        losses[("white",stage)] = loglossb(probsb)
        raw9 = d[[f"p_{x}_RF9_{stage}" for x in COLOURS]].to_numpy()
        rawb = d[f"p_white_RF_{stage}"].to_numpy()
        report.append({
            "outcome":"nine", "stage":stage,
            "fold_prior_baseline_logloss":float(baseline["nine"].mean()),
            "raw_class_weighted_logloss":float(logloss9(raw9).mean()),
            "prior_reweighted_logloss":float(losses[("nine",stage)].mean()),
            "prior_reweighted_balanced_accuracy":float(balanced_accuracy_score(
                labels,np.array(COLOURS)[probs9.argmax(axis=1)])),
        })
        report.append({
            "outcome":"white","stage":stage,
            "fold_prior_baseline_logloss":float(baseline["white"].mean()),
            "raw_class_weighted_logloss":float(loglossb(rawb).mean()),
            "prior_reweighted_logloss":float(losses[("white",stage)].mean()),
            "prior_reweighted_balanced_accuracy":float(balanced_accuracy_score(
                binary,probsb>=.5)),
            "prior_reweighted_ROC_AUC":float(roc_auc_score(binary,probsb)),
        })
    unique, indices = np.unique(d.genus.to_numpy(str),return_inverse=True)
    gsize = np.bincount(indices)
    rng = np.random.default_rng(seed)
    draws = np.array([rng.integers(0,len(unique),size=len(unique))
                      for _ in range(bootstrap)],dtype=int)
    contrasts = []
    for outcome in ("white","nine"):
        for old,new in [
            ("fold_prior_baseline",STAGES[0]),(STAGES[0],STAGES[1]),
            (STAGES[1],STAGES[2])
        ]:
            dloss = losses[(outcome,new)] - (
                baseline[outcome] if old=="fold_prior_baseline" else losses[(outcome,old)]
            )
            gs = np.bincount(indices,weights=dloss)
            vals = gs[draws].sum(axis=1)/gsize[draws].sum(axis=1)
            contrasts.append({
                "outcome":outcome, "baseline":old,"new":new,
                "new_minus_old_logloss":float(dloss.mean()),
                "CI95":list(map(float,np.quantile(vals,[.025,.975]))),
                "resamples":int(bootstrap), "seed":int(seed),
            })
    return pd.DataFrame(report),contrasts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--artifact",required=True)
    p.add_argument("--output",required=True)
    a = p.parse_args()
    d = read_source(a.artifact)
    rows,contrasts = analyze(d)
    out = Path(a.output)
    if out.exists():
        raise ValueError("Refuse to overwrite existing independent audit")
    out.mkdir(parents=True)
    rows.to_csv(out/"train_fold_prior_sensitivity.csv",index=False)
    (out/"train_fold_prior_genus_bootstrap.json").write_text(
        json.dumps(contrasts,indent=2)+"\n"
    )
    receipt = {"original_actions_artifact_sha256":SHA,
               "photo_colour_is_unverified_proxy":True,
               "method":"train-fold-prior-only reweighting, NOT trained calibration",
               "n_species":int(len(d)),"n_genera":int(d.genus.nunique()),
               "n_bootstrap":2000}
    receipt["outputs_sha256"] = {
        p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()
    }
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(rows.to_string(index=False))
    print(json.dumps(contrasts,indent=2))


if __name__=="__main__":
    main()
