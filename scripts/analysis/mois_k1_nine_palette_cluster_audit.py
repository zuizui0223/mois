#!/usr/bin/env python3
"""Source-locked exploratory nine-palette clustering of REAL MOIS K1 ROI evidence.

The nine frozen flower_fraction_* channels preserve much more information than
the legacy four discrete hues. The *original four classes remain unchanged*;
an unsupervised palette archetype is never used to reclassify failed ROI/photo
acquisition as white or pigmented.

Only source-classified photos enter model fitting. A genuinely positive and
finite nine-palette mixture from the explicit "ambiguous palette composition"
status can receive an EXPLORATORY cluster projection, not an original hue
label. ROI/flip-gate failures cannot be "rescued" by relabelling their pixels.

Distances use sqrt(normalized 9-palette proportions), a Hellinger embedding
for compositional colour. Cluster count is selected from prespecified k based
on holdout silhouette + seed stability. No geographic, species or climate
features enter clustering; result is descriptive, not causality or native
flower prevalence. This is NOT image CIELAB, which needs original ROI pixels.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

PALETTE = ("white","yellow","orange","red","pink","magenta",
           "purple","blue","bronze")
FRACTIONS = tuple("flower_fraction_"+name for name in PALETTE)
HUES = ("white","yellow_orange","red_pink","blue_purple")
ACQ_FAIL = "image_acquisition_failed"
CLASSIFIED = "classified_four_state_morph"
AMBIGUOUS = "not_evaluable_ambiguous_palette_composition"
SCHEMA = "mois_k1_true_original_roi9_unclassified_reason_and_cluster_pilot_v1"
K1_SOURCE_TAG = "mois-full93m-tier06-balanced-onephoto-K1-new-original-warm8-RGB-v1"
EXPECTED_SOURCE_SHA = "0d3de3374a06b1e7b9ad0281e2ab7844632ed139db871488524061191c8f34c1"
EXPECTED_RECEIPT_SHA = "e86463bc61f9eebd8cbdbe01269f767c140bf243f9e374027008e9f095eb9044"
EXPECTED_ATTEMPTS = 76753
EXPECTED_CLASSIFIED = 27689
EXPECTED_ACQUISITION_FAILED = 203
EXPECTED_ACQUIRED_UNKNOWN = 48861
SEED = 20261025


def require(ok, message):
    if not ok:
        raise ValueError("MOIS_K1_9PALETTE: "+message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def stable_sample_indices(ids, n, seed=SEED):
    """Photo-ID-only deterministic selection: never pick training IDs by colour."""
    import numpy as np
    ids=list(map(str, ids))
    keys=[(hashlib.sha256(f"{seed}|{pid}".encode()).digest()[:12],i)
          for i,pid in enumerate(ids)]
    keys.sort(key=lambda x:x[0])
    return np.asarray([i for _,i in keys[:min(n,len(ids))]],dtype=int)


def nine_palette_from_frame(df):
    import numpy as np
    import pandas as pd
    require(set(FRACTIONS) <= set(df.columns),
            "Source signed measured photos do not retain all nine palette channels")
    raw=np.column_stack([
        pd.to_numeric(df[col],errors="coerce").to_numpy(float)
        for col in FRACTIONS])
    finite=np.isfinite(raw).all(axis=1)
    nonnegative=(raw>=0).all(axis=1)
    finite_bound=(raw<=1+1e-6).all(axis=1)
    total=np.where(finite[:,None],raw,0).sum(axis=1)
    valid=finite & nonnegative & finite_bound & (total>0)
    mixture=np.full(raw.shape,np.nan,dtype=float)
    mixture[valid]=raw[valid]/total[valid,None]
    return raw,valid,mixture


def cluster_sqrt_palette(known_palette,photo_ids,
                         k_values=(3,4,5,6,7,8),
                         train_limit=16000,score_limit=1200,
                         seed=SEED):
    """Select k from compact deterministic subsets; no photo outcome rerun."""
    import numpy as np
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.metrics import adjusted_rand_score, silhouette_score

    P=np.asarray(known_palette,float)
    require(P.ndim==2 and P.shape[1]==len(PALETTE) and len(P)==len(photo_ids),
            "Expected one nine-palette vector per source-classified photo")
    require(np.isfinite(P).all() and (P>=0).all()
            and np.allclose(P.sum(axis=1),1,rtol=0,atol=1e-6),
            "Nine-palette composition must sum to one and remain finite")
    require(len(P)>=30,"Too few real nine-palette classified photos for a cluster pilot")
    k_values=tuple(k_values)
    require(k_values and len(set(k_values))==len(k_values)
            and all(type(k) is int and 2<=k<=32 for k in k_values),
            "Invalid prespecified cluster search")
    require(max(k_values)<min(train_limit,len(P)),
            "Cluster candidate exceeds real original photo training capacity")
    X=np.sqrt(P)
    chosen=stable_sample_indices(photo_ids,train_limit,seed)
    score_i=stable_sample_indices(photo_ids,score_limit,seed+1)
    X_train,X_score=X[chosen],X[score_i]
    require(len(X_score)>=max(k_values)*4,
            "Insufficient signed photo holdout for clustering stability")
    diagnostics=[]
    for k in k_values:
        def estimate(random_state):
            est=MiniBatchKMeans(n_clusters=k,random_state=random_state,
                                n_init=3,max_iter=120,batch_size=1024,
                                reassignment_ratio=0.01)
            est.fit(X_train)
            return est
        a,b=estimate(seed+k),estimate(seed+100+k)
        lab_a,lab_b=a.predict(X_score),b.predict(X_score)
        required_labels=len(set(lab_a))
        require(required_labels>=2,"Degenerate cluster result with fewer than two observed clusters")
        silhouette=float(silhouette_score(X_score,lab_a,metric="euclidean",
                                          sample_size=None))
        ari=float(adjusted_rand_score(lab_a,lab_b))
        require(math.isfinite(silhouette) and math.isfinite(ari),
                "Nonfinite original photo palette cluster diagnostic")
        # Penalize seeds whose assignments are not repeatable on the SAME fixed
        # photo-ID-only validation set. This is an exploratory, not inferential,
        # data-driven k; all candidate scores are retained for inspection.
        selection=silhouette-0.15*(1-max(0,ari))
        diagnostics.append(dict(k=k,silhouette=silhouette,seed_ari=ari,
                                penalized_score=selection,
                                n_train=len(X_train),n_holdout=len(X_score),
                                actual_occupied_clusters=required_labels))
    chosen_k=max(diagnostics,key=lambda d:(d["penalized_score"],-d["k"]))["k"]
    model=MiniBatchKMeans(n_clusters=chosen_k,random_state=seed+chosen_k,
                          n_init=3,max_iter=120,batch_size=1024,
                          reassignment_ratio=0.01).fit(X_train)
    centroid_X=model.cluster_centers_
    center=np.square(np.maximum(0,centroid_X))
    center=center/center.sum(axis=1,keepdims=True)
    # Deterministic canonical cluster names; original fitting indices remain
    # available separately for exact reproducibility.
    order=sorted(range(chosen_k),
                 key=lambda j:tuple(-round(float(x),12) for x in center[j]))
    lut={source:i+1 for i,source in enumerate(order)}
    centroids=[]
    for source in order:
        p=center[source]
        row=dict(cluster_id=f"C{lut[source]:02d}",
                 original_kmeans_cluster_index=int(source),
                 dominant_palette=PALETTE[int(np.argmax(p))])
        for n,value in zip(PALETTE,p):
            row[f"palette_{n}"]=float(value)
        centroids.append(row)
    # Return a safe projection method implemented with squared distances to
    # the frozen centers; sklearn .predict() is deterministic here.
    def assign(P_external):
        x=np.asarray(P_external,float)
        require(x.ndim==2 and x.shape[1]==len(PALETTE)
                and np.isfinite(x).all() and (x>=0).all()
                and np.allclose(x.sum(axis=1),1,atol=1e-6,rtol=0),
                "Cannot assign invalid or technical-unknown colour fractions")
        if len(x)==0:return []
        raw=model.predict(np.sqrt(x))
        return [f"C{lut[int(v)]:02d}" for v in raw]
    return dict(k=chosen_k,diagnostics=diagnostics,centroids=centroids,
                assign=assign,training_photo_hash=hashlib.sha256("\n".join(
                    map(str,np.asarray(photo_ids,dtype=object)[chosen])).encode()).hexdigest(),
                n_training_photos=len(X_train),
                model_fit_only_source_classified=True)


def run(source, receipt, out, *,
        expected_source_sha=EXPECTED_SOURCE_SHA,
        expected_receipt_sha=EXPECTED_RECEIPT_SHA,
        expected_attempts=EXPECTED_ATTEMPTS,
        expected_classified=EXPECTED_CLASSIFIED,
        expected_failed=EXPECTED_ACQUISITION_FAILED,
        expected_other=EXPECTED_ACQUIRED_UNKNOWN,
        k_values=(3,4,5,6,7,8),
        train_limit=16000,score_limit=1200):
    import numpy as np
    import pandas as pd
    source,receipt,out=Path(source),Path(receipt),Path(out)
    require(not out.exists(),"Refuse to overwrite a prior source-locked palette result")
    require(sha(source)==expected_source_sha and sha(receipt)==expected_receipt_sha,
            "Frozen original measured-photo parquet/receipt byte SHA changed")
    proof=json.loads(receipt.read_text(encoding="utf-8"))
    require(proof.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
            and proof.get("measured_photo_parquet_sha256")==expected_source_sha
            and proof.get("attempted")==expected_attempts
            and proof.get("classifiable")==expected_classified
            and proof.get("acquisition_failed")==expected_failed
            and proof.get("unclassifiable_acquired")==expected_other
            and proof.get("technical_unknowns_not_assumed_nonwhite") is True
            and proof.get("raw_licensed_photo_pixels_redistributed") is False,
            "Original signed FCP photo/ROI totals and unknown status not proved")
    # load only canonical actual-original ROI columns plus nine frozen colours
    df=pd.read_parquet(source)
    require({"photo_id","observation_id","inat_taxon_id","morph",
             "measurement_status","image_sha256"} | set(FRACTIONS) <= set(df.columns),
            "Signed K1 photo source lacks source ID/hue/evaluable features")
    pid=df.photo_id.fillna("").astype(str)
    obs=df.observation_id.fillna("").astype(str)
    states=df.measurement_status.fillna("").astype(str)
    morph=df.morph.fillna("").astype(str)
    valid_id=pid.str.fullmatch(r"[0-9]+")
    require(len(df)==expected_attempts
            and valid_id.all() and pid.is_unique
            and obs.ne("").all() and
            df.inat_taxon_id.fillna("").astype(str).str.fullmatch(r"[0-9]+").all(),
            "Actual image attempt/observation/taxon identities missing")
    known=morph.isin(HUES).to_numpy(bool)
    fail=states.eq(ACQ_FAIL).to_numpy(bool)
    require(np.array_equal(known,states.eq(CLASSIFIED).to_numpy(bool))
            and morph.isin((*HUES,"mixed_uncertain")).all()
            and known.sum()==expected_classified
            and fail.sum()==expected_failed
            and (~known & ~fail).sum()==expected_other,
            "K1 source four-hue status does not agree with real-photo receipts")
    image_sha=df.image_sha256.fillna("").astype(str)
    require(image_sha[~fail].str.fullmatch(r"[0-9a-f]{64}").all(),
            "Acquired image missing original source-pixel SHA")
    raw, palette_ok, palette=nine_palette_from_frame(df)
    statuses=states.to_numpy(str)
    amb=(statuses==AMBIGUOUS) & ~known & ~fail
    valid_train=known & palette_ok
    valid_amb=amb & palette_ok
    require(int(valid_train.sum())>=min(30,expected_classified),
            "Source-known original photo ROI has no usable nine-colour vectors")
    require(not np.any(palette_ok & fail),"Acquisition failure cannot carry real flower palette mass")
    # An unknown's positive palette does not imply a trustworthy flower ROI.
    # Only explicit ambiguity among the ORIGINAL terminal statuses is projected.
    candidate_nonclass=(~known & ~fail & palette_ok)
    cluster=cluster_sqrt_palette(palette[valid_train],
                                 pid.to_numpy(str)[valid_train],
                                 k_values=k_values,train_limit=train_limit,
                                 score_limit=score_limit)
    assignments=np.full(len(df),"",dtype=object)
    assignments[valid_train]=cluster["assign"](palette[valid_train])
    assignments[valid_amb]=cluster["assign"](palette[valid_amb])
    cluster_context=np.full(len(df),"no_evaluable_flower_colour_ROI",dtype=object)
    cluster_context[valid_train]="observed_original_four_hue_classified"
    cluster_context[valid_amb]="exploratory_palette_ambiguous_unclassified"
    cluster_context[known & ~palette_ok]="original_hue_known_but_palette_features_missing"
    require(np.all(assignments[~(valid_train | valid_amb)]==""),
            "Technically unknown/failed flower images were assigned cluster labels")
    counts=Counter(map(str,statuses))
    status_table=[]
    for status,n in sorted(counts.items()):
        flag=statuses==status
        status_table.append(dict(measurement_status=status,
            n_photo_attempts=n,
            source_classified_four_hue=int(known[flag].sum()),
            image_acquisition_failed=int(fail[flag].sum()),
            finite_nonzero_nine_palette_vectors=int(palette_ok[flag].sum()),
            exploratory_unclassified_palette_candidate=int(valid_amb[flag].sum()),
            positive_palette_but_failed_other_quality_gate=int(
                (palette_ok & ~known & ~fail & ~amb & flag).sum()),
            original_source_classification_unchanged=True))
    labelled=df[["photo_id","observation_id","inat_taxon_id","morph",
                 "measurement_status"]].copy()
    labelled["palette9_evaluable_numeric_vector"]=palette_ok
    labelled["nine_palette_exploratory_cluster"]=assignments
    labelled["cluster_assignment_provenance"]=cluster_context
    labelled["source_four_hue_classified_unchanged"]=known
    for j,color in enumerate(PALETTE):
        labelled["continuous9_"+color]=palette[:,j]
    out.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-nine-palette-source-",dir=out.parent))
    try:
        pd.DataFrame(status_table).to_csv(stage/"genuine_original_photo_missingness_by_status.csv",
                                           index=False,lineterminator="\n")
        pd.DataFrame(cluster["diagnostics"]).to_csv(stage/"cluster_k_stability_and_silhouette.csv",
                                                    index=False,lineterminator="\n")
        pd.DataFrame(cluster["centroids"]).to_csv(stage/"nine_palette_cluster_archetypes.csv",
                                                  index=False,lineterminator="\n")
        labelled.to_parquet(stage/"original_K1_photo_9palette_and_exploratory_clusters.parquet",
                            index=False,compression="zstd")
        selected=cluster["k"]
        cluster_counts=Counter(assignments[valid_train])
        report={
            "schema":SCHEMA,
            "source_original_photo_RGB_release":K1_SOURCE_TAG,
            "original_photo_parquet_sha256":expected_source_sha,
            "original_original_photo_receipt_sha256":expected_receipt_sha,
            "source_real_original_photo_attempts":len(df),
            "source_four_hue_original_classifications_UNCHANGED":int(known.sum()),
            "source_image_acquisition_failures_UNCHANGED":int(fail.sum()),
            "source_image_acquired_but_unclassifiable_UNCHANGED":int((~known & ~fail).sum()),
            "real_photo_9palette_finite_nonzero_vectors":int(palette_ok.sum()),
            "original_four_hue_classified_WITH_nine_palette_features":int(valid_train.sum()),
            "source_classified_missing_nine_palette_features":int((known & ~palette_ok).sum()),
            "source_ambiguous_palette_composition_total":int(amb.sum()),
            "source_ambiguous_palette_composition_with_usable_features":int(valid_amb.sum()),
            "all_unclassified_acquired_with_positive_palette_vectors":int(candidate_nonclass.sum()),
            "other_nonclassifiable_with_9palette_NOT_treated_as_flowers":int(
                (candidate_nonclass & ~amb).sum()),
            "cluster_fit_source_status":CLASSIFIED,
            "model_fitted_ONLY_on_source_four_hue_classified_original_ROIs":True,
            "cluster_assignments_on_known_original_photo_rois":int(valid_train.sum()),
            "exploratory_ambiguous_unknown_palette_clusters_NOT_original_hue_labels":int(valid_amb.sum()),
            "source_technical_unknown_photos_promoted_to_four_hue":0,
            "number_clusters_selected_exploratory":selected,
            "cluster_selection_objective":"holdout_Hellinger_silhouette_minus_0.15_times_seed_instability",
            "cluster_n_train_blind_photo_hash_sample":cluster["n_training_photos"],
            "cluster_training_photo_ID_only_sample_SHA256":cluster["training_photo_hash"],
            "source_classified_cluster_sizes":dict(sorted(cluster_counts.items())),
            "cluster_photo_fraction_features":list(FRACTIONS),
            "latent_colour_features":"nine_frozen_ROI_biological_palette_fractions_NOT_CIELAB",
            "cluster_reuses_original_photo_SHA_and_does_not_reinfer":True,
            "no_new_source_photos_downloaded":True,
            "no_new_original_photo_ROI_evaluation_or_altered_frozen_palette":True,
            "not_botanically_verified_blooming_petal_wild_native_samples":True,
            "no_globally_unbiased_white_flower_rate_estimated":True,
            "output_sha256":{p.name:sha(p) for p in sorted(stage.iterdir()) if p.is_file()},
        }
        (stage/"nine_palette_original_RGB_AUDIT_receipt.json").write_text(
            json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
        os.replace(stage,out)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_K1_NINE_PALETTE_REAL_SOURCE_DIAGNOSTIC",json.dumps({
        "original_attempts":len(df),
        "classified_original_four_hues":int(known.sum()),
        "downloaded_but_unclassifiable":int((~known & ~fail).sum()),
        "palette_valid_known":int(valid_train.sum()),
        "ambiguity_with_evaluable_palette_features":int(valid_amb.sum()),
        "other_failed_quality_gates_with_positive_palette":int(
            (candidate_nonclass & ~amb).sum()),
        "automatic_cluster_archetypes":selected,
        "original_four_hue_labels_modified":0,
    },sort_keys=True),flush=True)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",required=True,type=Path)
    parser.add_argument("--receipt",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    run(args.source,args.receipt,args.output)


if __name__=="__main__":
    main()
