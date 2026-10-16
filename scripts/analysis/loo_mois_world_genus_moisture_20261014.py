#!/usr/bin/env python3
"""Retrospective genus influence sensitivity on the source-frozen world moisture slopes.

This is a leave-ONE-genus-out point-estimate diagnostic, not a new CI, null test,
lineage contrast or biological replication. Private source genus moments must not
be published or passed to anonymous reviewers without licence/privacy approval.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

PRIVATE_MOMENT_SHA256 = "77dff4acc5756610e463bec3e5877c711f1695d8984be0b675e86f7d8afb4aae"
AXES = ("wc_bio12", "wc_vapr_annual_mean")
AXIS_COLUMNS = {"wc_bio12": 4, "wc_vapr_annual_mean": 8}
HUES = ("white", "yellow_orange", "red_pink", "blue_purple")


def replay(moments:Path, original_ledger:Path)->dict:
    sha=hashlib.sha256(moments.read_bytes()).hexdigest()
    if sha != PRIVATE_MOMENT_SHA256:
        raise ValueError("Private genus sufficient statistics do not match the exact original source")
    with np.load(moments,allow_pickle=False) as f:
        cross=np.asarray(f["xtx_genus"],dtype=np.float64)
        rhs=np.asarray(f["xty_genus"],dtype=np.float64)
        counts=np.asarray(f["n_species_per_genus"],dtype=np.int64)
    assert cross.shape==(1934,9,9) and rhs.shape==(1934,9,4)
    assert counts.shape==(1934,) and int(counts.sum())==20694 and counts.min()>=3
    assert np.isfinite(cross).all() and np.isfinite(rhs).all()
    with original_ledger.open(newline="",encoding="utf8") as f:
        original=list(csv.DictReader(f))
    assert len(original)==1280

    all_cross=cross.sum(axis=0)
    all_rhs=rhs.sum(axis=0)
    output={}
    total_sign_changes=0
    for axis in AXES:
        cols=[0,1,2,AXIS_COLUMNS[axis]]
        base=np.linalg.solve(all_cross[np.ix_(cols,cols)],all_rhs[cols])[-1,:]*100
        loo_cross=all_cross[np.ix_(cols,cols)][None,:,:]-cross[:,cols,:][:,:,cols]
        loo_rhs=all_rhs[cols][None,:,:]-rhs[:,cols,:]
        by_genus=np.linalg.solve(loo_cross,loo_rhs)[:,-1,:]*100
        assert by_genus.shape==(1934,4) and np.isfinite(by_genus).all()
        assert np.allclose(by_genus.sum(axis=1),0,atol=1e-10)
        differences=by_genus-base[None,:]
        summary={}
        for i,hue in enumerate(HUES):
            matched=[row for row in original if (
                row["mode"],row["quality"],row["region"],row["adjustment"],
                row["weight"],row["anchor"],row["hue"])==
                ("climate","all","world","genus_geo","genus_equal",axis,hue)]
            assert len(matched)==1
            assert abs(base[i]-float(matched[0]["direction_pp_per_env_global_sd"]))<1e-10
            flips=int(np.count_nonzero(np.sign(by_genus[:,i])!=np.sign(base[i])))
            total_sign_changes+=flips
            summary[hue]={
                "all_genera_original_point_pp":float(base[i]),
                "leave_one_genus_out_min_pp":float(by_genus[:,i].min()),
                "leave_one_genus_out_max_pp":float(by_genus[:,i].max()),
                "largest_absolute_point_change_pp":float(np.max(np.abs(differences[:,i]))),
                "p95_absolute_point_change_pp":float(np.quantile(np.abs(differences[:,i]),.95)),
                "point_sign_changes_over_1934_genus_exclusions":flips,
            }
        assert summary["white"]["point_sign_changes_over_1934_genus_exclusions"]==0
        output[axis]=summary
    assert total_sign_changes==0
    return {
        "schema":"mois_world_genus_equal_predeclared_two_moisture_axis_loo_20261014_v1",
        "status":"PASS_1934_INDEPENDENT_SINGLE_GENUS_DELETIONS_ON_SOURCE_LOCKED_TWO_MOISTURE_AXES",
        "source_private_genus_sufficient_moments_sha256":PRIVATE_MOMENT_SHA256,
        "source_original_measurement_denominator_species":20694,
        "source_original_nominal_genera":1934,
        "n_single_genus_delete_refits_per_axis":1934,
        "n_focal_original_environment_axes":2,
        "n_photographed_hues_per_axis":4,
        "n_total_leave_one_genus_out_slope_values":1934*2*4,
        "all_eight_original_full_model_points_match_frozen_source_ledger":True,
        "all_eight_slope_signs_unchanged_in_all_single_genus_deletions":True,
        "new_statistical_significance_test":False,
        "confidence_intervals_reestimated":False,
        "licence_and_geoprivacy_review_required_before_releasing_source_moments":True,
        "private_moments_uploaded_to_public_repository":False,
        "no_source_photo_ids_species_genus_names_or_site_locations_exported":True,
        "interpretation":"Local single-genus deletion stability within this source and original regression scope; does not establish geography, strict-photo, phylogeny, MNAR, causal water selection, or study-wide inference robustness",
        "axis_hue_results":output,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--private-moments",type=Path,required=True)
    p.add_argument("--original-ledger",type=Path,required=True)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    result=replay(a.private_moments,a.original_ledger)
    if a.output:
        a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
    compact={axis:{h:{key:v[key] for key in ("leave_one_genus_out_min_pp","leave_one_genus_out_max_pp","point_sign_changes_over_1934_genus_exclusions")} for h,v in rows.items()}
             for axis,rows in result["axis_hue_results"].items()}
    print(json.dumps({"status":result["status"],"partial_result":compact},indent=2))


if __name__=="__main__":
    main()
