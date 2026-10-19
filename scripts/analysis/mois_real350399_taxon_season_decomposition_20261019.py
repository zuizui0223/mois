#!/usr/bin/env python3
"""350,399 actually measured MOIS photographs: exact source-season taxon mixture.

The outcomes are classifier ROI four-hue labels, not botanist-verified flowers.
All new taxa-prioritised arms are source-selected, not random global flora.
This analysis REUSES exactly verified frozen photo pixels; it never infers a
flower colour from an original photo ID, a taxon name, or an unknown ROI.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    assemble, summary, require, sha256,
)
from scripts.analysis.mois_full149653_original_photo_season_taxon_composition import (
    decompose_photographic_taxa,
)

TOTAL_ATTEMPTED=350399
REPORT="mois_350399_actual_photo_season_taxon_four_term_decomposition.json"
FIGURE="mois_350399_hemisphere_season_taxon_components.csv"


def photographic_calendar(d):
    z=d.loc[d.classified & d.date.notna() &
            d.latitude.between(-90,90) &
            d.longitude.between(-180,180) &
            d.inat_taxon_id.ne("")].copy()
    m=z.date.dt.month
    northern=z.latitude.ge(23.5)
    southern=z.latitude.le(-23.5)
    valid=(northern | southern) & (
        m.isin((6,7,8)) | m.isin((12,1,2))
    )
    z=z.loc[valid].copy()
    m=z.date.dt.month
    z["zone"]=np.where(z.latitude.ge(0),
                       "north_extratropical","south_extratropical")
    z["season"]=np.where(
        (z.latitude.ge(0) & m.isin((6,7,8))) |
        (z.latitude.lt(0) & m.isin((12,1,2))),
        "summer", "winter")
    z["season_year"]=z.date.dt.year+m.eq(12).astype(int)
    require(z.photo_id.is_unique and z.classified.all() and
            z.season.isin(("summer","winter")).all(),
            "The classified original-photo seasonal frame is not source consistent")
    return z


def observation_equal_white(z):
    """A photographic-observation sensitivity, not census flower abundance."""
    s=z[z.observation_id.ne("")].copy()
    if s.empty:
        return {"observations":0, "white_label_fraction_mean_over_observations":None}
    # Some iNaturalist observations have more than one actual source image:
    # collapse these rather than calling every photo an independent visit.
    o=s.groupby(["observation_id","zone","season","inat_taxon_id"],
                observed=True).white.mean().reset_index()
    require(len(o)==s.observation_id.nunique(),
            "One source observation has inconsistent taxon, zone or season")
    return {
        "observations":len(o),
        "photos_with_observation_id":len(s),
        "max_images_per_observation":int(
            s.groupby("observation_id",observed=True).size().max()),
        "white_label_fraction_mean_over_observations":float(o.white.mean()),
        "seasons":{
            zone:{
                season:{"observations":int(len(q)),
                        "observation_equal_white_fraction":float(q.white.mean()) if len(q) else None}
                for season,q in ((season, o[(o.zone==zone)&(o.season==season)])
                                 for season in ("summer","winter"))
            }
            for zone in ("north_extratropical","south_extratropical")
        },
    }


def estimate(d, *, seed=20261019, draws=399):
    z=photographic_calendar(d)
    report={"photographed_original_RGB_attempts":len(d),
            "classified_extratropical_dated_summer_winter_photos":len(z),
            "hemispheres":{},
            "source_arm_seasonal_photo_counts":{}}
    for zone in ("north_extratropical","south_extratropical"):
        sub=z[z.zone.eq(zone)]
        result=decompose_photographic_taxa(sub,draws=draws,seed=seed)
        report["hemispheres"][zone]=result
    report["observation_equal_sensitivity"]=observation_equal_white(z)
    for name, group in z.groupby("arm",observed=True):
        report["source_arm_seasonal_photo_counts"][name]={
            zone:{
                ss: {"classifiable_dated_photo_ROIs":len(u),
                     "white_fraction":float(u.white.mean()) if len(u) else None,
                     "taxa":int(u.inat_taxon_id.nunique())}
                for ss,u in ((ss, group[(group.zone==zone)&(group.season==ss)])
                             for ss in ("summer","winter"))
            }
            for zone in ("north_extratropical","south_extratropical")
        }
    return report


def run(core,broad,wave2,cached8192,tail16384,output):
    out=Path(output)
    require(not out.exists() or not any(out.iterdir()),
            "Do not overwrite real ecological photographic source results")
    all_photo,arms=assemble(core,broad,wave2,cached8192,tail16384,strict=True)
    require(len(all_photo)==TOTAL_ATTEMPTED and all_photo.photo_id.is_unique,
            "Expected exactly 350399 separately measured unique original photo pixels")
    require(len(arms)==5 and len(arms[-1])==16384,
            "The complete authentic 16384-photo acquisition is absent")
    pre=pd.concat(arms[:4],ignore_index=True)
    require(len(pre)==334015 and pre.photo_id.is_unique,
            "Earlier 334015 authentic source checkpoint no longer identical")
    core_only=arms[0]
    result={
        "schema":"mois_350399_actual_original_RGB_photo_source_season_taxon_mix_v1",
        "source_sha256":{name:sha256(p) for name,p in
                         zip(("core317631","cumulative_broad317631",
                              "wave2_4096","cached8192","tail16384"),
                             (core,broad,wave2,cached8192,tail16384))},
        "all_photo_source_35e4":summary(all_photo),
        "pre_tail_334015_photo_source":summary(pre),
        "tail_16384_photo_source":summary(arms[-1]),
        "frame_comparison":{
            "core_317631":estimate(core_only,seed=20261019,draws=399),
            "before_tail_334015":estimate(pre,seed=20261020,draws=399),
            "completed_350399":estimate(all_photo,seed=20261021,draws=399)},
        "sampling_scope":"All actual measured photograph IDs, not original 93,209,351 source photo pixels. Original source-arm sampling probabilities of taxon-first photographs unknown.",
        "interpretation_limits":[
            "The four-term Kitagawa identity exactly reconstructs unpaired seasonal photographic-label differences and does not identify phenology, colour plasticity, thermal responses or pigment costs.",
            "A taxon absent from winter photographs is not assumed to be winter-absent, biologically nonflowering, or nonwhite.",
            "Taxon-within-season terms compare source photos at different sites and observer effort, not individual flowers.",
            "The image ROIs are machine-labelled four hues, not independent botanist-confirmed open petal colours.",
            "Whole-taxon resampling intervals are exploratory, not site-independent or familywise adjusted.",
            "The observation-weighted check is diagnostic, not correction for unsampled individuals or detection."
        ],
        "all_original_photo_IDs_source_unique":True,
        "core_first_broad_313535_duplicate_IDs_counted_once":True,
        "no_unmeasured_source_photo_IDs_assigned_hues":True,
    }
    figures=[]
    for frame,entry in result["frame_comparison"].items():
        for zone,r in entry["hemispheres"].items():
            for term,col in [
                ("within_common_taxon_photographic_white_rate_component_pp","within_shared_taxon_colour_rate"),
                ("changed_common_taxon_relative_photo_frequency_component_pp","shared_taxon_photo_frequency"),
                ("summer_photo_exclusive_taxon_support_component_pp","summer_only_source_taxon_support"),
                ("winter_photo_exclusive_taxon_support_component_pp","winter_only_source_taxon_support")]:
                figures.append(dict(frame=frame,hemisphere=zone,
                                    component=col,photo_white_difference_pp=r[term],
                                    raw_photo_white_summer_minus_winter_pp=r["raw_photo_white_summer_minus_winter_pp"],
                                    common_taxa=r["source_named_taxa_common_both_summer_winter"],
                                    source_photo_summer=r["summer_classified_original_photos"],
                                    source_photo_winter=r["winter_classified_original_photos"]))
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(figures).to_csv(out/FIGURE,index=False,lineterminator="\n")
    result["figure_source_sha256"]=sha256(out/FIGURE)
    (out/REPORT).write_text(json.dumps(result,indent=2,sort_keys=True,
                                      ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_REAL_350399_EXACT_TAXON_PHOTO_SEASON_DECOMPOSITION",
          json.dumps({frame:{
              z:{k:round(float(a["hemispheres"][z][k]),5) for k in
                 ("raw_photo_white_summer_minus_winter_pp",
                  "within_common_taxon_photographic_white_rate_component_pp",
                  "changed_common_taxon_relative_photo_frequency_component_pp",
                  "summer_photo_exclusive_taxon_support_component_pp",
                  "winter_photo_exclusive_taxon_support_component_pp")}
              for z in ("north_extratropical","south_extratropical")}
              for frame,a in result["frame_comparison"].items()},sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("core","broad","wave2","cached8192","tail16384","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    run(a.core,a.broad,a.wave2,a.cached8192,a.tail16384,a.output)


if __name__=="__main__":
    main()
