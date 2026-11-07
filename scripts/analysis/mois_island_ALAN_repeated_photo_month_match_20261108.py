#!/usr/bin/env python3
"""MOIS observed-season match for repeated island/mainland photo species sets.

Use source-dated, original strict9-success photographs of exactly accepted
WCVP taxa in genuine 2019/21/23/24/25 island-mainland photo-cell cohorts.
At the SAME physical source 0.05° cell+GSHHG landmass, compare two actually
photographed years. Distinguish taxa Jaccard for pairs with same majority
photograph month, for different majority months, and for only genuinely
overlapping photograph months within each year.

This tests observation-season comparability, NOT real change/turnover of field
plant communities, species population abundance, adaptation or ALAN mechanism.
No new photographs, no inferred anthesis/month, no missing ROI white label or
rejected MNAR/tipping runs. No undated photos silently assigned a month.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
import numpy as np
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis import mois_island_ALAN_repeated_sites_photo_year_20261108 as repeat
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import cell_id

SPECIES="wcvp_accepted_plant_name_id"
YEARS=(2019,2021,2023,2024,2025)
EXPECTED={"strict_5p3r3s":7982,"measured_3p2r2s":20543}

def original_photo_source_months(original,wcvp,eligible):
    """Actual original accepted, strict9 photo and observed calendar month."""
    x=original.loc[
        original.source_subarm.eq("historical_main") &
        original.strict9 &
        original.strict9_white.notna() &
        original.latitude.between(-90,90)&
        original.longitude.between(-180,180)].copy()
    dates=x.observed_on.astype(str)
    x["year"]=pd.to_numeric(dates.str.slice(0,4),errors="coerce")
    x["month"]=pd.to_numeric(dates.str.slice(5,7),errors="coerce")
    x=x.loc[x.year.isin(YEARS)].copy()
    x.year=x.year.astype(int)
    lookup=wcvp.loc[wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
        ["inat_taxon_id",SPECIES]].copy()
    panel.require(lookup.inat_taxon_id.is_unique,"Accepted taxon photo source duplicated")
    x=x.merge(lookup,on="inat_taxon_id",how="inner",validate="many_to_one")
    x=x.loc[x[SPECIES].notna()].copy()
    panel.require(x.photo_id.is_unique,"Source original ROI photo repeated across accepted taxa")
    x["cell_id"]=cell_id(x.longitude,x.latitude)
    z=eligible[["year","target_id","land_id","land_class"]].copy()
    split=z.target_id.astype(str).str.split("|",regex=False,expand=True)
    panel.require(split.shape[1]==3 and split[0].eq(z.year.astype(str)).all()
        and split[2].eq(z.land_id.astype(str)).all(),
        "Original source knownland-year target identity mismatch")
    z["cell_id"]=split[1].to_numpy()
    panel.require(not z.duplicated(["year","cell_id"]).any(),
        "Ambiguous original 0.05° cell was falsely assigned to island/mainland")
    x=x.merge(z,on=["year","cell_id"],how="inner",validate="many_to_one")
    panel.require(x.photo_id.is_unique,
        "Source photo ID duplicated by geographic/cohort join")
    count_all=len(x)
    valid=x.month.between(1,12)&x.month.notna()
    dated=x.loc[valid].copy()
    dated.month=dated.month.astype(int)
    dated["repeat_site_id"]=dated.cell_id.astype(str)+"|"+dated.land_id.astype(str)
    return dated,{
       "accepted_original_strict9_source_photos_joined":count_all,
       "actual_observed_calendar_month_photos":len(dated),
       "source_observation_date_month_unavailable_photos":int(count_all-len(dated)),
       "original_photograph_calendar_month_only_not_inferred_flowering_month":True,
       "photo_years_are_original_and_not_imputed":list(YEARS),
    }

def actual_photo_year_month_sets(d):
    """A real original photo is counted once. Tie in dominant month -> earliest month."""
    ret={}
    if not len(d):
        return ret
    grouped=d.groupby(["repeat_site_id","year"],sort=True)
    for (site,year),photos in grouped:
        counts=photos.groupby("month").size()
        maxcount=int(counts.max())
        dominant=int(counts.index[counts.eq(maxcount)].min())
        monthsets={
            int(m):frozenset(rows[SPECIES].astype(str))
            for m,rows in photos.groupby("month")
        }
        panel.require(all(1<=m<=12 for m in monthsets),
            "Source observation month outside actual Gregorian calendar")
        ret[(site,int(year))]={
            "dominant":dominant,
            "all_taxa":frozenset(photos[SPECIES].astype(str)),
            "monthsets":monthsets,
            "n_source_photos":len(photos),
        }
    return ret

def observed_photo_year_pairs(repeated,month_data):
    """True consecutive observation years of same site, not invented annual gaps."""
    mapping=actual_photo_year_month_sets(month_data)
    pairs=[];undated_pairs=0
    for site,g in repeated.groupby("repeat_site_id",sort=False):
        g=g.sort_values("year")
        records=list(g[["year","land_class"]].itertuples(index=False,name=None))
        for (y0,cls0),(y1,cls1) in zip(records[:-1],records[1:]):
            panel.require(y1>y0 and cls0==cls1,
                "Genuine repeated photo-site geography/year drift")
            a=mapping.get((site,int(y0)));b=mapping.get((site,int(y1)))
            if not a or not b:
                undated_pairs+=1
                continue
            union=a["all_taxa"]|b["all_taxa"]
            shared=a["all_taxa"]&b["all_taxa"]
            panel.require(len(union)>0,"Original photo species sets are empty")
            both_months=set(a["monthsets"])&set(b["monthsets"])
            matched0=frozenset().union(*(a["monthsets"][m] for m in both_months))
            matched1=frozenset().union(*(b["monthsets"][m] for m in both_months))
            matched_union=matched0|matched1
            matched_jaccard=(len(matched0&matched1)/len(matched_union)
                             if both_months and matched_union else None)
            dist=abs(a["dominant"]-b["dominant"])
            gap=min(dist,12-dist)
            pairs.append({
               "land_class":cls0,
               "calendar_year_gap":int(y1-y0),
               "same_dominant_photo_month":a["dominant"]==b["dominant"],
               "dominant_photo_month_circular_distance":gap,
               "at_least_one_shared_actual_photo_month":bool(both_months),
               "jaccard_all_actual_dated_photos":len(shared)/len(union),
               "jaccard_restrict_to_same_actual_calendar_months":matched_jaccard,
               "no_same_taxon_across_actual_photo_years":len(shared)==0,
               "shared_month_no_overlap":bool(both_months and matched0.isdisjoint(matched1)),
            })
    return pd.DataFrame(pairs),undated_pairs

def real_photo_month_support(d,repeat_frames):
    """Aggregate only, no original photo ID/location species rows are exported."""
    pairs,unusable=observed_photo_year_pairs(repeat_frames,d)
    out={
        "photo_month_year_pairs_with_valid_original_dates":len(pairs),
        "repeated_photo_year_pairs_without_date_comparison":unusable,
        "actual_photo_month_not_inferred_phenological_flowering_month":True,
        "all_year_pairs_contain_both_years_same_genuine_GSHHG_photo_site":True,
        "results":{},
    }
    for land in ("island","mainland"):
        land_pairs=pairs.loc[pairs.land_class.eq(land)]
        parts={
          "all_dated_pairs":land_pairs,
          "same_dominant_photo_month":land_pairs.loc[
              land_pairs.same_dominant_photo_month],
          "different_dominant_photo_month":land_pairs.loc[
              ~land_pairs.same_dominant_photo_month],
          "at_least_one_genuinely_shared_photographed_month":land_pairs.loc[
              land_pairs.at_least_one_shared_actual_photo_month],
          "different_dominant_month_2plus_month_apart":land_pairs.loc[
              land_pairs.dominant_photo_month_circular_distance>=2],
        }
        summary={}
        for name,part in parts.items():
            vals=part.jaccard_all_actual_dated_photos
            month_match=part.jaccard_restrict_to_same_actual_calendar_months.dropna()
            summary[name]={
               "observed_repeated_siteyear_pairs":len(part),
               "observed_pair_source_Jaccard_mean":float(vals.mean()) if len(vals) else None,
               "observed_pair_source_Jaccard_median":float(vals.median()) if len(vals) else None,
               "zero_shared_photographed_accepted_taxon_pairs":int(
                   part.no_same_taxon_across_actual_photo_years.sum()),
               "fraction_zero_observed_accepted_taxon_overlap":float(
                   part.no_same_taxon_across_actual_photo_years.mean()) if len(part) else None,
               "genuine_same_month_subcomparison_pairs":len(month_match),
               "genuine_same_month_source_Jaccard_median":float(
                   month_match.median()) if len(month_match) else None,
               "genuine_same_month_source_Jaccard_mean":float(
                   month_match.mean()) if len(month_match) else None,
            }
        out["results"][land]=summary
    return out

def run(paths,out):
    out=Path(out)
    panel.require(not out.exists(),"Original immutable evidence cannot be overwritten")
    frames=tier.build_frozen_five_years(paths)
    cohorts=tier.by_year_tier(frames)
    canonical,wcvp=panel.load_frozen_photos(paths)
    report={}
    for label,n in EXPECTED.items():
        data=cohorts[label]
        panel.require(len(data)==n and data.target_id.is_unique,
                      "Original measured photo year cohort changed")
        prepared,repeat_source=repeat.repeat_site_frame(data)
        dated,quality=original_photo_source_months(canonical,wcvp,data)
        panel.require(dated.target_id.isin(set(data.target_id)).all(),
                      "Source month photo joined an ineligible colour-cell-year")
        summary=real_photo_month_support(dated,repeat_source)
        report[label]={
          "source_photo_date_provenance":quality,
          "original_repeat_site_support":{
              "repeated_true_source_cell_years":len(repeat_source),
              "repeated_true_source_sites":int(repeat_source.repeat_site_id.nunique()),
              "repeated_island_siteyear":int(repeat_source.is_island.sum()),
          },
          "month_comparable_observed_photo_taxon_pairs":summary,
        }
    receipt={
        "schema":"MOIS_true_original_photo_month_overlap_repeat_island_sites_20261108_v1",
        "original_observation_years":list(YEARS),
        "source_original_strict9_reused_only":True,
        "actual_month_observed_on_original_photograph_not_inferred_bloom_month":True,
        "photo_no_record_not_plant_absence":True,
        "date_month_missing_excluded_not_filled":True,
        "no_new_original_photo_RGB_download_or_ROI_inference":True,
        "no_2022_2020_source_substitution":True,
        "no_MNAR_tipping_or_IPW":True,
        "original_input_SHA256":{k:panel.sha(v) for k,v in paths.items()},
        "results":report
    }
    out.mkdir(parents=True)
    (out/"photo_month_matched_source_taxon_overlap_receipt.json").write_text(
        json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_TRUE_SOURCE_MONTH_MATCH_REPEATED_PHOTO_SITES",
        json.dumps(report,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    keys=("canonical","wcvp",
      "geo2019","viirs2019","pop2019",
      "geo2021","viirs2021","pop2021",
      "geo2023","viirs2023","pop2023",
      "geo2024","cn2024","uc2024","uc2024receipt",
      "geo2025","geo2025receipt","viirs2025")
    for key in keys:
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in keys},a.out)
