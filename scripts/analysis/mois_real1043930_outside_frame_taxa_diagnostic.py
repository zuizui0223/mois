#!/usr/bin/env python3
"""Locate what 32,710 attempted but colour-less photo IDs outside frozen
42,111-taxon/93m source are, and where the 20,884 genuine upstream no-nine
taxa fail. Do not invent 9-palette or colour from old 4-class labels.

This is a metadata-only source diagnosis (no image downloads).
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_real1043930_data_structure_audit import (
    PHOTO_SHA,PHOTO_UNION_SHA,UPSTREAM_SHA,sha,require
)
from scripts.analysis.mois_real1043930_global_zonal_colour import latitude_zone

SCHEMA="mois_true_1043930_outside_93m_and_zero_colour_taxa_diagnostic_v1"

def examine(source,union,upstream,out):
    source,union,upstream,out=map(Path,(source,union,upstream,out))
    require(not out.exists(),"Cannot rewrite an immutable source diagnosis")
    require(sha(source)==PHOTO_SHA and sha(union)==PHOTO_UNION_SHA
            and sha(upstream)==UPSTREAM_SHA,"Original frozen photo or taxon register bytes changed")
    p=pd.read_parquet(source,columns=["photo_id","inat_taxon_id","species",
           "source_arm_original","original_four_hue","continuous9_strict",
           "continuous9_exploratory_ambiguous","latitude","longitude","date"])
    q=pd.read_parquet(union,columns=["photo_id","observation_id","arm",
                    "inat_taxon_id","measurement_status"])
    c=pd.read_csv(upstream,dtype={"inat_taxon_id":str})
    require(len(p)==len(q)==1043930 and len(c)==42111
            and p.photo_id.astype(str).tolist()==q.photo_id.astype(str).tolist()
            and p.inat_taxon_id.astype(str).tolist()==q.inat_taxon_id.astype(str).tolist(),
            "Original source-photo IDs not conserved for diagnosis")
    p["photo_id"]=p.photo_id.astype(str)
    p["taxon"]=p.inat_taxon_id.astype(str)
    p["status"]=q.measurement_status.fillna("").astype(str)
    p["observation_id"]=q.observation_id.fillna("").astype(str)
    p["genus"]=p.species.fillna("").astype(str).str.split().str[0]
    p["lat_zone"]=[latitude_zone(v) for v in pd.to_numeric(p.latitude,errors="coerce").to_numpy(float)]
    p["four"]=p.original_four_hue.isin((
             "white","yellow_orange","red_pink","blue_purple"))
    valid=p.continuous9_strict
    names=set(c.inat_taxon_id)
    p["in_93209351_source_taxon_frame"]=p.taxon.isin(names)
    outside=p[~p.in_93209351_source_taxon_frame]
    require(len(outside)==32710 and not outside.continuous9_strict.any(),
            "Original 1m source outside 42,111 photo frame changes detected")
    def freq(label,series):
        z=series.value_counts(dropna=False)
        return [{"dimension":label,"value":str(k),"photo_attempts":int(v)}
                for k,v in z.items()]
    detail=(freq("original_acquisition_arm",outside.source_arm_original)+
            freq("ROI_terminal_status",outside.status)+
            freq("latitude_zone",outside.lat_zone))
    by_frame=p[p.in_93209351_source_taxon_frame].copy()
    valid_taxa=set(p.loc[valid,"taxon"])
    missing=c[~c.inat_taxon_id.isin(valid_taxa)].copy()
    existing_ids=by_frame.groupby("taxon").agg(
       original_photo_attempts=("photo_id","size"),
       source_legacy4_photo_count=("four","sum"),
       original_source9_ambiguous=("continuous9_exploratory_ambiguous","sum"),
       unique_original_observations=("observation_id","nunique"))
    missing=missing.join(existing_ids,on="inat_taxon_id")
    for key in ("original_photo_attempts","source_legacy4_photo_count",
                "original_source9_ambiguous","unique_original_observations"):
        missing[key]=missing[key].fillna(0).astype(int)
    require(len(missing)==20884
            and int(missing.n_unique_source_photo_ids.sum())==2290570,
            "No-original-nine-taxonomic photo opportunities unexpectedly changed")
    has4=missing.source_legacy4_photo_count.gt(0)
    bins=[]
    for lower,upper in ((0,0),(1,9),(10,99),(100,999),(1000,9999),(10000,10**10)):
        z=missing.loc[missing.n_unique_source_photo_ids.between(lower,upper)]
        bins.append({
           "upstream_original_photo_opportunity_range":f"{lower}-{upper}",
           "missing_valid_nine_taxa":len(z),
           "source_photo_opportunities_of_missing_taxa":int(z.n_unique_source_photo_ids.sum()),
           "upstream_taxa_missing_strict9_but_historical4_exists":int(
               z.source_legacy4_photo_count.gt(0).sum()),
           "legacy4_photo_labels_for_missing9_taxa":int(z.source_legacy4_photo_count.sum()),
           "source_image_attempts_for_missing9_taxa":int(z.original_photo_attempts.sum()),
           "observations_with_source_image_attempt_for_missing9_taxa":int(
               z.unique_original_observations.sum()),
        })
    report={
        "schema":SCHEMA,
        "outside_42111_source_photo_ids":len(outside),
        "outside_42111_source_unique_taxa":outside.taxon.nunique(),
        "outside_42111_source_original_4hue_photo_labels":int(outside.four.sum()),
        "outside_42111_source_valid9_photo_ROIs":int(outside.continuous9_strict.sum()),
        "outside_42111_source_ambiguous9_photo_ROIs":int(
            outside.continuous9_exploratory_ambiguous.sum()),
        "outside_42111_taxa_NOT_in_93209351_eligible_source_taxon_register":True,
        "photo_ID_membership_in_93m_source_NOT_proven_from_taxon_lookups":True,
        "outside_42111_source_photo_characterization":detail,
        "missing_original9_with_eligible_source_taxa":len(missing),
        "missing_original9_taxa_with_at_least100_upstream_photos":int(
            missing.n_unique_source_photo_ids.ge(100).sum()),
        "missing_original9_taxa_having_old_four_hue_label_only":int(has4.sum()),
        "missing_original9_taxa_total_known_old_four_hue_photos":int(
            missing.source_legacy4_photo_count.sum()),
        "missing_original9_taxa_source_image_attempts":int(missing.original_photo_attempts.sum()),
        "missing_taxon_opportunity_bins":bins,
        "source_registry_SHA256":UPSTREAM_SHA,
        "original_photo_pixels_posthoc_downloads":0,
        "source9_SHA256":PHOTO_SHA,
        "source_photo_ID_union_SHA256":PHOTO_UNION_SHA,
    }
    out.mkdir(parents=True)
    (out/"mois_1043930_source_outside_frame_missing_colour_taxa.json").write_text(
        json.dumps(report,sort_keys=True,indent=2,allow_nan=False)+"\n")
    pd.DataFrame(detail).to_csv(out/"mois_32710_taxa_outside_frame_photo_failure_types.csv",index=False)
    pd.DataFrame(bins).to_csv(out/"mois_20884_source_taxa_no_valid9_colour_by_photo_opportunity.csv",index=False)
    print("MOIS_1043930_REAL_OUTSIDE_FRAME_AND_TAXON_GAP",json.dumps(report,sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("source","union","upstream","output"):
        p.add_argument("--"+name,required=True,type=Path)
    z=p.parse_args()
    examine(z.source,z.union,z.upstream,z.output)
if __name__=="__main__":main()
