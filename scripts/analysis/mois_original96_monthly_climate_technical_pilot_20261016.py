#!/usr/bin/env python3
"""96 ACTUALLY photographed FCP ROI-v4 technical-pilot flower hues × monthly climate.

This source is a FOUR-TAXON, NONREPRESENTATIVE image-processing repair PILOT,
not a new random sample of 42,111 species. Source photo IDs/coordinate rows
never enter the public output. The original ROI-classifiable denominator
is 20/96 with only three white photos; NO worldwide or causal regression.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import sys
import zipfile
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts.analysis.mois_join_monthly_dynamic_climate_20261016 import (
    DynamicNormals, COLOUR
)
from scripts.analysis.mois_seasonal_flower_colour_opportunity_20261016 import local_season

ORIGINAL_CSV_SHA="911ec9edd10a7126bee28cfa67738cef0ff4d275b1f0ca1c06e3f26e7e290039"
ORIGINAL_PILOT_N=96
ORIGINAL_CLASS_N=20
ORIGINAL_WHITE_N=3
SOURCE_CSV="real_original_roi_v4_new_photo_pilot_repaired.csv"
SOURCE_RECEIPT="real_original_roi_v4_new_photo_pilot_repaired_receipt.json"
SCHEMA="mois_real_fcp96_original_photo_monthly_climate_technical_pilot_v1"
STATUS="ACTUAL_96_IMAGE_ROI_TECHNICAL_PILOT_CLIMATE_JOIN_NOT_GLOBAL_SEASON_EFFECT"
DYNAMIC_COLS=("monthly_tavg_c","monthly_srad_kj_m2_day","monthly_prec_mm",
              "monthly_vapr_kpa","monthly_vpd_kpa")
SOURCE_COLS={
    "inat_taxon_id","species","observation_id","photo_id","measurement_id",
    "observed_on","latitude","longitude","morph","measurement_status",
}
REGIONS=("NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL","TROPICAL")


def audit(pilot_zip:Path, source:DynamicNormals, output:Path, *,
          synthetic:bool=False):
    if output.exists():
        raise ValueError("Refuse overwrite of already reported original photo source climatology")
    if source.synthetic is not synthetic:
        raise ValueError("Synthetic raster normals cannot be published as original climate")
    with zipfile.ZipFile(pilot_zip) as z:
        if set(z.namelist())!={SOURCE_CSV,SOURCE_RECEIPT} or z.testzip() is not None:
            raise ValueError("Not the admitted original 96 ROI-v4 source pilot ZIP")
        raw=z.read(SOURCE_CSV)
        source_receipt=json.loads(z.read(SOURCE_RECEIPT))
    observed_sha=hashlib.sha256(raw).hexdigest()
    if (not synthetic and observed_sha!=ORIGINAL_CSV_SHA
        or source_receipt.get("schema")!="fcp_real_original_roi_v4_case_sensitive_extension_repair_v1"
        or source_receipt.get("status")!="REPAIRED_REAL_PHOTO_TECHNICAL_PILOT_NOT_42111_CENSUS"
        or source_receipt.get("source_photo_id_count")!=ORIGINAL_PILOT_N
        or source_receipt.get("n_actual_original_image_measurements_total")!=ORIGINAL_PILOT_N
        or source_receipt.get("n_classified_four_state_photo_rows")!=ORIGINAL_CLASS_N
        or source_receipt.get("terminal_csv_sha256")!=observed_sha
        or source_receipt.get("original_method",{}).get("location_blind") is not True):
        raise ValueError("Not a source-attested remeasured original ROI-v4 96 photo pilot")
    reader=csv.DictReader(io.StringIO(raw.decode("utf8")))
    if not SOURCE_COLS.issubset(reader.fieldnames or ()):
        raise ValueError("The original FCP technical pilot lacked photo site/date and ROI")
    seen_photo=set()
    seen_observation=set()
    taxa=set()
    grouped=defaultdict(lambda:{
        "n":0,"class":0,"white":0,"classifiable_by_taxon":defaultdict(lambda:[0,0]),
        "sums":defaultdict(float),"valid_for_monthly_climate":0,"classifiable_with_climate":0})
    statuses=Counter()
    totals=Counter()
    for r in reader:
        pid,obs=r["photo_id"],r["observation_id"]
        if not pid or not obs or pid in seen_photo or obs in seen_observation:
            raise ValueError("Original 96 image/observation IDs duplicated or absent")
        seen_photo.add(pid);seen_observation.add(obs)
        taxa.add(r["inat_taxon_id"])
        valid=r["measurement_status"]=="classified_four_state_morph"
        if valid!=(r["morph"] in COLOUR):
            raise ValueError("Original ROI method flower hue classification was modified")
        domain,localseason,month,_=local_season(r["latitude"],r["observed_on"])
        if domain not in REGIONS or localseason=="UNKNOWN_SEASON":
            raise ValueError("This original 96 case pilot did not contain reliable geo/calendar season")
        climate=source.site_month(r["latitude"],r["longitude"],r["observed_on"])
        row=grouped[(domain,localseason)]
        row["n"]+=1
        row["class"]+=int(valid)
        row["white"]+=int(valid and r["morph"]=="white")
        if valid:
            t=row["classifiable_by_taxon"][r["inat_taxon_id"]]
            t[0]+=1;t[1]+=int(r["morph"]=="white")
        statuses[r["measurement_status"]]+=1
        totals["n"]+=1;totals["class"]+=int(valid)
        totals["white"]+=int(r["morph"]=="white")
        if climate["source_climate_join_status"]=="VALID":
            totals["source_monthly_climate_join"]+=1
            row["valid_for_monthly_climate"]+=1
            if valid:
                row["classifiable_with_climate"]+=1
                for key in DYNAMIC_COLS:
                    row["sums"][key]+=float(climate[key])
        else:
            totals["source_monthly_climate_unavailable"]+=1
    if not synthetic and (
        len(seen_photo)!=ORIGINAL_PILOT_N
        or len(taxa)!=4 or totals["class"]!=ORIGINAL_CLASS_N
        or totals["white"]!=ORIGINAL_WHITE_N
        or statuses.get("not_evaluable_roi_or_flip_gate")!=67
        or statuses.get("not_evaluable_ambiguous_palette_composition")!=9):
        raise ValueError("Original 96 photo technical cohort denominator/hue changed")
    result=[]
    for (domain,season),g in sorted(grouped.items()):
        # Classifiable source photos with missing climate are not silently
        # assigned a monthly source value. Keep explicit two denominators.
        nclass=g["class"]
        n_valid_climate=g["valid_for_monthly_climate"]
        n_valid_class=g["classifiable_with_climate"]
        q={
            "domain":domain,"local_season":season,
            "sampled_original_photo_opportunities":g["n"],
            "four_hue_classifiable":nclass,
            "classified_white":g["white"],
            "white_fraction_among_classifiable":g["white"]/nclass if nclass else None,
            "species_equal_white_fraction_in_same_season":
                sum(nw/nc for nc,nw in g["classifiable_by_taxon"].values())/len(
                    g["classifiable_by_taxon"]) if nclass else None,
            "n_taxa_with_at_least_one_classifiable_photo":len(g["classifiable_by_taxon"]),
            "n_original_monthly_climate_joined_of_all_source":n_valid_climate,
            "n_classifiable_with_known_monthly_climate":n_valid_class,
            "source_classified_monthly_climate_means":(
                {key:g["sums"][key]/n_valid_class for key in DYNAMIC_COLS}
                if n_valid_class else None
            ),
        }
        result.append(q)
    report={
        "schema":SCHEMA,"status":STATUS,
        "source_original_FCP_96_pilot_csv_sha256":observed_sha,
        "source_48_worldclim_monthly_manifest_sha256":source.manifest_sha,
        "n_original_source_photo_IDs":len(seen_photo),
        "n_original_distinct_observations":len(seen_observation),
        "n_original_taxa":len(taxa),
        "n_roi_original_four_colour_classifiable":totals["class"],
        "n_roi_original_classified_white":totals["white"],
        "n_monthly_source_site_climate_joined":totals["source_monthly_climate_join"],
        "n_monthly_source_site_climate_missing":totals["source_monthly_climate_unavailable"],
        "roi_status_counts":dict(statuses),
        "observed_season_rows":result,
        "original_photo_selection_is_nonrepresentative_four_taxon_technical_pilot":True,
        "original_source_population_42111_taxa_photo_IDs_joined":False,
        "prespecified_site_cell_season_FE_within_month_estimated":False,
        "seasonal_climate_white_association_or_fitness_estimated":False,
        "no_photo_identifiers_coordinates_or_individual_observer_IDs_published":True,
        "climatology_years":"1970-2000 NOT actual photographic weather",
        "synthetic_test_only":bool(synthetic),
        "note":"Twenty original ROI classifiable photos and three white across four deliberately chosen species; source photos had high nonclassification. Group climate means are simple descriptive original site-month averages and cannot estimate population seasonal effects."
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf8")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-96-source-zip",type=Path,required=True)
    p.add_argument("--worldclim-dynamic-root",type=Path,required=True)
    p.add_argument("--worldclim-dynamic-sha256",required=True)
    p.add_argument("--out-source-audit",type=Path,required=True)
    a=p.parse_args()
    climate=DynamicNormals(a.worldclim_dynamic_root,a.worldclim_dynamic_sha256)
    d=audit(a.original_96_source_zip,climate,a.out_source_audit)
    print("MOIS_REAL_ORIGINAL96_PHOTO_DYNAMIC_MONTHLY_CLIMATE_TECHNICAL_PILOT",
          json.dumps(d,sort_keys=True))


if __name__=="__main__":main()
