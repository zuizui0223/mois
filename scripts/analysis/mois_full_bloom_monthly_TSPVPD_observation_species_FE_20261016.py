#!/usr/bin/env python3
"""FULL original FCP 93M-photo dated BLOOM + four real monthly climate axes.

Reads ONLY FCP's completed 42,111-taxon ROI-v4 flower-positive source file
which has valid ORIGINAL observed_on, genus, taxon, photo and observer UUID
under SHA lock. NEVER uses the older 100543-photo cohort or an arbitrary
40-per-taxon seasonal subset.

Every positive original PHOTO is admitted. Multiple photographs of the
SAME source observation first become one observation-mean white vote; then
within each month/geographic climate cell, each named species contributes
one species-equal mean. This stops repeat multi-photo observations from
masquerading as independent blooming individuals. It does NOT undo biased
iNaturalist observer choice, source ROI nonresponse or hidden flowering taxa.

Reports full-snapshot photo and observation denominators, geographic cell
monthly climate join support, and exploratory cell/month/latitude-zone FE
models. Exact date/flower/ROI source membership is never inferred here.
No individual observation UUID/photo ID/observer/exact location is public.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sqlite3
from collections import Counter,defaultdict
from pathlib import Path

try:
    from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
        source_observation_day,climatic_season,sha,
        FULL_PHOTO_OPPORTUNITIES
    )
    from scripts.analysis.mois_original_full_fcp_dated_flower_ROI_20261016 import (
        SCHEMA as SOURCE_SCHEMA,STATUS as SOURCE_STATUS,MORPHS
    )
    from scripts.analysis.mois_join_monthly_dynamic_climate_20261016 import (
        DynamicNormals,REAL_MANIFEST_SHA
    )
    from scripts.analysis.mois_monthly_white_species_equal_cell_FE_20261016 import (
        fit_fe,DYNAMIC
    )
    from scripts.analysis.mois_complete_original_date_first_blooming_photo_census_20261016 import (
        SCHEMA as DATE_FIRST_COMPLETE_SCHEMA,
        STATUS as DATE_FIRST_COMPLETE_STATUS
    )
except ModuleNotFoundError:
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
        source_observation_day,climatic_season,sha,
        FULL_PHOTO_OPPORTUNITIES
    )
    from scripts.analysis.mois_original_full_fcp_dated_flower_ROI_20261016 import (
        SCHEMA as SOURCE_SCHEMA,STATUS as SOURCE_STATUS,MORPHS
    )
    from scripts.analysis.mois_join_monthly_dynamic_climate_20261016 import (
        DynamicNormals,REAL_MANIFEST_SHA
    )
    from scripts.analysis.mois_monthly_white_species_equal_cell_FE_20261016 import (
        fit_fe,DYNAMIC
    )
    from scripts.analysis.mois_complete_original_date_first_blooming_photo_census_20261016 import (
        SCHEMA as DATE_FIRST_COMPLETE_SCHEMA,
        STATUS as DATE_FIRST_COMPLETE_STATUS
    )

SCHEMA="mois_original_full_photo_observation_species_equal_monthly_climate_FE_v1"
STATUS="COMPLETE_SOURCE_DATED_BLOOMS_JOINED_TO_ORIGINAL_MONTHLY_CLIMATE"
TEMP_SCHEMA="mois_original_all_photo_bloom_observation_month_species_sqlite_v1"
PRIVATE_FIELDS=(
    "latitude_domain","climate_grid_row","climate_grid_col",
    "original_month","inat_taxon_id","genus",
    "n_original_independent_observations","n_original_classifiable_photographs",
    "n_original_white_photographs",
    "mean_observation_white_fraction",
    "site_month_tavg_c","site_month_solar_kj_m2_day",
    "site_month_prec_mm","site_month_derived_vpd_kpa",
)
MAX_SITE_CLIMATE_CACHE_CELLS=4096


def evaluate(positive_ndjson:Path,source_receipt:Path,
             climate:DynamicNormals,private_aggregates:Path,
             public_result:Path, *,
             synthetic=False,min_repeat_cells=30,min_public=10)->dict:
    for path in (private_aggregates,public_result):
        if path.exists() or path.with_name(path.name+".INCOMPLETE").exists():
            raise ValueError("Do not overwrite a previous complete monthly original photo-climate analysis")
    if min_public<10 or min_repeat_cells<2:
        raise ValueError("Cannot silently relax public privacy/source geography support")
    source=json.loads(source_receipt.read_text(encoding="utf8"))
    full_roi_complete=(source.get("schema")==SOURCE_SCHEMA and
                       source.get("status")==SOURCE_STATUS)
    original_date_first_complete=(
        source.get("schema")==DATE_FIRST_COMPLETE_SCHEMA and
        source.get("status")==DATE_FIRST_COMPLETE_STATUS)
    full_source_opportunity_count=(
        source.get("n_original_FCP_all_source_photo_ROI_results")
        if full_roi_complete else
        source.get("n_original_FCP_all_source_photo_opportunities"))
    if (not (full_roi_complete or original_date_first_complete)
        or source.get("private_actual_bloom_photo_month_hue_manifest_sha256")!=sha(positive_ndjson)
        or source.get("WorldClim_monthly_T_S_P_VPD_source_site_link_not_yet_run") is not True
        or source.get("full_42111_original_source_not_selectively_date_sampled") is not True
        or source.get("per_species_source_photo_cap","INVALID") is not None
        or source.get("synthetic_test_only") is not synthetic
        or climate.synthetic is not synthetic
        or original_date_first_complete and (
            source.get("n_original_date_known_photo_ROI_image_terminals")!=
            source.get("n_original_source_date_known_photo_candidates")
            or source.get("all_original_date_eligible_photo_pixels_attempted") is not True
            or source.get("all_missing_date_source_photos_were_excluded_before_roi") is not True)
        or not synthetic and (
            full_source_opportunity_count!=FULL_PHOTO_OPPORTUNITIES
            or climate.manifest_sha!=REAL_MANIFEST_SHA)):
        raise ValueError("Unsealed full FCP source photo ROI or complete DATE-FIRST original bloom/photo monthly source")
    private_aggregates.parent.mkdir(parents=True,exist_ok=True)
    public_result.parent.mkdir(parents=True,exist_ok=True)
    stage=private_aggregates.with_name(private_aggregates.name+".INCOMPLETE")
    public_stage=public_result.with_name(public_result.name+".INCOMPLETE")
    scratch=private_aggregates.with_name(private_aggregates.name+".SQLITE_INCOMPLETE")
    conn=None
    try:
        conn=sqlite3.connect(scratch)
        conn.execute("PRAGMA journal_mode=DELETE")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.execute("PRAGMA temp_store=FILE")
        conn.execute("""CREATE TABLE original_observations (
          uuid TEXT PRIMARY KEY,
          source_taxon INTEGER NOT NULL,
          genus TEXT NOT NULL,
          observer TEXT NOT NULL,
          domain TEXT NOT NULL,
          original_month INTEGER NOT NULL,
          climate_row INTEGER NOT NULL,
          climate_col INTEGER NOT NULL,
          n_bloom_classified_photos INTEGER NOT NULL DEFAULT 0,
          n_original_white_photos INTEGER NOT NULL DEFAULT 0
        )""")
        counted=Counter()
        lat_domain_photo=defaultdict(Counter)
        with positive_ndjson.open(encoding="utf8") as f:
            for line in f:
                r=json.loads(line)
                if (r.get("actual_photo_contains_classifiable_flower") is not True
                    or r.get("original_snapshot_20260927_date_proven") is not True
                    or r.get("original_ROI_v4_hue") not in MORPHS):
                    raise ValueError("Source photo has no original measured classified flower ROI")
                day,state=source_observation_day(r.get("original_observed_on"))
                if day is None or state!="ORIGINAL_SOURCE_DATE_KNOWN" or \
                   r.get("observed_month")!=day.month:
                    raise ValueError("Selected original flower lacks immutable observed_on date")
                lat,lon=float(r["latitude"]),float(r["longitude"])
                if not (math.isfinite(lat) and math.isfinite(lon)
                    and -90<=lat<=90 and -180<=lon<=180):
                    raise ValueError("Source classifiable flower has unusable geographic coordinates")
                domain,season=climatic_season(lat,day.month)
                if r.get("latitude_domain")!=domain or r.get("local_season")!=season:
                    raise ValueError("Original dated photo hemisphere/season source changed")
                if (type(r.get("original_taxon_id")) is not int
                    or r["original_taxon_id"]<=0
                    or not isinstance(r.get("original_genus"),str)
                    or not r["original_genus"]
                    or str(r.get("original_species_name","")).split()[0]!=r["original_genus"]
                    or not str(r.get("original_photo_id","")).isdecimal()
                    or not str(r.get("original_observation_uuid","")).strip()
                    or not str(r.get("original_observer","")).strip()):
                    raise ValueError("Incomplete original source photo/taxon/observer evidence")
                cell=climate.cell(lat,lon)
                row,col=cell if cell is not None else (-1,-1)
                white=int(r["original_ROI_v4_hue"]=="white")
                obs=str(r["original_observation_uuid"])
                cur=conn.execute("""INSERT INTO original_observations
                    (uuid,source_taxon,genus,observer,domain,original_month,
                     climate_row,climate_col,n_bloom_classified_photos,
                     n_original_white_photos)
                     VALUES(?,?,?,?,?,?,?,?,1,?)
                     ON CONFLICT(uuid) DO UPDATE SET
                       n_bloom_classified_photos=n_bloom_classified_photos+1,
                       n_original_white_photos=n_original_white_photos+excluded.n_original_white_photos
                     WHERE source_taxon=excluded.source_taxon
                       AND genus=excluded.genus
                       AND observer=excluded.observer
                       AND domain=excluded.domain
                       AND original_month=excluded.original_month
                       AND climate_row=excluded.climate_row
                       AND climate_col=excluded.climate_col""",
                    (obs,int(r["original_taxon_id"]),r["original_genus"],
                     str(r["original_observer"]),domain,day.month,
                     row,col,white))
                if cur.rowcount!=1:
                    raise ValueError("Duplicate FCP observation UUID changed source taxon/site/date")
                counted["n_positive_bloom_photos"]+=1
                counted["n_white_bloom_photos"]+=white
                lat_domain_photo[(domain,day.month)]["photo_n"]+=1
                lat_domain_photo[(domain,day.month)]["white_n"]+=white
                if counted["n_positive_bloom_photos"]%10000==0:conn.commit()
        conn.commit()
        if counted["n_positive_bloom_photos"]!=source.get(
           "n_original_date_and_site_known_ROI_bloom_hues") or \
           counted["n_white_bloom_photos"]!=source.get(
           "n_original_date_known_white_flower_photos"):
            raise ValueError("Actual full 93M original ROI positive flower totals drifted")
        n_observations=conn.execute("SELECT COUNT(*) FROM original_observations").fetchone()[0]
        if n_observations>counted["n_positive_bloom_photos"]:
            raise ValueError("Original source photo/independent observation denominator inverted")
        conn.execute("CREATE INDEX obs_zone_month_taxon ON original_observations(domain,original_month,source_taxon)")
        conn.execute("CREATE INDEX obs_cell_month_taxon ON original_observations(climate_row,climate_col,original_month,source_taxon)")
        conn.commit()

        # One UUID/observation vote for photos that came from the same
        # original source plant observation, then mean across independent
        # original observations within each species for seasonal community.
        domain_stats=defaultdict(lambda:{"taxa":0,"obs_equal_white_sum":0.,
                                        "observation_n":0})
        taxa_named=set()
        for domain,month,tid,genus,nobs,avg_obswhite in conn.execute(
             """SELECT domain,original_month,source_taxon,genus,
                       COUNT(*),AVG(1.0*n_original_white_photos/n_bloom_classified_photos)
                FROM original_observations
                GROUP BY domain,original_month,source_taxon,genus"""):
            bucket=domain_stats[(domain,month)]
            bucket["taxa"]+=1
            bucket["obs_equal_white_sum"]+=avg_obswhite
            bucket["observation_n"]+=nobs
            taxa_named.add(tid)

        climate_cells=[]
        n_photo_matched=n_obs_matched=n_photo_missing=n_obs_missing=0
        n_genus_source=set()
        with stage.open("w",newline="",encoding="utf8") as sink:
            writer=csv.DictWriter(sink,fieldnames=PRIVATE_FIELDS)
            writer.writeheader()
            current=None
            values=None

            def flush(current,values):
                if current is None:return
                zone,row,col,month=current
                p=climate.profile((row,col)) if row>=0 and col>=0 else None
                if p is None:
                    return
                normal=[
                    float(p["tavg"][month-1]-p["tavg"].mean()),
                    float(p["srad"][month-1]-p["srad"].mean()),
                    float(p["prec"][month-1]-p["prec"].mean()),
                    float(p["vpd"][month-1]-p["vpd"].mean())
                ]
                if not all(math.isfinite(v) for v in normal):
                    raise ValueError("Original WorldClim monthly profile invalid after source raster SHA")
                row_model={
                    "cell":f"{row}_{col}",
                    "domain":zone,"month":month,
                    "n_species":values["n_species"],
                    "white_species_equal":values["species_white_sum"]/values["n_species"],
                    **dict(zip(DYNAMIC,normal))
                }
                climate_cells.append(row_model)

            sql="""SELECT domain,climate_row,climate_col,original_month,
                           source_taxon,genus,COUNT(*) AS n_obs,
                           SUM(n_bloom_classified_photos) AS n_p,
                           SUM(n_original_white_photos) AS n_w,
                           AVG(1.0*n_original_white_photos/n_bloom_classified_photos) AS obs_white
                    FROM original_observations
                    GROUP BY domain,climate_row,climate_col,original_month,source_taxon,genus
                    ORDER BY domain,climate_row,climate_col,original_month,source_taxon,genus"""
            for zone,row,col,month,tid,genus,nobs,npics,nwhite,obs_white in conn.execute(sql):
                current_key=(zone,row,col,month)
                if current!=current_key:
                    flush(current,values)
                    current=current_key
                    values={"n_species":0,"species_white_sum":0.}
                values["n_species"]+=1
                values["species_white_sum"]+=obs_white
                n_genus_source.add(genus)
                p=climate.profile((row,col)) if row>=0 and col>=0 else None
                if p is None:
                    n_obs_missing+=nobs
                    n_photo_missing+=npics
                    continue
                n_photo_matched+=npics
                n_obs_matched+=nobs
                writer.writerow({
                    "latitude_domain":zone,
                    "climate_grid_row":row,"climate_grid_col":col,
                    "original_month":month,"inat_taxon_id":tid,"genus":genus,
                    "n_original_independent_observations":nobs,
                    "n_original_classifiable_photographs":npics,
                    "n_original_white_photographs":nwhite,
                    "mean_observation_white_fraction":obs_white,
                    "site_month_tavg_c":float(p["tavg"][month-1]),
                    "site_month_solar_kj_m2_day":float(p["srad"][month-1]),
                    "site_month_prec_mm":float(p["prec"][month-1]),
                    "site_month_derived_vpd_kpa":float(p["vpd"][month-1])
                })
            flush(current,values)
        if n_photo_matched+n_photo_missing!=counted["n_positive_bloom_photos"] \
           or n_obs_matched+n_obs_missing!=n_observations:
            raise ValueError("Full source photo bloom to month-specific climate loss not conserved")
        if not climate_cells:
            models=[{"fixed_effect_model":x,"status":"NO_ORIGINAL_SOURCE_SITE_MONTH_WITH_VALID_CLIMATE"}
                    for x in ("cell_only","cell_and_month","cell_and_latitude_zone_month")]
        else:
            models=[
                fit_fe(climate_cells,include_month=False,min_cells=min_repeat_cells),
                fit_fe(climate_cells,include_month=True,min_cells=min_repeat_cells),
                fit_fe(climate_cells,include_month=True,zone_month=True,min_cells=min_repeat_cells),
            ]
        present=defaultdict(set)
        for r in climate_cells:present[r["cell"]].add(r["month"])
        public_groups=[]
        for (zone,month),counts in sorted(lat_domain_photo.items()):
            n=counts["photo_n"]
            if n<min_public:
                public_groups.append({
                    "latitude_domain":zone,"original_month":month,
                    "public_suppressed_below":min_public
                })
                continue
            record=domain_stats[(zone,month)]
            public_groups.append({
                "latitude_domain":zone,"original_month":month,
                "n_all_original_image_verified_classifiable_flower_photos":n,
                "n_white_classifiable_original_photos":counts["white_n"],
                "photograph_white_rate_given_real_flower_ROI":counts["white_n"]/n,
                "n_independent_original_observations_with_bloom_photo":record["observation_n"],
                "n_source_species_with_original_bloom_photos":record["taxa"],
                "species_equal_observation_equal_white_fraction":(
                    record["obs_equal_white_sum"]/record["taxa"]
                    if record["taxa"] else None),
            })
        report={
            "schema":SCHEMA,"status":STATUS,
            "source_full_42111_original_FCP_bloom_receipt_sha256":sha(source_receipt),
            "source_positive_dated_bloom_private_photo_manifest_sha256":sha(positive_ndjson),
            "worldclim_original_48_source_manifest_sha256":climate.manifest_sha,
            "private_independent_observation_species_cell_month_climate_sha256":sha(stage),
            "n_original_FCP_all_photo_opportunities_processed_with_ROI":(
                full_source_opportunity_count if full_roi_complete else None),
            "n_original_FCP_all_photo_opportunities_reconciled":
                full_source_opportunity_count,
            "n_original_date_known_photo_ROI_terminals":(
                source.get("n_original_date_known_photo_ROI_image_terminals")
                if original_date_first_complete else None),
            "source_full_93m_image_ROI_or_original_date_first_strategy":(
                "FULL_ORIGINAL_IMAGE_ROI_ALL_SOURCE_PHOTOS" if full_roi_complete else
                "FULL_ORIGINAL_METADATA_DATE_FIRST_IMAGE_ROI_ONLY_DATE_VALID"),
            "source_date_unknown_photo_pixels_not_analyzed_by_design":(
                source.get("n_original_date_unknown_photo_opportunities_not_analysed_for_bloom")
                if original_date_first_complete else None),
            "n_original_bloom_and_source_date_qualified_photo_IDs":counted["n_positive_bloom_photos"],
            "n_original_white_bloom_date_known_photos":counted["n_white_bloom_photos"],
            "n_distinct_original_blooming_photo_observation_UUIDs":n_observations,
            "n_original_named_taxa_with_one_or_more_source_date_known_bloom_photos":len(taxa_named),
            "n_original_genera_with_date_known_classifiable_bloom":len(n_genus_source),
            "n_bloom_photos_monthly_worldclim_matched":n_photo_matched,
            "n_bloom_photos_no_monthly_worldclim":n_photo_missing,
            "n_original_independent_bloom_observations_climate_matched":n_obs_matched,
            "n_original_independent_bloom_observations_missing_monthly_worldclim":n_obs_missing,
            "n_source_climate_cell_months_with_bloom_species":len(climate_cells),
            "n_source_climate_cells_with_any_classified_original_bloom":len(present),
            "n_climate_cells_with_multiple_distinct_blooming_source_months":
                sum(len(v)>=2 for v in present.values()),
            "real_original_source_calendar_month_by_domain_photo_and_species_white":
                public_groups,
            "three_geographic_month_effect_model_support_and_exploratory_estimates":models,
            "source_photo_sampling_cap_per_species":None,
            "original_fcp_positive_93million_snapshot_all_photos_scope":True,
            "within_observation_multiple_photos_dependent_and_aggregated_before_species":True,
            "species_equal_photographed_bloom_rate_is_not_all_wild_flowering_plant_abundance":True,
            "photo_roi_missingness_from_full_source_must_be_read_separately":True,
            "observation_and_genus_cluster_dependence_not_fully_removed":True,
            "monthly_weather_is_1970_2000_climatological_not_photo_year":True,
            "synthetic_test_only":bool(synthetic),
        }
        public_stage.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf8")
        conn.close();conn=None
        scratch.unlink()
        os.replace(stage,private_aggregates)
        os.replace(public_stage,public_result)
        return report
    except Exception:
        if conn is not None:conn.close()
        for file in (scratch,stage,public_stage):
            if file.exists():file.unlink()
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-full-date-and-bloom-positive-photo-ndjson",required=True,type=Path)
    p.add_argument("--original-full-positive-bloom-receipt",required=True,type=Path)
    p.add_argument("--verified-original-worldclim48-root",required=True,type=Path)
    p.add_argument("--verified-original-worldclim48-manifest-sha256",required=True)
    p.add_argument("--out-private-observation-species-site-month-climate",required=True,type=Path)
    p.add_argument("--out-public-species-equal-month-white-climate",required=True,type=Path)
    args=p.parse_args()
    climate=DynamicNormals(args.verified_original_worldclim48_root,
                           args.verified_original_worldclim48_manifest_sha256,
                           max_cache_cells=MAX_SITE_CLIMATE_CACHE_CELLS)
    data=evaluate(args.original_full_date_and_bloom_positive_photo_ndjson,
                  args.original_full_positive_bloom_receipt,climate,
                  args.out_private_observation_species_site_month_climate,
                  args.out_public_species_equal_month_white_climate)
    print(json.dumps({k:data[k] for k in (
        "status","n_original_bloom_and_source_date_qualified_photo_IDs",
        "n_distinct_original_blooming_photo_observation_UUIDs",
        "n_climate_cells_with_multiple_distinct_blooming_source_months"
    )},sort_keys=True))


if __name__=="__main__":main()
