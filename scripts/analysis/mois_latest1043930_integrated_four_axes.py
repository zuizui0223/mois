#!/usr/bin/env python3
"""Reconcile the 2026-10-25 LATEST real-1,043,930 MOIS four-axis analyses.

No photo image inference. This receipt collects independently completed
photo SHA, per-photo environmental GSHHG, seasonal taxon-decomposition and
classification selection-bias results. Explicitly contrasts genuine
350,399-photo earlier estimands without merging old coefficients into newest.
"""
from __future__ import annotations
import argparse, hashlib, json,csv
from pathlib import Path

SCHEMA="mois_actual_1043930_original_photos_all_four_gradients_current_20261025_v1"
HUES=("white","yellow_orange","red_pink","blue_purple")
def need(value,msg):
 if not value:raise ValueError("MOIS1043930_FOUR_AXES: "+msg)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as p:
  for b in iter(lambda:p.read(4<<20),b""):h.update(b)
 return h.hexdigest()
def read(path):
 return json.loads(Path(path).read_text(encoding="utf-8"))
def main_latest(source,climate,water,elev_island,season,bias,prior,out):
 source,climate,water,elev_island,season,bias,prior,out=map(Path,
  (source,climate,water,elev_island,season,bias,prior,out))
 need(not out.exists(),"Do not overwrite already interpreted/latest analysis")
 s,c,w,g,t,b,old=map(read,(source,climate,water,elev_island,season,bias,prior))
 need(s["schema"]=="mois_fully_signed_latest_original_photo_colour_1043930_v1"
    and s["source_pixel_attempts"]==1043930
    and s["source_four_hue_classified"]==459235
    and s["latest_fcp5_genuine_complete_groups"]==35,
    "Not actual complete 35-group photo ROI source")
 need(c["schema"]=="mois_real_latest1043930_photo_id_aligned_worldclim_humidity_rain_altitude_v1"
    and c["original_exact_photo_attempts"]==1043930
    and c["four_hue_classified_climate_complete"]==458174
    and c["new_original_photo_ID_climate_sites_sampled"]==693531
    and c["new1043930_true_photo_source_sha256"]==s["photo_parquet_sha256"],
    "Climate table is not latest original-photo identity matched")
 need(w["schema"]=="mois_original1043930_four_hue_latest_climate_source_conditional_water_altitude_v1"
    and w["complete_original_actual_photo_attempts"]==1043930
    and w["original_automatic_four_hue_ROI_success"]==459235
    and w["source_four_hue_parquet_SHA256"]==s["photo_parquet_sha256"]
    and w["source_annual_climate_parquet_SHA256"]==c["new1043930_climate_parquet_sha256"],
    "Water altitude regression did not use same latest complete original photo IDs")
 need(g["schema"]=="mois_latest1043930_source_photo_GSHHG_island_and_nearby_altitude_v1"
    and g["real_original_source_photo_attempts"]==1043930
    and g["latest_four_hue_classified"]==459235
    and g["source_shas"]["colour"]==s["photo_parquet_sha256"]
    and g["source_shas"]["climate"]==c["new1043930_climate_parquet_sha256"],
    "Island landmass environmental response source SHA mismatch")
 need(t["schema"]=="mois_full1043930_real_photos_original_season_taxon_decomposition_v1"
    and t["full_photo_attempts"]==1043930
    and t["classified_four_hues"]==459235
    and t["source_sha256"]==s["photo_parquet_sha256"],
    "Seasonal taxonomy analysis source does not cover all latest photos")
 need(b["schema"]=="mois_real1043930_five_axis_photo_colour_classifiability_bias_v1"
    and b["overall"]["original_true_photo_attempts"]==1043930
    and b["overall"]["four_hue_classified"]==459235
    and b["source_photo_sha256"]==s["photo_parquet_sha256"]
    and b["source_climate_sha256"]==c["new1043930_climate_parquet_sha256"]
    and b["source_GSHHG_photo_parquet_SHA256"]==g["output_geography_photo_ID_parquet_SHA256"],
    "Colour classifiability selection audit is incomplete or conflates sources")
 need(old["source_actual_photo_attempts"]==350399
    and old["source_four_hue_classified"]==164297,
    "Not previous actual-photo 350399 four-axis reference")
 # Keep same old 5deg nominal genus by environment, even though the latest
 # additionally conditions on distinct K1/FCP arms (different sampling).
 previous=old["axes"]
 old_loc=previous["wet_dry"]["univariate_within_genus_region_and_source_arm_slopes"]
 new_loc=w["moisture_elevation_within_cohort_genus_geocell_univariate"]
 new_partial=w["simultaneous_climate_correlation_diagnostic"]
 def oldval(predictor):
  rows=[r for r in old_loc if r["grid_degrees"]==5
    and r["cohort"]=="combined_arm_conditioned"
    and r["predictor"]==predictor]
  need(len(rows)==1,"Missing frozen old humidity/altitude regression "+predictor)
  r=rows[0]
  return r["pp_per_full_source_taxon_cell_sd"],r["genus_bootstrap_nominal_95_pp"]
 def newval(predictor):
  v=[r for r in new_loc if r["grid_degrees"]==5 and r["predictor"]==predictor]
  need(len(v)==1,"Missing newest 5degree genus local slope "+predictor)
  r=v[0]
  return r["white_pp_per_original_taxon_cell_SD"],r["nominal_95_genus_cluster_bootstrap_pp"]
 def partial(predictor):
  r=[x for x in new_partial if x["grid_degrees"]==5 and x["predictor"]==predictor]
  need(len(r)==1,"Missing newer simultaneous correlation sensitivity")
  return r[0]
 summary=[]
 for label,predictor in [
   ("annual_precipitation","log1p_annual_rainfall"),
   ("air_vapour","mean_atmospheric_vapour"),
   ("elevation_source_within_genus","elevation")
 ]:
  older,oldci=oldval(predictor);new,newci=newval(predictor)
  detail=partial(predictor)
  summary.append({
   "axis":label,"method":"within source-cohort genus geographic 5degree cell",
   "old_original_attempts":350399,"latest_original_attempts":1043930,
   "old_effect_white_percentage_points_per_SD":older,
   "new_effect_white_percentage_points_per_SD":new,
   "new_nominal95_lower":newci[0],"new_nominal95_upper":newci[1],
   "simultaneous_other_climate_partial_new_pp":
      detail["cohort_genus_region_fully_partial_white_pp_per_SD"],
   "partial_new_nominal95":
      detail["cluster_genus_nominal95_partial_ci_pp"]
  })
 old_elev=previous["elevation"]["landmass_5deg_313535_baseline"]
 new_elev=g["elevation_same_landmass_nearby_5deg"]
 summary.append({
   "axis":"nearby_high_altitude",
   "method":"same landmass nearby 5degree source taxa-equal cells >=1500 minus 0-500m",
   "old_original_attempts":350399,"latest_original_attempts":1043930,
   "old_effect_white_percentage_points_per_SD":old_elev["high_minus_low_white_pp"],
   "new_effect_white_percentage_points_per_SD":new_elev["high_minus_low_pp"],
   "new_nominal95_lower":new_elev["nominal_neighbourhood_t95_pp"][0],
   "new_nominal95_upper":new_elev["nominal_neighbourhood_t95_pp"][1],
   "n_neighbourhoods_old":old_elev["paired_landmass_geographic_neighbourhoods"],
   "n_neighbourhoods_new":new_elev["paired_same_landmass_neighbourhoods"]
 })
 old_island=previous["islands"]
 new_island=g["island_mainland_environment_matched_and_spatial_regression"]
 summary.append({
   "axis":"island_matched",
   "method":"environment matching island mainland within same coarse geographic zone",
   "old_original_attempts":350399,"latest_original_attempts":1043930,
   "old_effect_white_percentage_points_per_SD":old_island["matched_pair_white_pp"],
   "new_effect_white_percentage_points_per_SD":new_island["matched_island_minus_mainland_white_pp"],
   "new_nominal95_lower":new_island["matched_nominal_region_cluster_95_pp"][0],
   "new_nominal95_upper":new_island["matched_nominal_region_cluster_95_pp"][1],
   "n_pairs_old":old_island["matched_pairs"],"n_pairs_new":new_island["optimal_match_pairs"],
   "new_raw_island_minus_mainland_pp":g["island_minus_mainland_raw_white_pp"],
   "new_spatial_climate_adjusted_island_minus_mainland_pp":
       new_island["spatial_climate_adjusted_island_minus_mainland_white_pp"],
   "new_adjusted_nominal95":
       new_island["adjusted_nominal_spatial_cluster_95_pp"]
 })
 for hemi in ("north_extratropical","south_extratropical"):
  oldrow=previous["season"]["hemispheres"][hemi]
  newrow=t["hemispheres"][hemi]
  summary.append({
   "axis":"season_"+hemi,
   "method":"unpaired original photo summer-minus-winter white rate",
   "old_original_attempts":350399,"latest_original_attempts":1043930,
   "old_effect_white_percentage_points_per_SD":oldrow["raw_photo_white_summer_minus_winter_pp"],
   "new_effect_white_percentage_points_per_SD":newrow["raw_photo_white_summer_minus_winter_pp"],
   "new_shared_taxon_photo_frequency_term_pp":
      newrow["changed_common_taxon_relative_photo_frequency_component_pp"],
   "summer_classified_photos_latest":newrow["summer_classified_original_photos"],
   "winter_classified_photos_latest":newrow["winter_classified_original_photos"]
  })
 report={
   "schema":SCHEMA,
   "global_original_photo_ID_attempts_latest":1043930,
   "actual_classifiable_four_hue_photographs_latest":459235,
   "source_hues":s["four_hue_counts"],
   "source_colour_unknowns":584695,
   "annual_environment_classified_complete":458174,
   "actual_new_sites_climatology_calculated":693531,
   "source_GSHHG_previous_plus_new_photo_IDs":1043930,
   "source_arm_FCP5_groups_sealed":35,
   "recomputed_four_axes":summary,
   "newest_detail_by_axis":{
      "wet_dry_elevation":w,
      "island_altitude":g,
      "summer_winter_taxon_photo":t,
      "classifiability_selection_by_source_photo":b
   },
   "individual_input_SHA256":{
      k:sha(p) for k,p in [
        ("latest_photo_union",source),("climate_receipt",climate),
        ("wet_dry",water),("island_altitude",elev_island),
        ("season_taxon",season),("observation_bias",bias),("previous350399",prior)
      ]},
   "photo_image_downloads_for_this_reanalysis":0,
   "ROI_four_hue_unknowns_relabelled":0,
   "ecological_mechanistic_causality_or_native_wild_status_verified":False,
   "multiple_exploratory_axes_not_familywise_corrected":True,
   "old_and_new_coefficients_differ_in_source_recruitment_and_conditioning":True,
 }
 out.mkdir(parents=True)
 (out/"mois_latest1043930_four_axis_comparison.json").write_text(
  json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
 pd_fields=sorted(set().union(*(set(x) for x in summary)))
 with (out/"mois_latest1043930_four_axis_comparison.csv").open("w",newline="") as f:
  wri=csv.DictWriter(f,fieldnames=pd_fields)
  wri.writeheader()
  for row in summary:
   wri.writerow({k:json.dumps(v) if isinstance(v,(list,dict)) else v
                 for k,v in row.items()})
 excerpt={
  "latest_source_attempts":1043930,"four_hue":459235,
  "complete_climate":458174,
  "latest_fcp5_groups":35,
  "comparison":[{k:v for k,v in row.items() if k in (
   "axis","old_effect_white_percentage_points_per_SD",
   "new_effect_white_percentage_points_per_SD","new_nominal95_lower",
   "new_nominal95_upper","n_pairs_new","n_neighbourhoods_new",
   "new_raw_island_minus_mainland_pp",
   "new_spatial_climate_adjusted_island_minus_mainland_pp")}
   for row in summary]
 }
 print("MOIS_ACTUAL_LATEST1043930_FOUR_AXIS_NEW_VS_350399",
       json.dumps(excerpt,sort_keys=True,allow_nan=False),flush=True)
 return report

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for key in ("source","climate","water","elev-island","season","bias","prior","output"):
  p.add_argument("--"+key,type=Path,required=True)
 a=p.parse_args()
 main_latest(a.source,a.climate,a.water,a.elev_island,a.season,a.bias,a.prior,a.output)
if __name__=="__main__":main()
