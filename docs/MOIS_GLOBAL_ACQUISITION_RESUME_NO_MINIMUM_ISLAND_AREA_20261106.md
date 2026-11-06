# MOIS global environmental acquisition resume — 2026-11-06

## Scope authorized by the user

Continue data acquisition alongside the analysis. MOIS includes continental land AND islands worldwide; it is not an island-only clone of the island repository. There is no 5 km2 minimum landmass area or +/-60-degree acquisition cutoff. Landmass class, size and latitude are retained as attributes, not selection thresholds. Preserve prior source-arm distinctions and original strict9/technical unknown states.

## Existing data recovered, not newly measured

The user-supplied reproducibility ZIP `MOIS_global_all_land_no_5km2_cutoff_20261106.zip` has SHA256 `5f340bd1e2ec932ad221a42059b21e7f5bb9866c24694a870b53449fc29676c5`. Its global queue `GLOBAL_no_area_cutoff_environment_queue.csv.gz` has the recorded SHA256 `d8a2bb2152a334706588a9434f777dc4d3310ae911e28bd1ac720d028c0f5702`. The previous independent verification reports 513,339 location-year-landmass targets reconciling 1,032,685 source observations with usable coordinates/year, plus 197 unresolved source observations retained outside the request queue.

The previous collector attempted one request and stopped on local NETWORK_UNAVAILABLE, with zero new numerical values and zero new photos. On this resume the local runtime still cannot resolve the two external ArcGIS hosts. GitHub Actions is being tested as the execution route; writing this note alone does not constitute successful data collection.

## Bounded first acquisition and checks

Use unused request keys from the existing queue, with deterministic target selection independent of flower colour and classification success. Start with 2020, for which the existing WorldPop service has same-year coverage. Each target requests its original January/April/July/October VIIRS values and 2020 population density, with original response bytes and SHA recorded. Four months are not an official annual VIIRS composite. A target is complete only when explicit quality checks pass; planned targets are never counted as acquired values.

Successful or terminal no-data requests must be cached by product, time and coordinates. On authentication errors, provider limits or network failures, checkpoint and stop instead of repeating large jobs. A failed image worker is not restarted: this scope acquires environmental numerical data only, no original photos or ROI inference. Later population years require a separately validated same-year product, never implicit substitution of 2020 values.

No main-branch edits, merges, source ROI relabelling or discarded MNAR sensitivity analyses are authorized by this scope. Keep PR #199 draft.
