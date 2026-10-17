# MOIS: increase the existing global floral-colour dataset by reusing FCP original photographs

**Active research question:** global *between-species* photographed flower-colour composition and its biogeographic, geological, static and time-varying environmental associations. Original photos may overlap with FCP. **Only FCP's within-species D/W/50-km estimand is excluded**, not the image sources. The frozen original MOIS GEB source analysis remains unchanged.

## Source-level increase without reclassifying existing photos

The historical MOIS photographic source has **100,543** unique photo IDs, with **45,953** classified four-state hues and 54,590 no definite colour. FCP has **149,900** *already measured* additional original photographs from three species-disjoint high-depth cohorts (500×100 + 500×100 + 499×100), with **3,411 literal photo IDs shared across both studies**. The merged original-photo-ID universe therefore contains **247,032** IDs, of which **146,489 are non-overlapping with the old MOIS atlas**. The number of added *classifiable* hue photographs is derived afresh from the real signed files in the GitHub Actions job; it is not assumed to equal 146,489.

The source-sharing **is intentional and allowed**: photo labels already measured in the fixed ROI system are reused, without computing FCP's intraspecific spatial diversity D, nine-colour W or RGFCA hotspot probability in MOIS. Same source means it is not an external independent biological replication.

The actual photo-level materialized output contains each original photo exactly once, identifies whether its measured colour came from the old MOIS archive or FCP, and carries original latitude/longitude, source photo date (when available), existing WorldClim annual temperature/rain/solar proxies and SoilGrids-derived soil pH, organic carbon, N, clay and available-water proxy **where those source columns were already present**. No missing topography or new monthly climate is filled with proxy averages. Unclassifiable and historical blank original colours remain missing, not discarded or called non-white.

The output also contains taxon-wise four-colour counts and the source-specific number of classifications so **high-depth FCP species do not silently dominate** the original species/genus-balanced MOIS estimates.

## Static versus time-dependent ecological predictors

| Scientific block | Real variables | Rule |
|---|---|---|
| Spatial/static coordinates | latitude, longitude, elevation | Control broad geographic sampling and altitudinal climate covariation; elevation does not change by flower season |
| Relatively stable substrate | soil pH, soil organic carbon, soil N, soil clay, available-water capacity | Keep fixed within an original mapped site; classify this as site background, not a monthly rainfall |
| Terrain | slope, aspect, terrain ruggedness, topographic position | Add only after source-pinned DEM raster derivation and an actual photo-site join; none currently in the original MOIS archive |
| Long-term climate | original `wc_bio1/bio5/bio12/bio15` | Climatic normals for comparisons among sites, **not** the weather in the observation month |
| Actual observation-month climate | monthly temperature/tmin/tmax, solar radiation, precipitation, deficit, VPD | Require an original observation date, site coordinates, real-year TerraClimate data and exact source CF units; compute within-site same-calendar-month 1991–2020 anomalies separately from site averages |
| True local clock hour | hour-of-day / daypart | **Not in the original 100,543 archive or the frozen 3,430,502-source `observed_on` snapshot.** Leave unknown unless an authenticated `time_observed_at` or equivalent original timestamp/timezone is recovered. Month, calendar season and daypart are different axes. |

### Nonoverlapping MOIS estimand

For photographed flower white versus chromatic composition, a planned enlarged MOIS analysis will compare:
1. Among-taxon/genus environmental gradients after taxonomic and geographic standardization (the main global macroecological question).
2. Geographically and taxonomically matched observational sampling across observed months/seasons, with changing climate separate from longitude/latitude/elevation/site soil and terrain.
3. Source-status controls (historical shallow sampling vs FCP species-concentrated high depth), observer sampling, missingness, and source-identity duplicates. This is *not* within-species polymorphism-maintenance, FCP D/W photo-pair diversity, or a fitness experiment.

The first runnable product in this PR is the **actual 247,032-source-photo materialization and source-faithful four-colour composition table**, not an unsupported claim that all new photos have independent biological flower-stage verification or month-specific climate. The subsequent signed photo-metadata and gridded covariate join runs only on dates and sources that genuinely exist.

## Reproduction

`.github/workflows/mois-measured-original-247032-photo-union.yml` fetches the three exact, original measured FCP CSVs by immutable Git commits and validates their SHA256, fetches the original MOIS 100,543-unique photo archive by immutable GitHub Actions artifact ID and SHA256, computes the union, and publishes the **real compressed photo-level hue data** together with taxon tables, photographic classifiability rates and complete source receipts. It rejects any change in the known 3,411 shared photo IDs or 1,693 concordant overlapping classifiable labels, and does **not** alter historical macro-models.
