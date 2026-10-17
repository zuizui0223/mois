# MOIS primary Article: environmental heterogeneity behind common white floral colour

**18 October 2026 — independently sourced exploratory ecological analysis specification.**

## Main question

Worldwide original flower-like photographic ROIs have a visibly white category
nearly as large as the three other visible hue categories combined. In the
actually verified 313,535 source-photo-pixel dataset, **71,427 / 151,239
classifiable ROIs = 47.23% photographic white**. This is not a fixed adaptive
50:50 equilibrium, a systematic compensation mechanism across places, or a
measured metabolic or reproductive pigment trade-off.

**Test:** Does the relative prevalence of photographed WHITE versus chromatic
hues differ between water-limited, high-light, hot and moist regimes? Do these
photographic gradients survive source geography, nominal genus, DEM terrain and
soil adjustments? Is the near-half global proportion a mixture of distinctly
more coloured and more white sampled ecological contexts rather than a
universal colour proportion?

A visible-white flower may contain UV-absorbing flavonoids and other pigments.
Yellow/orange does not equal red/pink or blue/purple; the global other-colour
pool is **not** an anthocyanin-specific outcome. Pigment production costs,
pollinator benefits, optical properties, herbivory, natural selection and
flower-specific fitness were not measured here, and are therefore competing
hypotheses rather than estimated biological pathways.

## Already physically measured source

Read only the source-preserving GitHub Release:
https://github.com/zuizui0223/mois/releases/tag/mois-primary-313535-fullfourth4096-featurelake-v2
containing independently signed three original photographed photo-ID-aligned
Zstd Parquet datasets: all actual ROI-v4 four-hue classified/unknown source
attempts; original site long-term BIO5/BIO12/BIO14/BIO15/annual-solar
environment; partially covered historical soil and DEM terrain data. Verify
the complete original photo, climate, soil/terrain SHA256 fingerprints and
actual 313,535 photo attempts /151,239 four-hue labels /151,230 dated
four-hue observations before treating the data as source-measured colours.
The **separately season-gap-targeted 256-photo cohort is explicitly excluded**.

**Time:** The annual WorldClim 2.1 variables correspond to 1970–2000 climate
normals, not weather observed during the source photograph's date. The
WorldClim solar field is already converted from kJ m^-2 day^-1 to W m^-2.
Seasonal fluctuations / actual observed-year months are Supplementary
Information (SI) and **do not require a plant to change its colour within
a flowering season**.

## Three genuinely informative analyses

1. **Rainfall × irradiance and rainfall × heat environmental partition.**
   Define LOW/MIDDLE/HIGH rainfall, annual radiation and warmest-month heat
   tertiles *from all attempted photograph sites*, regardless of whether the
   four-hue classifier succeeded. In two independent 3×3 comparisons
   (rain×solar and rain×heat), record actual attempted photographs, classifier
   successes/unknowns, observed white and each other hue, photo-weighted,
   equal nominal-taxon and equal-genus photographed white rates. No category
   is equated to biochemical pigment assay.
2. **Latitude/solar partially conditioned climate direction.** Within
   hemisphere ×10° absolute latitude ×solar tercile, compare WET and DRY
   photo-white fractions only when each has ≥80 classifiable original ROIs
   and ≥20 independently labelled source taxa. Equal-stratum mean and 999
   exploratory resamples describe geographic/photo mixture variability.
   This is **not longitude/site/species matching**, not a causal drought
   response, and not an out-of-source confidence family.
3. **Taxon/genus-centered environment and static site blocks.** Aggregate
   actual original white photographed fraction for named source taxa, retain
   nominal genera with ≥2 taxa, standardise geography and annual predictors
   on matched source taxa and fit equal-genus-weighted genus-centered WLS
   with genus-cluster nominal intervals. Compare geography+annual climate
   against the same support including standardised rain×solar and
   rain×heat terms (BH q across two exploratory interaction contrasts).
   On two **separate complete-case matched photo populations**, compare
   baseline geography+climate against added (a) DEM elevation/slope/TRI,
   (b) original soil pH/SOC/N/clay/available-water proxy. **Only subtract
   R² within one matched block**, because soil/terrain have unequal and
   incomplete support. This is in-sample association, not geographically
   transferable predictive validation or an experiment.

Reference freeze: prior global-original 309,439 source analysis estimated
photo-white and annual precipitation +1.99 percentage points per
environment SD, genus-wide maxT95 [+0.55,+3.42], and atmospheric vapour
+3.16 [+1.04,+5.28]; both lost sign confirmation under coarse geographic
conditioning. The new analyses are exploratory and **not** a second
independent bootstrap or a correction for all project-history comparisons.

## Reproducible outputs

Source code: scripts/analysis/mois_global313535_water_light_heat_habitat_partition.py

Workflow: .github/workflows/mois-global313535-white-water-light-heat-habitat.yml

The source-verified run produces four CSVs: rainfall×solar and rainfall×heat
composition, latitude/solar-stratified wet–dry contrasts, multivariate
nominal-genus model terms and within-cohort static-block ΔR². A JSON contains
original photo-source SHA256, sample denominators, 3×3 cutoffs, 999-stratum
bootstrap summaries, model-support counts and limitations. Archive data
derivatives as a new versioned immutable Release; leave original copyrighted
photo pixels and unmeasured labels out of this result.

**Publication priority:** foreground source-wide white/common visible-colour
prevalence as a surprising descriptive phenomenon and its heterogeneity
between environments. Water, irradiance, temperature, soil and terrain
are explanatory *candidates*. The energetic/optical/pollination
cost–benefit balance is a subsequent mechanistic hypothesis rather than
the tested and demonstrated explanation.

Seasonal/monthly global photo-colour environment contrasts belong in SI.
Main manuscript: docs/MOIS_PRIMARY_NEAR_HALF_WHITE_CHROMATIC_ENVIRONMENT_ARTICLE_20261017.md.
