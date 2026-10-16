# Mois: seasonal photographed-white fractions versus changing climate, controlling stable geography

**Prespecified clarification — 16 October 2026.** This is a statistical design for the ONE global interspecific flower-colour paper, **not a new fitted ecological result**. The primary phenomenon is **seasonal replacement of flowering/photographed taxa and their visible hue composition**. Do not turn the main analysis into a test of whether one individual plant changes its petal colour between seasons.

## Main estimand: species-equal photographed-white among flowers observed in each place/month

For a fixed geographic analysis cell `c` and observed **local calendar month** `m`, let `S[c,m]` be the set of original-source taxa with at least one original image in the independently selected, successfully four-hue ROI-classified **visible flower** class. For each species `s`, first estimate `w[s,c,m]`, its sampled original photographic probability of white (between 0 and 1), including the original observation/photo inclusion probabilities and classifiability denominator. Then define:

```
W_species_equal[c,m] = mean_over_observed_flowering_species(w[s,c,m])
W_photo_weighted[c,m] = proportion of classified selected source photographs that are white,
                        with their ORIGINAL PHOTO inclusion/design weights.
```

**PRIMARY** = species-equal photographed-white fraction, with a separate white-photo share as a sensitivity and **all ROI unclassifiable/no flower photo denominator counts reported**. This avoids a highly photographed single species dominating the composition. However, an iNaturalist photographed and classified flower is **not a random census of all truly blooming wild taxa**; the paper must say `photographed blooming-source composition`, not worldwide wild-flora blossoming prevalence.

The estimand intentionally changes when **the species being photographed with flowers change by month**. A main regression that conditions on exact species or genus fixed effects would remove part of that biological community-composition question; use a separate genus-controlled or same-species matched sensitivity to answer whether the association remains beyond taxon turnover.

## Partition environment by temporal behaviour and spatial confounding

| Source variable | Temporal role | Role in main analysis |
|---|---|---|
| Observation month and year | Different each observation | Calendar/local-season labels, month and year controls, not substituting upload date or the later annotation date |
| WorldClim monthly temperature (`tavg`), monthly solar radiation (`srad`), monthly precipitation (`prec`) | Periodic climatological cycle in each cell | **Main dynamic covariates** with months aligned to original `observed_on` |
| Monthly water demand (VPD if computed from correctly dimensioned monthly temperature/vapour pressure) | Periodic annual climatology | Dynamic evaporative-demand candidate; independent standardized units and formula validation required |
| Latitude and longitude | Fixed per point, related to all long-term climate | **Spatial confounders and photoperiod context**. Do not interpret coordinate coefficients as biological mechanisms |
| Site elevation; terrain slope/aspect/roughness | Stable site attributes | Stable location controls **if genuine source DEM is available**; original photo-site elevation alone is not evidence of slope/aspect/TPI rasters |
| Soil pH, organic carbon, nitrogen, clay, available-water proxy | Slowly varying site attributes at study scale | Stable background **in across-location models**, optional effect modifiers of seasonally varying water supply |
| Original photo/observer selection, Flowers annotations, ROI classifiability | Changes by user, date, site and photographed taxon | Observation-ascertainment sensitivity and selected-photo weights; lack of annotation is NOT evidence of no flower |

The original 13-variable macroecological historical model is a **spatial cross-species baseline**. These monthly values form a clearly prespecified **dynamic extension**, not a replacement or fabricated extension of the original 13 coefficients. Monthly historic WorldClim 1970–2000 normals reflect the **usual** conditions at the observation site in month `m`, **not the specific weather in its observation year**. Absolute latitude, longitude sine/cosine and stable soil/terrain values do not need to be repeated as separate causal season covariates in a within-cell specification.

## Model A — MAIN: same cell, across seasons (stable geography controlled)

For each dynamic climate axis `X` and original geographic cell `c`, compute `mean12(X[c,1..12])` across **all twelve source climatology months**, NEVER as a mean over observer-submitted photo dates (because the observation calendar is biased). Define `deltaX[c,m] = X[c,m] - mean12(X[c,1..12])`, independently for monthly temperature, solar radiation, precipitation and water demand. A possible species-/observation-weighted hierarchical model is:

```
logit Pr(photographed four-hue-classifiable flower is white at c,m)
   = alpha[cell=c]
   + delta_calendar_month[month=m, local hemisphere context]
   + beta_temp    * delta_tavg[c,m]
   + beta_solar   * delta_srad[c,m]
   + beta_water   * delta_prec[c,m]
   + beta_demand  * delta_vpd[c,m]
   + observation-year / source-observer uncertainty structure
```

**`alpha[cell]` absorbs stable within-cell geography**, including coarse-grained longitude, latitude, long-term soil, mean elevation, terrain and other time-invariant place differences. Accordingly **no separate coefficient for soil/elevation/latitude/longitude is identifiable simultaneously with exact cell fixed effects** if that variable is constant within the cell. Within larger 0.5°–1° analysis cells, fine-scale soil and terrain may remain heterogeneous; coarser cell-fixed effects only control the *cell-level* average and need a sensitivity to the grid resolution. Use two-level/original-site as possible; **site/date support must be checked before selecting grid size**. Varying observer effort, seasonal classifier quality and unobserved blooming taxa are *not* removed by cell fixed effects.

An additional month indicator absorbs synchronized seasonal effects (global phenological and camera/calendar changes) and makes climate coefficients depend on **different seasonality profiles across cells**; compare **with and without** hemisphere-specific month indicators, rather than treating the climate slope of either specification as causal. If a cell has only one sampled month, it supplies no within-cell seasonal information and must not be counted as a matched-cell seasonal estimate. Report the number of shared cells and independently photographed species supporting each month comparison.

## Model B — secondary: across geography and between-cell climate gradient

For the same source photo sample (no new hue filters), separately examine differences **between cells** using `mean12(X[c,*])`, the site-average annual precipitation and temperature/radiation, genuinely available terrain/soil background and a **regularized smooth geographic field** `f(latitude, longitude)` / spatial random effect. This asks which places have high observed white-flower composition, not what happens when one place moves from its dry to wet time of year. Do not claim an independently identified climatic `beta` if broad annual rainfall and latitude/spatial smoothness are strongly confounded; quantify predictor collinearity, available support and held-out geographic/cross-genus predictive gain. Report the original static 13-variable model in its old retrospective cohort with its own denominator, **never directly compare R² across different samples**.

**Do NOT put `alpha[cell]`, latitude/longitude, local soil/elevation/terrain and spatial smooth terms indiscriminately into one enormous model**, then attribute any coefficient that survives to mechanism. Exact cell fixed effects already consume cell-level static predictors; a competing between-cell association needs the separate model.

## Dynamic climate is also mutually correlated

Hot months may be sunny and dry; rainfall, vapour pressure and temperature frequently covary by climate regime. Report (a) separate one-climate-axis + fixed geographic/control models, (b) a regularized joint dynamic T/S/P/VPD model with predeclared scale and collinearity/rank diagnostics, and (c) held-out **geographic cells** and **genera** where adequate data permit. If temperature/radiation change together everywhere, their partial effects may be indistinguishable, even with cell fixed effects. A predeclared water-demand/vapour-pressure proxy is not a free fourth independent treatment.

Soil moisture storage or topographic aspect can **moderate** the seasonally variable climate exposure (e.g., soil available-water × dry month; slope aspect × sunlight): an interaction in an otherwise within-cell analysis may be estimable even if the stable main effect is absorbed, but it needs strong cell×season support and is **secondary**, not scanned to rescue weak effects.

## Two required answers, NOT one conflated number

1. **Observed ecological composition:** Across different local months, does the **species-equal fraction of classified photographed blooming taxa that are white** vary with same-cell monthly sunlight, temperature and water supply? The community's constituent flowering taxa MAY turn over by month; that is the object of study.
2. **Controlled taxonomic sensitivity:** Do corresponding photographic colour associations survive within identical genus/species sets in cell/month matches, with separate photo-classification and observer biases? This is conditional taxonomic association, **not individual-flower plasticity**.

Both analyses need the true (not yet complete) original FCP 42,111-taxon 93,209,351-photo-ID frozen source, an independently selected Mois v1/v2 sample with ORIGINAL observation dates/weights, unchanged ROI-v4 source photo-hue measurement, and original-location monthly temperature, solar, precipitation and vapour-demand profiles. Existing 12-month WorldClim **precipitation** archive has been source-hash verified; production monthly tavg/srad/vapr original ZIPs and their raster unit/resolution calibration have **not** been shown to be acquired or reconciled. Observation-only current iNaturalist Flowers annotations are incomplete (~6.20% of photographed Angiospermae observed during 2025 carried a positive Flowers annotation in an independently measured API audit); they may prioritize processing but are not a complete flowering record. Hence no new seasonal white–solar/temperature/water empirical coefficient is claimed here.


## Original 48-raster T/S/P/vapour admission and photo-month proof (16 October 2026)

**Source availability has advanced from hypothetical to verified originals.** The official WorldClim 2.1 historical monthly **10-minute, 1970–2000** raster ZIPs for temperature, solar radiation, precipitation and water vapour pressure were acquired and **all 48 original monthly GeoTIFFs verified together on exactly the same WGS84 north-up grid**, including source units. The full [independent archive hash/units source receipt](../results/mois_worldclim21_dynamic_48_real_source_receipt_20261016.json) is grounded in the [successful real source Actions run](https://github.com/zuizui0223/mois/actions/runs/38848205285): ZIP SHA256 tavg 5e567dcfe94379b94229492849ce91078b1c6e5210aaf435fba449fae6b95405, srad c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026, prec 1090f578f6402672a2b6438e252156f6cf867a25002a8ad6ffcbfa19e516ddc6, vapr 0501b5b5cbbc6e2fc2a3801efa9a126a5f9669f0f9a00642d2519c2324ce726d. Combined 48-layer GeoTIFF manifest SHA256 **9d0cf37a64534542d7788debef28c13f106a15c3bf7bf91be0653f3ce106baae**. These are **real external climate normal pixels**, not flower pixels.

[Executable source-bound photo-site/month join](../scripts/analysis/mois_join_monthly_dynamic_climate_20261016.py) reads the original sample's selected observation date and protected original source photo coordinates, extracts the corresponding four original month-grid climate values, and derives **site-specific deviation from the twelve-month climate-normal mean**. The original source units are °C, kJ m⁻² day⁻¹, mm/month and kPa; the solar values can be converted to MJ m⁻² day⁻¹ by dividing by 1,000 *after* source validation. A diagnostic approximate vapour-pressure deficit is calculated as

\`\`\`
sat_vapour_pressure_kPa = 0.6108 * exp(17.27*Tavg_C / (Tavg_C + 237.3))
monthly_VPD_kPa ≈ max(0, sat_vapour_pressure_kPa - monthly_vapr_kPa)
\`\`\`

This is **VPD estimated from MONTHLY mean air temperature and vapour pressure**; by nonlinear averaging it is not identical to the daily/instantaneous VPD actually experienced by the petals or observers. Months without exact observation dates, source climate cells with missing values, site identity conflicts and source raster SHA/unit drift stop or remain explicitly missing; nonclassifiable photos remain nonclassifiable.

### Real original FCP 96-photo technical pilot, not a new worldwide floral conclusion

As a **technical validation only**, a separate 96-original-image FCP ROI-v4 artifact (four selected species, SHA256 of already remeasured ROI terminal CSV **911ec9edd10a7126bee28cfa67738cef0ff4d275b1f0ca1c06e3f26e7e290039**) was recovered directly from [its FCP original image repair artifact](https://github.com/zuizui0223/fcp/actions/runs/38673995478). The independently [source-validated Mois photo-month×48-climate join run](https://github.com/zuizui0223/mois/actions/runs/38848886635) extracted all four environmental variables for **all 96 original source photographs**, with **20 classifiable four-hue photographs and 3 white**. Northern extratropical local spring had 7 classified/1 white (14.3%); summer 9/1 (11.1%); winter 3/1 (33.3%); the remaining one classified image occurred in a tropical calendar Q4, with no white. Original photo coordinates, observation IDs and photographer identities were **not deposited** in the Mois result. [Derived site-season climatology aggregate receipt](../results/mois_actual_original96_photo_monthly_tavg_srad_prec_vpd_technical_pilot_20261016.json) contains only the small climatic group summaries.

The North seasonal **MEANS among classified-photo source locations**, in order spring/summer/winter, were mean air T 13.54/17.58/10.17 °C; daily solar input 22.40/26.03/10.83 MJ m⁻² day⁻¹; monthly precipitation 28.43/9.00/87.00 mm/month; and approximate monthly VPD 0.731/1.146/0.522 kPa. This provides a **real test of month+location climate extraction from previously measured original photo-pixel hue classes**, but there are just 20 colours, with different sampled species/locations and severe photo classification nonresponse. **No worldwide original 42,111-taxon season×climate coefficient or statistically supported effect exists on these 96 photos.** The technical pilot is *not* a probability sample of the complete original-source FCP photo population and must not be combined as an independent ecological replication with the old historical 100,543-photograph analysis.

### Primary source-ready species-equal cell-fixed effect estimator

The independent [species-first geographic-cell × observation-month climate FE estimator](../scripts/analysis/mois_monthly_white_species_equal_cell_FE_20261016.py) is now executable. It uses **one species vote per source climate pixel and observation month** after within-species original-observation-probability weighting; a separately identified pixel-weighted variant will be a sensitivity. It requires the fully **SHA-sealed** Mois sampled photo source and exact 48-raster original climate-site month join, retains all failed photo classifier results as missing outcomes (no white/coloured reclassification), and aggregates true original site/month climate values. It fits **cell-only** and **cell PLUS shared calendar-month fixed effects** using sparse weighted Frisch–Waugh–Lovell demeaning, never loading dense millions-of-cells dummy matrices.

A test deliberately injects a common calendar-month effect into synthetic photographed-white fractions, demonstrating why a spatial-only model cannot isolate conditional T/S/P/VPD partial associations and why **adding month fixed effects is a different question**. If too few geographic cells have classified photographs in multiple months, or the four dynamic axes are linearly dependent, the tool returns explicit **NOT_ESTIMABLE** status instead of zero coefficients. It does **not** currently compute genus/observer/cell-cluster confidence intervals, heldout-region predictive comparison, or other causal identification: any finite point coefficient is marked EXPLORATORY and needs further robustness before a manuscript result. In the actual technical 96-photo pilot, **60 different 10-minute climate grid cells** contain photographed source opportunities and **11 of them** were photographed in two or more observation months. After the original image ROI-v4 gate, only **16 cells** contain a classifiable flower colour and just **3 cells** have classifiable flowers in two or more months. Consequently the joint four-dynamic-axis geographic-cell plus calendar-month fixed-effects white-colour model is **NOT ESTIMABLE** on this tiny technical pilot; that does not imply an absence of a climate relation. These are the verified original [source-measured monthly pilot receipt](../results/mois_actual_original96_photo_monthly_tavg_srad_prec_vpd_technical_pilot_20261016.json) and [successful source calculation run](https://github.com/zuizui0223/mois/actions/runs/38849479497), not a claim about seasonal support in the future full 42,111-taxon source.

The full Mois original study still requires FCP's durable all-42,111-taxon source-ID freezer and substantially more classifiable source original photographs before this fully source-bound seasonal estimator can be run for the main study. Only climatological source raster inputs, the small real technical photo-month join and synthetic estimator tests have actually been executed.
