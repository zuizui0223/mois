# Seasonal floral photographic colour: taxon turnover, apparent white enrichment, and the lack of local summer–winter contrasts

**MOIS research result, 17 October 2026.** This is the source-authenticated seasonality follow-up to the interspecific global moisture study, not the within-species photo-colour polymorphism estimand of `fcp`. The [fully measured 267,425-photo MOIS source artifact](https://github.com/zuizui0223/mois/actions/runs/39002957961) contains 129,795 automatic four-hue labels; 81,867 of these have an original observation date. **Neither automated hue classification nor an observation date establishes botanically verified open flowers.** The frozen source photo-CSV SHA256 is `fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9`.

## 1. The summer-white photographic excess is real as a photograph count, but is not fixed-species composition

White-labelled photographs among classifiable source photographs in the extra-tropical **Northern hemisphere**: summer 13,312/27,210 (48.9%), winter 1,481/3,246 (45.6%): **+3.30 percentage points**. **Southern hemisphere**: summer 1,627/3,030 (53.7%), winter 765/2,098 (36.5%): **+17.23 percentage points**. These are *photographic ROI labels*, subject to time-variable detection and source-selected taxa.

The first between-taxon finding is that the original-photo taxa appearing in one sampled season are not the same set in the other. Across all geographically available dated/classifiable Southern photos, **433 named taxa** are photographed in summer and **322 in winter**, but only **212** in both seasons (Jaccard **0.3904**). The corresponding Northern overlap is **450** of 1,622 union taxa (Jaccard **0.2774**).

Only **48.5%** of classified Southern *summer* photographs and **57.2%** of *winter* photographs come from those 212 common taxa. Thus the unpaired full-atlas summer/winter difference mixes taxon-frequency and season-exclusive photographed-taxon-support changes. These are **not proofs that a taxon biologically ceases flowering**; absence from citizen-science photographs is not an absence observation.

### Exact inclusive-source decomposition

Using the symmetric Kitagawa identity on named taxa actually photographed in *both* seasons, then adding algebraically exact corrections for the season-only photographed taxa, reconstructs the whole source-photo difference without assigning unobserved taxa an imaginary opposite-season white frequency. Report all parts so a common-species decomposition is not falsely treated as a decomposition of the full atlas.

| Source-photo summer − winter (percentage points) | North | South |
|---|---:|---:|
| Within common photographed taxon, white-label-rate component | +1.99 | +3.27 |
| Shift in *relative photo frequencies* of common taxa | −1.95 | +11.76 |
| Taxa observed in summer only: support correction | +3.70 | +4.00 |
| Taxa observed in winter only: support correction | −0.44 | −1.80 |
| **Reconstructed actual all-photo season contrast** | **+3.30** | **+17.23** |

This is *exact accounting in the photo census*, **not** unique causal attribution. The component for a taxon observed only in a season is a missing-support correction, not an estimate of its hue in the season when no photo exists. The [separate >=2 photo/taxon whole-taxon bootstrap decomposition](MOIS_REAL_267425_SEASONAL_WHITE_CLIMATE_NORMAL_DIAGNOSTIC_20261017.md) found the Southern common-taxon photo-frequency term +11.91 pp (exploratory **95% bootstrap percentile interval [+0.19,+23.64]**); it is modestly positive but sensitive to photo-count threshold.

## 2. Geographic and chronological positivity collapses

The same taxon across *the world* is not the same taxon within one co-located community or one observed year. A **0.25° lat/lon source cell** is already on the order of several tens of kilometres, not a known biological population. We independently enumerate complete original `taxon_id × coarse geography × local meteorological season-year` photographic support for both **summer and winter**. December is assigned to the *following* winter/summer season-year and Southern seasons are hemisphere-reversed.

| Matched original-photo design | North matched groups | South matched groups |
|---|---:|---:|
| Same named taxon globally, any site/year | 450 taxa | 212 taxa |
| Same taxon and same 0.25° cell, **any year**, at least one photo per season | 131 | 70 |
| Same taxon, same 0.25° cell **and season-year**, at least one photo per season | 27 | 17 |
| Same taxon, same 0.25° cell **and season-year**, ≥2 independent observation IDs per season | **1** (4 photos) | **0** |
| Same taxon and coarser **5° cell/year**, ≥2 observations per season | 26 (171 photos) | 13 (95 photos) |

The coarse 5° rescue is merely a **sensitivity of source opportunity**, not valid within-local-population inference: the grid can span hundreds of kilometres. The paired site × taxon × year photo sample is essentially absent at ecologically defensible geographic support. The separate model of *any two local seasons* rather than summer–winter uses 61 groups/47 taxa/384 source photos and remains too sparse and floral-stage-unverified to establish a universal within-site seasonal water/solar response.

## 3. Original-photo data generation and measurement regimes matter

The original *high-depth FCP photo source reused for MOIS* (`new_fcp`) shows Southern summer–winter automatically white photo difference **+18.80 pp** among 2,790 summer and 1,913 winter classified photographs. In the **additional small ROI-wave sources** (`fcp_additional_ROI`), the corresponding difference is **+0.52 pp**, with 240 summer and 185 winter labels. Both photo groups are **related source programmes**, not independent ecological replicas; the second has a different photo/taxon composition and incompletely executed FCP wave4 shards. The source heterogeneity itself argues against treating the broad global photographic summer gradient as a universally maintained pigment effect.

## 4. Ecological conclusion and publication direction

The evidence supports **seasonal turnover in the set and relative photographic observation frequency of named plant taxa**, associated with large changes in the *observed photographic* abundance of white labels, particularly in Southern extratropical source photographs. It does **not** yet distinguish biological flowering phenology from seasonal observer effort or photograph selection: those require photo-specific open-flower validation, effort data, credible absence windows and/or field floral surveys. The within-taxon label-rate component cannot be equated to plastic pigmentation or evolutionary loss/gain of floral pigments.

**Do not claim seasonal white flowers increase because of sunlight or water** from the current WorldClim 1970–2000 month-normal analyses. A true test would compare actual observation-year/month temperature, precipitation, radiation, atmospheric VPD and climatic-water-deficit anomalies among botanically verified flower photographs where *both* seasonal values are actually observed in the *same biological locality*. Static soil, elevation, terrain and lat/lon represent among-site baseline differences; they cannot act as manipulated seasonal exposures.

This is a productive negative identifiability result: the strong photographic seasonal pattern is quantifiable, its taxonomic turnover can be partitioned, and its missing local temporal counterfactual can be demonstrated without overinterpreting a species composition map.

### Reproducibility

- [Frozen numeric decision](../results/mois_267425_seasonal_taxon_spatial_support_decision_20261017.json) — source SHA, every key original-photo numerator, two source strata, exact photo-mixture component and spatial-support cutoffs.
- `scripts/analysis/mois_seasonal_spatial_identifiability_20261017.py` recreates an authenticated complete original photo/time/site × season support grid for 0.25°, 0.5°, 1°, 2° and 5° cells, both original season-years and pooled years, with minimum one or two distinct observations in both seasons.
- `.github/workflows/mois-267425-seasonal-spatial-identifiability.yml` downloads original *already measured photographic ROI* artifact, runs unit tests, checks actual source and frozen key values and uploads the entire site-scale table. It never calls FCP's within-species `D/W` spatial phenotype estimator.
