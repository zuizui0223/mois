# MOIS × island ALAN: ROI-free source first, pollinator-redundancy interaction second (2026-11-06)

## Primary question, not a conclusion

**Does island setting modify the association of artificial night-light radiance with the photographed plant assemblage, especially among species with independent evidence of moth visitation?**

Island pollinator redundancy is a hypothesis. Not all islands lack bees, not all moth visitation is effective nocturnal pollination, and white flowers are not automatically moth dependent. Mainland non-association in PR #199 neither proves nor rules out an island interaction. The old positive 47-pair direct-colour contrast is retained as a genuine but selected mainland-dominated observational result.

## Executed technical advance and scope

- Source-complete **1,032,882** canonical original observation rows, original taxonomy and coordinates; use original **historical_main** as a source-consistent first cohort without selecting on strict9 automatic ROI success.
- One accepted WCVP name per source taxon via frozen exact-match source; unresolved taxonomy stays unknown.
- Independently documented **High/Medium, resolved, single-state** original Island species flower colour; unresolved/mixed colour remains missing, not nonwhite.
- **61,967** GloBI interaction records are flower-visit *positive* evidence; no GloBI record is **not** negative pollinator evidence. Shared moth/butterfly visitors must not be counted as independent comparison species.
- 20+ km2 source GSHHG-style island polygon support; only exact geographic containment; overlapping island boundaries excluded, source geometry provenance SHA recorded. Original Island geometry source preferred; a new Natural Earth fallback must be explicitly labeled as a different candidate universe.
- The photo and year denominator includes **2012–2024** where available, without attributing dates/colour to rows lacking evidence.
- Existing true same-year 2020 and 2024 **four-month VIIRS** cells are support/contrast diagnostics only; they are not an EOG annual VNL panel, and earlier/later photo years are not assigned 2020 values.
- A future all-year annual VIIRS CSV may join only by identical photo cell and year and must include coverage and annual-product version; EOG VNL V2.2 is documented for 2012–2024, with 2012 partial-year coverage.

Code: `scripts/analysis/mois_alan_island_roi_free_support_20261106.py`. Tests: `tests/test_mois_alan_island_roi_free_support_20261106.py`. Authentic-source run: `.github/workflows/mois-ALAN-island-ROI-free-source-support.yml`.

The public result is island×year and calendar-year **aggregate support counts** only (original photos, distinct observation/species, documented flower colour, positive moth evidence, intersection, 2020/2024 light cell range). Suppress small direct-colour strata. Never redistribute third-party species-specific Island colour labels, photo coordinates or article source data.

## Native-status gate

The official WCVP TDWG3 native designation is **too coarse** to declare an observation of a plant native on a specific oceanic island. A source-reviewed `(island_id, accepted_species, establishment_status, review_status=accepted)` data product is required for a native-only PRIMARY effect. Until supplied, island-native and island-introduced count endpoints are **unresolved**, not zero and not globally estimated. The optional exact-status CLI joins only reviewed accepted records.

## What a support-passing follow-up would need

1. Native flora and introduced flora modelled **separately**; within-island ALAN range and independently matched elevation, land cover, population, photo effort and year. Use island random effects or matched within-island gradients; no pretend island-level independent photo replicates.
2. Real EOG same-year annual VIIRS and cloud-free coverage for geographic photo cells; do not splice prior four-month values. VIIRS is upward radiance, not measured flower or moth spectral illumination.
3. Additional iNaturalist **insect** observation data are a different sampling frame from the 93m original angiosperm-photo pool. Compare moth observation rates only with observation effort/season and insect recording changes; absent photos are NOT evidence of moth decline or floral extirpation.
4. Independent night-opening/effective pollen deposition/seed-set evidence would be needed before upgrading from recorded moth *visits* to biological **dependence** or a mechanistic mediation claim.

### Explicit estimands

- **Now:** sample support and availability. What fraction of exact island original plant photo observations can be linked to a documented accepted-species flower colour, GloBI positive moth flower visit, native establishment records and same-year satellite radiance?
- **Later if supported:** island–mainland effect modification in observational photo-assemblage composition after geography, sampling and time adjustment, with separately reported native and introduced strata.
- **NOT inferred:** species loss, extinction, effective nocturnal pollination, causal ALAN effects, or a globally representative island native-flora white rate.

No new original photo pixels, original ROI inference, flower colour imputation, MNAR scenarios, or main-branch merging in this stage.

## Actual island-source audit — 2026-11-06 (completed)

**Original successful run:** https://github.com/zuizui0223/mois/actions/runs/42171173575
(aggregate-only artifact ID 13481850056). All nine source-island unit tests, exact historical numeric source SHA checks, GSHHG geometry validation, and the final two CSV byte-SHA checks passed. The original Island alpha1 polygon artifact was no longer retrievable, so **GSHHG 2.3.7 source polygons were rebuilt independently**: 8,260 raw island units, 20 km2 analysis threshold. Geometry SHA is in the private-to-public aggregate receipt. This rebuild is NOT assumed byte-identical to a prior 8,264-unit corrected Island analysis.

- Historical-main 2012–2024 eligible original observed photo frame (before island polygon): **524,661**. Unlike the old strict9 outcome this selection does not condition on pixel ROI success.
- Exact 20+ km2 island polygons contain **61,734** genuine photo/observation IDs across **889** island units and **3,575 island×year** rows.
- Original strict9 successes within those island images: **21,954** (diagnostic ONLY). All-photo eligibility is **2.812×** this ROI-success subset, but independent documented single-state colour is not available for every added photo.
- Independently documented High/Medium single-state flower-colour matching: **28,024 photos**. Positive GloBI moth family flower-visit evidence: **14,041 photos**. **5,328 photos** have both separate species-level evidence sources.
- Islands with >=20 total original images: **212**; >=100: **66**; >=5 independent documented-colour photos: **288**; >=5 positive moth-visit AND colour photos: **103**. Islands with 3+ known-colour accepted species in at least one year: **261**. Islands with 3+ moth-visited accepted species in at least one year: **206**. At least two photo years: **622 islands**.

**Independent source-only power gate**: https://github.com/zuizui0223/mois/actions/runs/42171461128
(artifact ID 13481685619), PASS with all original aggregate output SHA checks.

- Existing **quarter-month**, not annual, same-year VIIRS support covers only **19 islands / 25 island×year** rows. At least three distinct source VIIRS cell exposures: **10 island×year** (all ten have >=3 independently coloured species and >=5 positive-moth+colour photos).
- Within those 10 viable photo island-years, median observed within-island **log1p(radiance) 90th–10th percentile contrast = 2.2404**, maximum **3.2807**. These are raw quarterly-proxy spatial exposure ranges, **not annual EOG VNL brightness**, not population/land-cover matched and not an ALAN effect.
- Because 2020/2024 four-month VIIRS source positions were previously selected mainly for mainland 0.05° cells, **the light-coverage bottleneck does not imply island botanical support is absent**. It is primarily an incomplete spatial light exposure join.
- There is currently NO independently accepted island-specific native/introduced source table in this run: all such local establishment statuses are unknown. The present counts are not a native-only analysis, and neither the 2020/2024 light model nor the proposed moth-mediated island interaction was fitted.

**Decision:** continue island-first rather than stop. The immediate next data acquisition bottleneck is versioned **2012–2024 EOG VNL V2.2** same-year island-cell annual radiance with valid coverage; that should be followed by audited island-specific establishment status and a separate effort-controlled iNaturalist moth observation frame. Do not mislabel positive GloBI moth visits as nocturnal dependence or interpret missing photo records as plant extinction. No original pixels acquired and prior PR #199 results retained.
