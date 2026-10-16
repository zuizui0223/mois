# MOIS | Seasonality of white flowers and moisture-related environment

**Project home: `zuizui0223/mois`, 16 October 2026.** This is the active owner of seasonal photographed flower-colour × climatic seasonality, separate from the high-depth polymorphism paper in `fcp`.

## Focused ecological question

Across the world's georeferenced flowering photographs, does the share of *white versus chromatic flower* differ between seasons at **the same taxon and sampled locality**, and do these changes coincide with actual observation-month temperature, radiation and water availability? Does the direction persist after accounting for geography, photographer, source image classifiability and phylogenetic/taxonomic structure?

A photographed colour state is not evidence of a genetically segregating morph, within-individual colour change, pollinator selection or pigment energetic cost. The study is observational and explicitly not causal.

## Data cohorts must remain separate

1. **Prior global Mois interspecific manuscript:** 100,543 original photo opportunities across 42,111 taxa; 45,953 classifiable visible hue labels. Historical macro association source results are frozen; see `docs/MOIS_GLOBAL_MOISTURE_MANUSCRIPT_20261014.md`. This is an *old measured cohort*, not the current season-specific image work.
2. **New frozen 2026-09-27 source frame acquired by FCP:** 42,111 outcome-blind candidates, 93,209,351 eligible source photograph opportunities; a first species-stratified 1,000-taxon sample contains **3,430,502** globally deduplicated authentic photo IDs with verified original observed-on dates and 2,416,225 photos having independently recorded location uncertainty at most 5km. The first up-to-128 photo/taxon RGFCA *pre-colour* balanced inclusion comprises 106,376 IDs. These are not yet flower phenotypes.
3. **Actual first exact-original photo pixels:** FCP [run #38901816890](https://github.com/zuizui0223/fcp/actions/runs/38901816890) scored 64 original photos across *Cyperus glomeratus* and *Leptomeria drupacea*: 27 four-state ROI classifiable, but 4 include ≥10% white/yellow regions within the same ROI. Excluding such organ mixtures leaves 8 white and 15 chromatic **automated photo labels**; 37 are unclassifiable. These small-flowered/fruit-like images need independent verification of *open flower organs* before any primary white-rate denominator. They establish pipeline function, **not a seasonal ecological effect**.
4. **Separate MOIS temporal-source supplemental candidate:** source-only matching yielded 268 site×season-years in 39 taxa with ≥24 independent observations in each of two different local seasons, climate years through 2025. A further original-photo/observer-disjoint, pre-colour reservation selected *Asclepias syriaca* at one 0.25° geocell in season-year 2025, with 64 original photo IDs per season (32 first-wave + 32 reserve per season), **zero new seasonal-wave pixels opened**. Neither observation opportunity nor floral annotation alone proves the ROI shows open flowers.

Do **not** add cohort counts, compare 8/15 against the global historic white fraction, or report apparent month effects as if 106,376 or 3.43M photos were measured.

## Measurement, date and climate pipeline

- FCP supplies immutable original photo-ID/observed_on/observer/position records and SHA-sealed blind ROI-v4/EfficientSAM flower-colour measurements. Source-specific date, location and observer data never enter the pixel classifier.
- MOIS post-ROI first preserves *all attempts*: source photo, date status, independent observation, confidence of flowering-stage identification, four-state ROI classifiability, white/yellow within-organ mixture, and contradictions between multiple photos of one observation. Failed ROIs and absent iNaturalist phenology annotations remain **unknown**, not nonflowering.
- `scripts/acquisition/audit_fcp_original_rgfca_monthly_photo_opportunities.py`: source-only complete observed-date census. No colour output.
- `scripts/analysis/fcp_original_dated_flowering_season_qc.py`: separate real flower ROI attempts and source-date candidates; observation-weighted white fraction by month and coarse source geocell.
- `scripts/analysis/fcp_match_original_season_site_year.py`: co-located, same-species, same meteorological season-year matching; report naive unmatched fraction alongside matched estimates to expose season/geography mixing. **Five-degree cells are far too broad for a population or causal inference.**
- `scripts/analysis/fcp_prepare_observed_month_terraclimate.py`: independently SHA-linked post-ROI observed-month climate-cell request keys; only source-recorded fine geographic accuracy, no false 2026 weather imputation. The matched climate values themselves are **not yet retrieved or used in a fitted effect**.
- Next preregistration-in-project (not external preregistration): retrieve independently verifiable exact observation-month solar radiation, air temperature, precipitation, TerraClimate deficit/VPD, derive *within-site temporal anomalies* and stable soil/elevation/terrain adjustments, identify season variation within each species/cell/year and handle repeated observer/source-photo effort. Run observer-balanced, flower-stage verified controls before any ecological claim.

Climate-normal mean is not the actual weather of an observation month. Long-term soil/elevation/topography are between-site adjustments, not time-varying seasonal effects. No 2026 month can borrow a 2025 monthly field and masquerade as real 2026 values.

## Inference boundaries and authorship

This stream does not merge into the frozen GEB historical coefficient ledger until source identity, independent flower-organ classification and independent seasonal ecological estimates exist. The code/data copy is recorded in `docs/MOIS_SEASONAL_FCP_SOURCE_TRANSFER_20261016.json`. The original FCP historical files remain accessible for verification only; **new interpretation/analysis belongs to MOIS**.

The main existing `mois` interspecific GEB manuscript and its original mixed-source limitations remain unchanged. FCP's New Phytologist within-species W/D spatial tests also remain independent.
