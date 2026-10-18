# MOIS follow-up: does "temperature is strongest" survive incomplete soil sampling?

The published source-verified 313,535-original-photo 33-variable PCA/block
ablation report selected **temperature** as a *tentatively strongest predictive
group* when both five-fold held-out-genera and held-out-10-degree-geography
rankings agreed. However, full 33-variable complete cases comprise only
**30,109 measured classifiable photographs**, a small nonrandom part of all
**151,239 classified photographic four-hue ROIs**. Soil was imported for only
a portion of the original source; no missing pixels or soil values were invented.

This follow-up tests **whether the apparent temperature prominence is
a property of the soil-complete source subset** rather than something
reproducible across the original full photographic atlas.

## Cohort-conditional planned/implemented comparison

On the identical frozen source photo IDs and their full original actual
WorldClim v2.1 19 BIO + annual VAPR + annual solar and DEM data, compare:

1. All actual **climate-complete** original classified photos:
   11 temperature, 9 moisture, 1 annual solar variable.
2. **Terrain-complete** original photos: the same 21 climate variables
   plus elevation, slope, roughness, TRI, TPI, eastness, northness.
3. The original **all-33 complete-case photo set** (soil and terrain
   available), but fit ONLY the same 21 climatic predictors, with
   identical grouped genus/source geographic treatment.
4. The already source-signed **all-33** five-block model on the exact
   **same** case subset as (3), published here:
   https://github.com/zuizui0223/mois/releases/tag/mois-313535-all19-BIO-vapr-PCA-ablation-v1

**Important:** compare the climate ranking between (1) and (3) to
diagnose soil-related selection, and (3) against (4) to see whether
conditioning on soil/terrain changes climatic predictive contributions.

The model target is each named source taxon's actual original
photographed WHITE ROI fraction among classifiable visible-hue ROIs,
not actual population white-flower abundance, chemical pigmentation,
anthesis, reproductive fitness or biological balancing selection.
Reject imputation for unmeasured soil; repeat PCA and StandardScaler
fit separately on each training fold only; use fivefold grouped
held-out nominal genera and approximate 10° taxon centroid geographic
blocks, with equal genus weights. The score is held-out weighted R².
For each conditioning group, compare the full selected block model
to the same model with that group removed on the **same source taxa**.

The source climate normals are **1970–2000** long-term annual WorldClim,
not weather measured in photo month/year. Monthly variation stays SI.
The result cannot indicate that a flowering plant changes petal colour
within its lifetime or that visible-white petals have no UV pigments.
No photo IDs from targeted supplementary seasonal enrichment enter the main
global sampling estimator.

## Selection and inference limits

- Source soil measured/retrieved only for some original photos can
  select a different taxonomic/geographic population. Testing all
  predictors on the soil-complete subset is a valid matched comparison
  but is **not** a universal world flora result.
- PCA component count is chosen for **environmental variance**, not
  white-hue significance. The train-fold-only PCA prevents leakage but
  cannot identify a causal abiotic driver.
- BIO indices are derived from climate and correlated; dropping one
  group or predictor measures *unique conditional predictive
  contribution*, which can be small or negative when information is
  redundant. Rankings may differ by cohort and generalization target.
- Holdout taxa and geographic 10° blocks are approximations, not
  exhaustive independent field sampling of flowering populations.
- Larger data collection remains desirable but must increase genuine
  original photo RGB pixel classifications and environmental coverage,
  not merely append provider's 3.43 million metadata photo IDs.

Source code:
scripts/analysis/mois_global313535_climate_subset_sampling_robustness.py

Workflow:
.github/workflows/mois-313535-climate-only-vs-soil-sampled-importance.yml

Outputs: SHA-checked three source-cohort support summaries, held-out
genus/geographic climatic group ablation rankings, same-photo all33
versus climate-only comparison, interpretable JSON and source hashes,
archived as a separate GitHub Release.
