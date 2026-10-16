# Original-source environmental blocks and taxonomic-scale transfer

**15 October 2026 — retrospective source-only ecological analysis.** Exact original FCP measured-flower-photo ZIP SHA256 2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0 and expanded environment ZIP SHA256 6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca were independently verified and reopened. The analysis retains **16,801 original soil36-complete classifiable taxa, 1,665 nominal genera and 26,163 originally classified photographs** within those taxa. Source colour labels, photo identities and environmental joins were not changed; no sensitive original photo IDs, taxon names, coordinates or attributed images were placed on GitHub.

## Why this is an ecological question

Previous checks found a modest positive world congeneric photographed-white association with annual rainfall and atmospheric vapour, and a positive conditional in-sample soil available-water *proxy* coefficient. Yet soil blocks did not improve absolute photographed-white species-fraction prediction for previously unseen genera. These results need not conflict because an **in-sample congeneric explanation**, **out-of-genus absolute species-white forecast** and **out-of-genus transfer of a congeneric white/environment slope** have different targets.

## 1. Actual within-genus explained white-photo variance

The original source species-level photographic-white fraction and every environmental predictor were centred on their sampled genus. Each genus contributes equal total weighted least squares (WLS) weight. The 13 ecological covariates were divided into five fixed, overlapping-information predictor blocks, but there are no duplicated predictors: (i) absolute latitude and periodic longitude sine/cosine, (ii) annual precipitation BIO12 and atmospheric vapour, (iii) warmest-month temperature BIO5 and mean solar radiation, (iv) point elevation, and (v) soil pH, organic carbon, nitrogen, clay and available-water proxy.

The entire 13-predictor genus-centred WLS explained **0.6202% (R²=0.0062024)** of the sampled within-genus variation in originally classifiable-photo white species fractions. This is not a species-tree variance component, not the variance of plant abundance and not pollen/seed fitness.

The exact **five-block Shapley decomposition** enumerates all 31 nonempty possible block subsets and allocates each block's average increase in descriptive in-sample R² over all 5! addition orders. This sums exactly to the original R². It does *not* causally allocate plant trait evolution to mutually exclusive physical drivers:

| Environmental block | Shapley share of total within-genus photo-white variance, percentage points | Unique partial R² conditional on other four blocks, percentage points |
|---|---:|---:|
| Broad geographical position | **0.1792** | 0.1608 |
| Soil5 including water-retention proxy | **0.1728** | 0.1256 |
| Annual rainfall + atmospheric vapour | **0.1433** | 0.1188 |
| Point altitude | **0.0773** | 0.0045 |
| Warm-month temperature + mean solar radiation | **0.0477** | 0.0111 |
| **Total of five Shapley shares** | **0.6202** | *Partial values overlap and do not sum* |

The in-sample 999-whole-genus bootstrap of descriptive Shapley shares was completed with frozen seed 202610151610. All R² Shapley shares are **nonnegative by mathematical construction**; do **not** read a bootstrap positive lower endpoint as evidence of significant effects, water-mediated selection, an independently tested soil mechanism or genuine partial identifiability. In particular, point altitude added after all four other blocks uniquely explained just **0.0045 percentage points** of the within-genus white-photo variance, consistent with climatic and geographic redundancy, not proof that selection cannot operate along terrain gradients.

## 2. Holdout evaluation after matching the phenotype target

The original soil36 heldout five-fold absolute-white test used the exact same 16,801 source taxa and 1,665 genera and held out entire genera with SHA256-deterministic fold assignment. It showed water-related climate additions improved predicting **absolute species photographic-white fractions** relative to broad geography by **+0.002355** (ordinary whole-genus fixed-loss bootstrap95 +0.001125 to +0.003608), while soil5 additions to the full climate/elevation setup worsened absolute prediction by **−0.000435** (−0.000700 to −0.000171).

To check genuine *transfer of a within-genus ecological species difference*, the new test uses the **same five folds and same taxa**, but its test response is a species photographic-white fraction minus its **own heldout sampled genus's mean**, and its environmental inputs are likewise centred on that heldout genus. Training genera are disjoint. Features are training-only standardized using equal-genus weights and fitted with fixed ridge penalty 0.01 without residual-response clipping; uncertainty uses 1,999 genus resamples (seed 202610151710) of fixed heldout losses.

| Source 16,801-species heldout target: within-genus white-photo fraction deviation | Heldout squared error | Gain relative to geographic predictor baseline |
|---|---:|---:|
| Geographic coordinates only | 0.128554 | Reference |
| Geographic coordinates + annual rainfall/vapour | **0.128322** | **+0.000232 [−0.000104,+0.000564]** |
| Geographic coordinates + water + warm-month temperature + solar + altitude | 0.128488 | +0.000066 |
| Previous complete climate/elevation block + soil5 | 0.128570 | −0.000016 |

The relevant **incremental soil5 versus identical full climate/elevation** contrast on the genus-residual response is **−0.000082 [−0.000310,+0.000127]**, not a conclusive transfer gain or loss. The new within-genus water-only gain interval likewise includes zero. These are **small** compared with the heldout genus-centred baseline loss about 0.128.

**Critical interpretation:** These heldout data use the *observed mean of the test genus itself* to define the species-relative response. They test whether previously learned **within-genus slopes** transfer across newly sampled genera, **not** whether one could prospectively predict the absolute white frequency of an entirely unknown genus without observing its colour. The previous absolute-white model is the latter different task. Treating the two losses and their units as identical scientific outcomes is incorrect.

## 3. Ecological implications and negative gates

- **Supported:** More than one environmental block has descriptive within-genus photographed-hue information; soil5 is not universally irrelevant, yet all 13 variables explain only **0.62%** of the sampled within-genus species variation. In-sample soil and climate associations are insufficient to guarantee heldout intergeneric prediction on either target.
- **Biogeography:** Stronger world genus-mean and broad white–moisture results can reflect where environmental niches and different photographed taxonomic assemblages are sampled. This is compatible with floristic sorting and clade/observer structure, not a proven historical species-tree signal.
- **Environmental dependence:** Soil36 predictors are nonrandomly available along geographic and altitude gradients; many source photographs have no reliably classifiable white/nonwhite floral phenotype. High-dimensional soil texture terms are partially redundant; the selected soil5 avoids direct sand-silt-clay closure.
- **Unmeasured:** Genetic white/coloured morph frequencies, floral reflectance and pigment chemistry, pollinator pollen export or receipt, survival and seed fitness, and any selected global 50:50 equilibrium.
- **True terrain:** The original dataset has individual photo-site point altitude only. Slope, aspect, topographic position (TPI), ruggedness (TRI), geomorphic shelter and catchment wetness have **not** been fitted or measured at the source flower sites.

For an independent future biogeographical/terrain test, [RESOLVE Ecoregions 2017](https://ecoregions.world/) (Dinerstein et al. 2017; DOI 10.1093/biosci/bix014) has real realm/biome/ecoregion polygons and [EarthEnv global multi-scale topography](https://www.earthenv.org/topography) (Amatulli et al. 2018; DOI 10.1038/sdata.2018.40) provides slope, TPI, TRI and related landforms at different resolutions. Source licences and versions, realistic island and coastline matching, sample sizes, spatial support and a local private coordinate overlay need checking. **Do not submit original sensitive photo positions to a third-party public geocoding API or report actual terrain as tested until those rasters have been obtained and joined.**

## Reproducibility

Companion original-ZIP script: scripts/analysis/mois_environment_block_shapley_residual_cv_20261015.py, requiring the owner-held original photo ZIP and environmental ZIP with complete SHA256 verification. Public aggregate-only source receipt: results/mois_env_shapley_and_genus_residual_transfer_20261015.json. Independent original-source calculations and bootstrap/heldout comparisons were run locally on the verified ZIPs. Hosted GitHub CI has **no private photo or site ZIPs** and checks syntax, aggregate source cohort invariants, synthetic WLS/Shapley arithmetic and publication wording; it must never be described as a raw-source replicated ecological experiment.

**Manuscript decision:** Keep the positive original 20,694-species congeneric white–water-climate association as the working paper's main result; supplement it with the more scale-specific interpretation here (**SI S2.16**) without claiming soil absence, causal pigmentation ecology, or actual terrain morphology.
