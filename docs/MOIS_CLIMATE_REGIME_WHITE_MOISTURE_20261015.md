# Does the congeneric photographed-white moisture gradient depend on climatic exposure?

**Original-source exploratory analysis — 15 October 2026.** This is a retrospective photo-colour ecology comparison using the exact original source measurements, not an assay of pigment cost, reproductive fitness, balancing selection, or a 50:50 genetic polymorphism.

## Ecological question

The original world genus-centred photographed-white coefficients were **+2.254545 percentage points/source SD** for annual precipitation BIO12 and **+3.951850** for atmospheric vapour. The original source uses 20,694 qualified species and 1,934 genera, but each genus need not share the same direction. Might the positive species-level rain–white relationship be confined to genera inhabiting sunny, warm or wet climates—an observable but noncausal expectation motivated by some pigment-protection hypotheses?

## Exact original photo and environmental source

Two source archives passed full SHA256 and ZIP member integrity checks: measured photo tables SHA256 2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0 and environmental tables SHA256 6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca. The historical classification and photo-ID deduplication retain 100,543 distinct source photo opportunities and 45,836 classifiable photographs with all climate26 predictors, including 25,028 nominal species. Requiring at least three species per genus yields **the identical 20,694 species, 1,934 genera and 38,498 photo outcomes** in the main macroecology.

Using exactly these species, calculate the genus-equal species-average of each environmental modifier: archived annual mean solar radiation, warmest-month maximum temperature (BIO5) and annual precipitation (BIO12). Split the 1,934 original genera into **967 lower-mean and 967 higher-mean** genera at the observed median for each modifier. Within each half, independently refit the archived genus-centred photographic-white versus one moisture-axis slope adjusted for the same broad absolute-latitude and circular-longitude covariates. Each genus receives one equal total regression weight and the original photo-site environmental SD remains unchanged. The entire unsplit original worldwide white/BIO12 and white/vapour coefficients are reproduced to less than 1e-9 percentage points before subgrouping.

Uncertainty: 1,999 frozen-seed whole-genus bootstrap draws independently within each climate half; the difference between half-specific slopes has a percentile 95% interval using those draws. These are **ordinary observational within-design bootstrap intervals**. They are not the historical original world 24-maxT simultaneous intervals, not study-wide post-selection adjustment, and not a causal interaction estimate. All group divisions and cutpoints remain fixed during resampling. Each group differs in its original species and climatic support; the three grouping variables are examined post hoc on the **same** original photo cohort, not three independent biological replications.

## Original-source annual precipitation results

| Original genus-mean context | Lower context: BIO12 slope [95%] (pp/source SD) | Higher context: BIO12 slope [95%] | Higher minus lower [95%] |
|---|---:|---:|---:|
| Solar radiation | **+2.076 [+0.581,+3.617]** | **+2.481 [+0.567,+4.432]** | +0.405 [−1.996,+2.930] |
| Peak temperature BIO5 | **+1.781 [+0.180,+3.338]** | **+2.485 [+0.704,+4.323]** | +0.704 [−1.583,+3.060] |
| Annual rainfall BIO12 | **+3.381 [+1.059,+5.640]** | **+2.021 [+0.623,+3.497]** | −1.359 [−4.074,+1.314] |

All six original-source **point slopes remain positive** with ordinary positive-side percentile intervals. All three *higher-minus-lower* intervals include zero; these results do not establish an increase or decrease of the moisture–white gradient with sunlight, maximum temperature or rainfall context. An observed high-solar coefficient exceeding the low-solar coefficient is not itself evidence that anthocyanin provides an adaptive thermal or ultraviolet benefit.

## Atmospheric water vapour comparison

| Original genus-mean context | Lower context: vapour slope [95%] (pp/source SD) | Higher context: vapour slope [95%] | Higher minus lower [95%] |
|---|---:|---:|---:|
| Solar radiation | +4.192 [+1.438,+6.685] | +4.159 [+1.485,+6.615] | −0.032 [−3.519,+3.584] |
| Peak temperature BIO5 | +3.879 [+1.284,+6.505] | +4.031 [+1.360,+6.776] | +0.152 [−3.406,+4.056] |
| Annual rainfall BIO12 | +4.455 [+1.511,+7.342] | +3.749 [+1.351,+6.185] | −0.706 [−4.360,+3.210] |

The descriptive positive direction occurs on both sides of each sampled climate grouping for both moisture axes. This does **not** mean a climate response is universal across genera or stable within genuine populations, and cannot repair unresolved photo hue missingness or spatial sensitivity.

## Biological interpretation

**Supported observation:** In the genus-equal broad historical photographic sample, the positive *congeneric* water–white association is present among both lower- and higher-solar, lower- and higher-peak-temperature, and relatively drier and wetter groups of sampled genera. It is not confined to the high-light, hot or wet half of the sampled genera. This gives ecological context for the main global source signal without treating the original strict-photo and local-geography nonconfirmations as evidence that the full result disappears.

**Not identified:** The dataset contains no calibrated floral VIS/UV reflectance, floral or vegetative pigment concentration, pigment-energy expenditure, pollinator reward, pollen export, seed/fruit relative fitness, heritable morph identity or morph-frequency dynamics. Thus no selective maintenance of a global half-white ratio, compensation of pigment protection against metabolic cost, or climate-caused relative survival advantage can be claimed.

**Next biological test:** Sample populations with independently confirmed white and coloured floral morphs, validate flower-organ reflectance and photoprotective compounds, manipulate water and irradiance alongside pollinator access, and track both pollen and seed fitness plus survival across morph frequencies and years. Test an environmental-regime × colour-morph interaction in measured *relative fitness*, with non-floral traits and shared genetic background controlled. A photo-only global colour–climate slope is not that interaction.

## Reproducibility and publication scope

- Fully reproducible authorised original-source script: scripts/analysis/mois_replay_climate_regime_moisture_20261015.py. It reads the owner-held original photo/environment ZIPs, checks both SHA256 hashes, reconstructs the sample and both original slopes, and accepts --write and --check modes.
- Public aggregate-only result: results/mois_climate_regime_moisture_20261015.json. Original source ZIPs, site coordinates, taxon/photo identifiers, genus data rows and photo credits are **not** present in public GitHub.
- Manuscript placement: Supporting Information S2.8 and scientific readiness record, while the main abstract and central positive within-/between-genus scientific inference stay unchanged.
- Draft PR is unmerged and not submitted. Original photo rights/attribution and location privacy, stable anonymised reviewer data, DOI, author approval and journal-format figure proofs are pending.

**Scientific decision:** The broad photo-white versus moisture association persists across sampled genus climate contexts, but the source does not identify a direct pigment-protection or balancing-selection fitness trade-off.
