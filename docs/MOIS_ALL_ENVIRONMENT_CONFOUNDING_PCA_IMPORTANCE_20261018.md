# MOIS: all environmental variables, collinearity, PCA and most influential factor

**Main Article environment-confounding analysis. 18 October 2026.**

## Ecological question

The verified global original-photo source contains **313,535 actually processed
RGB photo attempts** and **151,239 visible four-hue labels**, **71,427 white**.
The worldwide near-half photo-white share is a descriptive photograph/taxon
sample, **not** an exact 50:50 evolutionary equilibrium or a fitness optimum.

Find the strongest *source-adjusted association/predictive environmental
group* for white versus chromatic photographed hues after accounting for
covariate correlation. Maintain global broad flora as the primary outcome;
observed month/season climate is a separate Supplementary Information (SI)
question and not a same-individual flower recolouring design.

PCA is useful to describe gradients but PC1 maximises *environmental*
variance, not variance in observed floral colour. Comparing univariate
rainfall, vapour pressure, light and heat regressions or their raw
standardised coefficients would assign redundant correlations to arbitrary
predictors. Reporting one "winner" without matched source support and
out-of-sample validation is misleading.

## Complete, source-backed predictor acquisition

Use the exact previously released [313,535 photo-ID source feature lake](https://github.com/zuizui0223/mois/releases/tag/mois-primary-313535-fullfourth4096-featurelake-v2),
do not download raw licensed photos again and do not add the separate 256
purposefully season-enriched source into global frequency estimates.

Expand site-level long-term WorldClim **1970–2000** climate to all 19 BIO
indices, the 12-month mean vapour normal and the already verified physically
correct annual solar normal (21 variables). Point-sample a **2.5 arcminute
coarse** DEM for elevation, slope, roughness, TRI, TPI and aspect components
(7 variables). Retain 5 historical SoilGrids-derived source variables
(pH, organic carbon, nitrogen, clay, soil available-water proxy) only where
they truly exist. This gives **33 measured abiotic candidate fields**, in
addition to nuisance latitude, absolute latitude, longitude sin/cos and
source-recruitment cohort.

Source and model checks:
- All original photo IDs and measured white/other four-hue ROI classifications
  must exactly match the authenticated 313,535-source Parquet snapshot.
- Original WorldClim 19-BIO and solar ZIPs are verified by pinned SHA256 and
  the vapour archive is compared with the **independently completed** original
  309,439-photo WorldClim/vapour model receipt.
- Existing BIO5/12/14/15 site values are checked against freshly downloaded
  original rasters for changed location, period, or units.
- All environmental source photos are retained, including hue-unknown
  photographs; unknown hue is not silently coded as coloured.
- Soil/terrain gaps are displayed, not treated as zero, not extrapolated to
  unmeasured photo sites and not patched with unverified global point values.
- The original Wm2 source-unit error is already corrected in the frozen
  313,535-photo feature lake; do not divide solar again.
- Same latitude does not imply the same moisture/humidity. Annual rainfall,
  dry-month rain, vapour, precipitation seasonality, heat, solar and terrain
  are distinct contextual descriptors. **Distance to coastline is not yet a
  source-verified predictor** and is explicitly absent from this ranking.
- WorldClim BIO indices derive from temperature/precipitation and are often
  strongly correlated. A 33-variable model is not 33 independent tests.

## Matched complete-case predictive comparison

Build three supporting photo-denominator panels (climate only,
climate+terrain, climate+soil) to show what was observed, but make the primary
**all-five-block ranking on exactly the same complete-case original photos**,
before computing each named-taxon photographed-white fraction. Restrict to
qualifying source genera with >=2 nominal taxa. A single environmental
predictor cannot gain strength simply by admitting 100,000 photographs
that lack soil observations.

Five blocks:
1. Moisture: BIO12–BIO19 and mean annual VAPR.
2. Temperature: BIO1–BIO11.
3. Solar: converted annual srad normal.
4. Terrain: coarse elevation, slope, roughness, TRI/TPI, eastness/northness.
5. Soil: pH, SOC, N, clay and measured available-water proxy.

Within each cross-validation **training fold only**, standardise all block
variables and reduce collinearity with a within-block PCA retaining 90%
variance (up to six principal components per group). Standardise the five
nuisance geography/source variables in training folds also. Train
**equal-nominal-genus-weight Ridge regression** on source named-taxon
photo-white fraction, and evaluate bounded photo-white predictions on
(a) five whole-genus holdout folds and (b) five approximately 10-degree
photographed taxon-centroid geographic folds.

For each holdout design compare geography/source-only, full model,
full-minus-one-environmental-block and geography plus one block on the
**same held-out taxa**. Use reduction in held-out weighted prediction
error (R² gain) as environmental group *conditional predictive importance*.
Also fit a raw-standardised all-33-predictor model and remove one original
variable at a time, retaining **all other source predictors in the fold**:
this is an exploratory marginal conditional ranking and may be nearly zero
under strong collinearity.

In parallel, compute full-panel *unsupervised* 33-variable global PCA and
per-block PCA loadings, axis variance and strong Pearson correlations as
environmental gradient interpretation, **never** as evidence of the most
important flower colour trait determinant.

A group is named **tentatively strongest** only when its full-minus-drop
held-out R² is positive and ranks first in both whole-genus and broad
geographic holdout designs. If directions diverge or are nonpositive,
report "no unique strongest source-transportable environmental group."
Even a robust "winner" is not a causal selection factor or measured
pigment biosynthetic fitness benefit.

## Interpreting any outcome

The previous global genus-equal 309,439-photo six-source maxT/999 water
analysis found +1.99 white percentage points per environmental SD annual
precipitation and +3.16 for atmospheric vapour globally; neither maintained
a signed simultaneous interval under genus×coarse geographic conditioning.
The newer 313,535-photo exploratory models had annual rain +1.92 nominal
cluster 95% [+0.15,+3.69] but climate+geography explained only ~0.4% of
within-genus taxon photographic-white variation. **The new predictive
group ranking is different from these regression coefficients** and must
not be described as independent confirmation.

Outcomes are photo-label frequencies and image-detection-success rates, not
census counts of flowering plants. Source taxa are not a dated phylogeny;
source photographs are not a common-garden or measured pollination fitness
trial; WHITE flowers may still have UV-absorbing pigmentation; "coloured"
is not equivalent to an anthocyanin assay. Missing-colour MNAR,
observer effort, source taxon turnover and geographic confounding remain.
No meaningful ecological finding should be selected solely for being
positive.

## Durable outputs and continued photo acquisition

Scripts:
- scripts/analysis/mois_original313535_full19_bio_vapr_terrain_lake.py
- scripts/analysis/mois_global313535_all_environment_pca_ablation.py

Workflow:
- .github/workflows/mois-original313535-all19-bioclimate-pca-ablation.yml

The workflow emits a SHA-signed extended climate+DEM photo-level Parquet,
environmental missingness CSV, descriptive 33-variable PCA axes/loadings,
correlation-pair diagnostics, exact-population 5-block genus/geography heldout
ablation R², raw individual predictor conditional ablations, and a
machine-readable explanation of any strongest factor. A trusted main push
creates a versioned GitHub Release. Incremental subsequent 256-image RGB
workers may proceed separately and only their new, verified photo IDs and
environment sites should be appended. Refit PCA and held-out model snapshots
after sufficiently large **genuine image-pixel** batches rather than
reprocessing every photo when a 256-photo shard finishes.
