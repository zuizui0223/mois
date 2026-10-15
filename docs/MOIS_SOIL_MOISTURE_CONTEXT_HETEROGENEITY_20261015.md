# Are soil prediction failures confined to particular macroclimatic environments?

**Original-photo SHA256 independently checked · 15 October 2026 · supplemental ecological analysis.** This is a paired generalisation diagnostic on the existing photograph-colour dataset, **not** a new dataset, ecological selection experiment, independently validated pigment trait, soil water observation, or measured geomorphology.

## Scientific question

The existing *mois* global soil36-complete photo cohort (16,801 sampled species, 1,665 genera) shows both modest soil association and a small negative change in unseen-genus prediction when adding soil to the full 26 climate variables. **Could this apparent global soil prediction disadvantage hide a positive soil contribution in dry, hot, low-elevation or high-elevation taxa?** To distinguish descriptive context dependence from a post hoc search, we compare *identical frozen held-out predictions* stratified by source environmental medians and explicitly test the differences among strata with **paired whole-genus resampling**.

## Independent reproduction of the fixed source model

The two real original FCP photo and environmental ZIPs have been read and checked against their **full original SHA256** values; source photographic hues, unique photo identities, site environments and classifiability were not changed. The source sample was explicitly verified as **31,925 soil36-complete classifiable images**, with **26,163** underlying images across **16,801 species in 1,665 sampled genera** after the ≥3 species/genus eligibility gate.

An **independently implemented** 5-fold, SHA256-hashed whole-genus ridge procedure recreated prior repository prediction errors without using prior fitted predictions, and each fold standardized all predictors using training species **with genus-equal weights only**. The target is original taxon-average classifiable photographic white fraction, and the loss is **equal-genus mean squared prediction error on unseen genera**, not a test of seed production, plant fitness or flower morph equilibrium.

| Predictors accompanying geographic latitude/circular longitude | Independently reconstructed unseen-genus squared error |
|---|---:|
| Geography alone | 0.229810064 |
| Annual rainfall BIO12 + annual atmospheric vapour | **0.227455352** |
| Soil5 alone | 0.228778485 |
| Other 24 climate features (excluding BIO12 + annual vapour) | 0.228173205 |
| Full climate26 | 0.228132885 |
| Full climate26 + soil5 | 0.228537347 |
| Full climate26 + all soil10 | 0.228661099 |

The exact independently reconstructed soil5 addition to full26 climate is **−0.000404461** (1,999-genus bootstrap ordinary95 **[−0.000674388,−0.000129229]**) in heldout improvement. The historical source full26 comparison has the same fitted point value; changing bootstrap seed changes percentile endpoints slightly. This agreement materially reduces the risk that earlier block comparisons were transcription errors or artifacts of a shared result-reading function.

## Environmental subgroup check: same heldout predictions

For each species in the fixed soil36 population, split at the original *species-level* median, a different estimator than a priori dividing continents or genera. Each half is evaluated with **one vote per genus represented in that half**, so any genus with source sampled species on opposite sides of the median can occur in **both halves**. The same pre-existing genus folds remain fixed in prediction fitting; environmental groupings are applied **after prediction**. The median thresholds are: annual rainfall BIO12 **834 mm**, BIO5 warmest-month maximum temperature **27.528** in source native degrees, elevation **547.5 m**, absolute latitude **32.277°**. These source medians are NOT universal arid/humid biome definitions or climate treatments.

| Original species-mean grouping | Low-source group: soil5 additional gain vs full26 | High-source group: soil5 additional gain vs full26 |
|---|---:|---:|
| Annual rainfall | **−0.000278** [−0.000697,+0.000126] | **−0.000480** [−0.000918,−0.000005] |
| Peak temperature BIO5 | **−0.000509** [−0.000976,−0.000018] | **−0.000240** [−0.000619,+0.000160] |
| Elevation | **−0.000219** [−0.000617,+0.000164] | **−0.000360** [−0.000786,+0.000073] |
| Absolute latitude | **−0.000322** [−0.000794,+0.000128] | **−0.000589** [−0.001023,−0.000157] |

The soil-block addition has a **negative point estimate for every tested original-photo climate/elevation/latitude median half**. Nevertheless, many of these ordinary source-bootstrap 95% intervals span zero, and eight reuses of one opportunistic dataset are not independent demonstrations. The result does not show a soil field has no causal fitness or ecological mechanism.

The **water2-versus-geography** gain points were positive in all eight source halves, including cooler and warmer source species (**+0.001765** vs **+0.003102**), lower and higher elevation (**+0.002184** vs **+0.001998**), lower and higher latitude (**+0.001625** vs **+0.002782**), and lower versus higher rainfall (**+0.000741** vs **+0.001586**). These climate subgroup contrasts similarly should not be promoted into a causal effect of temperature, radiation, elevation or latitude on pigmentation.

## Test context heterogeneity, not just individual strata

The *same global 1,665-source-genus multinomial bootstrap vector* was used to compare the high-minus-low effect within each of the four source environmental splits for **two** intended estimands: water2 gain relative to geography and soil5 addition relative to all26. A simultaneous maxT 95% critical value of **2.7490** was obtained across all eight selected *post hoc* high-minus-low contrast differences. No selected high-minus-low interval excluded zero. Illustrations:

| Source grouping | High minus low: water2 gain [eight-contrast maxT95] | High minus low: soil5 extra gain [eight-contrast maxT95] |
|---|---:|---:|
| Annual rainfall | +0.000845 [−0.002493,+0.004183] | −0.000202 [−0.001035,+0.000631] |
| Peak temperature | +0.001338 [−0.001761,+0.004436] | +0.000269 [−0.000565,+0.001104] |
| Elevation | −0.000186 [−0.003384,+0.003012] | −0.000141 [−0.000965,+0.000682] |
| Absolute latitude | +0.001157 [−0.001940,+0.004255] | −0.000267 [−0.001168,+0.000635] |

This selected-eight-family maxT only adjusts the eight deliberately recorded exploratory comparisons, **not the entire preceding *mois* search**, and its error bars are calculated from heldout predictions **without retraining models within bootstrap replicates**. The site/photo classifier and species/environment sample remain nonrandom. Thus there is **no supported interaction-style conclusion** that the incremental soil disadvantage or the water2 predictability is biologically larger in particular climate/elevation halves, even if some within-half intervals exclude zero.

## Explicit source-topography path rather than inferring relief from altitude

One original photo-site altitude does **not** identify hillside slope, uphill exposure, ridge/valley, TPI, TRI, thermal shelter, cold-air drainage or orographic water concentration. None of the source ZIPs includes the DEM neighborhood needed to derive those independently.

A citable **actual** external terrain resource has been verified: Amatulli *et al.* (2018), *Scientific Data*, DOI [10.1038/sdata.2018.40](https://doi.org/10.1038/sdata.2018.40), with downloadable derivatives including slope, aspect, TPI, TRI, vector ruggedness, curvatures and classified landforms at **1, 5, 10, 50 and 100 km** spatial grains. Data series PANGAEA DOI [10.1594/PANGAEA.867115](https://doi.org/10.1594/PANGAEA.867115) reports **CC BY 3.0**. See [EarthEnv topography](https://www.earthenv.org/topography). **Those raster layers have not been downloaded, extracted against the private source photo-site locations or linked to photographed floral colour.** They provide a defensible next measurement source, not a positive or negative terrain result. The existing projected DEM mathematical unit tests at `scripts/analysis/mois_dem_terrain_metrics_gate_20261015.py` verify only synthetic slope/aspect/TPI/TRI functions. For a genuine global terrain test, first prepare rights-cleared 1/5/10km terrain grids and source-site matching/QC, then add them to the **identical source eligible species/genus heldout prediction designs** with no geographic leakage.

## Deliverables, evidence status and editorial decision

- Original source ZIP-resolved independent climate/soil heldout replay: `scripts/analysis/mois_soil_moisture_context_20261015.py`; requires owner-held original photo+environment ZIPs and supports `--write`/`--check`.
- Aggregate-only numerical receipt with exact sample support, subgroup uncertainty and eight maxT differences: `results/mois_source_soil_moisture_context_eight_contrast_20261015.json`.
- Frozen original source full26 comparison: `results/mois_full26_climate_soil10_transfer_20261015.json`.
- Existing matched original soil36 full-audit source: `docs/MOIS_BIOGEOGRAPHY_CLIMATE_SOIL_ELEVATION_TERRAIN_AUDIT_20261015.md`.

**Ecological editorial decision:** Keep the primary worldwide photographed-white/moisture association and the descriptive importance of water-related climate. Treat the lack of additional full26+soil predictive gain as a **model/sample-specific transfer result**, not a rejection of edaphic selection. Do not add a new headline based on unverified subgroup interactions, on genus-free continent fractions, or on terrain geometry that the source did not measure. The main GEB manuscript remains Draft; human authorship/permissions and licensed independent reviewer source archive are still unresolved.
