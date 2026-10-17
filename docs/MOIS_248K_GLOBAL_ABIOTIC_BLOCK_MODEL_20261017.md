# MOIS: Global photographed four-hue responses to geography, climate, landscape and soil

**Active study:** `zuizui0223/mois`, broad species- and genus-level macroecology. Images shared with `fcp` are allowed, but **FCP's within-species spatial polymorphism D/W/50-km local depletion is not an outcome of this study**.

## Data and estimands

The real source consists of **248,312 unique original photo IDs** (original MOIS historical + FCP high-depth + six independently measured ROI technical groups) and **120,807 actual four-state automated photographed-hue classifications**. Fully unclassifiable images remain in source denominators but cannot be coded as white or coloured petals. Public photo source does not prove genotype, reproductive fitness, flower tissue confirmed by human botanists, or balancing selection.

The enriched site fields were **actually extracted** from original WorldClim 2.1 raster source in successful [WorldClim action 38992341847](https://github.com/zuizui0223/mois/actions/runs/38992341847): 19 BIO long-term variables at ~18 km 10 arcminute grid; original photo-site elevation and ~4.6 km terrain slope/roughness at 2.5 arcminute; year-round mean solar radiation and source-month temperature/precipitation/radiation *normals*. Archived soil variables exist only for historically covered sites and are **not imputed** for all new photos.

### Matched source-photo block modelling

The original photo source/lat-lon/elevation and genus are included in each analysis; `new_source_fraction` controls the shallow historical MOIS versus high-depth/new ROI sampling identity. Nominal genus is the first Latin genus word in source scientific name; species are keyed by iNaturalist taxon ID. Four outcomes are species-level photo shares of **white, yellow/orange, red/pink, blue/purple**.

Five independent descriptive sensitivity analyses are computed:

| Analysis | Source-photo support | Nested comparison |
|---|---|---|
| Climate over geography | All complete global long-term BIO+radiation+geo photo sites | geographic baseline → +BIO1/BIO5/BIO12/BIO15/mean solar |
| Terrain | Same original photographs with genuine ~4.6 km terrain values | geo → geo+climate → +slope/TRI |
| Soil | Same original photos with all five original soil fields | geo → geo+climate → +pH/SOC/N/clay/water-capacity |
| Minimum two classifiable original photos per named species | Climate-complete photos under stricter species measurement depth | geo → geo+climate |
| Historical MOIS source alone | Original global shallow photo source, not high-depth FCP | geo → geo+climate |

Every nested block uses **identical original photo support**, not merely a new dataset with different named taxa. At the named-species level, each species has equal weight, all predictors and photo shares are centred within genus, and each genus is weighted equally. Additional covariates with no within-genus contrast are dropped explicitly. Coefficient units are **percentage points in photographed hue fraction per between-taxon standard deviation** of the covariate. The nominal 95% intervals use a genus-cluster sandwich covariance and do *not* provide the 999 genus bootstrap / maxT confirmatory uncertainty of the original historical GEB analysis. Exploratory Benjamini–Hochberg q-values are computed over the four-hue × climatic/new-block predictors, with exact tested family specified in output.

Independent predictive diagnostics exclude entire genera in five hash-deterministic CV folds; CV R² is evaluated for within-genus deviations about the genus mean on genuinely heldout genera. In-sample variance fractions, CV R², and incremental change for each covariate block are all output. A tiny within-genus R² or negative heldout R² must be interpreted as **poor practical ecological prediction**, even if individual exploratory coefficients have a small nominal p-value.

### Temporal environmental variables remain separate

The source-recorded months are used to describe photographed hue opportunities by calendar month, hemisphere and source. **Annual/monthly WorldClim 1970–2000 normals are not actual observation-year weather.** The corresponding TerraClimate actual-year monthly temperature, precipitation, solar radiation, vapour deficit and water deficit require a separate full original-date/coordinate join; there is no evidence of a universal within-individual seasonal flower-colour switch. Valid dayparts such as dawn/noon/dusk need actual UTC/local timestamps with timezone, which are absent from this photograph archive.

## Actual executed results (verified 17 October 2026)

Source-verified real run: [Actions #38992961375](https://github.com/zuizui0223/mois/actions/runs/38992961375), artifact `12086155833`. Immutable decision record: [`results/mois_248k_global_abiotic_block_decision_20261017.json`](../results/mois_248k_global_abiotic_block_decision_20261017.json).

- Full global genus comparison: **120,381** original photos with complete climate/geography; **22,705 named taxa across 2,858 nominal genera** contribute after requiring multiple species/genus. Geographic baseline explained **R²=0.003111**, adding 5 BIO/radiation climate terms raised this to **0.004498** (ΔR²=0.001387). On five-fold **whole-genus-heldout** comparisons, predictive R² actually moved **0.001844 → 0.001379**, not an improvement in practical generalization.
- Exploratory photo-white vs annual BIO12 rainfall slope in the joint genus-centred, source/geographic/elevation/temperature/solar adjusted model = **+2.1381 percentage points per species-covariate SD**, genus-cluster nominal 95% **[+0.7497,+3.5265]**, *p*=0.002542. Crucially, across the **four hue × five climate predictor** family, Benjamini–Hochberg *q*=**0.050840**: narrowly **not** below 0.05, and still not the original simultaneous maxT inferential framework. The direction remains aligned with the frozen historical MOIS global +2.255 pp/SD but is **not a new independent-source replication**.
- Terrain block: identical original-photo support **104,809** classified photos; **20,639 taxa/2,691 genera**. Incremental in-sample ΔR²=**0.000321**, but whole-genus heldout ΔR²=**−0.000302**. White photographic fraction effects: slope **+1.49 pp** (nominal 95% [−0.35,+3.33], *q*=0.397); TRI **−1.41 pp** (nominal 95% [−3.58,+0.77], *q*=0.621). These are **~4.6 km grain** terrain values, not field-scale landscape adaptation.
- Soil block: identical original-photo support **31,943** classified photos, **18,411 taxa/2,464 genera**. Incremental in-sample ΔR²=**0.001451**; genus-heldout ΔR²=**+0.000490**, both small. Original site-map available-water proxy vs photo-white **+1.80 pp** (nominal 95% [+0.56,+3.05]) but BH *q*=**0.184**, **not supported after the declared exploratory multiplicity correction**. No source-independent soil/fitness cost claim.
- Stricter **≥2 classified original photos per species**: **7,604 taxa/1,333 genera**; white/BIO12 +2.86 pp/SD nominal 95% [+0.83,+4.88], BH *q*=**0.0568**, again above the threshold.
- The original historical-source-only diagnostic returns white/BIO12 **+2.08 pp** nominal 95% [+0.68,+3.48], BH *q*=0.0708. The broader photo union is not an independently sourced positive result.

**Biological conclusion at this stage:** worldwide genus-level photographed white composition retains a weak positive rainfall-associated direction, but geographic, source, taxonomic, image-quality and multiple-testing qualifications are material. The additional terrain/soil features do not yet provide a robust new floral-adaptation mechanism, and the minute within-genus/heldout explained variance cautions against calling these globally predictive explanatory models.

### Reproduction

`scripts/analysis/mois_248k_global_genus_abiotic_blocks.py` and its [Actions workflow](../.github/workflows/mois-248k-global-four-hue-abiotic-block-model.yml) read the exact photo × static environment archive from original GH Actions without creating new data or importing FCP polymorphism metrics, then output the descriptive four-hue coefficient ledger, same-photo nested variance comparisons, date/hemisphere colour-photo summaries, per-model support denominators and SHA256 receipts.

**This is exploratory macroenvironmental source analysis, not a new frozen independent-source replication or a natural-selection experiment.** The original historical 100,543-photo global MOIS water gradients remain frozen. The enlarged 248,312-photo model is a candidate next-stage GEB environmental paper, for which independent source-photo quality, observation effort, taxonomic and regional matching, and multiplicity/transport robustness must accompany any positive claim.
