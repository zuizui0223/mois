# Worldwide photographed floral white and 13 environmental variables: pooled versus within genus

**16 October 2026 · direct reconstruction from exactly the original FCP source ZIPs.** This is a new **analysis of the existing measured floral photo-colour and photo-site environment tables**, not new plant sampling, an independent flower-colour dataset or a field experiment. Both original ZIPs passed exact SHA256 and member CRC checks. Photo SHA256: `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`; environment SHA256: `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`.

## Ecological question

The earlier **0.6202%** explained variance was for **differences among species WITHIN their sampled nominal genera**, after subtracting each genus's sampled photographed-white mean. It **did not** describe the full cross-taxon/whole-world white-colour variation. The appropriate follow-up is to fit exactly the same **13 environmental variables** to photographed-white species proportions **without genus centring**; compare equal-species and equal-genus weights; quantify how much within- versus between-genus colour structure exists *in this same complete-case source population*.

## One exact source cohort and two weighting conventions

- Soil36-complete source: **31,925** four-hue-classifiable photo records among **20,786 named taxa/4,853 nominal genera**. Photo/site climate and all ten original soil terms must be present. Main matched genus-qualified subset (≥3 classified-photo species per genus): **26,163 source photographic hue observations, 16,801 named taxa, 1,665 nominal genera**. **Both original photo IDs and soil/geographic join are unchanged.** The 13 predictors are complete for both subsets.
- Outcome: each sampled species' proportion of *original classifiable* photographs classified `white`, with `yellow_orange`, `red_pink`, and `blue_purple` forming the other three categories. This is NOT a genetic white morph/individual fitness variable or the fraction of all (including unclassifiable) photos that would be white.
- Environmental variables, 5 nonoverlapping named blocks of 13 columns: **Geography** (absolute latitude, sine longitude, cosine longitude; 3); **Moisture** (annual precipitation BIO12, mean atmospheric vapour; 2); **Temperature/light** (BIO5, annual solar radiation; 2); **Elevation** (original photo-site point altitude; 1); **Soil5** (pH, organic carbon, nitrogen, clay, modelled available-water *proxy*; 5). Predictor values are species averages of original photograph-site extractions. Coefficient scales are original *soil36 photo-site* predictor standard deviations, not SDs in this reduced species subset.
- Primary model: ordinary *in-sample* weighted least squares (WLS) of full species-level white-photo proportions on intercept plus all 13 variables **without genus fixed effects**, weighting each of 1,665 source genera equally in total (`w_species=1/n_genus`). The **secondary equal-species** model weights each eligible species once. Also refit the 13-variable regression **after genus centring** to reproduce the old original same-source R². Neither in-sample R² is a cross-validated measure of forecasting previously unseen species.

## New source-replayed environmental explanation

| Identical 16,801 taxa / 1,665 genera | All-white species variance explained (R²) | Interpretation |
| --- | ---: | --- |
| All 13 predictors, **pooled white (genus-equal weights)** | **2.000%**; 999 whole-genus bootstrap95 **[1.509,3.006]%** | Captures differences between and within sampled genera; no genus identity term fitted |
| All 13 predictors, **pooled white (species-equal weights)** | **1.282%** | More heavily sampled genera contribute more species |
| **Genus indicators only**, genus-equal | **44.315%** | Description of in-sample variation assigned to 1,665 fitted genus averages; not 44% inherited trait variance or validated phylogenetic signal |
| **Genus indicators + all 13 variables**, genus-equal | **44.661%** | Environmental terms beyond genus means add **0.3454 percentage points** of *total* observed white variation |
| All 13 predictors **within-genus white deviations**, genus-equal | **0.6202% of within-genus variance** | Previously reported source result independently reproduced to 1e-10; this is a distinct smaller denominator |
| All 13 predictors **within-genus white deviations**, species-equal | **0.3399% of within-genus variance** | Different weighting of sampled species, same physical rows |

The contrast between **44.315%** for genus averages and **2.000%** for observed environment-only covariates cannot be interpreted as an additive attribution of causation. The statistically valid *same-sample* identity is `R2(genus indicators + within-genus 13) = R2(genus indicators) + [1 − R2(genus indicators)] × R2(within genus 13) = .446607`. The genus-only number is intrinsically boosted by fitting **1,665 category means**; holding original sizes fixed and randomly reassigning the observed species photo-white proportions still gives a genus-equal average **19.063%** explained variance (1,999 permutations, central95 **17.782–20.318%**). Nevertheless, observed 44.315% exceeds that specific conditional null (`P=0.0005`). Under equal-species weighting, in-sample genus means explained **34.259%** against shuffle null **9.902%**. This establishes *sampling-frame taxonomic colour clustering*, not molecular phylogeny or taxon-specific natural selection.

**All soil-complete source taxa without the ≥3-species genus gate:** in the full **20,786 sampled species/4,853 nominal genera**, a **species-equal** pooled 13-environment WLS explained **1.2807%** of observed white-photo-fraction variance. This closely matches **1.2823%** for the 16,801-genus-qualified taxa, showing the pooled species-equal environmental R² is not just created by excluding one-/two-species genera. The 20,786 species have no uniform congeneric sampling and are not suitable for naively carrying over a genus-centred source result. The full-soil source is not all 42,111 unfiltered photograph opportunities: nonclassifiable white labels and incomplete soil layers remain missing.

## Which of the 13 predictors is associated with pooled white?

Conditional point coefficients below are **percentage points of species-level photographed-white proportion per 1 original soil36 source photo-site SD**, in the *same* pooled **genus-equal 13-variable** fitted model. Uncertainty uses **999 resampled whole genera**, refitting all 13 with an intercept; reported maxT 95% endpoints are simultaneous **within this 13-parameter family**, not across all original FCP macroecology models or spatial trials.

| Predictor | Pooled WLS point pp/source SD | Genus-bootstrap 13-family maxT95 |
| --- | ---: | ---: |
| Absolute latitude | **+4.793** | [−0.701,+10.287] |
| Longitude sine | +0.063 | [−2.090,+2.216] |
| Longitude cosine | **−3.132** | **[−5.209,−1.056]** |
| Annual precipitation BIO12 | **+3.770** | **[+0.817,+6.723]** |
| Mean water vapour | +2.789 | [−2.912,+8.490] |
| Warmest-month temperature BIO5 | −1.679 | [−5.330,+1.972] |
| Annual solar radiation | +0.344 | [−3.071,+3.760] |
| Source photo-site elevation | −1.505 | [−5.508,+2.497] |
| Soil pH | +0.596 | [−2.199,+3.390] |
| Soil organic carbon | +0.157 | [−2.857,+3.170] |
| Soil nitrogen | +0.983 | [−1.945,+3.912] |
| Soil clay | +0.061 | [−2.091,+2.213] |
| Soil available-water proxy | +0.835 | [−1.289,+2.960] |

**Biological interpretation:** Annual precipitation still has a **positive independent-model partial point coefficient** for pooled sampled white flowers after conditioning on the 12 other selected covariates, with an interval excluding zero in this one 13-variable exploratory family. Atmospheric vapour shares much variation with rainfall, latitude, altitude and other climate layers, so its **individual partial** interval includes zero when all terms enter together. Sine and cosine are a **paired cyclic longitude encoding**; their individual signed values must not be interpreted as independent ecological East–West effects. Each soil's individual partial interval spans zero, even though the soil block has descriptive shared information. None establishes that high-rainfall environments *caused* changes in floral pigments.

## Environmental block accounting on pooled white (not causal explained-variance partition)

The 13 variables fall into five disjoint *column* blocks with potentially overlapping geographic/climate/soil information. Exhaustively fit all 31 nonempty block subsets and split the **2.000%** descriptive in-sample R² via a five-block Shapley attribution. The Shapley shares are mathematically nonnegative and sum exactly, and their genus bootstrap intervals should NOT be read as separately significant ecological effects.

| Block | Pooled global white R² Shapley allocation (percentage points of total sampled white variance) | Increment in pooled R² given all other four blocks |
| --- | ---: | ---: |
| Geographic position (3) | **0.6490** | 0.4386 |
| Rainfall + vapour (2) | **0.5574** | 0.2728 |
| Soil5 (5) | **0.4022** | 0.0388 |
| Elevation (1) | **0.2373** | 0.0162 |
| Temperature + solar radiation (2) | **0.1540** | 0.0254 |
| **All13** | **2.0000** | Partial increments overlap; **do not add** |

Same species, but **within genus** the five-block allocations were (geography **0.1792**, moisture **0.1433**, soil5 **0.1728**, elevation **0.0773**, temperature/light **0.0477**) percentage points of **within-genus** variance, totalling **0.6202**. Different denominators and statistical targets make the pooled and within Shapley shares non-subtractable as a direct causal partition.

## Verification, robustness and stop conditions

- **Original-source replayer:** `scripts/analysis/mois_global_pooled_white_13env_20261016.py`; reproducible with the exact owner-held photo and environment ZIPs, source identity checks, same 16,801/1,665/26,163 support, cross-checked original 0.6202% within-genus R², 999 whole-genus resamples, 1,999 genus-size-preserving source-label permutations, and additional all-soil-eligible 20,786 species sensitivity.
- **Source-public aggregate result:** `results/mois_global_overall_white_13env_20261016.json`, containing cohort denominators, same-source R², all 13 slopes and simultaneous intervals, five Shapley contributions, genus-shuffle null, and precise source hashes **without taxon identifiers, protected coordinates, or original photo IDs**. The public hosted CI cannot refit from the rights-controlled ZIPs and tests their aggregate values and synthetic linear algebra, not new raw-source observations.
- **Additional independent numerical cross-check:** direct weight-sqrt least-squares of the original sampled species and all 13 columns gave exactly the same pooled R² and slopes as the centred normal equations for both species and genus weighting (absolute numerical differences under 1e-10).
- **Limitations:** photo-hue success is geographically nonrandom, missing colour exceeds half of all opportunities, soil36 data are nonrandom complete cases, clades here are nominal genera not independently verified molecular phylogeny, soil available water is a numerical source proxy, and **elevation is not terrain slope/aspect/TPI/TRI**. None of the model R² values establishes biochemical pigment cost, pollinator benefits, stabilizing selection, heritable white-morph equilibrium or generalizable local causality.

**Ecological bottom line:** In the world-scale soil-complete source, **moisture and geography are directionally associated with which species' photos are white**, but the full available environmental predictors together account for only about **1–2%** of worldwide photographed-white species differences. **Taxonomic composition is far more structured in this sampled species panel**, and the difference between pooled and congeneric gradients is biologically worth explaining, without treating taxonomic labels as an evolutionary mechanism.
