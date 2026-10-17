# MOIS: exact photographic taxon turnover explains much of global seasonal white-label composition

**Actually authenticated original source (17 October 2026).** The [successful original photo source/date run #39020466312](https://github.com/zuizui0223/mois/actions/runs/39020466312) contains **302,242 genuine original image-pixel processing attempts**, **146,061 actual four-hue automatic ROI outputs** and **146,052 classifiable original photo IDs with an exact historical source or separately verified public iNaturalist observation date**, from **25,638 original named taxa**. This *global among-taxon photographic composition* estimand belongs to MOIS, NOT the separate FCP within-species D/W or 50 km morphological polymorphism analysis.

## Current full source photo proportions, not flower abundances

| Photo-labelled extra-tropical original source | Summer white / four-hue photos | Winter white / four-hue photos | Unpaired photographic-white difference |
|---|---:|---:|---:|
| Northern | **21,448/43,717 = 49.06%** | **2,164/4,662 = 46.42%** | **+2.643 pp** |
| Southern | **3,404/6,721 = 50.65%** | **1,843/4,685 = 39.34%** | **+11.309 pp** |

All input photons/pixels and four-hue ROI labels were previously measured; photo-date restoration (original 3,430,502-row iNaturalist/FCP source observation UUIDs + independent 1,732 numeric observation ID dates) did **not** reclassify any pixel or presume photograph-specific anthesis.

Of the **7,207 Northern** source taxa photographed with one or more classifiable images in summer and/or winter, only **647** were photographed in **both** seasons (Jaccard **0.08977**); the South shares only **391 of 4,364** source taxa (Jaccard **0.08960**). These are photograph-set overlaps, **not** estimates of field flowering species presence, and original taxon IDs are not fully resolved time-calibrated phylogenetic lineages.

## An exactly additive decomposition of the *entire* photographed source

A common-taxonomic-support symmetric Kitagawa identity partitions summer–winter source photo-white fractions into relative photo frequencies and conditional photo-white rates for taxa observed in BOTH seasons. To reconstruct the difference across the *entire* actual photo dataset, add separately the exact seasonal corrections for taxa photographed only during one season; never assign those species any opposite-season hypothetical colour.

| Additive term, source photo-white summer minus winter | North (pp) | South (pp) |
|---|---:|---:|
| Within shared photographed taxa: difference in photo-white rate | **+2.761** | **+3.516** |
| Within shared photographed taxa: change in relative *photographic frequency* | **−0.965** | **+7.689** |
| Taxa with only source summer hue photos | +3.004 | +3.832 |
| Taxa with only source winter hue photos | −2.158 | −3.728 |
| **Entire source summer−winter photographed-white difference** | **+2.643** | **+11.309** |

The identity is mathematically exact on the classified original photo census, not a causal share of biological evolutionary pigmentation or flowering community abundance. For the Southern difference, **+7.69 of the +11.31 pp** reflects changes in photographic representation of taxa that appear in both seasons, whereas **+3.52** is the shared-taxon white-photograph-rate component and source-exclusive taxa provide a net **+0.104 pp**. This pattern cannot be inferred solely from a claim that warmer/drier/brighter months physiologically turn flowers white.

A separate 999-resample **whole nominal taxon** exploratory percentile bootstrap for the Southern common-taxon component yields:
- Within-common-taxon original photo white-rate component: **+3.516 pp [~+0.072,+7.032]**.
- Common-taxon relative photographed-frequency component: **+7.689 pp [~−0.423,+14.999]**.

Both are PHOTO-SAMPLING rather than independently verified biological mechanisms. These intervals are **not** geographically or observer-disjoint, not multiple-testing adjusted, and do not cover possible photo-specific organ-classifier error or the optional choice of analytic geographic grain. They should not be represented as identifying a causal taxonomic turnover effect.

## Same photographed locality/year is a separate, stricter ecological estimand

With actual source observation day restored on almost every classifiable ROI, the corresponding **0.25° source cell × identical meteorological season-year** photographed *community* summer–winter comparison requires at least five distinct original iNaturalist observation IDs and two named taxa in **both** seasons. This provides **24 Northern cell-years /347 classified photo ROIs** with mean source photo-white **+5.39 pp [−5.17,+15.05]** and **27 Southern cell-years /535 classified photos** with **+3.57 pp [−7.90,+12.84]**, geographically clustered exploratory 999-bootstrap percentile 95% intervals, all spanning zero.

The raw whole-Southern-source **+11.31** photographed white-label summer excess therefore is not independently confirmed as a single **positive local, same-year** phenomenon; a 0.25° geographic cell still spans tens of kilometres and can include more than one real ecological population. Taxa shared globally in summer/winter may be represented by different individuals, year/site, observer and flower-specific organ status. All-source original photos missing in winter are NOT evidence that winter flowers were biologically absent.

## Audit, figure and scientific boundaries

The new source-signature-locked `scripts/analysis/mois_full146052_original_photo_season_taxon_composition.py` reads the **exact final authentic 302,242 original-image/photo-date archive** (signed GZIP SHA256 `9c9692238f58cc61561270f7c07b1d782ce711adeade9c58550866f668fb4690`), applies correct hemisphere meteorological seasons and checks taxon and source-photo denominator additivity to less than 10−8 pp. It uses the *same* preexisting original MOIS symmetric photographic Kitagawa method and 999 resamples, not the previously sampled 127,504/140,247/143,770 photographic dates. The accompanying GitHub Actions test joins only the already completed authenticated original photo source — no new image downloads or provider API calls — and exports the full-source source/photo-taxon decomposition receipt and a four-term publication-figure CSV.

**Research claim allowed:** global source season effects in photographed white-hue composition depend strongly on which named taxa are photographed and the relative frequency of their images, and differ from the geographically constrained local-photo-assemblage signal. **Research claim NOT allowed:** universal seasonal pigment plasticity, common garden fitness, pigment energy costs, pollinator choice, actual observation-year moisture/VPD/solar effects or 50:50 evolutionary colour balance. Source photo-specific botanically confirmed open-flower stage, observational effort and independent phenology remain unmeasured.
