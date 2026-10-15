# Ecological scale separation: global water–white gradients are not a transferable congeneric colour rule

**15 October 2026 — cross-receipt synthesis for the Mois/GEB Draft PR #1.** This note independently recomputes arithmetic from previously source-verified, **aggregate-only** receipts. It is **not** a new source-photo model, additional independent replication, preregistered hypothesis, phylogeny, or measured fitness result. Runnable audit: `python scripts/analysis/mois_ecological_scale_synthesis_20261015.py --check`.

## A single scientific question at four observational scales

**Why can wetter photographed source environments have more white-flower labels globally, although water climate is of little use in predicting which of two new congeneric species will be whiter?** An association at one trait-composition scale need not be a transportable individual-species contrast at another scale.

| Level; source population | Reconstructed finding | What it estimates |
| --- | --- | --- |
| Photographic opportunity; 100,543 unique photos | 45,953 classifiable; 46.95% white **among classifiable**; 54,590 without a hue | Observer/classifier-conditioned photo composition; not a stable 50:50 biological equilibrium |
| Genus clustering; 20,725 classified taxa /1,935 nominal genera | **34.38%** of source taxon photo-white variation described by genus means, compared with **9.33%** under the original genus-size-preserving permutation reference | Strong association with nominal genus membership, *not* independently demonstrated phylogenetic conservatism |
| Moisture between versus within; **same 20,694 taxa/1,934 genera** | Geographically adjusted BIO12 slope of genus means **+8.555** versus congeneric species **+2.255 pp/source SD**; between-minus-within **+6.300** (ordinary paired genus bootstrap 95% **+3.666 to +8.779**) | Two observational gradients on identical taxa; difference is *not* a causal partition of adaptation |
| Congeneric sign mosaic; same 1,934 genera | Rainfall-positive covariance contributions **830**, negative **719**, photo-white invariant **385**; **53.58%** positive among the 1,549 variable genera | Slight positive directional imbalance but many counterexamples, not 830 significant genus-level selection coefficients |
| Environmental explanatory fit; **different 16,801 soil36 taxa/1,665 genera** | Thirteen geography, moisture, temperature, radiation, elevation and soil5 predictors jointly describe **0.6202%** of in-sample genus-centred photo-white variance | Small descriptive within-genus fit under coarse photographed source environment and selected complete cases |
| Prospective difference for two held-out congeneric taxa; same 16,801/1,665 | Extra rainfall/vapour reduces geography-only *pairwise difference* MSE by **0.21–0.26%** over four reused five-fold genus partitions; every fixed-prediction bootstrap95 includes zero | No reliable transferable new-genus species-*difference* rule established |

**The statistical non-equivalence is part of the result.** The genus-clustering panel has 20,725 taxa, the climate26 slope panel has 20,694, and the soil36 residual/pairwise panel has 16,801. The 34.38% genus-composition fraction and 0.6202% within-genus environmental fit are **not a common-sample variance decomposition**; neither should be divided by the other as a biological effect size. The +8.555 and +2.255 slopes *are* matched on the same taxa and environmental scaling, but target different covariance axes.

## New normalized check of the existing actual-source pairwise losses

The previous source replayer reports four whole-genus-heldout, five-fold predictions with all within-genus unordered species pairs counted equally **within genus**, then genera equally. Its pairwise predicted difference never uses the held-out genus's observed white-photo mean. The frozen source outcomes are consulted only to score those predictions.

| Genus assignment | Geography-only pair MSE | Geography + rain/vapour pair MSE | Relative loss reduction over geography |
| --- | ---: | ---: | ---: |
| Original | 0.316285 | 0.315620 | **0.210%** |
| Salt 202610151 | 0.316278 | 0.315454 | **0.260%** |
| Salt 202610152 | 0.316273 | 0.315477 | **0.252%** |
| Salt 202610153 | 0.316179 | 0.315450 | **0.230%** |

For the original heldout split the incremental pairwise MSE gain is **+0.000665**, ordinary whole-genus fixed-loss bootstrap95 **[−0.000227,+0.001506]**. The null-pairwise MSE is **0.316750**, so the combined geography + water model's reduction over that null is only **0.357%**, although the incremental water fraction above uses geography as its denominator. These are **two different baselines**. All four water-gain intervals cross zero; four salted splits re-use the same underlying taxa, and cannot be pooled as four independent ecosystem replicates. The soil5 increment over selected full climate is sign-unstable across the four splits. This analysis adds no new P-value and does not turn uncertainty into proof of no effect.

## Biological interpretation worth retaining in the Research Article

**Supported description:** Photographically white colour covaries with wetter conditions at both nominal genus-average and congeneric species-average scales across a broad sampled world flora. Genus identity is associated with substantially structured photographed-white composition, and the overall congeneric positive slope averages both positive and negative genus contributions. Coarse environmental measurements explain little of the within-genus photographic trait contrast and scarcely generalize to new-genus species-pair colour differences.

**Mechanisms remain distinct competing hypotheses:** (a) historical lineage and geographic species-pool assembly; (b) local pollinator or unmeasured microhabitat gradients; (c) pigment-linked water/light trade-offs under different reproductive and survival regimes; and (d) colour-classifier and photographic observation selection. This source tests **none** of these four as an identified causal mechanism. The separate 54.3% unclassified photographic opportunities and unequal soil36 completeness, especially at low elevations, prevent treating a near-half photo-label ratio as balancing selection or the observed photo-white fraction as an unselected global flora census. One site elevation does not measure hillside slope, aspect, TPI, TRI or water flow accumulation; no real linked DEM terrain effect has been estimated.

**Falsifiable next-data discriminator:** First obtain independently verified organ-level VIS/UV reflectance and occurrence localities on phylogenetically sampled, geographically replicated species. Compare lineage/region-predicted *composition* with paired congeneric *trait differences* in prespecified climate strata. Only subsequent common-garden and field fitness measurements of authenticated segregating flower-colour phenotypes—pollen export/receipt, seed quality and survival under water, radiation and pollinator manipulations—can test a pigment benefit/cost crossover. Additional analysis of the same photographic labels cannot establish that trade-off.

**Manuscript decision:** Retain the main positive global moisture–white pattern, but frame its distinctive ecological contribution as the **cross-scale mismatch between genus composition, congeneric covariance and poor pairwise transfer**. Keep published maxT intervals for the original 24-contrast family as the primary source statistics. The new MSE normalization is only a descriptive restatement of the already source-verified pairwise losses. The paper remains **Draft, unmerged and unsubmitted** pending licensed, stable anonymous reviewer data, source-location privacy, coauthor approval and production checks.

## Frozen evidence used

- `results/mois_original_genus_clustering_20261015.json`
- `results/mois_genus_between_within_moisture_20261015.json`
- `results/mois_congeneric_moisture_genus_sign_mosaic_20261015.json`
- `results/mois_env_shapley_and_genus_residual_transfer_20261015.json`
- `results/mois_pairwise_congeneric_photo_white_transfer_20261015.json`

All original source photo, photographer, species and site-identifying records remain outside public GitHub.
