# A positive global congeneric moisture slope is not a universal genus-level direction

**Original-source analysis checkpoint — 15 October 2026.** This report continues the interspecific photographic flower-colour paper and interprets the source-verified positive world genus-equal white–moisture association in terms of **heterogeneity among nominal genera**. It is a post hoc source-conditional sign/contribution analysis. It is not another source photo classifier, a phylogeny, independent ecological data collection, a pigment-metabolic measurement, or a balancing-selection experiment.

## Why test this?

The global species-level photographed-white/BIO12 coefficient is **+2.255 pp per original source environmental SD** among 20,694 qualified species and 1,934 nominal genera, with an original within-design 24-maxT 95% interval [+0.319,+4.190]. The genus-mean white/BIO12 coefficient on precisely those species was independently measured as +8.555 pp/SD. These positives could be mistakenly read as *most genera individually evolve or respond in the same way*. Genus membership averages conceal directionally divergent congeneric photo-colour/environment associations. The question is whether the globally positive genus-centred regression reflects predominantly positive genus contributions or the net balance of heterogeneous positive, negative, and zero source contributions.

## Exact original source and observational estimand

The two original source ZIPs were recovered and verified by full SHA256 and ZIP-integrity checks, with the photographic and climate26 join on original taxon, observation and photo IDs. The analysis holds the **same** eligibility support as the primary genus-equal model: 45,836 photo-climate26 complete original images, 20,694 qualified species, 1,934 genera and 38,498 classifiable photo outcomes in the qualified taxa. The original FCP four-state visible hue labels are unchanged.

For each named taxon, let `y_i` be its *proportion of classifiable original flower photos labelled white*. Both focal moisture axes (BIO12 annual precipitation and mean annual atmospheric vapour pressure) are normalized by the **original photo-site source SD**. Within each nominal genus, centre `y_i`, moisture and broad site geometry (absolute latitude and circular longitude sine/cosine). Each genus receives one equal total model weight. By the Frisch–Waugh–Lovell identity, partial the focal centred moisture against the globally pooled, genus-centred geography variables and allocate to genus `g` its weighted cross-product:

```text
c_g = sum_{species i in g} (1/n_g) * x_residual_i * (y_i - mean_g(y))
d_g = sum_{species i in g} (1/n_g) * x_residual_i^2

global slope = sum_g(c_g) / sum_g(d_g)
```

Here `c_g > 0` means a **positive contribution to the globally adjusted covariance numerator**, **not** a statistically significant genus-specific selection coefficient or even a free-standing fully conditioned genus-only regression. A genus whose classified photo-white fractions do not vary across its species has `c_g=0` by construction. This avoids spuriously assigning geographically adjusted contributions to groups with constant source colour.

The exploratory conditional null **independently permutes the observed species-photo-white fractions within each genus**, preserving genus sample sizes, the complete genus-level distribution of source photo-white fractions, the original species environment and coordinates, source photo depth and original hue-classifier output. It asks whether the number of positive genus contributions exceeds what that restricted *exchangeability model* predicts. Both moisture axes share the same 1,999 permutations for each species/genus qualification threshold. Original genus identities and site/photo rows remain private.

## Principal results on the identical main source species

| Outcome | Annual precipitation (BIO12) | Annual vapour pressure |
|---|---:|---:|
| World genus-equal original WLS slope (pp/source SD) | **+2.255** | **+3.952** |
| Positive contributing genera / all 1,934 | **830** | **841** |
| Negative contributing genera | 719 | 708 |
| No within-genus photographic white variation | 385 | 385 |
| Positive share among 1,549 *colour-variable* genera | **53.58%** | **54.29%** |
| Genus-preserving permutation mean positive share | 49.96% | 49.66% |
| Source-permutation 95% envelope | 47.58–52.49% | 47.19–52.23% |
| Exploratory upper-tail one-sided permutation *P* | **0.0040** | **0.0005** |

Thus the positive worldwide coefficient reflects a **modest excess of positive signs among genera with measurable within-genus photo-white variation**, not uniform positive direction: 719 (BIO12) and 708 (vapour) genera contribute in the opposite direction, while 385 have no colour variation in the selected photo measurements. Among *all* 1,934 eligible genera, 42.92% contribute positively to the BIO12 numerator, 37.18% negatively and 19.91% zero; the denominators are not interchangeable.

The **largest 10% of genera by positive FWL covariance magnitude** account for **60.28% (BIO12)** or **58.00% (vapour)** of the *positive covariance mass*. This indicates pronounced concentration of mathematical positive evidence. It is **not** the share of the overall regression coefficient causally generated by those genera, and positive and negative genus numerator terms cancel in the global weighted sum. No genus identities or coefficients are disclosed.

## Denominator sensitivity: more thoroughly photographed genera

| Minimum eligible sampled taxa per genus | Taxa / genera | Genera with variable photo-white | Positive among variable genera, BIO12 | Positive among variable genera, vapour |
|---|---:|---:|---:|---:|
| ≥3 (main source) | 20,694 /1,934 | 1,549 | 53.58% | 54.29% |
| ≥5 | 18,164 /1,184 | 1,058 | 53.78% | 53.02% |
| ≥10 | 14,130 /555 | 531 | 55.93% | 53.48% |
| ≥20 | 9,675 /225 | 223 | 60.99% | 60.09% |

For the ≥20 subset, the original source has 136 positive versus 87 negative variable-genus contributions for BIO12 (permutation mean positive fraction 49.84%; source-conditional exploratory *P*=0.0005). These larger-sample subsets also change the taxonomic and ecological source population, so their trend in positive fractions is **not** a causal sample-depth effect or independent biological replication. There are eight examined sign tests (two moisture axes × four thresholds), along with older exploratory analysis families; the displayed *P* values are **not corrected for the entire scientific search**.

## Ecological result and interpretation

**Supported in the studied photographic source:** The worldwide photographed-white/moisture slope is positive under the original broad genus/geographic adjustment and is a **net of heterogeneous genus-level positive, negative, and zero contributions**. A small direction imbalance persists under shuffling that preserves taxonomic photo-white composition, source climate and geographic support. More extensively sampled genera show a higher observed positive fraction but remain mixed. This is an ecological-scale fact about *photographically classified species-level traits*, not evidence that every genus follows a common universal rain-to-colour rule.

**Not identified:** The 385 genera with no classified-white variation may be biologically colour-conserved, sparsely or selectively photographed, uniformly mislabelled, or simply lacking observed hue contrasts. Opposing genus covariance signs could reflect species sorting, historical trait conservation, shared pollinators, correlated stress traits or systematic photo capture; their relative fitness effects are **not** separated. Frequency-dependent white–coloured polymorphism, stabilizing 50:50 equilibrium, pollen-reward advantages, and flavonoid pigment metabolic costs remain unmeasured.

**Biological follow-up:** Reconstruct reliable phylogenetic relationships and validated petal reflectance in a stratified sample of positive/negative photo-colour genera; test local community co-occurrence and pollinator guild variation at harmonized sites. A direct **fitness crossing** test requires actual genetically verified white/colour phenotypes, manipulations of environmental regime and morph frequency, and pollen/seed/survival outcomes. Additional post hoc photo-only filters cannot manufacture this causal evidence.

## Source, rights and manuscript placement

- Original-photo [SHA256] `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0` and environment [SHA256] `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca` were verified in the owner execution environment. Original licensed photo/site/taxon records are **not public GitHub content** and are not anonymous reviewer data.
- Complete source replayer `scripts/analysis/replay_mois_genus_sign_mosaic_20261015.py` uses the previously checked shared source loader and can be rerun with both owner-held ZIPs and `--write`/`--check`. The source was recomputed and aggregate-sign/shuffle results checked independently in the owner runtime, without changing original phenotypes.
- Public **aggregate-only numerical receipt** `results/mois_congeneric_moisture_genus_sign_mosaic_20261015.json`, plus synthetic numerical invariants in `tests/test_mois_genus_sign_mosaic_20261015.py`, carry the denominator and interpretation boundary.
- Manuscript: a concise observation in Results 3.2/Discussion 4.1; full methods and null model in Supporting Information S2.7. The original signed overall result and its scope-defined maxT intervals are **not replaced** by these exploratory genus sign statistics.
- Submission remains **Draft/not submitted** pending licensed anonymous source data/DOI, author approvals and final figure/reference production.

**Decision:** The novel ecological interpretation is *heterogeneous congeneric rain–white response combined into a positive macroecological average*. There is evidence of a small source-conditional excess of positive genus-level covariance contributions, not the same water-colour association in every lineage, and not a verified pigment cost/benefit trade-off.
