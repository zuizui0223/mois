# Why geography coefficients change when atmospheric vapour is included

**Source-verified quantitative follow-up — 15 October 2026.** This document derives a weighted-least-squares omitted-variable identity for the original flower-photo white response. Its purpose is **clarity about covariate sharing**, not a search for new significant climate predictors or a claim that atmospheric water causally mediates altitude and flower-colour adaptation.

## Ecological question

The recent original-source `mois` geography analysis observed relatively fewer photographed-white flowers at high elevations in the *raw source photo* comparison, and a negative **genus-centred** white/elevation coefficient in a broad geographic model. After accounting for atmospheric vapour pressure, the elevation point coefficient became positive and its within-model confidence interval included zero. Absolute latitude also changed considerably.

**What mathematically explains those changing signs if the original 20,694 species, 1,934 genera, photo hue labels and model weights are exactly fixed?**

## Source and original population

The complete historical FCP measured-photo ZIP (SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`) and photo-site environment ZIP (SHA256 `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`) were independently SHA256-checked and reopened. Preserving the unique original photo ID intersection and original four-state classifier gate yielded **100,543** unique photographic opportunities and **45,836** climate26-qualified classifiable images. Requiring ≥3 taxa per nominal genus yields **20,694 distinct species /1,934 genera and 38,498 original classified photos**, precisely the source population of the main climate26 genus-centred result. The original white/BIO12 slope was rederived as **+2.254545 pp/source SD**.

For every eligible species, compute its original classified-photo white fraction and average original site elevation, absolute latitude, longitude sine/cosine, BIO12 rainfall and annual atmospheric vapour. Centre all predictors and the photographed-white response **within nominal genus**. Each sampled species receives weight `1/n_g`, so each genus contributes equal total regression weight. All explanatory columns use the original full photo-climate26 source SD. There are no new photo traits or georeferenced field observations.

## Exact coefficient-change algebra

Let `Y` be the centred white-photo fraction, `G` the genus-centred geography design matrix (absolute latitude, circular longitude sine/cosine, elevation), and `V` the genus-centred atmospheric vapour. Let `W` implement one equal vote per genus. When fitting models on the **same original species**, standard weighted-regression algebra gives:

```text
beta_geography_only = beta_geography_given_vapour
                      + beta_vapour_given_geography * gamma_vapour_on_geography

gamma_vapour_on_geography = (G' W G)^(-1) G' W V
```

The same relationship holds with BIO12 added to both versions of `G`. This is a numerical **omitted-variable identity**. It is *not* a counterfactual intervention, path-specific selection coefficient, climate mediation estimate or plant-physiology test. Source numerical consistency is better than 1×10⁻⁹ pp before rounding.

### Identical-species source results

| Fixed original-source model pair | White coefficient excluding vapour | White coefficient with vapour | Change by including vapour [paired genus-bootstrap 95%] |
|---|---:|---:|---:|
| Geography-only → plus vapour, altitude (+1,000 m) | **−1.903** | **+1.011** | **+2.914 [+1.110,+4.759]** pp |
| Geography-only → plus vapour, absolute latitude (+10°) | +0.483 | **+3.083** | **+2.600 [+0.995,+4.233]** pp |
| Geography+BIO12 → geography+BIO12+vapour, altitude (+1,000 m) | **−1.564** | +0.483 | **+2.047 [+0.148,+3.986]** pp |
| Geography+BIO12 → geography+BIO12+vapour, absolute latitude (+10°) | +1.101 | **+2.830** | **+1.729 [+0.127,+3.370]** pp |

In the geography-only version, vapour itself has a **+5.106 pp/source vapour SD** partial photo-white coefficient. Because standardized vapour decreases with both genus-centred latitude and elevation conditional on the other broad geographic covariates, leaving it out algebraically moves the latitude and elevation photo-white coefficients downward. The sign changes thus have an exact **conditional-covariate accounting explanation**, without showing which variable is biologically causal.

The weighted genus-centred source correlations of vapour with absolute latitude, elevation and annual precipitation are **−0.517, −0.489 and +0.445**, respectively. Geographic covariates jointly explain **74.53%** of the **weighted, within-genus vapour variation** in this exact source species panel (`VIF = 3.927`). Including BIO12 in that prediction increases shared explanation to **78.57%** (`VIF = 4.667`). These are diagnostics of overlapping source covariates, not estimates of physical photosynthetic moisture limitation, pigment budget or direct climatic selection.

The paired 95% intervals were calculated from **1,999 common genus multinomial bootstrap resamples** (frozen seed 202610151216) for each pair of nested models, so the changed coefficient is evaluated on aligned observations and resampling units. Intervals are **ordinary post hoc within-design percentiles**; they are *not* corrected for the multitude of geographic adjustments previously inspected or the study-wide source model search. The original 1999-bootstrap **geo4 maxT** intervals for each *individual fitted model* reported in S2.9 remain unchanged and should not be conflated with these paired shift intervals.

## What this adds ecologically

**Supported:** Absolute latitude, elevation and moisture are entangled across photographed plant species, including when species and geography are centred within sampled genera. Hence raw high-altitude coloured-flower representation and model-conditional high-latitude positive white fractions cannot be interpreted as universal independently acting physiological effects. In particular, changing the conditioning set can change a signed coefficient even though neither the plant photographs nor the studied species change.

**Not established:** That altitude *causes* flower pigmentation, atmospheric water pressure *mediates* an adaptive high-altitude effect, high-elevation flowers gain measured UV/thermal protection from anthocyanins, unpigmented flowers enjoy a known resource-saving advantage, or either colour is maintained at an evolved global 50:50 fitness equilibrium. Photographic classifiability is missing for >54% of the original photo opportunities; original photo sites and climatic averaging are non-probability sampled, and the source remains vulnerable to clade, location and observer structure.

**The paper's main biological result remains:** A broad congeneric, photographed-white association with source-site annual precipitation (+2.255 pp/SD; original 24-contrast within-family maxT interval +0.319 to +4.190) and atmospheric vapour (+3.952 [+1.199,+6.704]), alongside positive genus-mean gradients and heterogeneous per-genus contributions. The new coefficient identity is a supplementary *interpretive check*, not a reason to abandon that supported worldwide within-design association or imply stronger generality than the original image/geography robustness studies justify.

## Reproducibility and responsible source governance

The portable self-contained source replayer is `scripts/analysis/mois_geography_vapour_ovb_20261015.py`. It requires the rights-controlled **original photo and environment ZIPs**, validates both full SHA256 hashes, independently recreates the taxon-mean, genus-equal source WLS and its original white/BIO12 coefficient, then computes the exact identity and 1,999 paired source genus-bootstrap intervals. The `--check` command compares the owner-held source against committed **aggregate-only** receipt `results/mois_geography_vapour_ovb_20261015.json`.

The original owner-side source execution and independent numerical check passed using the historical ZIPs in the authorised runtime. Hosted CI validates syntax, source-aggregate invariants and synthetic algebraic identities **without holding those ZIPs**. The public receipt contains no named taxa, photo IDs, source coordinates, raw image data or identified genus-specific results.

The detailed existing raw-geography and conditional-geo4 comparisons are in [Supporting Information S2.9](MOIS_GEB_SUPPORTING_INFORMATION_20261014.md); the new exact identity is in **S2.10**. Manuscript remains GEB Draft PR #1 and **has not been submitted**. A stable licensed anonymous-reviewer dataset/DOI, original photographer rights, source geoprivacy and final coauthor/figure production remain open.
