# Original-source worldwide moisture slopes: leave-one-genus-out influence sensitivity

**Scope:** retrospective *descriptive influence audit*, 14 October 2026. This note verifies whether a **single named source genus**, rather than many nominal genera in aggregate, dominates the two original worldwide genus-equal photo-colour/moisture point coefficients. It is *not* a new confidence interval or a new ecological source of flowers. Original model outcomes and all historical significance decisions remain frozen.

## Computation

The original historical source (100,543 deduplicated photo IDs/42,111 named taxa, 45,836 climate-complete classifiable photo rows) supplied **20,694** eligible named species in **1,934** qualifying nominal genera. Original genus-equal regression used within-genus centred species-mean photo-colour fractions, broad `abs(lat), sin(lon), cos(lon)` adjustments and the source-global environmental SD. Two pre-existing central moisture axes were considered: annual precipitation `wc_bio12` and mean atmospheric water vapour `wc_vapr_annual_mean`.

A private original-source-derived, *anonymous by nominal genus* `1934×9×9` crossproduct table and its `1934×9×4` response products enabled deleting one entire genus at a time and re-solving the originally fitted geographic-plus-single-environment normal equations, **with the original source scaling held fixed**. Each deletion removes an entire qualifying nominal genus while keeping the original species-colour/outcome sampling and model definition otherwise unchanged.

This is a local statistical influence test with **1,934 deletion refits for each of the two original axes**, each returning four photographed hue coefficients. Exactly **15,472** leave-one-genus-out hue estimates were calculated. Full-score original coefficients match the source 1,280-row model file. The original source input/group moment SHA256 is `77dff4acc5756610e463bec3e5877c711f1695d8984be0b675e86f7d8afb4aae`, but the private group moment matrix has **not** been publicly deposited or approved for external redistribution.

## Results

| Photo-colour outcome | Full-sample rain/BIO12 slope (pp/source SD) | Across 1,934 single-genus deletions | Full-sample vapour slope | Across 1,934 deletions |
|---|---:|---:|---:|---:|
| White | +2.255 | **+2.174 to +2.339** | +3.952 | **+3.843 to +4.068** |
| Yellow/orange | −1.221 | −1.320 to −1.133 | −1.385 | −1.546 to −1.226 |
| Red/pink | −0.446 | −0.506 to −0.380 | −1.509 | −1.587 to −1.412 |
| Blue/purple | −0.588 | −0.663 to −0.520 | −1.057 | −1.167 to −0.976 |

**All 1,934 deletion fits preserved the sign of all four colour-category points on both original moisture axes.** Maximum absolute white photo-colour coefficient change induced by deleting any one nominal genus was **0.08481 percentage points** for rain and **0.11643** for vapour. Across all 8 photographed hue coefficients, the single largest absolute deletion change was **0.16075 pp**, for yellow/orange with vapour.

The original source-scope **maxT 95% intervals remain those of the full model and were not re-estimated in the single-genus-deletion experiment**. A sign-stable point estimate can still be statistically uncertain, especially for individual yellow/orange, red/pink and blue/purple hues (whose original simultaneous intervals include zero). No new post hoc positive-significance finding is introduced.

## Ecological interpretation and reproducibility boundary

Within the specified broad, genus-equal photographic model, the two signed white–moisture point estimates are not explained *solely* by the removal/inclusion of one especially influential genus. This does **not** mean the ecological response is phylogenetically conserved across genera, that the photographed white–rainfall slope transfers to a different geographical flora, or that no collection of related genera or sampling artifacts can distort the result. The original strict-image and same-genus×coarse-cell source sensitivity models, missing-colour MNAR, botanically sourced hue nonconfirmation and distant geographic prediction remain material limitations.

- Public **aggregate-only** result receipt: `results/mois_world_primary_moisture_single_genus_influence_20261014.json`, retaining only whole-source point estimates, min/max deletion point ranges, 95th percentile absolute deletion differences and number of sign reversals. It contains **no original genus identities, photo IDs, original detailed coordinates or per-genus coefficient rows**.
- Privately runnable source programme: `scripts/analysis/loo_mois_world_genus_moisture_20261014.py --private-moments [authorized NPZ] --original-ledger results/fcp_macro_direction_20261013/direction_ledger.csv`.
- Public CI verifies the source aggregate receipt against the exact original frozen CSV and science nonclaim; it **cannot independently recreate the 1,934 group exclusions without the privately held sufficient-statistic moment NPZ**.
- The anonymous genus moments encode some original geographical covariance despite stripped source IDs, and were not uploaded to GitHub or placed in a stable reviewer archive. External sharing requires **documented human licence/geoprivacy authorization**.

**Next scientific question unaffected:** verifying the broad slope cannot replace field pigment and individual/reproductive phenotype evidence needed to infer water-dependent selection. No assertion of genetic white-coloured morph maintenance is supported.
