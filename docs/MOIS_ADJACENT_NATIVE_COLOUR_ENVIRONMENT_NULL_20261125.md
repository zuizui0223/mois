# MOIS — Environmental contrast and native adjacent floral colour beta diversity (2026-11-25)

## Actual completed analysis, MOIS inputs only

Extends the 6,962 adjacent 100-km grid pairs of the preceding native-colour beta-diversity analysis. No Island trait data, new species occurrence sources, new photos, or ROI inference. Existing integrated GBIF/WCVP and environment inputs only. The local reconstruction used the preceding archived decoder and region map, including its documented three geographic corrections; all 2,271,852 native species-cell records, all 18,577 x 9 signed grid counts, 34,781 donor species, and 6,962 edges were reproduced. This is a count-matched replay, not a new full GEOS native-region overlay. The previous full9 and white-fixed global JS means for **all 999 draws each** were replayed with maximum numerical difference exactly 0.

### Sources

- Native reference artifact 14710329666, ZIP SHA256 c031cf8fe3f80d6cf868ee20477fdb4c550b836435b79ae16d75e771c564f951; member native_photo9_georegion_cells.csv.gz, SHA256 e07c5906c10bc6944189053ee11ef35da6b76f0f93365b74ac0d13f1d71b327e.
- Native elevation geography artifact 14400177816, ZIP SHA256 cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd; member global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz, SHA256 449521d176ef8b1e029b5236ed464224d20b7538de1fb8d1f1552127d40c34ec.
- Existing MOIS environment artifact [14397855889](https://github.com/zuizui0223/mois/actions/runs/44236695589/artifacts/14397855889), ZIP SHA256 4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6; member all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz, SHA256 8bb6b408b76771bdeb3682799e1ff748da833c7a8362906b665480405d65d0df. The environment member is the original preselected 12,295 complete-cell table, not all 18,577 cells. All 6,962 primary candidate edges have both endpoint temperature, precipitation and elevation.
- Reconstructed local native NPZ SHA256 b290f13acdca2264b17c378cdafc7af5d56464952788dbad8f0e710da9d18693 (raw native pair matrix not redistributed).

All three ZIP and member hashes checked against fixed pins and receipt. Coordinate joins one-to-one.

## Design

Primary outcome: total nine-colour adjacent-cell JS divergence (natural logs), multiplied by 1000. Main sample: endpoints each >=100 native species, horizontally/vertically adjacent, both in the same projected 500-km block, retain blocks with at least two edges: **5,492 edges, 462 blocks**. The covariates of interest are absolute endpoint differences in temperature, log1p(annual precipitation), and median land elevation; all three enter jointly, along with mean values of these environmental measures, mean log richness, absolute difference in log richness, species Jaccard dissimilarity and edge orientation. Standardize covariates by retained edge sample SD before block centering. The primary SDs are 1.40093 C, 0.161460 log1p(mm), 315.958 m, respectively. Block fixed effects are implemented by within-block centering. Same design and edge masks for every permutation. The two-sided rank test compares observed regression slopes against slopes from 999 globally synchronous species-colour label permutations, preserving genus, native footprint size bands, native six-sector affinity (>=65%), occurrence membership and richness. Seed 20261125. Report observed coefficient minus mean null coefficient, not an effect estimated as independent edge-level observations. Holm3 for the primary full9 total outcome across three environmental contrasts; Holm6 for two secondary components x three contrasts. White-fixed sensitivity holds all white identities constant; total and conditional nonwhite slope deviations there are identical, not independent tests. These choices were specified before this local numerical run, not an externally preregistered confirmatory test; the preceding results had already been inspected.

### Primary estimates: full9

| Environmental difference | Observed slope | Null mean slope | Excess slope | Two-sided p | Holm3 |
|---|---:|---:|---:|---:|---:|
| Temperature | .063993 | .085588 | -.021595 | .732 | .732 |
| log1p annual precipitation | .132644 | .053333 | +.079310 | .040 | .120 |
| Elevation | .134864 | .010918 | +.123947 | .064 | .128 |

No primary environment variable passes Holm3 at .05. The precipitation direction is exploratory. The null interval for precipitation slope is [-.0201333,.1280971], not an effect CI.

### Secondary components

| Component | Exposure | Excess slope | p | Holm6 |
|---|---|---:|---:|---:|
| White/nonwhite | Temperature | -.022188 | .340 | .680 |
| White/nonwhite | Precipitation | +.025923 | .044 | .220 |
| White/nonwhite | Elevation | +.045637 | .034 | .204 |
| Conditional nonwhite | Temperature | +.000593 | .998 | .998 |
| Conditional nonwhite | Precipitation | +.053387 | .142 | .568 |
| Conditional nonwhite | Elevation | +.078310 | .176 | .568 |

White-fixed conditional nonwhite sensitivity: temperature -.027548, p=.588; precipitation +.076042, p=.024 (Holm3=.072); elevation +.076182, p=.118 (Holm3=.236). Thus the precipitation tendency is not limited to white/nonwhite composition; still no corrected support.

### Design sensitivity (precipitation, excess total JS slope)

| Design | Edges | Blocks | Full9 excess | p | White-fixed excess | p |
|---|---:|---:|---:|---:|---:|---:|
| Primary 500 km | 5492 | 462 | .079310 | .040 | .076042 | .024 |
| Same edges, omit species-turnover covariate | 5492 | 462 | .083081 | .028 | .078438 | .018 |
| Within 1000-km block FE | 6235 | 180 | .081859 | .022 | .065900 | .018 |
| All edges, no spatial FE | 6962 | 1 global intercept | .060562 | .060 | .051621 | .040 |

The sensitivity tests are overlapping and not separately confirmatory. Elevation weakens under 1000-km FE (full9 p=.380; whitefixed p=.342). Do not infer altitude adaptation. The Jaccard turnover control changes the conditional estimand, not a causal mediator test.

## Additional conditional spatial holdout

Same primary 5,492 edges, 177 nonoverlapping projected 1000-km blocks assigned to five balanced folds (1098,1098,1098,1099,1099 edges); no grid endpoint appears in both train and test. Fit ordinary linear models with training-only scaling, no block FE at prediction (new block intercept unavailable). Target is observed minus null mean JS x1000. Baseline has mean climate/elevation, richness, Jaccard and orientation; add rainfall difference alone or all three differences. The null reference was calculated conditional on the full known catalogue: this is a conditional spatial residual prediction diagnostic, not independent species-level prediction or an independent validation of botanical colour. Added post-result.

| Outcome | Baseline MSE | +rain difference MSE | +all three MSE |
|---|---:|---:|---:|
| Full9 total residual | 3.947427384 | 3.945244459 | 3.947299122 |
| White-fixed nonwhite residual | 2.450840169 | 2.449576998 | 2.453171225 |

Full9 rainfall-only MSE reduction .002182925 (0.05530%), 20,000-block-bootstrap interval [-.004899155,+.009717982] conditional on OOF predictions; improvement in 3/5 folds. White-fixed corresponding .001263171 (0.05154%), interval [-.005114193,+.007734742], 3/5 folds. The full9 all-three model improves by only .000128261, interval [-.009465094,+.010162165]. White-fixed all-three model is worse. No demonstrated independent geographic predictive benefit. No claim of absence of an effect or equivalence.

## Verification and reproducibility

Local execution completed for 1,998 synchronous species-label draws, 72 coefficient rows, 2 x 6,962 edge aggregate rows and 10,984 held-out prediction rows. Exact replay of previous global draw means, independent rank p and mean checks for all 72 coefficients, additive JS component closure in slope arrays <=1e-15, independent MSE recomputation, source SHA and coordinate/count gates passed. The holdout output initially failed JSON serialization for numpy.int64; cast the block count to int and reran successfully, with numeric checks after output. The entire local code, source pins, numeric summaries, arrays, edge aggregates, predictions and SHA manifest are in conversation artifact MOIS_adjacent_environment_20261125.zip (2,581,612 bytes, SHA256 3aa14f1fc9865cc6567aa4c9dcc17dd1182cd80b9b39588b81210ba20eef25b1). The existing reconstruction source and previous draws are in MOIS_native_beta_null_20261125.zip. No raw native species-cell rows or photos exported. No GitHub Actions rerun claimed. Main unchanged; results recorded in draft PR #206.

## Ecological interpretation

The main global adjacent composition result remains near the genus/native footprint/affinity expectation. Local environmental dissimilarity may modify colour turnover, with a consistent exploratory precipitation slope across sensitivity designs, including when white species are fixed. But multiple-test correction and spatial prediction do not support a robust additional environment effect. This is a species-composition association conditional on the frozen photographic proxy, not direct evidence of pigment physiology, pollinator selection, local flower colour, or LGM refugia. A useful next MOIS-only diagnostic is colour-specific directional turnover along the moisture gradient, using wetter-minus-drier endpoint species proportions and synchronized nulls (not absolute JS alone), with a fixed multiplicity family and shared/exclusive membership checks.
