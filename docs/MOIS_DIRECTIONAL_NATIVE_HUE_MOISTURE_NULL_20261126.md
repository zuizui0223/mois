# MOIS — Directional moisture-associated change in native plant colour composition (2026-11-26 JST)

## Executed result: MOIS only

Follow-up to the November 25 adjacent-cell environmental JS analysis. Only existing MOIS data were used; no Island trait records or new photos, occurrence or environment sources. Source reconstruction via the previous archived decoder and region map again yielded 5,493,931 GBIF species-cell pairs and 396,142,825 records, 2,271,852 strict-native pairs, 34,781 native donor species and exact agreement for the 18,577 x 9 reference counts. The original GEOS overlay was not independently rerun, so equality of the aggregate counts is not a claim of byte-identical original memberships. Reconstructed native NPZ SHA256 b290f13acdca2264b17c378cdafc7af5d56464952788dbad8f0e710da9d18693 (same as previous).

Original adjacent cells: 6,962. Two exact precipitation ties were excluded before examining outcomes, giving 6,960 oriented edges. Primary within-500-km-block sample: 5,491 edges, 462 blocks, one fewer edge than the previous absolute-JS analysis because of a precipitation tie. All source ZIPs and selected members checked against the same fixed SHA256 and original receipts as the preceding report:
- native reference 14710329666: c031cf8fe3f80d6cf868ee20477fdb4c550b836435b79ae16d75e771c564f951
- geography/elevation 14400177816: cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd
- contemporary climate 14397855889: 4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6
- previous local reconstruction and environment bundles: aa1867ea665a6e95268c6308b29a63a03afe8623a4ccb4ae6bda95e565b0bc23 and 3aa14f1fc9865cc6567aa4c9dcc17dd1182cd80b9b39588b81210ba20eef25b1.

## Estimand and design

Response for each hue = 100*(species proportion in wetter endpoint minus species proportion in drier endpoint), not abundance, within-species colour or JS divergence. Primary coefficient: slope on positive log1p annual rainfall difference, per sample SD (.161465), while adjusting for signed wet-minus-dry temperature/elevation contrasts, mean temperature/log-rain/elevation, mean log richness, signed log-richness difference, species Jaccard dissimilarity and edge orientation; 500-km spatial FE. Covariates standardized before block centering. Pre-specified within this execution, not independent preregistration because previous results were known. Signed differences in confounders differ from the earlier absolute-JS design intentionally.

Strict null: synchronous species labels permuted within genus x strict-native occupied-grid breadth bin x six-zone native affinity (65% majority, otherwise mixed); 19,111 mobile species in 3,638 strata. 999 draws, seed 20261125; all fixed species memberships and source richness retained. A sparse species-weight operator makes the linear slope and mean contrasts fast; checked against direct hue count arithmetic and direct regression on the observed vectors and three sampled draws. Primary Holm9. Separate exploratory sensitivity families: no turnover adjustment, 1000-km blocks, no spatial FE, mean wet-minus-dry, same-support full vs exclusive species (>=30 exclusive species on each side), and white-fixed conditional composition with >=30 nonwhite species per endpoint (Holm8). Each uses the same fixed conditioning structure, not independent replicates.

## Primary numerical results

| Hue | Observed slope pp/SD | Null mean | Excess | Rank p | Holm9 |
|---|---:|---:|---:|---:|---:|
| white | +.139321 | +.119340 | +.019981 | .584 | 1 |
| yellow | -.106403 | -.107008 | +.000605 | .952 | 1 |
| orange | +.005651 | -.003239 | +.008890 | .222 | 1 |
| red | +.065772 | +.054921 | +.010851 | .432 | 1 |
| pink | -.025057 | -.019768 | -.005289 | .718 | 1 |
| magenta | +.002593 | -.000129 | +.002723 | .750 | 1 |
| purple | +.013387 | +.010938 | +.002448 | .858 | 1 |
| blue | -.024799 | -.016166 | -.008634 | .518 | 1 |
| bronze | -.070465 | -.038890 | -.031575 | .208 | 1 |

Observed positive white/red and negative yellow slopes largely reproduced by the conditional null. Yellow is nearly identical. White +.139 vs +.119 and red +.066 vs +.055 are not independent support for colour-specific sorting. The raw mean wet-minus-dry (same edges) white +.196792 pp, yellow -.228849 pp, red +.043149 pp; none of the nine mean contrasts passes the correction. No causal explained fraction inferred from point-estimate ratios.

## Sensitivities and joint vector test

| Design | Edges | Blocks | Minimum hue p | Minimum within-family Holm | Joint exploratory p |
|---|---:|---:|---:|---:|---:|
| Primary slope | 5491 | 462 | .208 | 1.000 | .834 |
| Omit turnover | 5491 | 462 | .108 | .972 | .698 |
| 1000 km | 6234 | 180 | .162 | 1.000 | .647 |
| No spatial FE | 6960 | 1 | .124 | 1.000 | .485 |
| Mean wet-minus-dry | 5491 | 462 | .212 | 1.000 | .723 |
| Exclusive-support full assemblage | 4723 | 439 | .170 | 1.000 | .705 |
| Exclusive species only | 4723 | 439 | .038 | .342 | .504 |
| White-fixed conditional nonwhite8 | 5473 | 460 | .132 | 1.000 | .595 |

Exclusive-species bronze slope -.241671 vs null -.118162, excess -.123509, p=.038, Holm9=.342; bronze is not a verified brown-petal class. This is not a robust independent hue effect. White-fixed donor support: 7,124 mobile species, 1,632 variable strata. All white identities and nonwhite denominators fixed, no independent hue supported.

The joint test was a post-result exploratory diagnostic for distributed weak directional change. For each design, use the observed vector and 999 null vectors together to compute a common mean and covariance; calculate the squared Mahalanobis distance for each, with a pseudoinverse to handle compositional closure. Symmetric across observed/null entries; upper rank p includes observed. Rank 8 for 9 hues, rank 7 for fixed-white conditional 8. It is conditional on the chosen exchangeability model and not an independently confirmed test. None indicates a clear directional deviation. No multiple-testing claim based on joint tests.

## Interpretation and next step

The previous exploratory precipitation effect on unsigned JS divergence (unadjusted p=.040, Holm3=.120) does not imply a globally coherent wet-versus-dry hue direction. This directional result cannot confirm a universal moisture-induced shift toward particular colours. It is compatible with taxonomic/geographic assembly generating most of the observed slope, but does not establish causal mediation, absence of selection or equivalence. A signed mean can cancel across regions; local contrasts might vary in direction. A useful next bounded analysis is to test regional heterogeneity of directional slopes using the existing six geographic sectors and joint synchronous permutation statistic, followed by leave-one-sector-out prediction; do not search individual hue x region combinations for uncorrected positives. Such analysis is proposed, not yet executed. No Island input.

## Verified outputs and limitations

1,998 actual synchronous draws, 71 coefficient rows, 8 joint statistics. Three representative draws per null reproduced the previous global JS vector exactly (maximum absolute error 0). Independently recalculated all coefficients from saved arrays, rank p, Holm and means (maximum numerical discrepancy 4.31e-16), maximum additive compositional closure error 4.83e-15. This does not independently revalidate the botanical photo hue proxy.

Conversation reproducibility archive MOIS_directional_moisture_20261126.zip, 572,086 bytes, SHA256 39728c13c7644b512628df14a2623be5f0f5ada56c853edee8688b19da705619. Includes actual executed directional.py (SHA 76ab939344d7c230404fe2c5c0d123a7314d01d4514255f22bea5ae3dc9890e0), directional_verify.py, reused env_beta.py, all full9/fixed-white draw arrays, 71-row CSV, summary and verification JSON, README and manifest. Both output NPZ hashes c18aae24c5b11a4f019661a5d54016c681992c68bdb3ff1fe54ea79ec214d415 and eef7cbee2258aa071a9270dd8518323297c72b4159f1af01be36f7a7b33c37cb. Code and arrays are in the local downloadable archive rather than the GitHub commit; this commit records the scientific report. ZIP members verified after packaging. No raw species-cell rows exported. No GitHub Actions rerun or main merge; existing PR #206 remains draft.
