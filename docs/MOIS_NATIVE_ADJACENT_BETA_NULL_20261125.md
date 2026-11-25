# MOIS — Native adjacent colour beta diversity under synchronous species-label nulls (2026-11-25)

## Completed empirical analysis; MOIS data only

This extends the adjacent-cell Jensen–Shannon (JS) decomposition in PR #202, not the proposed Island trait join. No Island records or new trait/occurrence/image acquisition were used. The frozen MOIS GBIF catalogue and WCVP native-range layers already used for PR #206 are retained.

- [GBIF species/cell source #14166701658](https://github.com/zuizui0223/mois/actions/artifacts/14166701658): ZIP SHA256 `422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2`; species-cell Parquet SHA256 `e00b0408055984abfe4ec0ff0a24833fc148888afd1ea0e0c9b270a077d08b11`; coverage catalogue CSV SHA256 `6ce4ac5289dccfdfce3dae946de3b1617bbc646eba76de75086c59795e5a4603`.
- [Native reference #14710329666](https://github.com/zuizui0223/mois/actions/runs/44914103213/artifacts/14710329666): ZIP SHA256 `c031cf8fe3f80d6cf868ee20477fdb4c550b836435b79ae16d75e771c564f951`; native nine-colour count/coordinate-sector CSV SHA256 `e07c5906c10bc6944189053ee11ef35da6b76f0f93365b74ac0d13f1d71b327e`.
- Previously used MOIS WCVP source #12174290403 ZIP SHA256 `c10d8fc7bf39b0cb404d5eaa09e3e520596b8d1a9284e643dd007d5b485a93cf`; both name and regional-status member SHA values verified against the original receipt.
- Original TDWG3 geometry is the same previously pinned source commit 52da7828aba9d461dd133c27b3bd7a4407161f54 and Git blob 91104e5159e31f88154833a51d3b0d1c9271083f. The exact source blob was retrieved via GitHub, not a new occurrence/trait source.

### Reconstruction gates and actual execution limitations

PyArrow was unavailable locally. A narrow Parquet decoder for the original dictionary/plain, ZSTD, non-null columns recovered **5,493,931 distinct species-cell rows**, 396,142,825 filtered occurrence records, and exactly the original per-species footprint and occurrence totals for all 37,104 accepted source species. Reconstructed the native status using original WCVP exact accepted name crosswalk (36,762 matches) and the previous TDWG3 25-subpoint threshold (>=20 unambiguous same-region hits). The alternative scanline overlay initially differed at 3 cells; rechecking their 25 projected points gave MAG 23/25 and 25/25 at the longitude wrap cells (-174,67) and (173,67), and MDG 18/25 at (47,-23), below the strict threshold. Corrected those specific projection/boundary assignments by the independent geometry check, not by choosing a colour match. Final assigned strict cells: 9,288. Final native pairs: **2,271,852**, 34,781 native donor species. All **18,577 × 9** native counts exactly match the signed reference, with zero discrepancies. This is strong aggregate replay validation, not a claim of byte-identical original species-cell membership or a fresh full GEOS make_valid overlay; the local geometry reconstruction used the pinned Git blob, without independent SHA256 revalidation of the entire geometry file.

The actual numerical run used six models, 999 global simultaneous species-label draws each, seed 20261125 (independent random keys sorted within each fixed stratum). No re-estimation of hue, location or native status in a permutation. Genus, fixed native footprint size bins (1,2–5,6–20,21–100,101+) and coordinate-sector affinity (>=65% of native cells in one of six sectors, otherwise mixed) define nested donor strata. Labels are permuted once globally across species, not independently per cell. For white-fixed model, all white species identities and their cell counts are invariant. Scope remains photo-colour proxies.

## 1. Adjacent whole assemblages

Horizontal/vertical adjacency in EPSG:6933 projected cells, both endpoints >=100 native species. The same 6,962 edges as the descriptive beta analysis. Observed mean JS = 0.004095388025525331 (natural logs); binary white/nonwhite = 0.000508548352937373; weighted conditional nonwhite = 0.003586839672587940.

| Null model | Exchangeable species | Variable strata | Mean total JS ×1000 | Observed minus null / null | Two-sided rank p |
|---|---:|---:|---:|---:|---:|
| Native genus only | 28,461 | 2,482 | 4.007195 | +2.20% | .176 |
| Native genus × footprint | 23,211 | 3,414 | 4.023302 | +1.79% | .206 |
| Native genus × footprint × affinity (primary) | 19,111 | 3,638 | 4.025834 | +1.73% | .178 |
| Same strict null, white species fixed | 7,124 | 1,632 | 4.046850 | +1.20% | .232 |

The primary total null 95% reference interval is [0.0039316608673, 0.0041279706544]. All total deviations are positive, not lower as predicted by excess local colour similarity. None resolves a difference from the conditional null. Do not claim equivalence or that colour cannot matter. Regional primary total tests: northern 4,525 edges p=.208; southern 2,421 p=.488; within 500-km blocks 5,549 p=.188. These overlap the global sample and are exploratory.

Binary white/nonwhite JS excess is +8.61% (genus p=.020), +6.73% (genus×footprint p=.028), then +3.84% (strict p=.142). Weighted conditional nonwhite excess under strict is +1.43% (p=.320); with white identities fixed +1.37% (p=.232). No stable strict whole-assemblage binary or nonwhite effect.

The same strict draws give mean within-cell H9 (>=30 species) 1.2152685262 observed versus 1.2275551271 null, difference -0.01228660, close to the previously reported native-only pool result -0.012388. Local alpha concentration and adjacent beta dissimilarity are different estimands and do not logically imply one another.

## 2. Exclusive-species diagnostic

Shared species necessarily share the fixed catalogue label; to test whether shared species conceal a colour effect, remove all common species from both endpoints. Across the original 6,962 edges mean Jaccard species similarity is 45.45% (median 44.25%). Restrict to **6,024** edges with >=30 exclusive species per endpoint; mean exclusive left/right and shared species counts 147.55/143.25/286.66. This is a different selected estimand, not a direct same-sample ratio with the full-edge result. Maintain all species identities and native footprint strata in the permutation.

| Component | Observed JS ×1000 | Strict null mean ×1000 | Relative difference | Rank p |
|---|---:|---:|---:|---:|
| Total 9 | 25.066200 | 24.939435 | +0.51% | .482 |
| Binary white/nonwhite | 3.360192 | 3.198474 | +5.06% | .020 |
| Conditional nonwhite | 21.706008 | 21.740961 | -0.16% | .836 |
| Conditional nonwhite, white-fixed null | 21.706008 | 21.694850 | +0.05% | .920 |

White-only excess is an exploratory signal: it is not the primary test, whole-assemblage strict binary p=.142, and even a 3-component Holm adjustment in the exclusive diagnostic would give p=.060 for the smallest p=.020. No independent claim of white selection. In particular, the descriptive 87.6% nonwhite contribution is not evidence of stronger biological filtering. Exchanging exclusive species also remains conditional on the source catalogue and taxonomic strata; no local adaptation, phylogenetic causation or botanical colour calibration is identified.

## 3. Evidence and reproducibility

Six actual permutation arrays (4 × 999 × 4 groups × 3 metrics and 2 × 999 × 3 metrics) were saved locally and independently checked finite. Maximum total-minus-binary-minus-conditional closure error 4.43e-17. The corresponding aggregate array NPZ SHA256 values:
- genus: 95ea916822b510c4812e5f3d5b77e1ba112312185b48f4ee0c9792e88ba0d2bc
- genus_native_width: fd953b6d05b4b65d4e5aabba771f16dc5064f69de07ea88de0a84fb8ea1d6ece
- genus_native_width_zone: 246984ace1c9978bb0d281e977c23e87716f9b3360ea0f549b62ea6c1dac4894
- genus_native_width_zone_whitefixed: 89ade147ebd2af1f31ad4f1986f121bd15fa55f4a98be45e85d409445168281f
- exclusive_full9: efa807d01b8a4c4dd958d2a7b378d6c4866772997aefbfb98768649c04f0de64
- exclusive_whitefixed: 1a655223ccf293edfc31af00e01031c1f18266e696caa1820237e26395c0e2d5

Combined repository analysis script: `scripts/analysis/mois_native_adjacent_beta_null_20261125.py` (requires native NPZ reconstructed from the original source, with the signed reference CSV and exact source count gates). The actual local execution scripts, alternative decoder, region map, numeric summaries and all six draw arrays were bundled in the conversation artifact `MOIS_native_beta_null_20261125.zip`, SHA256 `aa1867ea665a6e95268c6308b29a63a03afe8623a4ccb4ae6bda95e565b0bc23`. No raw species-by-cell rows were exported. Repository script consolidation was not rerun in GitHub Actions; do not claim CI success or byte-identical consolidated rerun. Main remains unchanged; this extends the existing ecology draft PR #206.

## Scientific next step

Keep primary beta result (no extra local similarity under defined null) and alpha result distinct. A direct test of environmental gradients in beta residuals is a logical MOIS-only extension: retain the same adjacency and shared/exclusive membership, calculate signed null-centered per-edge JS residuals, compare environmental contrast (existing elevation, precipitation and temperature only) against a geographically blocked baseline, and use a synchronized global permutation to test the slope rather than treating edges as independent. This is a proposed next step, not a completed result.
