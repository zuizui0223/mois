# MOIS-only spatial colour composition beta diversity (2026-11-25)

## Scope and actual execution

Following the explicit instruction to use this database only, the Island external trait join was stopped. No Island records were read for this analysis. Both inputs are existing frozen MOIS numerical products, not new external downloads of biodiversity/environment observations. The WCVP regional status and GBIF distributions embedded in the MOIS products remain the previously integrated data.

Local empirical execution completed on the two signed ZIPs:
- MOIS artifact [14710329666](https://github.com/zuizui0223/mois/actions/runs/44914103213/artifacts/14710329666), SHA256 c031cf8fe3f80d6cf868ee20477fdb4c550b836435b79ae16d75e771c564f951; member native_photo9_georegion_cells.csv.gz, SHA256 e07c5906c10bc6944189053ee11ef35da6b76f0f93365b74ac0d13f1d71b327e.
- MOIS artifact [14400177816](https://github.com/zuizui0223/mois/actions/runs/44242421559/artifacts/14400177816), SHA256 cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd; member global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz, SHA256 449521d176ef8b1e029b5236ed464224d20b7538de1fb8d1f1552127d40c34ec.

Both outer ZIP and inner members matched the source receipts. Exact grid_x/grid_y one-to-one join of 18,577 cells; all nine hue species counts reproduce the supplied richness, strict-native pair count 2,271,852 and moderate-native count 2,604,870. No restricted taxon-level trait rows or images were used. No new GBIF, WCVP or Island source was joined.

## Methods

On the existing EPSG:6933 100-km grid, use unique horizontally or vertically adjacent pairs, not diagonal pairs, with minimum documented native richness 100 in both cells (30 as sensitivity). Each species is counted once per cell in the frozen photo-hue catalogue, not weighted by abundance. For probability vectors p and q, JSD = H((p+q)/2) − (H(p)+H(q))/2 (natural logs). Exact chain rule splits this into binary white/nonwhite JSD plus weighted conditional JSD within the nonwhite categories. Let w and v be nonwhite fractions of p and q and u=(w+v)/2; conditional component = u H(m_nonwhite/u) − (w H(p_nonwhite/w)+v H(q_nonwhite/v))/2, with zero-mass terms zero. All pair identities checked at tolerance 1e−12. Independent KL category-sum check maximum error 5.83e−16. Relative contribution = sum(conditional component)/sum(total JSD), not mean edge fractions. Resample 500-km projected blocks of edge midpoints 20,000 times; report descriptive percentile interval conditional on observed data, not a species-independence test. Compare equal-block weighting and within-block edges. Coordinate adjacency is not identical to geodesic separation.

## Results (actual local calculation)

| Cohort | Minimum richness per endpoint | Edges | Mean JSD | Nonwhite conditional share | Spatial resampling 95% range |
|---|---:|---:|---:|---:|---:|
| Strict native | 100 | 6,962 | 0.0040953880 | 87.5824% | 86.9343–88.2128% |
| Strict native | 30 | 10,153 | 0.0097794613 | 87.9598% | 87.3146–88.5607% |
| Moderate native | 100 | 7,939 | 0.0040057006 | 87.5827% | 86.9797–88.1668% |
| Strict native, bronze excluded, ≥100 retained species | 100 | 6,603 | 0.0035887732 | 84.3142% | 83.4948–85.1162% |
| Strict native, north | 100 | 4,525 | 0.0037286369 | 88.0645% | 87.2969–88.7992% |
| Strict native, south | 100 | 2,421 | 0.0047711379 | 86.9479% | 85.8695–88.0062% |
| Strict native, within 500-km block | 100 | 5,549 | 0.0041043895 | 87.6032% | 86.9139–88.2805% |

Main strict n≥100 estimate (548 midpoint blocks): binary mean JSD 0.00050854835, conditional mean JSD 0.00358683967, equal-block share 87.3456%. Same 6,603 edges retained after bronze exclusion: full9 share 87.3876%, no-bronze share 84.3142%. Thus bronze removal itself does not eliminate the observed dominance of the conditional term. These cohorts overlap, not independent replications.

### Dimension-matched diagnostic

Each of the nine categories versus all remaining categories is separately projected to a binary composition on the identical 6,962 edges. Mean binary JSD (x 1,000):
white 0.50855; yellow 0.51020; orange 0.50072; red 0.49588; pink 0.45190; magenta 0.58501; purple 0.53157; blue 0.48863; bronze 0.52061.

White-vs-rest is not an exceptional spatial contrast on the same-dimensional scale. These projections overlap and cannot be summed. The 87.6% result compares eight-category conditional variation to a binary component, so it is an information decomposition, not evidence that nonwhite biological selection is seven times stronger. Absolute JSD is small and sensitive to richness threshold. No richness-matched species-label null, measurement error calibration, phylogenetic test or between-hemisphere causal inference was carried out. Repeated species across cells mean the block interval is a sensitivity description, not a universal species-population confidence interval.

## Ecological interpretation and next step

A white/nonwhite-only map loses most of the observed nine-category local compositional information. This supports making full nine-hue community composition a central analysis in the MOIS paper, rather than treating nonwhite species as interchangeable. It does not yet identify climate adaptation, local floral colour evolution, species replacement versus shared species, or pollinator effects. The true floral-organ colour remains unvalidated and bronze is not verified brown petals.

The next MOIS-only analysis should hold the same adjacent cell pair set fixed while comparing observed JSD with synchronous genus × native footprint × geographic-affinity colour-label permutations on the existing MOIS species-by-cell source. This would separate descriptive dimensionality from excess spatial colour similarity relative to a defined taxonomic null. No new external floral trait source is required.

Code: `scripts/analysis/mois_internal_native_colour_beta_js_20261125.py`. Source ZIPs are passed as --native999 and --geo; --output must not exist. The empirical results were computed locally in two scripts and their logic consolidated into the repository script; no GitHub Actions execution or post-commit consolidated-script rerun is claimed here. Main unchanged; evidence draft PR #202 only.
