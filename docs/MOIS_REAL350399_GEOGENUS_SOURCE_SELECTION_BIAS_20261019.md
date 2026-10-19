# MOIS: taxon–geography matching changes the apparent WHITE excess in real photographed angiosperms

**Data freeze:** 19 Oct 2026; original-photo ROI pixels measured (not only photo-ID opportunities). [Executed authenticated GH Actions run](https://github.com/zuizui0223/mois/actions/runs/39216663652); [frozen complete result JSON and three CSV tables](https://github.com/zuizui0223/mois/releases/tag/mois-real350399-geogenus-matched-source-selection-v1).

## Research question and answer

Among 350,399 unique, genuinely original-image-measured, global plant photographs, is the near-half automatically WHITE-labelled composition insensitive to changing the *taxa and geographic sites* that are preferentially photographed?

**No: the method of recruiting photo taxa matters.** The full actual-pixel source has 164,297 classifiable four-hue ROIs, of which 78,026 are WHITE. The original deeply sampled core has 317,631 photo attempts and 152,827 classifiable hues, of which exactly 72,175 WHITE (47.2266%). The deliberately species/genus-expansion arms together add 32,768 genuine photo attempts with 11,470 classifiable hues, of which 5,851 WHITE (51.0113%). **The raw new-minus-core difference is +3.7847 percentage points**, under **different source-recruitment and missingness regimes**; it is not a randomized experiment.

The population of photographs is not a population of naturally blooming individuals. The 350,399 attempts contain 186,102 technically unknown/failed or unclassified hue outcomes, with separate technical acquisition failures retained. No unknown image ROI was treated as chromatic or counted as a flower without a classification. White is a visible ROI hue, not automatically no pigment or UV nonabsorbing.

## Taxon–geography matched results

We predeclared geographic 5°, 10° and 20° latitude/longitude grid cells (not equal-area) and required named taxa from BOTH photo-acquisition arms for a named genus × geographic cell to enter a match. Within each matched stratum, mean WHITE frequency is calculated per **nominal named source taxon** first, then averaged across taxa, to prevent a single highly photographed species from dominating. Match-specific new-minus-old differences are averaged by **nominal genus**. The 95% intervals use 999 exploratory genus-resampling iterations and do not remove phylogenetic or geographic dependence.

| Grid size | Threshold: named taxa per arm/stratum | Matched genus × cell strata | Matched genera | Equal-genus WHITE difference, pp | Exploratory genus bootstrap 95%, pp |
|---|---:|---:|---:|---:|---:|
| 5° | ≥1 | 1,513 | 1,193 | +0.924 | [−2.046, +3.860] |
| 10° | ≥1 | 2,009 | 1,556 | +0.635 | [−1.847, +3.070] |
| 20° | ≥1 | 2,933 | 2,134 | +1.823 | [−0.242, +3.887] |
| 5° | ≥2 | 144 | 143 | −1.316 | [−7.324, +4.945] |
| 10° | ≥2 | 241 | 237 | +1.565 | [−3.152, +6.470] |
| 20° | ≥2 | 572 | 551 | +1.419 | [−1.481, +4.488] |

The apparently strong **+3.78 pp raw contrast does not remain a well-identified signed effect** after matching genus and geography. This is a comparison of **different, much smaller matched-support source subpopulations**, not an algebraic proof that geography and genus causally account for exactly 3.15 points. Moreover, the taxon-equal within matched groups is a *different estimand* from raw photo weighting even before matching.

### Geography-only versus genus×geography: a scientifically important diagnostic

This *separate* actual-photograph comparison matches 5°, 10°, and 20° cells
**without requiring the same named genus**, still giving equal weight to
nominal taxa within each arm and then an equal weight to every matched
geographic cell. It is not the same estimand as the genus-equal test.

| Matching design | 5° cells | 10° cells | 20° cells |
|---|---:|---:|---:|
| Geographic cells with both source arms, ≥1 nominal taxon/arm | 630 | 259 | 99 |
| Geographic-only new-minus-core WHITE pp, ≥1 | **+5.503** | **+6.199** | **+4.630** |
| Geographic-only new-minus-core WHITE pp, ≥2 taxa/arm | +5.256 | +5.987 | +6.030 |
| Matched *genus* × geography, genus-equal WHITE pp, ≥1 | +0.924 | +0.635 | +1.823 |

The geographic-only 10° comparator retains **98.01% of core and 99.78% of
new taxon×cell occurrences**, whereas 10° genus×cell matching admits
just 6.05% / 24.90%. Geography-only equal-cell matching makes the new
source WHITE excess *larger*, not smaller. Requiring taxonomic support
changes both the average and the comparison population. This directly
supports the importance of interspecific photographic assemblage
composition and warns against a universal abiotic explanation.

Do **NOT** subtract the 10° +0.635 from +6.199 and claim that 5.56 pp is
the causal contribution of taxon turnover: geography-only versus
genus-and-geography estimators have different weighting and overlap.
A fully identified causal decomposition would need shared source
selection probabilities and adequate same-taxon/time/site replication,
which these species-first acquisitions do not contain.

### Why this is not proof of a universal unchanged ratio

- Of all classifiable named source taxon IDs, only **14** are observed in both deeply sampled and new species-first photo cohorts; **3,285 nominal genus names** overlap globally. Virtually no directly comparable species-level color observations were collected in the two arms. The same genus in a cell may still contain entirely different species.
- In the 10° genus×cell ≥1 matched design, **only 6.05% of core taxon–cell occurrences and 24.90% of new taxon–cell occurrences** qualify for the matched-support estimator. At ≥2 taxa per arm these shrink to **1.46%** and **5.53%**. This is a serious positivity/transport limitation, NOT evidence of a null biological colour effect.
- Even 5°–20° cells may include distinct populations, habitats, elevation, solar/rain conditions and observation years. The arm difference could also reflect original per-image classifiability, introduced plants, taxonomic naming, species turnover, observers and open-flower stage.
- Actual flower petal/anthesis has not been independently botanically verified at this full scale. UV appearance, anthocyanin chemistry, metabolic production costs, pollination and reproductive fitness are unmeasured.

## Consequence for the ecological article

**Supported descriptive insight:** although WHITE has appeared near half of automatically classifiable global flower-like ROIs, the measured distribution is sensitive to the definition and weighting of the photographed species/genera. A +3.78 pp shift from focusing acquisition on new taxa becomes between roughly +0.6 and +1.8 pp under genus-and-geography restricted standardizations at ≥1 taxon/arm/cell, with exploratory intervals crossing zero. Thus equal photo frequency, equal taxon frequency and genus/geocell-weighted patterns must remain distinct across macroecological interpretations.

**Not identified:** a 50:50 maintained evolutionary strategy, pigment production cost versus pollination benefit, lineage-independent selection, or an invariant worldwide WHITE frequency. The incoming entire 93,209,351-photo-ID source is an observation-photo sampling frame only; actual photo pixels are needed before color outcomes can be estimated. Old source arms lack a defensible retrospective global photo inclusion probability and cannot be directly inverse-probability weighted to the world flora.

**Next prioritized empirical test:** annual moisture/heat/solar effects on *true classified photo hues* with geographically and taxonomically overlapping source strata, evaluating omitted/missing strata explicitly rather than selecting one positive global direction. Environmental covariates must be source-matched, with coordinate uncertainty and technical missingness kept distinct. Seasonal observed-year responses remain supplemental and must not be substituted with long-term monthly normals.

## Source-specific reproducibility

- Analysis: `scripts/analysis/mois_real350399_geogenus_arm_white_balance_20261019.py`
- Contracts: `tests/test_mois_real350399_geogenus_arm_white_balance_20261019.py`
- Signed CI pipeline: `.github/workflows/mois-350399-geogenus-source-bias.yml`
- Durable machine-readable result: [Geo/genus matched v1](https://github.com/zuizui0223/mois/releases/tag/mois-real350399-geogenus-matched-source-selection-v1). JSON includes all five original source photo Parquet SHA256 values, matched-support denominators, four-hue contrasts, geographic-only comparisons and original photo quality attempts.
- Source photo/cohort exact duplicates audited in [previous 350,399-photo paper checkpoint](MOIS_REAL350399_ACTUAL_RGB_FOURTERM_SEASON_AND_TAXON_20261019.md).
