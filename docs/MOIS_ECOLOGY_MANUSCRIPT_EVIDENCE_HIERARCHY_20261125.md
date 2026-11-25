# MOIS — 全球花色組成の論文化・証拠階層（2026-11-25 JST）

## Status and study target

This is a cross-PR synthesis and arithmetic sensitivity audit of reported aggregate results, not a new empirical model fit, image measurement, independent source-byte verification, or botanical validation. Existing acquisition PR #200 and other analysis branches are unchanged. The primary estimand is the composition of **photo-derived visible-colour proxies of species recorded in regional plant assemblages**, not verified perianth colours, local population abundance, evolutionary transitions, or pollinator-mediated selection. Keep photographic observations, species-range summaries and species-by-cell assemblages separate.

The November 25 PR #200 audit reports 1,231,357 distinct formal photo IDs, 39,914 numerical taxa, 39,174 clear representatives and 37,105 strict accepted taxa. The ecological GBIF cube remains the frozen 37,104 accepted-species catalogue: 36,758 with valid geographic data, 5,493,931 species-by-100-km-cell pairs, 18,577 cells. Do not substitute the updated denominator for the frozen analysis or infer independent observations from photo IDs. No authenticated all-species GBIF SQL download is claimed (PR #201).

## Working title and central question

**Taxonomic and geographic context shapes global photographic flower-colour composition**

How much of the colour structure of recorded plant assemblages depends on the available regional species pool and its taxonomic composition, and does any association with environment remain when these are held more nearly constant?

Do not claim that environmental effects are absent or that composition causally mediates selection. Use a graded comparison across defined estimands and null pools.

## Evidence matrix (reported estimates; frozen denominators)

| Question / level | Result | Manuscript role |
|---|---|---|
| Observed photo assemblage, rainfall | 893,253 climate-complete observations, 350,081 strict ROI9 successes; technical MAR-IPW slope +1.314 pp/SD [0.163, 2.500], unweighted +1.595; supported target 871,767 observations | Positive descriptive result; not native flora or same-species response |
| Species-range rainfall | 36,581 species in full environmental panel; same-genus white contrast 29,381 species/2,365 genera; species-equal +2.294 pp/SD [0.090,4.499], genus-equal −0.078 [−3.514,3.358]; 18 hue-by-moisture comparisons none BH-supported | Context/weight sensitivity, not a universal moisture rule |
| Elevation, observed photo assemblage | Old near-neighbour −5.09 pp; on identical 350,081 observations taxon-equal white ROI slope −0.916 [−1.514,−0.319] to −0.145 [−1.006,+0.717] pp/SD after BIO5/log-rain/solar adjustment | Describe composition gradient; independent elevation association weak |
| Seasonal photo assemblage | South summer−winter +10.50 pp; common-species relative photographic frequency component +8.04 pp | Historical observation/phenology mixture; not individual plasticity |
| All GBIF assemblages, H9 | 12,728 cells with ≥30 species: observed H=1.199973; genus-only null H=1.268812, difference −0.068838 (p=.005); genus×occupancy-width null H=1.230773, difference −0.030799 (p=.010), 199 permutations | Main structural result, explicitly conditional on exchangeability of those nulls |
| Strict native assemblage, H9 | 7,231 cells: genus×all-GBIF width×geographic affinity deviation −0.024029 (p=.002, 999 permutations); native-only donor, width, geographic affinity deviation −0.012388 (p=.032, Holm4=.128) | Essential sensitivity: direction persists, confirmatory support does not under native-only pool |
| Island / high elevation, H9 | All-recorded flora island−continent residual +0.00138 (p=.935), high−low −0.01311 (p=.235); strict native island−continent +0.01008 (p=.555), high−low −0.00209 (p=.815) | No robust universal shift in entropy; distinct from white fractions; islands mostly large landmasses |
| Native vs introduced white, matched cells | 3,106 strict cells with ≥30 species per status, native 62.570%, introduced 65.547%; native−intro −2.977 pp [−3.700,−2.282]; minority-status richness weights −1.467 [−2.328,−0.620]; common-genus 2,516-cell subset weighted +0.140 [−0.779,+1.015] | Positive descriptive contrast, sensitive to pool and support |
| Native vs introduced, elevation | 337 high cells (≥1500 m) native−intro white −9.54 pp; 1,863 low cells (<500 m) −0.85 pp approximately; high−low −8.686 [−10.432,−6.722], minority richness weighted −9.128 [−10.948,−7.016]. High−low H9 residual contrast +0.13894 [+0.10733,+0.16687]; local 35 shared 500-km blocks +0.07315 [+0.02297,+0.12221]; adjusted altitude coefficient +0.00559 H/km [−0.02197,+0.03315] | Strong exploratory context-dependent contrast, not causal altitude effect |
| Same-cell same-genus colour matching | 1,363 cells, 236 common genera; high-elevation introduced−native matching +0.11133, spatial interval [.03870,.16401], crossed genus×space sensitivity [−.00373,.22693]; southern hemisphere +.09677, crossed [−.01589,.21323] | Dependence and repeated taxa limit generality |
| NPP/AGB/NDVI and LGM | Same 36,581 species, unseen-genus RF white AUC .562278→.561666→.562574; incremental intervals span zero; Hmsc GBIF record predictions improve with LGM in EU/NA, not EA, but not a consistent colour-specific benefit | Secondary negative/uncertain mechanism tests; never infer refugia or investment in flowers |
| Pollinators | GIFT independent modes 3,196 species, animal white 62.08%, wind white 71.36%; both modes only 11 genera/54 species. DoPI 10/63 same-publication comparisons BH-supported, 0/63 after same-publication common-genus restriction; grade1/2 only 0 | Independent ecological plausibility, not causal mechanism; positive visitation only |
| Introduction status & timing | WCVP red photo proxy within-genus +5.241 pp, q=.002362, but only 27 high-evidence red species and held-out genus logloss improvement uncertain. FirstRecords v4.0 2,799 species: historical baseline LL .543602, +white .543910, +9 hues .545110 | Exploratory supplement, not invasion advantage or post-introduction selection |

### New arithmetic synthesis of published aggregates

- In the same all-recorded 12,728-cell sample, adding occupancy-width constraints to genus-only label exchange reduces the magnitude of the H9 deficit by **55.3%** (1−0.030799/0.068838).
- Within the 7,231-cell strict-native analysis, changing from all-GBIF donor/width/affinity to strict-native donor/width/affinity reduces the magnitude by **48.4%** (1−0.012388/0.024029); corrected support is lost.
- Within the same historical observation regression, climate adjustment reduces the magnitude of the elevation-white ROI slope by **84.2%** (1−0.145/0.916).
- Within the same 3,106 matched strict cells, minority-status richness weighting reduces the magnitude of the native–introduced white difference by **50.7%** (1−1.467/2.977).

These are four separate, descriptive changes in point-estimate magnitudes under specified changes to baseline, adjustment or weighting. They are not four independent effect sizes, causal explained fractions, pooled meta-analysis, or new p-values. In particular the first two nulls condition on different exchangeability structures and the fourth changes the estimand.

## Manuscript structure (proposed; do not cherry-pick)

1. Global dataset and measurement scope: distinguish nominal taxon, accepted species, photo, observation and GBIF cell; show colour composition and evidence-status coverage (including uncertain/nonflower).
2. Environmental gradients at photo and species-range level: same-sample rainfall/elevation comparisons, genus weights, all nine classes, seasonal frequency decomposition.
3. Regional community structure: mapped H9 and white prevalence, explicitly paired null ladder and native-only sensitivity; positive all-recorded effect and conditional native result both visible.
4. Native–introduced regional contrast: same-cell estimates, elevation strata, common-genus support and crossed taxon/space sensitivity; avoid assigning high-elevation pattern to selection.
5. Discussion: independent pollination modes and visitation, modern productivity and LGM as measured candidate correlates with limited generalization. FirstRecords time analysis and exploratory red result in supplementary material.

### Draft abstract (not submission-ready)

The geographic distribution of flower colours is often interpreted through environmental and pollinator-mediated selection, yet regional plant assemblages differ in taxonomic composition and available species pools. We combined a large photographic colour-proxy catalogue with accepted plant names and global occurrence records to compare composition across photographic observations, species ranges and 100-km assemblages. In the historical climate-complete photographic dataset, white labels were more frequent under wetter conditions after adjustment for technical image-measurement response. At the species level, however, this relationship depended on the weighting of genera. Across 12,728 recorded assemblages, nine-colour entropy was lower than expected from within-genus and range-width-constrained label permutations, but the deficit weakened under native-only exchange pools and did not pass the designated multiplicity correction. Within 3,106 cells containing both recorded native and introduced plants, introduced species sets had higher white-label proportions, especially at high elevations; the size and interpretation of these contrasts depended on species-pool support, taxonomic replication and environmental adjustment. Productivity, historical climate and independently documented visitor guilds did not establish a general colour-specific mechanism. These results characterize context-dependent patterns in photographic colour proxies and emphasize the importance of species-pool definition in macroecological inference. Independent perianth-level validation is required before interpreting these patterns as botanical flower-colour composition.

## Next analysis gate (not yet executed)

Primary follow-up: repeat the **same native-only species-pool null ladder** and same-cell native–introduced elevation contrast on a prespecified, independently documented floral-colour subset (e.g. the Island direct high/medium, single-state colour records, exact accepted-name join; PR #197 describes the source). Do not use the photograph model to fill independent labels. Harmonize a priori only exact semantic colour classes, keep ambiguous/multicolour/bronze separate; report number of independently documented species, shared genera, status-richness thresholds, geography and selection coverage before testing. Compare photographic labels and independent labels on the identical retained taxa/cells, with the same null and synchronized species labels, before attributing any change to measurement. If joint coverage cannot sustain a native-only test, report insufficiency; do not relax thresholds after looking at effects. This is a validation/replication subset, not a population-representative truth dataset.

The existing blinded QA cohort in PR #204 is 180 species, 538 photos, **534** independent observations; 54-image first-look triage demonstrated nonflower risk, but expert perianth annotation remains 0/538 in that record. PR #200 November 25 audit further establishes all 113 later appended source records refer to /large.jpg, with 32 numerical strict9 successes, 29 clear representatives, 81 unknowns and **0 verified full-original-resolution successes within these 113**. Use "provider-rendered large images" for those rows; recorded payload hash is not full-original provenance. Historical full-image hash census is incomplete. No new labels or counts are asserted here.

## Evidence pointers (PR heads observed 2026-11-25)

- [PR #196](https://github.com/zuizui0223/mois/pull/196), [pinned moisture results](https://github.com/zuizui0223/mois/blob/8b80f8106ebfd8887b40e1fe58143de85a93c626/docs/MOIS_NOHUMAN_MOISTURE_REFIT_20261030.md).
- [PR #197](https://github.com/zuizui0223/mois/pull/197), [pinned composition results](https://github.com/zuizui0223/mois/blob/02cecc69f45ecdcaaae78fbddb47b58f3bb81ba6/docs/MOIS_COMPOSITION_ENVIRONMENT_RESULTS_20261030.md).
- [PR #200 audit](https://github.com/zuizui0223/mois/pull/200#issuecomment-7136714715), head a807659e828681ce2323c1d73df1fb747d498568.
- [PR #203](https://github.com/zuizui0223/mois/pull/203), [PR #204](https://github.com/zuizui0223/mois/pull/204), [PR #205](https://github.com/zuizui0223/mois/pull/205).
- [PR #206](https://github.com/zuizui0223/mois/pull/206), head 4e931f5cf6af8bdb85f0f77663e5d73c4f5dce30; [native-only pool result](https://github.com/zuizui0223/mois/blob/4e931f5cf6af8bdb85f0f77663e5d73c4f5dce30/docs/MOIS_NATIVE_FOOTPRINT_AND_INTRODUCED_RESULTS_20261122.md).
- [PR #207](https://github.com/zuizui0223/mois/pull/207). Earlier historical values retained in the November 17 evidence ledger in this branch.

No main merge, no new acquisition dispatch, no raw data redistribution, and no independent label or model rerun in this synthesis.
