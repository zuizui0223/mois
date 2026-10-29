# MOIS — 人手不要・元画像ROI9欠測バイアスの自動解析（2026-10-29）

**Execution policy: NO manual annotation, NO human review gates, NO new photograph acquisition and NO re-inference of original image ROIs.** This replaces the earlier manual 1,200+240 review proposal on draft PR #196. Those proposed human-review scripts/workflow were removed from the branch, not made mandatory elsewhere. Existing signed October 27 GEB manuscript data and external sampling arms remain unchanged.

## 1. What is frozen and already known

- Signed source: **93,209,351** photo ID opportunities, **53,652,241** source observations and **42,111 nominal taxa**. This is an iNaturalist opportunity frame, not a census of worldwide blooming native species or flowering abundance.
- Actual SHA-certified four-arm genuine original-photo union: **1,181,146** image attempts and **451,643** strict nine-palette ROI9 successes, leaving **729,503** technical unknowns. Use the corrected lat/long/date v2 only:
  https://github.com/zuizui0223/mois/releases/tag/mois-1181146-originalROI9-128bulk-georestored-seasons-and-512QC-20261027-v2
- Outcome-blind canonical one photo per original observation UUID: **1,032,882 independent represented observations**, **401,703** strict9 successes, **631,179** strict9 unknown, **148,264** extra photos of already represented observations. The original full photo QC frame is not deleted.
- Actual measured nine-component colour vectors were not modified. For the automatic **white-photo-label** estimator, success is original \`strict9==True\` *and* original \`original_roi_four_hue=="white"\`. The other three recorded hue classes are nonwhite. Any colour text from a **strict9 failure is always ignored**, even if another non-strict pipeline had a four-colour prediction.
- The 500px flower-box Lab pilot does not become flower-petal ground truth, and no manual-reviewed corolla/open flower or wild/native identity is assumed. Never pool targeted and source-reference arms unweighted and describe the result as worldwide natural white-flower prevalence.

## 2. Distinguish estimands instead of fixing every bias in one weight

**T0. Frozen automatic photograph labels, R=1**. Exactly measured, but conditional on strict9 classifiability. Within each original source arm, report raw fraction among success photos and taxon-equal among taxa with success.

**T1. Counterfactual photograph-label composition of the *selected original observation frame* if classifiability were not restrictive**. Not identified without assumptions because >60% of canonical observations are technical unknown. Apply area/season/arm poststratification under **explicit conditional MAR** only; complement it with assumption-free mathematical bounds and arbitrary MNAR odds shifts. Because some records may lack a visible open corolla, even these are **photographic label conventions**, NOT the fraction of natural blooming individuals.

**T2. Natural worldwide native plant community colour composition**. Not identified with iNaturalist photo opportunities: independent verified native species/cell inventories, source coverage and natural abundance/occurrence information would be needed. No auto ROI-success reweighting can manufacture this estimand. Do not rename T1 as T2.

## 3. Automatic method now implemented

Source implementation: \`scripts/analysis/mois_roi9_observation_frame_20261029.py\`, \`scripts/analysis/mois_observation_roi9_auto_bias_20261029.py\`.

1. Read the **SHA-matched source v2 parquet**. Fail closed if canonical photo count / ROI9-success / receipt source checksum differs. Select one photo per observation by pre-ROI, pre-hue SHA hash of observation UUID/photo ID. Keep all real records in an untouched source QC ledger.
2. Partition selected observations into **source arm × broad latitude region × local recorded season**. Regions: north/south tropical, north/south extratropical and pooled polar/coordinate-unknown. Seasons derive from actual observation date and hemisphere; preserve unknown.
3. Count each stratum: \`N_h\` all canonical observations, \`K_h\` strict original9 successes, \`U_h=N_h-K_h\` strict unknown, \`W_h\` actually recorded strict9 WHITE label among successes.
4. Report actual complete-case \`sum(W_h)/sum(K_h)\`; a conditional-MAR poststandardized *photographic* rate \`sum(N_h W_h/K_h)/sum(N_h)\` **only if every included stratum has K_h>0**. If a stratum has K_h=0, report an interval and the missing-support count instead of making up a point. Flag thin strata with K_h<30. Also report arms, zones, seasons and supported-taxon-equal summaries.
5. For all observations without a measured colour, compute worst-case mathematical bounds:
   \`[sum(W_h)/sum(N_h), (sum(W_h)+sum(U_h))/sum(N_h)]\`.
   These bounds deliberately make no claim about which photo really shows an open flower.
6. MNAR sensitivity: for technical unknowns in strata with observed support, vary the hypothetical unknown-white odds versus success odds at fixed log-odds shifts \`delta ∈ {-2,-1,-0.5,0,+0.5,+1,+2}\`. A Jeffreys-smoothed success reference \`(W_h+0.5)/(K_h+1)\` is used **only for this scenario grid** so tiny strata do not imply infinite odds; the MAR point itself uses unsmoothed measured fractions. Unobserved-support strata remain intervals even after assigning delta. These are scenario outcomes, not recovered empirical pigment measurements.

Both \`T1\` and the MNAR bounds apply **only to the four-arm source-selected observation frame**. Missing-at-random is untestable from the strict9 success-only labels; the scenarios show where conclusions can change. Technical flower absence, underdetected petals, exposure and genuine non-white flower are not distinguished by statistical correction.

## 4. What is auto-verified and what remains a separate scientific task

The workflow \`.github/workflows/mois-auto-original-roi9-bias-20261029.yml\` downloads only pre-existing SHA-audited Parquet + receipt, runs invariant tests and frozen source calculations, and emits machine-readable JSON and CSV. It does **not** download any photo RGB and does **not** require one human label.

Guardrails:
- Complete source/photo/observation identity and exact known successful nine-colour counts; preserve all unknowns as unknown.
- Support count, source arm and geography/season strata, plus explicit gap coverage.
- Outcome-blind canonical photo selector checked against shuffled input and colour-label inversion.
- MAR correction explicitly conditional on the chosen metadata strata; source observer selection/flowering opportunity bias **not** removed.
- No bootstrap or confidence interval is assigned to the partially identified MNAR scenario as though it were a known true observation. The older climate/moisture slope from the original GEB paper is not silently replaced by these descriptive photographic calculations.
- White photographic label is not UV/pigment chemistry, fitness, or an independently verified open wild-native flower.

**Stop criterion:** When the automatic estimator and sensitivity report have succeeded on all 1,032,882 represented source observations, do not hold MOIS hostage to additional manual review or 90m full-image inference. Move to an independently sourced environmental join and spatial-block analysis only when those inputs/estimands are ready.

## 5. Evidence vs assumptions

| Claim | Status |
| --- | --- |
| 1,181,146 real original images, 451,643 ROI9 success | Observed |
| Canonical 1,032,882 observations, 401,703 ROI9 success | Observed via deterministic source-only identity |
| Four-hue WHITE label among strict9 successes | Observed technical/photographic label |
| Missing success probability by known arm/area/season | Observed |
| Strict9 unknown WHITE distribution from those strata | Conditional assumption (MAR or explicitly stated MNAR sensitivity) |
| Open flower, corolla identity, native status, full 93m/world flora colour prevalence | Not identified |

The previous **1,200+240 human review** is explicitly canceled as a dependency. Any earlier downloaded annotation UI is legacy and not used in this program.
