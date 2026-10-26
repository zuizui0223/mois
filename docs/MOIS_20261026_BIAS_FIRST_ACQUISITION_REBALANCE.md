# MOIS bias-first acquisition: evidence, gaps, and selection gates (2026-10-26)

**Status:** source-backed audit and next-source *decision*, NOT a claim that new photo-ID targets or new 9-component floral colour pixels were acquired by this document. This proposal does not change or interrupt any running image workers.

## 1. Frozen evidence / incomparable denominators

- [Original 1,043,930-photo / 42,111 nominal source-taxon audit](https://github.com/zuizui0223/mois/releases/tag/mois-original1043930-source-structure-42111-taxon-colour-gaps-20261026-v1): 1,043,930 distinct attempted original photo IDs; 390,541 STRICT nine-component ROI colour outcomes; 42,111 nominal taxa in the **93,209,351-photo upstream opportunity frame**. A source-frame taxon ID is not necessarily a WCVP accepted species, a wild native individual, or a botanically verified flowering organ. Keep the 32,710 older attempts outside the 42,111-taxon source frame flagged rather than declaring them invalid observations or filling frame bins.
- [Original latitude ROI audit CSV](https://github.com/zuizui0223/mois/releases/download/mois-original1043930-source-structure-42111-taxon-colour-gaps-20261026-v1/mois_1043930_latitude_zone_photo_ROI_coverage.csv): northern temperate 670,525 attempts / 279,988 strict9 success (41.76%); northern tropical 85,111 / 21,431 (25.18%); southern temperate 194,446 / 62,505 (32.15%); southern tropical 46,312 / 7,354 (15.88%); southern subpolar just 29 / 4. Northern temperate photographs account for 71.69% of all strict9 successes. These are observation-photo opportunities, not expected global flora occurrence shares.
- [Island/mainland original photo audit CSV](https://github.com/zuizui0223/mois/releases/download/mois-original1043930-source-structure-42111-taxon-colour-gaps-20261026-v1/mois_1043930_island_mainland_photo_ROI.csv): island 137,076 / 45,101 (32.90% success), mainland 886,151 / 338,539 (38.20%). This is a technical/ecological mixture; island labels do not establish native island status.
- [Source-taxon opportunity bins](https://github.com/zuizui0223/mois/releases/download/mois-original1043930-source-structure-42111-taxon-colour-gaps-20261026-v1/mois_1043930_upstream_taxon_opportunity_coverage.csv) and [priority groups](https://github.com/zuizui0223/mois/releases/download/mois-original1043930-source-structure-42111-taxon-colour-gaps-20261026-v1/mois_42111_source_taxon_priority_group_summary.csv): 20,868 taxa with original-frame photo opportunities have **zero** strict9 measured outcomes (2,290,570 potential source photos; 2.46% of 93m). Another 16 nominal taxa have zero upstream photo opportunities. Only 5,278 taxa already with 10+ strict9 successes supply **343,296 / 390,541 = 87.9%** of all success measurements. Strict9 counts are not identical to the legacy four-hue counts.
- [Source-arm composition audit CSV](https://github.com/zuizui0223/mois/releases/download/mois-original1043930-source-structure-42111-taxon-colour-gaps-20261026-v1/mois_1043930_photo_source_composition.csv): historical FCP fixed-100-taxon cohort contributes 227,107 strict9 successes from 447,500 attempts, or **58.15% of all strict9**. The historical non-FCP main cohort has 114,522 strict9 successes from 464,267 attempts (24.67%). Do not pool heterogeneous source arms as if they were one probability sample.
- [Earlier disjoint new-wave seasonal audit](https://github.com/zuizui0223/mois/blob/main/docs/MOIS_FOUR_DISJOINT_TIER_SOURCE_COVERAGE_20261020.md): April–July = 55.06% of **276,437 pre-colour targets**, NOT proven 55.06% of the latest 1,043,930 attempts or strict9 successes. An updated source+ROI month/region cross-tab is still missing.
- [Actual 64-new-photo bounded geographically balanced pilot](https://github.com/zuizui0223/mois/releases/tag/mois-real20261026-novel64-source-geobalanced-original500px-floralLab-v1): 64 new disjoint source photo/observation IDs, chosen from 64 taxa with >=100 source photo opportunities and no previous strict9 outcomes, eight each in four geographic sectors × two previous-attempt-effort strata. New 500px YOLO flower-box Lab descriptors **are not** validated petal pixels, strict9, or original frozen YOLO11/EfficientSAM ROI-v4 outcomes. Do not promote or combine these fields.

## 2. Priority for acquisition based on *missing effective evidence*, not total images

| Source-frame stratum | Nominal taxa | Source photo opportunities | Previous strict9 | Acquisition decision |
|---|---:|---:|---:|---|
| P0: zero source opportunities | 16 | 0 | 0 | Cannot solve with current 93m photo roster; flag coverage limitation |
| P1: rare (1–9 source photos), no strict9 | 3,050 | 14,878 | 0 | Preserve priority / census feasible existing photo opportunities; record irreducible sampling sparsity |
| P2: >=10 source photos, no strict9 | **17,818** | **2,275,692** | 0 | **Primary taxon gap**; prefer taxa with 0–1 prior attempts; for >=5 prior ROI failures, audit flower organ/ROI before redownloading |
| P3: 1–4 strict9 photos | 12,923 | 10,431,283 | 28,283 | High priority for independent observation-level replication, especially missing geographic/season support |
| P4: 5–9 strict9 photos | 3,026 | 6,884,970 | 18,962 | Limited support-matched reinforcement; avoid inflation of already-supported regions |
| P5: >=10 strict9 photos | 5,278 | 73,602,528 | 343,296 | **No automatic additional enrichment**; retain in independent probability-sample control and paired comparisons |

These strata reflect **existing outcomes** and thus define an outcome-conditioned targeted *coverage* arm, not an unbiased 93m global prevalence probability sample. Never invent an unconditional sampling probability for historic photo cohorts.

## 3. Next source-only gap matrix before selecting a single additional image

Create a non-destructive availability table from **all original source candidates** and **all existing or reserved photos/observations**, at the grain:

`accepted_taxon_or_unresolved_source_ID × equal-area land grid × local season / month × observation_year × source region × source-arm × technical ROI state`

Retain hierarchy `photo_id → observation UUID → observer → taxon → genus → family`. Geographic latitude-zone totals are screening diagnostics; they **must not** become a 50:50 geographic flora target. Extract native/introduced/wild statuses separately; WCVP TDWG3 “native in region” cannot certify island nativeness or an individual wild plant.

Include in every joint stratum: `frame_photo_opportunities`, `unattempted_original_photo_ids`, `unattempted_observation_ids`, `previous_attempts`, `strict9_valid`, `four_hue_valid`, `failed_download`, `ROI_gate_unknown`, `ambiguous_palette`, `distinct_observers`, `distinct_accepted_taxa`, `classified_local_season_years`, `reserved_in_active_workers`, and `source_selection_probability_known`. Flag missing joint intersections: marginal region and taxon totals cannot establish that eligible undermeasured taxon–region–season cells are actually available. Do not create empty synthetic cells or presume a missing photo represents a biologically absent species.

Generate a before-colour **candidate-only** queue that is:
1. First anchored on P2 with genuinely unattempted observation IDs and P3 sites/seasons lacking independent replication, then P1 within their finite availability and P4 where joint support is weak. Hard-stop P5 targeted enrichment except a **separate pre-colour randomized control**.
2. Prioritized within available joints toward southern / northern tropics and underrepresented island, realm, biome, winter/local-season and fine-grid/accepted-taxon combinations; quantify the support first, and do not equalize photos to falsely assumed biological-area or phylogenetic baselines. Use original *local* season, not raw month without hemisphere.
3. Restricted to one photo from each new observation and no IDs already attempted, reserved by other active jobs, or sourced by the 64-image pilot. Stable hash within eligible observations; source rules frozen before any new hue result; include inclusion probability if derivable or mark targeted probability unidentified.
4. Logged with full upstream frame denominator, source release asset SHA256, per-row photo/observation/taxon/region/month source, original 9-colour measurement status, and a ledger of unavailable joints. No acquired photo without intact provenance.

## 4. Correct ROI missingness simultaneously

K1 original ROI-v4: [76,753 source attempts](https://github.com/zuizui0223/mois/releases/tag/mois-real-K1-76753-original-photo-nine-palette-and-unknown-classification-pilot-v1) include **46,288 original ROI/flip gate unknowns**, dwarfing the ambiguous positive palette subgroup (2,510). The southern tropical strict9 success rate (15.88%) and northern temperate rate (41.76%) differ substantially, but causal reason is unknown.

Audit a *pre-specified outcome-blind sample* of original K1 and historical failed images split by latitude, island/mainland, taxon, month, and source arm. Human review classes: no open flower visible; open flower visible but ROI missed; background or other organ detected; many/small flowers; rotation/flip or quality gate; photo acquisition failure. Only if errors are recoverable evaluate alternate full-resolution flower segmentation in a pilot against unmodified ROI-v4. Where no flower exists, a new source photo may be needed; where a flower is visible but segmentation fails, **reuse existing photo pixels** rather than spend acquisition effort on more of the same sampling bias. Keep the 64-image YOLO box Lab descriptors distinct from strict9 data.

## 5. Two cohorts / three report denominators, no indiscriminate pooling

- **Reference probabilistic cohort:** source-selected independently of measured hue, with observation and photo inclusion probabilities, including machine-classification failures. Design-corrected estimates only when actual inclusion probabilities/positivity are valid. Do not retrofit them onto the historic FCP taxa-first cohort.
- **Targeted gap-filling cohort:** prior strict9 failures and joint-support deficits determine selection; report improvement in coverage and ROI recoverability, not worldwide white-flower prevalence.
- **Botanically verified validation subset:** independent petal/anthesis annotations stratified by region/taxon/ROI status. This estimates detection and hue/organ measurement error, and prevents green leaf/brown tissue being interpreted as petal trait.

Publish before/after at frozen comparable grain: **fraction of 42,111 nominal taxa with >=1 valid strict9**; fraction with >=3 independent observations; number of supported equal-area taxon/genus×local-season-year cells; strict9 success per attempted original photo by geographic region and island flag; Gini/top-taxa success concentration; ROI failure reason by stratum; duplicates/positivity; robustness of existing white/water/elevation/island/season results to classification and weighting. For 9-colour trends retain continuous colour vectors and all missing statuses. No claim that a larger raw sample alone fixes bias.

## 6. Release/CI gate for an executable acquisition queue

Do not dispatch more full-res workers until a candidate queue is signed with:
- exact new source photo/observation IDs not appearing in 1,043,930 audited attempts, currently running and reserved photo jobs, or 64-photo pilot (full overlap proof, including legacy numeric UUID crosswalk uncertainty);
- authenticated all-frame source and 42,111-taxon opportunity table SHA, and exact counts before/after dedup and geographic/taxon/season sampleability;
- known photo/observation selection probability **where identifiable**, otherwise an explicit `targeted_nonprobability` estimator flag; phase-inference classified versus unclassified status remains observed, never imputed;
- separate fixed-source primary reference and targeted sample interpretations; holdout for confirmation; cost/throughput budget and quality-stop criteria;
- immutable machine-readable `source_gap_matrix.csv`, `targeted_original_photo_queue.csv`, `selection_receipt.json`, `coverage_before_after.csv` and explicit `pending_not_measured` states.

**Current deliverable:** evidence-based acquisition redirection and fail-closed queue contract. **Not yet executed:** full cross-tab source availability / new photo-ID roster / new strict9 ROI-v4 measurement. The 64-image pilot was already genuinely measured, but only with its explicitly different 500px YOLO flower-box method.
