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


## 7. Executed independent source/photo bias-repair trial, with actual postmeasurement effectiveness

**New signed evidence available 26 October 2026:** [real 2,048 genuine original-source photo ROI9 outcomes](https://github.com/zuizui0223/mois/releases/tag/mois-real2048-new-source-photo-geo-season-stratified-ROI9-20261026-v1) and [42,111-source-taxon BEFORE/AFTER audit](https://github.com/zuizui0223/mois/releases/tag/mois-original1045978-real401392-ROI9-42111taxon-bias-gap-delta-20261026-v1). Original 256-reference and 1,792-gap-targeted source observations were selected before opening image colours; eight 256-photo terminal units were processed by two successful original YOLO11/EfficientSAM worker jobs. The original collector failed solely on mismatched signed `design_cohort` versus `sampling_design_cohort` field names. [Data-only fix and production rerun](https://github.com/zuizui0223/mois/actions/runs/40320974849) recovered the original image terminals without any extra photo download or ROI inferencing. The original expensive image worker workflow's auto-push trigger has been disabled on its feature branch.

### Measured outcome

| Distinct original cohort | Attempts | Strict original ROI9 | ROI9 success | Prior zero-colour P2 attempts | Newly strict9-supported P2 source taxa |
|---|---:|---:|---:|---:|---:|
| Independent reference | 256 | 88 | 34.38% | 2 | 0 |
| Geo-season targeted | 1,792 | 594 | 33.15% | 258 | **85** |
| Combined two designs (NOT a combined probability sample) | **2,048** | **682** | **33.30%** | 260 | **85** |

Technical original photo outcomes: eight download failures, 1,358 acquired-but-4-hue unclassifiable, and 682 genuine strict9. Exactly zero 2,048-photo overlaps with the prior separate 64 genuinely new 500px-Lab pilot at both photo and observation UUID levels. The pilot method remains distinct and does **not** contribute an ROI9.

The signed **1,043,930 original attempts / 400,710 strict9** baseline plus the 2,048 source-disjoint real attempts gives **1,045,978 original photo attempts / 401,392 strict ROI9** (excluding the method-distinct 64). Within the fixed upstream 42,111 nominal source taxa, strict9-covered taxa rose **21,235 → 21,320 (+85)**; still missing strict9 fell **20,876 → 20,791**. A further **43 taxa** newly reached five strict9 images and **nine taxa** newly reached ten. 85 of the gains occur in taxa with **>=10 original source-photo opportunities and previously zero strict9** (P2).

### Important acquisition failure of the FIRST geographic balancing design

The geographic balance worked in the targeting sense: the targeted acquisition includes **347 N-tropical and 410 S-tropical** original attempts, delivering **122 and 144** new strict9 respectively. But taxonomic gap balancing was far weaker: among the 1,792 targeted photos, **788 (44.0%)** went to **P5 taxa already holding >=10 strict9 images**. P2, the genuinely unmeasured but eligible long tail, received just **258 (14.4%)**. Those 258 attempts alone produced **85 newly covered taxa (32.9% per attempted P2 photo)**. P3 received 485 and P4 260 attempts; P1 received one photo with no strict9. This is a **measured deficiency in the previous sampler objective**, not a generic suggestion to gather more photos.

### Next acquisition-policy gate: TAXON gap first, then region/season

* Keep the original frozen 256 reference and 1,792 targeted pilot distinct; do not retrofit universal frame probabilities. Any new reference arm remains separately independently randomized before colour selection.
* From the current **17,727 P2 still-uncoloured taxa**, rank *never-before-attempted observation UUIDs* within authenticated 93m source tier(s), stratified by **latitude zone × local season × source-taxon**; seek geographic and seasonal diversity among the P2 targets. Set **at least 60% of new targeted slots to P2** wherever independently sourced eligible P2 observations really exist, explicitly recording any shortfall from verified source shortages. Prefer taxa with 0–1 previous images attempted, and audit photo-specific ROI/flower-organ failures first for P2 taxa with >=5 historical failures.
* Do **not** allocate additional P5 *targeted* photos merely for geographic balance; allow P5 only in the separate reference cohort or a documented matching-support design. Use P3/P4 remainder for independent taxon × geography × local-season evidence; P1 limited by its rare source opportunity set, do not manufacture candidate IDs.
* Freeze one new independent observation per selected photo, stable before-image hash and source inclusion rules. Deduplicate against the **1,045,978 frozen image-attempt IDs**, the method-distinct 64 pilot, active FCP/K1 queues and any new reservations before downloading. Source-ID and observation-UUID disjointness require actual ID-level certificates, not simple arithmetic on prior group totals.
* Optimize **newly strict9-covered taxa and previously absent joint region/season support per compute-hour**, not raw photo count or the nominal fraction of northern/southern images. Keep all unclassifiable outcomes in each arm, separately estimate technical missingness, and publish before/after receipts.

**Execution boundary:** these are real executed 2,048-photo ROI and independently audited source-taxon improvements. The next >=60%-P2 source-only **photo-ID queue has NOT yet been selected or processed**, and no future ROI9 outcome for the next sample is claimed.
