# MOIS: 画像枚数を増やさずに観察・ROI欠測バイアスを識別するための設計（2026-10-29）

Status: **probabilistic observation-level sample and its blinded reviewer interface actually executed 29 Oct 2026**; human botanical review and ecological bias correction **NOT yet performed**. Keep the existing 27 Oct GEB single-paper frozen evidence intact. No replacement of strict9 unknown with a proxy.

### Executed source-only freeze (not biological estimates)

- [Authenticated 1.18m-source sampling run](https://github.com/zuizui0223/mois/actions/runs/40839008744), [blinded 1200+240 viewer run](https://github.com/zuizui0223/mois/actions/runs/40839599535); implemented in PR #196 without changing `main` or original colours.
- Exactly **1,181,146** original attempted photo IDs, all with observation UUID; **1,032,882** distinct observations after choosing a canonical original photo by preregistered SHA rank (independent of any hue or technical success); **148,264** extra photos of already represented observations suppressed from inference *without deleting the original photo QC ledger*.
- Among the canonical photos **401,703 original strict9 successes** and **631,179 strict9 unknown**. These differ from photo-level 451,643/729,503 precisely because photos are not independent observations.
- Exact review draw: **400 strict9 successes and 800 strict9 unknown**, **160 nonempty probability strata**, stratified on **source arm × broad latitude zone × recorded local season × original strict9 status**. 79 success strata, 81 unknown strata. Every nonempty audit frame stratum has `q_h = n_h/N_h > 0`; the observed frame is a mixed selected acquisition population, not globally probability-representative flora.
- Output **source-locked** roster SHA256 `01d882f9e7964a2573baea722523bf361495a44c0e8d73485e58053c94288d79`, original-photo-v2 Parquet SHA256 `faf6951505f1002a1327f64a2ced13431b562cf13893707125fafff918160f62`. Audit original labels are BLANK.
- An offline HTML review interface exists with **1,200 primary blinded images and 240 stratified reviewer-overlap images** (80 ROI successes, 160 unknown), loading only public medium images when displayed. It hides ROI status, estimated colours, source arm, region and season. NO human-verified flower presence, colour distribution or bias-corrected moisture slope has been obtained by building this interface.
- Estimator implementation is a separate **fail-closed** stage: refuse the point estimate if any human labels remain unfinished, use partial-identification bounds for uncertain flower/colour records, only estimate a source-observation-frame E1 ratio if completed labels allow it. No correction to the full 93m source or natural wild/native flora is identified.



## 0. Frozen source and what is actually identified

- Upstream: 93,209,351 source photo IDs / 53,652,241 observations / 42,111 nominal taxa. These are iNaturalist photographic opportunities, **not a randomized census of the world's wild flowering plants**.
- Frozen exact four-arm measured-photo union: 1,181,146 genuine attempted original photo IDs; strict original nine-component ROI = 451,643 (38.2%); strict9 unavailable = 729,503 (61.8%). **Use the latitude/date-restored v2** in DOI/release: https://github.com/zuizui0223/mois/releases/tag/mois-1181146-originalROI9-128bulk-georestored-seasons-and-512QC-20261027-v2
- N temperate strict9 success 42.0% vs S tropical 23.8%; technical missingness is nonuniform.
- Of the historical baseline, 148,264 photos are additional images from observations already represented by other photos. A photo ID is not an independent observation.
- A deliberately selected 512-photo subset drawn **only** from strict9 unknown photos has 495/512 YOLO500 flower boxes, but no botanical ground truth. The 128-image nested human review panel is diagnosis-oriented and was published with **zero annotations**; neither panel identifies worldwide missing-label frequencies. https://github.com/zuizui0223/mois/blob/main/docs/MOIS_ROI9_512_REAL_EXISTING_PHOTO_DETECTOR_FAILURE_DIAGNOSIS_20261027.md
- Independent 131,072 paired old ROI9 / 500px BOX photos contain 95 instances of 500px detector negative with strict original ROI9 success. Never drop all detector negatives without a prospective probability-based negative audit.
- Acquisition arms have different and in older historical arms sometimes unrecoverable first-stage selection probabilities. Photo-source selection is not natural plant abundance; **never report all-arm inverse-probability-weighted worldwide wild-native white-flora prevalence**.

## 1. Freeze THREE estimands before additional analysis

**E0 (descriptive; supportable now):** White/other visible-colour distribution **among successfully automated colour-labelled photographs**, reported per acquisition arm with own denominators. This is explicitly conditional on technical selection R=1 and is not a population estimate.

**E1 (primary correction target):** Among all selected **distinct observations** in a defensible source-reference sampling frame, and conditional on a *human-confirmed visible open flower*, the proportion with human-identified white corolla. Target population is the **source-covered iNaturalist observation frame**, not world flora. Compare by climate zone, place, year/local season. Source selection probability must actually be known for a source-wide estimate, otherwise report within-cohort only.

**E2 (macroecological secondary target):** For each supported equal-area spatial cell × time period, average **one taxon × cell × period** flower-colour estimate per taxon rather than letting photographed observations define taxon abundance. Report species-equal and, separately, genus-equal macro composition. Supported taxa and cells are those actually identified and estimable; unknown source taxa and lack of native-status support must be shown. This is *not* automatically the native flora of a cell.

A claim about **wild native community composition** requires an *independent native species-occurrence/checklist sampling frame* for each cell and coverage/occurrence probabilities; WCVP/native annotations or GBIF photo opportunities alone do not provide true local native abundance. Until externally identified, label E2 as "sampled taxon–cell assemblage" not natural floral community composition.

## 2. Reuse 1.18m images; one small probability-based MANUAL calibration audit

### 2.1 Audit frame (no new photo-ID recruitment, no inference rerun)

Use the frozen v2 source photo ledger. Before looking at any colour outcome or manual label, collapse photos to **one deterministic SHA-hash-selected original photo ID per observation UUID**, including original source arm, from all represented observations. Preserve a photo-level ledger for existing outcomes separately. If an observation straddles arms, record its multi-arm provenance; for probability inference use a single **reference-arm eligible** frame with known joint observation/photo inclusion or otherwise a source-arm-conditional descriptive target. Preserve the original 1.18m denominator in technical QC; audit inferential denominator is the independently selected **observation-level audit frame**, which must be counted, never guessed.

A probabilistic audit design must be **separate** from the existing 512/128 purposefully enriched failure-diagnosis samples. The older samples can inform failure taxonomy and interface design, not a population estimate; do not silently pool them as random.

### 2.2 Fixed first calibration sample: 1,200 already-attempted observations

- Select **400 strict9-success** and **800 strict9-unknown** original images, conditional on the observation-level frame above, with outcome-blind random within-stratum sampling. The deliberate failure oversample is permitted **only because its sampling probabilities are recorded**.
- **Implemented strata:** original source-design × broad latitude zone (N/S tropical, N/S extratropical, or a predeclared pooled polar/unknown category) × recorded local season (including an explicit unknown) × strict9 status. For the 160 empirically nonempty strata, allocate at least two cases when frame size and target budget permit; one case is only permitted in a one-observation census cell. Assign remaining cases by capped proportional allocation to reach exact 400/800 totals. Record source counts N_h and review counts n_h; the conditional **q_i=n_h/N_h** is nonzero for every represented source stratum. No new stratum allocation can be chosen after human outcomes are inspected.
- Do **not** force one per taxon/observer within a stratum after the draw: that breaks simple q_h. Instead record taxon, genus and observer clustering for subsequent uncertainty and coverage diagnostics. Avoid drawing more than one image from any original observation by construction.
- Retain unretrievable photos as a distinct inaccessible status; redraw only according to a recorded, probability-preserving replacement protocol (ideally no redraw; report access loss). Report accessible-denominator coverage.
- 1,200 is a **starting design budget**, not a claimed power guarantee for 32+ interactions. The reference arm alone needs enough eligible observations and label variation for its estimates. Release cell counts and uncertainty; borrow strength or aggregate where sparse rather than pretending fine geographic precision.
- Dedicate the frozen 128-image unknown-only sample to **mechanistic failure-mode diagnosis**, *not* to overall frequency estimation. Avoid re-annotating an existing audit ID unless the overlap and inclusion probabilities are explicitly resolved.

### 2.3 Blind human rubric

Two independently coded channels, masked to climate and model/ROI predicted hue:

1. Botanical: confirmed open floral corolla visible / no open flower / ambiguous; native/cultivated/unknown is a **separate** status based on external provenance, never assumed from a flower box.
2. Measurement: corolla colour visible (white, yellow/orange, red/pink, blue/purple, multicolour fractions or uncertain), clear floral organ vs bract/leaf/fruit/background; cause of original failure = flower absent, flower too small/occluded, wrong detector box, segmentation/flip gate, exposure/lighting, access failure, other.

Annotate a prespecified 20% overlap with a second blinded reviewer and adjudicate disagreements without concealing initial ratings. Store photo ID, annotation version, reviewer agreement, cause and exact missing/uncertain statuses. "YOLO flower box" != true blooming flower; "strict9 valid" != manually verified corolla.

## 3. Correction that is identifiable (and what cannot be corrected)

Let F_i=1 if a botanically confirmed open corolla is visible, and Y_i=1 if that corolla is visually white. Let A_i=1 if included in the random manual audit with known conditional q_i. For a fixed source/region/period g, calculate the finite-frame **direct two-phase audited ratio**:

    p_hat(g) = [sum_{audited i in g} (F_i * Y_i / q_i)] /
               [sum_{audited i in g} (F_i / q_i)]

If source-stage observation/photo inclusion pi_i is verified *for the specified upstream photographic frame*, replace 1/q_i with 1/(pi_i q_i). Otherwise the estimator only applies to the actual selected/cohort frame. Never invent retrospective pi for the pooled historical/targeted arms. Human-uncertain labels and inaccessible photos must remain separate; prespecify bounded alternatives, not zero imputation.

For precision, use **cross-fit calibration/GREG** combining frozen ROI9 and other universally recorded metadata as auxiliary predictors. At the audit-frame total level:

    T_hat(YF) = sum_{frame i} m_i + sum_{audited i} (F_i Y_i - m_i) / q_i

Do the same for F, and take ratio of totals. Region, season, source arm, taxonomy (with partial pooling), exposure/box flags when available and technical status can enter m_i; do not use the post-outcome "taxon gap priority" as a biological covariate. For a population not covered by the source frame, none of this removes volunteer plant/photographer selection bias. The direct audit ratio is the mandatory non-model comparator.

**Missing not at random:** An ROI-success propensity model fit only on successes cannot identify the white fraction among failures. The human audit observes actual colours in a probability sample of failures, making that portion testable. Where petals are still ambiguous or imagery inaccessible, produce explicit worst-case/tipping-point bounds. No universal missingness-at-random assumption.

## 4. Macroecological analysis: separate composition from colour response

Construct a frozen species × equal-area cell (initially nominal ~100 km resolution, to be optimized only by prespecified minimum support/coverage) × calendar-year/local flowering period panel. Use all observations via the observation-level calibration (do not throw away source opportunities solely to equalize photograph counts); aggregate *estimates* with one weight per supported taxon–cell–period. Report cell and taxon/observer coverage and external native-status confidence, with no extrapolation to cells without support.

Two distinct contrasts, each with the same source/provenance and calibrated outcomes:

- **Total sampled-assemblage contrast:** Compare macro-cell flower-colour composition along moisture/temperature/elevation/island gradients. Species composition is allowed to change; that is the biogeographic object of interest. This association is descriptive, not a climate causal effect.
- **Within-taxonomic/geographic support contrast:** Compare matched species or genus across settings with actual overlap, keeping the matched fraction and taxon/cell denominators visible. Species identity, within-region spatial structure, observer effects and source arm can be adjusted where estimable. Unlike the first contrast, this answers a *conditional* question and may have little support (do not treat an over-controlled zero coefficient as evidence that the ecological pattern is fictitious).

For seasonal hypotheses, use original observation year/month and local phenology, not a 30-year climate monthly mean masquerading as event weather. Separate species turnover, species photo-frequency changes and conditional photo colour; a photographed composition change is not within-individual pigment plasticity.

Precommit a primary climate contrast (moisture ↔ white) and treat elevation, island and local-season tests as distinct secondary contrasts. Spatially blocked geographic holdout (and genus-held-out checks for generalization), source/observer/observation-aware bootstrap or hierarchical uncertainty. A spatial cell fixed effect completely absorbs environmental variables constant within that cell: choose blocks appropriate to each estimand, rather than adding every geography predictor as "control."

## 5. Four-analysis table and credibility gates

| Analysis | Denominator | Main job |
| --- | --- | --- |
| A. Existing strict9-only | 451,643 auto-success of frozen 1.18m, with arm-specific breakdown | Source-dependent descriptive baseline **only** |
| B. Random manual audit | Prespecified independent observation-level 1,200, weighted by q_h, source frame specified | Flower presence, classification error, actual failure colour and adjusted observation-level rate |
| C. Calibrated corrected source contrast | B plus all observation-level auxiliary data; verified pi only where known | Higher-precision E1 and calibrated E2 for supported taxon–cells |
| D. Standardized and blocked macro contrast | E2 supported cell/taxon panel; overlap and native-status disclosed | Distinguish broad biogeographic composition from taxon-/region-conditional association |

Publish, **before ecological interpretation**, A→B→C→D shifts in white rate and climate slope; reference vs targeted arm composition, classifiability by region/season and the range of unclassifiable-true-flower fractions. If bias correction changes sign, interval spans null or meaningful no-overlap strata dominate, **downgrade to photo-conditional association**; never silently pick the positive model. A remaining contrast is not selection/fitness/anthocyanin mechanism evidence.

## 6. Execution boundaries and minimal next tasks

1. Freeze an **observation-level probabilistic audit-frame receipt**, with SHA-pinned corrected v2 photo source, one-photo-per-observation deterministic algorithm, stratum N_h/n_h, q_h, source-arm boundaries and overlap with old 512/128 (do not generate new colours).
2. Generate the 1,200-photo audit IDs from the **already existing** original-photo roster, reproducibly and blinded to colour; reuse the existing photo review interface with separate annotation files. No new 93m scan, new photo IDs, large-scale ROI reruns or full image re-download pipeline.
3. Annotate and independently adjudicate flower corolla truth in the random audit; use the existing 128 only for technical diagnosis. Until labels are actually collected, do not claim calibration achieved.
4. Compute direct design-weighted audit ratio and then cross-fit calibration; publish coverage and bootstrap/bound receipts, followed by a single moisture macro-model refit if it is meaningfully identified.
5. Keep the Oct 27 GEB existing-data manuscript on its frozen timeline. This is a prospective **bias-identification design**; do not replace its measured counts or make the existing published 128 "manual truth" by default.

**Stop criterion:** Once the primary association's corrected estimate, its source/region support and technical/botanical missingness bounds are clear, do not launch an image-count expansion just to narrow a misleading complete-case confidence interval. Inability to identify a native-world flora estimand is reported, not "fixed" by more iNaturalist photos.
