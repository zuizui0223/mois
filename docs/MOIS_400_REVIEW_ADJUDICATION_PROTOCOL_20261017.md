# MOIS: independent botanical two-reviewer reconciliation and third adjudication

**Scientific programme:** seasonal, global photographic floral colour ecology in `zuizui0223/mois`, **not** the `fcp` within-species `D`, `W`, or 50-km clustering paper.

**Status (17 October 2026):** All 400 original photo IDs are reserved in five species and independent-observer sampling is verified. Reviewer A and B source sheets each contain 400 entries, but **all `reviewer_id`, `open_flower_stage` and `roi_organ` outcomes are unfilled**. No open flower was biologically confirmed in this 400-photo wave. The source photo opportunities must not be relabelled as photos of blooming flowers.

## Frozen sources and the review firewall

The original photo sample is SHA-pinned by [MOIS source selection run 38929402279](https://github.com/zuizui0223/mois/actions/runs/38929402279). A and B blank review forms are independently randomized and SHA-pinned by [run 38930414968](https://github.com/zuizui0223/mois/actions/runs/38930414968). A coordinator-only 400-token→original photo/date/site mapping is never committed to Git or distributed to independent reviewers.

| File | Exact SHA256 |
|---|---|
| Reviewer A blank original | `576a55e7e613faf8c93652b4ac2d72b79f8a38bec4bfe8267c8d67bbc3fed322` |
| Reviewer B blank original | `bfbbfc92106585ca260fbd5a33573ad2598209d19907a52bbce2b941b1164203` |
| Coordinator-only original photo mapping | `6119a9606b4d90eb80a8a90e7fd2576631d14063de320e788e7c274fd050a3c4` |

The [new validator](../scripts/analysis/mois_validate_400_blinded_botanical_reviews.py) compares **complete original 400-token membership and immutable token/species/photograph URL/licence/attribution** across the original and returned A/B sheets. Changing a source photo or silently replacing a failed photograph is forbidden. Reviewer IDs must differ and a packet may have only one reviewer identity. An assigned reviewer with blank stage and organ remains **unreviewed**; partial or unsupported codes are rejected.

The approved stage codes are `open_flower`, `buds_only`, `fruit_only`, `no_flower_visible`, and `uncertain`. Organ codes are `petal_tepal`, `other_flower_organ`, `spikelet_glume_bract`, `fruit_seed`, `leaf_stem`, `background`, and `uncertain`.

## Workflow once the two forms are independently completed

1. Keep separate completed A and B files; do not let the second reviewer see the first reviewer's answer. Collect no RGB colour, environmental exposures or seasons in these review sheets.
2. Run the source-preserving offline audit using the original pinned blank forms and the **private** coordinator mapping.
3. For two complete agreeing votes, classify as verified open-petaloid (future ROI eligible), open other organ (not primary), verified non-open, or unresolved (uncertain). For partial or missing votes, retain pending, **not** negative flowering.
4. For two complete disagreeing votes, generate a separate **season/colour/first-review-decision-blind third-adjudication queue** with just the source photo, attribution and a token. A third reviewer must be different from both A/B reviewers. With a filled adjudication queue, re-run the audit with `--adjudication`; do not silently resolve contradictions by majority or by automated ROI colour.
5. Store `COORDINATOR_ONLY_400_original_photo_review_decisions.csv` and any filled review sheets outside the public repository, and publish only a redacted aggregate receipt as needed.

Example local execution, where `blank_A.csv` and `blank_B.csv` are the **original** files and `answered_A.csv` / `answered_B.csv` are independently filled copies:

```bash
python scripts/analysis/mois_validate_400_blinded_botanical_reviews.py \
  --blank-a blank_A.csv \
  --blank-b blank_B.csv \
  --coordinator-map coordinator_only_original_photo_mapping.csv \
  --filled-a answered_A.csv \
  --filled-b answered_B.csv \
  --outdir private_biological_review_qc
```

A second invocation can specify `--adjudication independent_third_answered.csv` to resolve conflicts, preserving all original A/B decisions and the source IDs.

**Crucial distinction:** `eligible_for_future_original_photo_ROI_measurement` means independent botanists agree that **an open flower and a petal/tepal are visible**. It does **not** say the ROI has been segmented, its hue classified, white counted, or that selection, plasticity, pigmentation costs or weather causation were demonstrated.

### Scientific analysis gate after human decisions

Only after botanical stages are confirmed should the unchanged photographic ROI method process source photos while preserving the original 400-photo selected denominator and all missingness. The ecological model must condition on the actual observed year and calendar month, site, species, image classifiability/quality, and observer/source effort. With only five taxa in largely temperate North America, global extrapolation is inappropriate. Even if an image gets two independent botanical approvals, the original iNaturalist source and automated colour measurement are still shared observational systems.

For automated regression tests, the PR-only CI reconstitutes the original **unanswered** 400 photos from the immutable source workflow, runs the validator, requires all 400 to stay pending, and uploads **only the aggregate redacted JSON**. There is deliberately no pathway in CI that uploads filled human reviews or the coordinator mapping to public workflow artifacts.
