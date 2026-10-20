# MOIS: one-pass full 93,209,351-photo target register before costly RGB

## Decision

**First register all original source photo IDs once; do not hard-delete photographs based on the current native/introduced metadata.** The verified 2026-09-27 source already has **93,209,351 uniquely indexed original photo IDs from 42,111 original source taxa** and corresponding observations in 128 source partitions. These IDs are metadata opportunities, **not** 93m photographs with measured floral colour. Re-reading 128 source pairs whenever the nested observation-hash probability changes from 0.0002 to 0.0004, 0.0008, etc., is wasteful.

The deterministic one-pass candidate registry (`scripts/acquisition/mois_onepass_full93m_target_tiers.py`) creates **14 mutually exclusive observation-hash score intervals**:

`(0, 0.0002), [0.0002, 0.0004), [0.0004, 0.0008), ..., [0.4096, 0.8192), [0.8192, 1)`.

All original photo IDs from an original observation remain in its same tier. The same seed `mois-global-observation-v1` and unchanged original hash method are used, so earlier nested selection remains a subset of later cumulative thresholds. A photo is **never discarded by taxon, genus, region or existing WCVP status**. Each archived row records its original photo ID, observation UUID, taxon ID, photographer/observer IDs, original observation date and site, license, original image format, date/coordinate availability and frozen model's supported-format flag. Previously measured photo IDs must be excluded only when **constructing a new pixel-processing queue**, not deleted from the complete source-opportunity register.

The source tiers do **not** assert an image depicts an open flower, classify hue, acquire image bytes, certify individual provenance, or give a photograph a native species status. Releasing 93m *target records* is not measuring 93m *flower colours*.

## Why not a native-only pre-filter?

The existing WCVP × coordinate overlay on **350,399 previously measured image attempts** describes 268,670 photographs with species recorded as native in their broad TDWG-3 botanical-country unit and 30,960 recorded as introduced. Many records are unmatched or geographically ambiguous. These proportions are descriptive of that historically selected image cohort; they **must not** be extrapolated as an authenticated native fraction for all 93m candidates.

- Broad TDWG-3 native status means only **species listed as native in a botanical country**. It does not identify a *wild individual*, a garden photograph, or nativeness on a **specific offshore island**. The last distinction matters for the proposed native-island white-colour hypothesis.
- Excluding introduced/uncertain records before measuring colour would destroy information about photographic sampling and potential native–introduced controls, disproportionately remove incompletely resolved regions, and could fabricate a geographic white-frequency pattern.
- Native/invasive status should be an **append-only, versioned, confidence-labelled annotation** on the full photo/observation manifest, following independent Kew name resolution and geospatial overlay. Keep `native_at_TDWG3`, `introduced_at_TDWG3`, `unknown`, `conflicting`, `unmatched`, and island-specific nativeness separate. Only the primary native-only *analysis stratum* uses any native restriction; the raw photo targets stay intact.

## Efficient work after registering the target pool

1. Build and SHA-verify the one-pass target register once from the existing full source Release. Source index is not re-downloaded separately for each probability threshold. Store each disjoint tier independently so later acquisitions read only newly admitted target IDs.
2. Use the original photo ID as the permanent deduplication key and maintain the original observation ID to track within-observation photo dependence and independent evaluation holdouts. Split target processing into resumable bounded **compute batches**, not species or geography inclusion caps. Freeze row-selection probabilities and target receipt before image colour inference.
3. Distinguish a missing/invalid date, missing coordinates, missing flower annotation, unsupported original image extension, image download failure, visually non-flowering photo, failure of floral ROI, and uncertain hue. Only unsupported image formats have an upfront mechanical compatibility flag; the absence of an iNaturalist **observation-level** flowering annotation is not evidence that a given photo has no open flower.
4. On genuinely downloaded source photos, first use a **cheaper low-resolution image/flower-stage candidate screen**, calibrated with blinded independent positive and negative photographic audit (including screen-negatives) before using it to avoid heavier YOLO11/EfficientSAM ROI-v4 computation. Preserve screening probability/confusion accounting; do not call screened-out photos chromatic or use them as observed nonflower absence. The unchanged frozen ROI-v4 yields the existing comparable four-hue labels only on processed photos.
5. Report progress using (a) metadata photo IDs enumerated, (b) truly attempted original image downloads, (c) botanical open-flower review, (d) 4-hue classifiable ROIs, (e) technical failures/unknowns, and (f) distinct taxa/observations/geographic cells. Do not treat metadata target counts as coloured flower pixels.

## Coverage and inference safeguards

The source photo catalogue itself was species-selected (42,111 original source taxa), geographically opportunistic and observer-biased, not a complete census of all wild plants worldwide. Compute-batch cutoffs do not equal biological selection probabilities. Preserve a probability-based observation cohort for the main descriptive estimate, and use separate, explicitly labelled coverage-enrichment batches for undersampled realms, seasons, altitudes, native island floras or taxa. Do not mix them into one unweighted world-flora prevalence estimate. Keep observation and photo-disjoint genuine future holdout before claims of confirmatory testing.

**Status of this new pipeline:** metadata-target registry code and source identity tests. A Release can only be reported complete after all 128 source partitions pass their original SHA checks, all 14 tier counts sum to **93,209,351**, and the immutable manifest is published. This work does **not** claim new image colour measurements or an inferred native/introduced fraction for the whole 93m frame.

## Source-tier-to-ROI execution bridge

`scripts/acquisition/mois_queue_from_full93m_target_tiers.py` reads **only one** completed disjoint tier from the signed registry, verifies all 128 original photo-part checksums and the original inclusion hash, removes photo IDs from **already completed** true-image measurements, preserves unsupported formats in a no-colour sidecar, and writes the existing frozen FCP ROI-v4 worker-shard schema. Its optional exhaustive observation compute buckets (default 64) divide **processing jobs**, not botanical target inclusion; observations' photos remain in the same compute bucket. Every bucket must eventually be attempted before interpreting an unweighted complete-tier colour sample. The `mois_full_index_real_original_photo_measure` executor remains unchanged; running a pre-colour queue alone has **not** downloaded any images or classified petals.