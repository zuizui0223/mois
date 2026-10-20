# MOIS: new original photos plus a no-exclusion low-resolution screen comparison

## Purpose and decision

Reduce costly image processing only when the reduction is supported by image-level evidence, without deleting native-status-unknown taxa, inventing colours, or redownloading previously measured photographs. The complete 93,209,351-photo target registry is already published. The original ROI-v4 runtime already avoids the EfficientSAM encoder/decoder when its detector supplies no boxes. A new low-resolution screen therefore does not automatically produce a speedup.

This implementation is **shadow measurement only**. Every acquired photograph receives the original full ROI-v4 + fixed-palette measurement regardless of the 320px screen result. The screen runs after each baseline call, using a separate detector instance with the same SHA-pinned detector weights. It does not feed boxes, masks, thresholds, species, coordinates or colour labels back into the baseline. Screen errors are recorded and do not remove baseline results. Thresholds 0.05, 0.10 and 0.25 are fixed before the new image outcomes are viewed.

## New image source, not an old-photo replay

Use the published source register's disjoint tier03: observation-hash scores in [0.0008, 0.0016), 74,883 source photo candidates. Process **bucket 0 of 128** of a separate deterministic observation hash. This is a bounded comparison cohort, not 74,883 completed pixels or a global flowering census. The exact selected/new image count is emitted by the real source queue.

All previous 612,483 genuine FCP/MOIS/p0002 photo IDs are removed using the byte-verified union. Active p0004 and planned p0008 observations have scores below 0.0008 and cannot overlap this tier. Each new image is fetched once using the same bounded, paced original loader; the screen uses an in-memory resized copy of those bytes. Native/introduced/unknown metadata does not hard-filter targets.

Do not interrupt or change the existing p0004 or p0008 jobs. This study has separate artifact names, release tag and concurrency group. When the rest of tier03 is processed, its prior-photo set MUST include the photo IDs in `mois-tier03-actual-RGB-lowres-shadow-bucket0-v1` once that release is complete. Do not redownload this comparison cohort as though it were unprocessed.

## Outputs and exact interpretation

The ordinary `measured_photos.parquet`, original method receipt and source manifest are preserved. A separate `screen_shadow.csv` associates the original measurement ID and image SHA with screen confidence, screen error, baseline ROI time, screen time and the original four-hue outcome. Acquisition failures remain in the ordinary source denominator and are not mislabeled screen negatives.

`screen_shadow_summary.json` reports the number of baseline-classifiable photographs a hypothetical screen would lose, separately for white, yellow-orange, red-pink and blue-purple. A hue with no baseline support has null retention, not a fabricated perfect score. Unknown original ROIs are counted separately and are not biological non-flowering negatives.

Projected savings equal summed baseline ROI call times on hypothetical rejects minus screen call times on all acquired images. These are **not realized production speedups**. They exclude download time, model initialization and palette classification, and include the cold first inference in each work unit. This shadow pilot itself does extra work. Its timing cannot be extrapolated directly to 93 million images.

The reference is the automated baseline, **not independent botanical truth**. No automatic exclusion is enabled by this pipeline, even if the shadow report looks favourable. Independent flower/petal review and a separate untouched validation cohort remain necessary before claiming botanical recall or deploying a phenotype-sensitive filter at scale.

## Verification and entry points

- `scripts/acquisition/mois_lowres_screen_shadow.py`: calls the exact frozen source main through an observation wrapper and collects sidecars without modifying the baseline result frame.
- `tests/test_mois_lowres_screen_shadow.py`: negative-screen baseline preservation, screen errors, baseline exceptions, non-finite confidences/timings, hue-specific accounting, duplicate IDs, empty acquired cohort, and the original-source evaluator seam. Synthetic test values are not ecological results.
- `.github/workflows/mois-tier03-lowres-screen-shadow.yml`: source-SHA verification, actual new image acquisition, full baseline plus shadow, complete ID reconciliation and data-only release.

Status must be reported separately as code/tests passed, source queue frozen, images actively being processed, and released measured results. The final screen retention, exact new four-hue count and timing outcome remain unknown until actual workers complete and outputs reconcile.
