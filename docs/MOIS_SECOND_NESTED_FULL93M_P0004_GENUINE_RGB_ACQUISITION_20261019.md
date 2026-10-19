# MOIS acquisition-only: second nested observation-photo RGB wave (p = 0.0004)

**Status:** conditional GitHub Actions acquisition pipeline, not yet photographed/mapped/analysed results.

The first [all-93,209,351-photo original-ID source](https://github.com/zuizui0223/mois/releases/tag/mois-full-93209351-indexed-photo-observation-source-v1) is fully source-signed. The first pre-colour observation-hash probability p=0.0002 selected **18,533** source photo IDs (71 previously measured, 18,462 genuinely new original-image attempts, zero unsupported source formats). Those p0002 photo-image workers must **all finish successfully** and publish their original YOLO11/EfficientSAM ROI-v4 `measured_photos.parquet` and `receipt.json` before this second pipeline launches.

## Current action

The new [nested p0004 acquisition workflow](../.github/workflows/mois-full93m-second-nested-p0004-RGB.yml) is intentionally triggered **only after successful completion of** the previously deployed [first true RGB acquisition](https://github.com/zuizui0223/mois/actions/runs/39232152549), not on a PR or an unrelated main push. A missing or failed first release is a hard gate, never a guessed source sample.

It downloads the full 128 original source-photo/observation shards from the original completed Release with checksums and the same 42,111 original named-taxon counts. It also fetches the **five previously measured MOIS original pixel Parquets and the completed p0002 real pixel Parquet**, validating that the combined previous source has **368,861 unique photographed original IDs** before selecting anything new.

A nested hash threshold p=0.0004 is then applied to the **same original observation IDs and original fixed seed**. Source observations already selected at p=0.0002 remain within the larger p0004 threshold, but their already measured photo IDs and the five historical cohorts are excluded from remeasurement. Thus the second pipeline genuinely adds image-source opportunities from previously unselected observations, still preserving all original 93,209,351 upstream photo candidates without taxon, genus or total census cap.

Every newly selected legitimate original photo is attempted using the same frozen image-specific source model, original provider rate gate and bounded original-image loader, with permanent image failure, unsupported source format (as a separate source-only ledger), and ROI technical colour uncertainty retained. The new source probability is known within the historical iNaturalist photo catalogue, **not** a probability that a living flower was present, observable or photographed in nature.

The second wave's final real original RGB attempt count, classifiable four-hue count and added taxa are **unknown until the p0004 source-index rerun and all new image processing finish**. The workflow publishes only a new immutable `mois-full93m-observation-second-original-RGB-p0004-v1` original-pixel Release if all source identity, new photo terminal and method signatures pass. The p0004 incremental result is not another copy of the p0002 photos.

**Scientific restriction:** complete actual pixel acquisition and source technical QA first. Native/introduced WCVP, UV, pollinator, island/ecological hypothesis and manuscript tests remain deferred until more genuine photo measurements exist. No source-only IDs are ever counted as coloured flowers.
