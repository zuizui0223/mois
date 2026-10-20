# MOIS actual original photo RGB checkpoint — 20 October 2026

## Result: new source-authenticated image outcomes added

Production [checkpoint run 39380302979](https://github.com/zuizui0223/mois/actions/runs/39380302979) completed, and its [immutable prerelease](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39380302979-v1) is public.

**636,604 distinct actual original-photo terminal IDs, including acquisition failures; 305,139 with a genuine four-hue automated classification**. Compared with the previous [checkpoint 39376435276](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39376435276-v1), this adds **8,192 photo attempts** and **3,113 newly classifiable photos**.

| Original-photo source | Completed attempted photo IDs | Four-hue classified | Pending originals in source roster |
| --- | ---: | ---: | ---: |
| Previous FCP/MOIS/p0002 immutable baseline | 612,483 | 295,801 | n/a |
| p0004 run 39353343403 | 16,384 | 6,425 | 1,735 |
| p0008 run 39371199127 | 5,120 | 1,936 | 31,255 |
| p0016 run 39374811084 | 2,048 | 778 | 71,824 |
| Completed tier03 320px shadow study 39369517433 | 569 | 199 | 0 |
| **Combined unique original-photo terminal IDs** | **636,604** | **305,139** | **104,814 remaining planned attempts** |

The current snapshot includes p0004 group 0–15 (16 groups), p0008 5 completed groups, p0016 2 completed groups and the entire 569-photo shadow comparison cohort. Remaining photo counts are computed from the exact source roster, **not** from a world-flora representative sample.

### Increment, original image failures and uncertain RGB states

Since the previous checkpoint:

| New source | Newly completed photos | Newly classified hues |
| --- | ---: | ---: |
| p0004 | 3,072 | 1,179 |
| p0008 | 3,072 | 1,156 |
| p0016 | 2,048 | 778 |
| **Total** | **8,192** | **3,113** |

Across **all 24,121 new actual original-image terminal outcomes** beyond the immutable 612,483-photo baseline: 24,056 images acquired, 65 technical acquisition failures, 9,338 four-hue classifications and 14,718 acquired-but-unclassifiable ROIs. These failures and unknowns are neither assumed chromatic nor interpreted as confirmed flower absence.

The new incremental 8,192-photo subset corresponds to 8,177 acquired originals and 15 acquisition failures; 3,113 image hues classifiable and 5,064 acquired-but-unclassifiable.

### Source authentication

- GitHub Actions checkpoint run: **39380302979**, passed both tests and recover jobs.
- GitHub Actions ZIP artifact: **12244677813**, ZIP SHA-256 **11f7af4a5aff9af914025097fff9aa6f80b06f92fb01dd40731c8699ecef40c2**, 24,665,801 bytes.
- Output files: `combined_actual_photo_checkpoint.parquet`, `new_completed_original_photo_ROIs.parquet`, `completed_photo_hue_counts.csv`, `checkpoint_receipt.json`.
- Independently extracted published ZIP and checked all three output file SHA256 hashes against `checkpoint_receipt.json`, recalculated the photo and classified counts from the CSV, and matched them to the JSON source-wave counters.
- This local container does not include `pyarrow` or `fastparquet`, so it did **not** independently reread all Parquet rows; source/photo ID and Parquet validation ran successfully in the original production recovery job.
- The original full 93,209,351-photo target registry remains metadata, **not 93 million measured flower hues**. Existing raw licensed photos are not redistributed.
- The data-only collector performed **zero new image downloads and zero colour reinference**. Active p0004, p0008 and p0016 image-acquisition jobs were not stopped or restarted.

### Methodological guardrails

The latest 320px shadow pilot would miss 7/199 baseline-classifiable images, including 5/93 labelled white, at only 5.20% projected ROI-call savings: **do not enable low-resolution image exclusion**. No hard native/introduced-only species filter. All original unknown outcomes and legacy labels remain immutable.

This release is a **partial computational-completion checkpoint**, not a completed probability stratum and not an unbiased global estimate of white-flower frequency. Use complete nested observation strata, taxonomic/geographic and seasonal weighting, and photo classification checks before ecological inference.
