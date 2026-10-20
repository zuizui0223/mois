# MOIS source-authenticated real image-colour checkpoint — 20 October 2026

## Current completion and next action

[GitHub Actions run 39388432954](https://github.com/zuizui0223/mois/actions/runs/39388432954) completed both test and recovery jobs successfully. Its immutable [real pixel release](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39388432954-v1) contains **645,507 unique terminal original-photo IDs** with **308,544 frozen four-hue automatic colour outcomes** (white, yellow-orange, red-pink, blue-purple). Terminal attempts include acquisition failures and unclassifiable originals, not just usable petals.

This recovers **3,072 additional original-photo terminals** and **1,162 additional four-hue classifications** since [checkpoint 39385699414](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39385699414-v1) (642,435 / 307,382). No original photograph was re-downloaded or re-inferred by the collector.

| Original source wave | Completed photo attempts | Four-hue classified | Acquisitions failed | Pending original-photo attempts |
| --- | ---: | ---: | ---: | ---: |
| Immutable FCP/MOIS/first93m baseline | 612,483 | 295,801 | Historical categories retained unchanged | – |
| Full p0004 (source 39353343403) | 18,119 | 7,085 | 54 | **0** |
| Partial p0008 (source 39371199127) | 10,240 | 3,890 | 14 | **26,135** |
| Partial p0016 (source 39374811084) | 4,096 | 1,569 | 9 | **69,776** |
| Completed tier03 320px shadow arm | 569 | 199 | 2 | **0** |
| **Combined distinct genuine photo terminal IDs** | **645,507** | **308,544** | — | **95,911** |

Among **33,024** new actual source-photo terminal outcomes (four waves beyond the 612,483-photo baseline), **32,945** images were acquired and **79** acquisitions failed; the acquired originals include 12,743 classifiable photos and 20,202 acquired-but-unclassifiable ROIs. The last 3,072 additional attempts include +2,048 p0008 photos (+760 classified) and +1,024 p0016 photos (+402 classified).

## Full GitHub job counts versus API pagination

The workflow-jobs connector by default returns only the **first 30 jobs** and must not determine the denominator in a multi-matrix workflow. A full GitHub REST `jobs?per_page=100` read returns **36 p0008 measurement groups** and **73 p0016 groups**. At the exact worker snapshot used here, p0008 was **10/36 complete, 4 active, 22 queued**, and p0016 was **4/73 complete, 2 active, 67 queued**; no failed measurement jobs were observed. The p0004 source is all **18/18 groups complete**. These group counts are computational progress, never inferred number of newly acquired flower pictures.

Future audits must paginate the GitHub Jobs API, or use immutable original source `group_ids` receipts, rather than treating the first results page's size as all scheduled jobs.

## Confirmed frozen RGB outcome counts

The newly completed p0008 10,240 original-photo terminal outcomes include white 1,912, yellow-orange 1,254, red-pink 369, blue-purple 355, and 6,350 mixed/unknown terminal labels (including 14 failed acquisitions). p0016 4,096 terminal originals include white 773, yellow-orange 502, red-pink 163, blue-purple 131 and 2,527 mixed/unknown labels (including 9 acquisition failures).

The p0004 complete cohort contains white 3,504, yellow-orange 2,316, red-pink 637 and blue-purple 628 automatically classifiable image ROIs. Across all 645,507 attempted original images, **147,058 were automatically labelled white**, **308,544 were classifiable into one of four hues**, while all technical missingness remained as observed.

The classifiable-photograph white share is descriptive **conditional on the model's automated ROI admission and the source photo-selection process**. Neither a roughly 50/50 white fraction in new cohorts nor the original 93,209,351 indexed photo-ID catalogue establishes an evolutionary equilibrium, pollen selection, species-level white frequency, confirmed anthesis, or native-only flowering composition.

## Data authenticity checks

The immutable GitHub Actions artifact [12248736241](https://github.com/zuizui0223/mois/actions/runs/39388432954/artifacts/12248736241) is **26,732,924 bytes**, ZIP SHA-256 **`bd14e431c511f57fcb7b6fe310024346626007b789b5051ab8d9843247fcd7b4`**.

Its contents are `checkpoint_receipt.json`, `combined_actual_photo_checkpoint.parquet`, `new_completed_original_photo_ROIs.parquet`, and `completed_photo_hue_counts.csv`. Independently reading the ZIP confirms that all three data-file SHA-256 hashes match the receipt's `output_sha256` object, and summing the CSV by wave and hue reproduces **645,507 / 308,544** and all source-wave counts.

The production collector independently verifies every completed photo-shard original source SHA, terminal CSV SHA and classifier/method receipt, checks for duplicate original photo IDs and preserves signed historical records. The local container has no `pyarrow` or `fastparquet`, so the local check does **not** independently iterate the Parquet rows; the successful original production job performs the Parquet write and source validation. Exact source and field validation at every image-pixel row has not been independently repeated locally.

## Next work on the main acquisition path

Continue the **95,911 uncompleted selected original photo attempts** in the existing p0008/p0016 pipelines. Do not start another copy of the same queue, do not reduce taxonomic coverage by excluding introduced/unknown native status, and do not enable the 320px photo deletion screen (which missed some white-labelled original ROIs). When new worker groups seal, source-verify them and append via another **data-only** checkpoint.

These checkpoints are partial computational-completion subsets, **not a probability-complete independent global colour sample**. The original 93,209,351 candidate photo records are metadata rather than 93m images with measured flower colours.