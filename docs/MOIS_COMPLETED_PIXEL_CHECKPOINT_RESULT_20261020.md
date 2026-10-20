# MOIS completed-pixel checkpoint: 20 October 2026

## Published, verified current data

Recovery run [39373763787](https://github.com/zuizui0223/mois/actions/runs/39373763787) completed its source checks, output upload and Release publication successfully. The immutable partial-data Release is [mois-real-pixel-checkpoint-39373763787-v1](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39373763787-v1).

| Actual original-photo source | Terminal photo attempts retained | Four-hue automatically classified |
|---|---:|---:|
| Previously published FCP/MOIS/p0002 baseline | 612,483 | 295,801 |
| Completed p0004 units recovered in this snapshot | 12,288 | 4,828 |
| Completed tier03 shadow cohort | 569 | 199 |
| **Combined unique photo IDs** | **625,340** | **300,828** |

New completed outcomes added: **12,857**, including **12,816 acquired images** and **41 acquisition failures**. Of the acquired images, **5,027** have automatic four-hue classifications and **7,789** remain unclassifiable. No acquisition failure or unknown ROI was imputed chromatic or non-flowering.

The p0004 increment includes 12,249 successful acquisitions, 39 failed acquisitions and 7,421 acquired but unclassifiable photos. The shadow increment includes 567 acquisitions, 2 failures and 368 acquired but unclassifiable photos. The shadow cohort is fully complete; the p0004 and p0008 source waves are not fully complete in this snapshot.

## Pending work is not counted as measured

At this checkpoint's artifact capture, **5,831 of 18,119 p0004 targets** remained outside recovered completed units. The **36,375 p0008 targets** were not yet present in completed worker artifacts, so zero p0008 rows are included here. This does NOT mean no p0008 processing has occurred: worker 118459467998 was separately confirmed in the actual image acquisition/ROI step after checkpoint publication. Acquisition run [39371199127](https://github.com/zuizui0223/mois/actions/runs/39371199127) remains separate and untouched.

The completed-unit subset is selected by computation completion and is **not a representative complete probability cohort**. Do not use this partial checkpoint alone as a global white-flower frequency estimate. Future captures use original photo IDs to avoid double-counting.

## Files and source checks

- `combined_actual_photo_checkpoint.parquet`: compact source-metadata/hue projection of all 625,340 unique photo IDs, with source arm and checkpoint wave.
- `new_completed_original_photo_ROIs.parquet`: full new measured rows, including preserved original palette fractions and source metadata.
- `completed_photo_hue_counts.csv`: actual outcome counts by source wave and original morph label.
- `checkpoint_receipt.json`: exact source manifests, worker artifact snapshots, completed/missing units, original terminal hashes and output file SHA256 values.

The original rich baseline Parquet and full-wave releases are not replaced. No original image download or model inference is performed during checkpoint recovery. Each included new unit was checked against its source roster, frozen model signature, image/terminal SHA, original hue fields and counters, then joined with the byte-verified baseline after requiring zero duplicate/overlapping photo IDs.

The first recovery attempt stopped before publication because its new-ROI unknown-label restriction was mistakenly applied to a legacy baseline with **one empty original morph label**. The fix preserves legacy labels and unknown flags rather than relabelling that row; new ROI states remain strictly checked. The regression failed before the fix and all 13 focused tests passed afterward. Production recovery then succeeded.

Published artifact ID: **12242935043**, `mois-actual-sealed-pixel-checkpoint-39373763787`, 21,994,982 bytes. ZIP SHA256: `1deac4924e5f465c719eaa527641388ba68b6d1a52a9529651d30a0bb02a1c5e`. The downloaded ZIP and every output checksum in its receipt were independently checked after publication.

## Efficiency decision

The independent 569-photo shadow comparison is complete: see [actual low-resolution result](MOIS_569_LOWRES_SCREEN_REAL_RESULT_20261020.md). Do not enable the current 320px detector as an exclusion filter: at threshold 0.05 it would discard 7/199 baseline-classifiable photos, including 5/93 white-labelled ROIs, for 5.20% projected ROI-call time savings. Those are pilot projections against automated baseline labels, not realized end-to-end speedups or independently validated botanical recall. Continue original full-resolution measurement without native-status hard filtering or low-resolution image exclusion.
