# MOIS: verified completed-pixel checkpoint request (2026-10-20)

This source-backed milestone intentionally triggers the existing **data-only** checkpoint workflow on main. It reads **previously completed and SHA-verified original-photo work units**; it neither downloads raw iNaturalist images nor reruns flower-colour inference. The ongoing p0004, p0008 and p0016 acquisition jobs are left untouched.

## Source snapshot that justifies a new immutable checkpoint

At the source check around 06:20 UTC / 15:20 JST on 2026-10-20:

| Fixed original-colour acquisition run | Completed sealed RGB worker groups | Total planned photos |
| --- | ---: | ---: |
| p0004 `39353343403` | 16 of 18 | 18,119 |
| p0008 `39371199127` | 5 completed | 36,375 |
| p0016 `39374811084` | 1 completed | 73,872 |
| shadow `39369517433` | 1 completed, full cohort | 569 |

The preceding immutable checkpoint is [39376435276](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39376435276-v1), with **628,412 unique original-photo attempts** and **302,026 machine-classified four-hue outcomes**. More worker artifacts have arrived since that checkpoint, but their new photo counts/classifiable labels MUST be recomputed by the validated production collector, not projected from the scheduled source rows.

The checkpoint source registry pins these exact four run IDs, their original photo roster SHA256 and the immutable FCP/MOIS prior baseline. New and previous photo IDs must remain disjoint, missing units must stay *pending*, acquired-but-unclassifiable images cannot become coloured flowers, and legacy unknown labels are preserved exactly.

## Operational result criterion

A successful [checkpoint workflow](https://github.com/zuizui0223/mois/actions/workflows/mois-completed-pixel-checkpoint.yml) will publish `mois-real-pixel-checkpoint-<run-id>-v1` with `checkpoint_receipt.json`, both Parquet outputs and photo/hue counts. Only after the run is successful and the new receipt matches all file digests should its totals be described as confirmed.

This milestone path is an explicit source-backed way to request an additional checkpoint during active long-running image acquisition without modifying or restarting the measurement workers. It does not schedule repeated artificial snapshots.

## Refresh request: complete p0004 source + new p0008/p0016 group terminals

This new source-verified snapshot was requested after the full p0004 terminal union was separately published as [mois-genuine630602-original-pixels-first-two-93m-waves-v1](https://github.com/zuizui0223/mois/releases/tag/mois-genuine630602-original-pixels-first-two-93m-waves-v1). The **previous partial pixel checkpoint** was [39381982167](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39381982167-v1), with 639,363 attempted photo IDs and 306,186 machine-classified four-hue images. The source-authenticated p0004 complete increment was independently verified at 18,119 attempts / 7,085 classifiable.

New archived image-worker groups are available beyond that checkpoint:
- p0004 39353343403: all 18/18 complete (unchanged since previous checkpoint)
- p0008 39371199127: 8 completed source-pinned groups; four acquisition workers active, none failed
- p0016 39374811084: 3 completed source-pinned groups; two acquisition workers active, none failed
- full shadow study 39369517433: 569 attempted and 199 classified (unchanged)

**Do not extrapolate a new classified count from group counts.** The existing collector must SHA-check and reconcile each original source photo / observation ID, byte-immutable model method and terminal outcome, then publish the exact newly classified and still-pending counts as an immutable checkpoint. Do not rerun colour inference, re-fetch images or add any low-resolution filter.

This documentation update intentionally triggers the existing one-off data-only `main` checkpoint. It does not schedule a recurring task or change image-acquisition workers.

## Further verified acquisition milestone — 2026-10-20, around 16:42 JST

After the previous immutable [checkpoint 39385699414](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39385699414-v1) published 642,435 terminal photo IDs with 307,382 frozen automatic four-hue classifications, the following additional worker files reached terminal-success status:

* p0004 39353343403: complete, all 18 of 18 archived; do NOT re-download or re-run
* p0008 39371199127: 10 of 28 completed original-photo worker groups; four currently processing, zero failures
* p0016 39374811084: 4 of 28 completed worker groups; two currently processing, zero failures
* tier03 shadow 39369517433: 569 original-photo attempts, already complete and permanently source-deduplicated

Use the pinned production collector on only the stored source SHA and actual image terminal receipts, without running image acquisition or colour inference. Publish an immutable intermediate checkpoint after verifying all newly completed units. The total new photo IDs and usable hue counts must come from exact original-photo terminal reconciliation, never assume one full group's pixels based on queued count.

This source-documented milestone is a one-off re-capture request through the existing push-triggered, data-only checkpoint workflow. It does not add an automatic recurring task.


## Source-verified checkpoint request — 18:21 JST, 2026-10-20

A fresh paginated GitHub Actions jobs read shows new completed **original pixel measurement** worker artifacts since checkpoint [39388432954](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39388432954-v1) (645,507 attempted photo IDs / 308,544 classifiable four-hue ROIs):

- p0004 original source run 39353343403: 18/18 full measurement workers successfully sealed, already in baseline checkpoint.
- p0008 new source run 39371199127: **16/36** sealed successful groups, 4 active, 16 waiting, **0 failures**; each completed group retains all original per-photo ROI outcomes.
- p0016 new source run 39374811084: **6/73** sealed successful groups, 2 active, 65 waiting, **0 failures**.
- Previous shadow pilot 39369517433: complete; its 569 original photo IDs are preserved and must not be counted twice.

The previous snapshot only contained p0008 **10** groups and p0016 **4** groups. This source-backed request triggers one more **data-only** recovery on the main branch. It must checksum-reconcile original photo/observation/source-shard IDs, immutable classifier metadata and terminal CSVs; never download the photos again, run new colour inference, or impute missing image colours.

**Outcome counts must be read from the new production recovery receipt, not inferred from these computational group counts.** Acquisition continues independently. The incremental snapshot is not a probability-complete sample or a global white-hue prevalence estimate.
