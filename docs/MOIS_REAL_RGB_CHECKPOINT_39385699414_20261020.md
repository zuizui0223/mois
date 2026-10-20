# MOIS: completed original-photo RGB snapshot 39385699414 — 2026-10-20

## Current location and next step

**642,435 unique original-photo terminal attempts, including 76 newly acquired-source image download failures; 307,382 photographs assigned one of the frozen four automated colour labels.** Validated as the latest immutable [MOIS intermediate data release](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39385699414-v1) created by completed [GitHub run 39385699414](https://github.com/zuizui0223/mois/actions/runs/39385699414).

| Source image wave | Actual terminal photo attempts | Frozen 4-colour classifications | Yet to complete |
| --- | ---: | ---: | ---: |
| Original previous MOIS/FCP + first 93m source draw | 612,483 | 295,801 | n/a |
| Complete p0004 original photo wave | 18,119 | 7,085 | **0** |
| p0008 disjoint new tier | 8,192 | 3,130 | **28,183** |
| p0016 disjoint new tier | 3,072 | 1,167 | **70,800** |
| Completed 320px shadow validation arm | 569 | 199 | **0** |
| **All unique original-photo terminal IDs** | **642,435** | **307,382** | **98,983** |

Relative to the preceding [checkpoint 39381982167](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39381982167-v1) with 639,363 attempted original photos and 306,186 classified, this adds **3,072 actual terminal outcomes and 1,196 classified hues**. The two increments are p0008 +2,048 attempts/+807 classified and p0016 +1,024/+389. p0004 was already complete in the previous checkpoint.

Among the **29,952** new terminal attempts beyond the original 612,483 baseline, **29,876** images were actually acquired, **76** acquisitions failed, **11,581** acquired photographs had four-hue classifications and **18,295** acquired photographs remained unclassifiable. Do not call the latter flowers of another colour, nonflowers, or native/non-native without independent evidence.

## Verification and immutable release

* Source artifact **12247083294**: \`mois-actual-sealed-pixel-checkpoint-39385699414\`, ZIP size **26,018,393 bytes**, ZIP SHA256 **\`2becaf55630b337f34674c5b512fd1844009dd6328a151b59f8e0efb5bf32978\`**.
* Independently downloaded and checked the ZIP against its Actions artifact SHA256. All four original published release assets matched GitHub's advertised SHA256, including the receipt and both Parquet outputs.
* Three output-file hashes inside \`checkpoint_receipt.json\` matched the ZIP members. Independent summation of \`completed_photo_hue_counts.csv\` reconciled to **642,435** terminal rows and **307,382** four-hue labels, broken out by all five arms above.
* The full byte/source/photo-ID/ROI terminal/model-method validation and Parquet roundtrip took place in the successful production recovery job. The local runtime did not independently re-read every Parquet row because it lacks \`pyarrow\` and \`fastparquet\`.
* Recovery is **data-only**: no original image redownload, no flower-petal inference or hue relabeling. Existing p0008 and p0016 worker groups continue in parallel.
* This is **a partial computational completion**. Unprocessed source photos remain pending; their colours are unknown and must not be called white/nonwhite. It is **not a complete global probability cohort** or a global wild-flora flower-colour estimate.

## Coverage and bias implication

The separately completed genuine p0004 image cohort was fully checked in [p0004 original-photo coverage audit](MOIS_P0004_FULL_WAVE_COVERAGE_QA_20261020.md). It contains 18,119 originals from 10,416 observations, with classifications across 2,867 nominal source taxon IDs; 326 IDs have at least five colour-labelled photographs. Only 1,087 p0004 attempted photos (6.0%) were in tropical latitudes, and automated classifiability is 34.2% there versus 40.1% in northern extratropical latitudes. The 4-colour proportions of this particular source arm cannot be interpreted as unbiased world-flora frequencies.

**Next acquisition priority**: finish the existing zero-overlap p0008 and p0016 probability tiers without losing technical-unknown statuses. Preserve the random-observation arm and consider a **separately labelled** coverage-enrichment source for neglected tropical/southern locations and under-observed taxa, with its own explicit selection design. The previously tested 320px detector must not exclude photos; it lost 7/199 baseline-classifiable photo colours including 5 white-labelled ROIs. Keep native/introduced range status as an optional covariate, never a photo-deleting prefilter.

The upstream **93,209,351** originals comprise an immutable photo-ID/observation metadata frame, not 93 million measured flower colours.
