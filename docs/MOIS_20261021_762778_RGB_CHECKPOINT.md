# MOIS 2026-10-21: 762,778 verified original-photo outcomes / 353,039 four-hue classifications

## Completed operational milestone

[GitHub Actions run 39544272359](https://github.com/zuizui0223/mois/actions/runs/39544272359) passed **tests, source-artifact freshness, and original-photo recovery**. The new immutable [six-source RGB release](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39544272359-v1) is published and its ZIP was downloaded and independently checked. It contains **762,778 distinct original photo terminal attempts**, including original-image download failures, of which **353,039** carry one of the same four frozen FCP full-resolution ROI labels.

| Signed original-image source arm | Genuine attempted photo IDs | Classifiable frozen hues | Pending original photo IDs |
| --- | ---: | ---: | ---: |
| Prior immutable FCP/MOIS/first93m source union | 612,483 | 295,801 | n/a |
| p0004, fully completed | 18,119 | 7,085 | **0** |
| p0008, fully completed | 36,375 | 13,874 | **0** |
| p0016 (41/73 original image work groups) | 41,984 | 16,163 | **31,888** |
| p0032 (35/145 work groups) | 35,840 | 13,810 | **112,231** |
| Separate tropics/southern-latitude geo05 enrichment (17/22 groups) | 17,408 | 6,107 | **4,870** |
| Independent 320px shadow evaluation, all complete | 569 | 199 | **0** |
| **Six-arm combined photo-ID deduplicated checkpoint** | **762,778** | **353,039** | **148,989** |

Relative to [previous immutable checkpoint 39534159275](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39534159275-v1) (**753,562 photo attempts, 349,606 classifiable hues**), actual completed photo attempts increased by **9,216** and four-hue labels by **3,433**, reconciled as:

| New original photo completion | New attempted original IDs | Newly four-hue classifiable |
| --- | ---: | ---: |
| p0016 new sealed group artifacts 38 → 41 | +3,072 | **+1,160** |
| p0032 new sealed groups 31 → 35 | +4,096 | **+1,552** |
| geo05 regional source new sealed groups 15 → 17 | +2,048 | **+721** |
| **Total new actual terminal attempts** | **+9,216** | **+3,433** |

Across the new source arms beyond the historical 612,483 terminal-photo baseline, the completed **150,295** additional original-photo attempts comprise **299 acquisition failures**, **149,996 acquired photos**, **57,238 classifiable frozen four-hue ROI labels**, and **92,758 acquired-but-unclassifiable** image outcomes. Do not silently treat technical unknowns as nonwhite or nonflowering, nor try to repair by inventing labels.

## Source and exact byte validation

- GitHub Actions data-only checkpoint **39544272359**, all three jobs success, no original photo download, no new frozen ROI inference.
- Actual checkpoint artifact: **12320906836**, ZIP length **53,486,206 bytes**, ZIP SHA-256 `471fff587060498b1dbf0f5eff3a2d94bd0073f1d7cf9f69af35d5075fbe0516`.
- Independently opened the exact downloaded GitHub ZIP in the active container, recalculated its SHA-256, and independently verified **every** output SHA-256 in `checkpoint_receipt.json` against its three output members. All matched.
- `completed_photo_hue_counts.csv` and its source-wave counters reconcile to **762,778 attempted original photo IDs / 353,039 classifiable outcomes**. Four-colour counts by wave were checked as part of the audit.
- The production code also validated each frozen photo-source shard, original photo identity, saved terminal ROI method and SHA, full 612,483-original historical baseline identity, and absence of duplicate source photo IDs. Local verification checked source/receipt hashes and CSV/JSON arithmetic; it did not independently reread every row of both Parquet outputs.
- Output assets: `checkpoint_receipt.json`, `completed_photo_hue_counts.csv`, `combined_actual_photo_checkpoint.parquet`, and `new_completed_original_photo_ROIs.parquet`. The photo table retains `source_arm`, original observation identity and any original `selection_probability` fields; pending photos are not fabricated as Parquet outcomes.

## Descriptive photo outcomes, not world-flora estimates

Of the **353,039** machine-classifiable photo ROIs, automatic hue labels are: **white 168,346**, **yellow-orange 111,586**, **red-pink 36,935**, **blue-purple 36,172** (the source dataset itself must control any updates to these counts). These labels describe the selected photographs and the frozen model's admitted ROIs only. They are **not** the proportion of native flowering plant species with white corollas, nor independent open-flower botanical verification.

The regionally enriched geo05 arm has **17,408 original-photo attempts**, **6,107 four-hue ROI outcomes** and original `selection_probability` values that vary by observation and latitude. It must not be pooled unweighted with the distinct, nested observation-hash probability source arms. Original photo observation IDs can have multiple photographs; hierarchical, flora/taxon, geographic, season and missingness adjustments are necessary before ecological inference.

The upstream **93,209,351 source photo IDs** are not 93 million actually measured pixel/hue outcomes; p0004 and p0008 have completed their separate observational photo cohorts, while p0016, p0032 and geo05 retain **148,989 selected original photo IDs** awaiting actual worker completion.

## Mainline continuation

Already running frozen original-image workers: **p0016 two**, **p0032 two**, **geo05 one**, without confirmed failures. The three-hour GitHub Actions source-aware watcher and source-SHA freshness gate can recover newly sealed workers into immutable *data-only* releases. Preserve original photo/observation IDs, regional enrollment propensities, true technical unknowns and the unchanged full-resolution ROI-v4 model. Do not redownload completed images, do not enable a loss-inducing low-resolution filter, and do not hard-delete introduced or unknown native-status targets.
