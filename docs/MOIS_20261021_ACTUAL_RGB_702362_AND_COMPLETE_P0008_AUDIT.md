# MOIS — 702,362 authentic original-photo outcomes and complete p0008 region/season audit

**2026-10-21 JST — source verification and next acquisition action**

## Latest immutable real-photo checkpoint

The new artifact-triggered data-only checkpoint [39452358860](https://github.com/zuizui0223/mois/actions/runs/39452358860), from [the source-audited Release](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39452358860-v1), finished tests and full original-photo recovery successfully. It has **702,362 distinct actual original-photo terminal IDs** and **330,221 genuinely automated full-resolution four-hue classifications**.

The previous [39450705160 checkpoint](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39450705160-v1) had 701,338 attempted / 329,812 classified, so this artifact-triggered update adds **1,024 p0032 original photo attempts and 409 valid four-hue ROIs**, not any photo-target metadata rows. This increment has 1 acquisition failure and 614 acquired-but-unclassifiable photos.

| Source | Actual finished original photo attempts | 4-hue classified outcomes | Still pending |
| --- | ---: | ---: | ---: |
| Source-authenticated historical FCP/MOIS/first p0002 photo baseline | 612,483 | 295,801 | n/a |
| Complete p0004 | 18,119 | 7,085 | 0 |
| Complete p0008 | 36,375 | 13,874 | 0 |
| p0016 unfinished original tier03 | 18,432 | 7,080 | 55,440 |
| p0032 unfinished original tier04 | 12,288 | 4,767 | 135,783 |
| Geo05 independent tropical/southern source | 4,096 | 1,415 | 18,182 |
| 569-image full-resolution shadow pilot | 569 | 199 | 0 |
| **Unique actually attempted original-photo IDs** | **702,362** | **330,221** | **209,405** |

All current known worker batches for p0016, p0032 and geo05 continue working; none of their originally scheduled pending photo IDs is imputed measured. The mainline keeps existing workers unchanged and uses the source-SHA artifact watcher instead of initiating duplicate original photo downloads.

**Independent output integrity:** artifact 12278318088, ZIP size 39,765,618 bytes, SHA-256 \`dfeb3064a2068abea3fed9af5fef6dbc0bb46b1cd1c92baced2f7cfc3b3fbbac\`. All three published output data-file hashes match the recovery receipt, and independent CSV summation confirms 702,362 photo terminals and 330,221 classifiable images. The source/ROI photo-identity validation ran in production. The local runtime verified ZIP/CSV/JSON SHA but lacks a local Parquet engine for re-reading the combined 702k photo rows.

## Complete p0008 actual photo-specific latitude × original observation month result

A separate, reproducible real-photo analyzer [PR #140](https://github.com/zuizui0223/mois/pull/140) passed CI and ran successfully as [production 39452760985](https://github.com/zuizui0223/mois/actions/runs/39452760985). Unlike source metadata checks, this **opened and validated every one of 36,375 p0008 original actual-image Parquet rows** using the runner's pyarrow installation. Results were published immutably at [mois-complete-p0008-original-photo-latitude-season-coverage-v1](https://github.com/zuizui0223/mois/releases/tag/mois-complete-p0008-original-photo-latitude-season-coverage-v1).

| Source latitude band | Source photo attempts | Actual 4-colour ROI outcomes | Measurable share | White among measurable |
| --- | ---: | ---: | ---: | ---: |
| Northern extratropics | 29,708 | 11,678 | 39.31% | 49.05% |
| Southern extratropics | 4,564 | 1,479 | 32.41% | 46.32% |
| Tropics | 2,103 | 717 | 34.09% | 48.68% |

Real botanical coverage remains extremely uneven: **20,773 distinct observations**, 8,191 observations with multiple photos; **6,549 nominal source taxon IDs attempted**, but only 4,453 with at least one successfully classifiable ROI and only 18 with at least 30 labelled photos. That is a **photo-level dataset**, not thousands of taxa each independently estimated for species/within-taxon seasonal flower whiteness.

Descriptively, among classifiable p0008 photographs, the white fraction was **53.14% for observations recorded in June** and **43.55% in October**. But taxon/observer/geography/phenology/photo selection and the seasonally unequal 30–42% ROI classifiability rate all confound this contrast. It is not evidence that individual taxa lose pigment in response to seasonal temperature or water. The original observed dates span 1974–2026, so no single year is held constant.

The direct p0008 output includes source-SHA-stamped per-latitude and per-original-month CSVs, the 36-region-month photo-table cells, per-latitude nominal taxon/observation coverage, and full JSON receipts. Independent ZIP/CSV hashes of [artifact 12277924138](https://github.com/zuizui0223/mois/actions/runs/39452760985/artifacts/12277924138) match SHA-256 \`142bceb023d3429c927bc9aaf988d5bdc9beec9b5648bdfc3a44009b0291e30e\`.

## Next acquisition and research line

1. Finish the already-running **p0016/p0032/geo05 original-image ROI-v4 workers** without changing the source selection and without redownloading previously attempted photo IDs.
2. Maintain the globally comparable observation-probability sample separate from the independently enriched geo05 tropical/southern original photos. Both source arm and real per-original-observation inclusion probabilities are saved in the combined checkpoint metadata projection.
3. Treat region/time white shares as photo-level descriptive results and model image-classification missingness. A within-species seasonal/abiotic claim needs within-taxon geographic and observed-date matching and more repeated genuine classifiable flower photos.

The upstream 93,209,351 original photo IDs remain a complete one-pass acquisition *target register*, not 93m measured flower hues. No original raw licensed photo bytes were redistributed in this audit.
