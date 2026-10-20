# MOIS: fifth tier04 p0032 verified original-photo source and five-wave checkpoint

Verified on 20 October 2026 (JST). This is an audit of source-target identity and genuine pixel outcomes, **not** a claim that all 93,209,351 source photo-ID metadata records have acquired image-colour labels.

## Acquisition mainline

The [fifth-wave original RGB run 39401484641](https://github.com/zuizui0223/mois/actions/runs/39401484641) was launched after merging [PR #132](https://github.com/zuizui0223/mois/pull/132). Its source-shard integrity and model contracts passed; model workers 0 and 1 entered original-image acquisition and full FCP frozen ROI-v4 evaluation, without 320px screening exclusion and without hard native/introduced photo filtering.

| Full-source nested interval and fifth-wave source | Count |
| --- | ---: |
| Source-signed original target photos in disjoint observation-hash tier04 [0.0016, 0.0032) | 148,983 |
| Previously attempted original photo IDs removed from new image queue | 912 |
| **Fifth-wave genuinely new original photos to attempt** | **148,071** |
| Distinct original observation IDs | **83,835** |
| Original FCP ROI source-photo work shards (at most 256 photos each) | **579** |
| Jobs of up to four work shards | **145** |
| Unsupported image-format targets | 0 |

The full upstream 93,209,351-photo target register was **not** re-enumerated. Only the signed tier04 archives were loaded. The source queue is independently SHA-verified: artifact [12254939208](https://github.com/zuizui0223/mois/actions/runs/39401484641/artifacts/12254939208), ZIP SHA-256 `b9b044f2c7b7b53a7c7cc34671de66a3985d463aee9282a0f6be55a19e15df37`; `manifest.json` SHA-256 `3628794fbb70f9e1515dabea6db731862265c0ad6c251732cfd9ab4df4badce8`.

Independent ZIP inspection verified **all 579 source CSV checksums**, exactly 148,071 distinct original photo IDs, 83,835 distinct observation IDs, and the frozen model input and source-tier fields. Pinned full-reservation checks recorded **zero overlapping photo IDs and zero overlapping observation IDs** against:
- p0004 signed source run 39353343403: 18,119 previously scheduled photos;
- p0008 source run 39371199127: 36,375 scheduled originals;
- p0016 source run 39374811084: 73,872 scheduled originals;
- completed tier03 shadow source run 39369517433: 569 attempted originals.

**Only 148,071 scheduled new photo attempts were frozen here, not 148,071 measured photo colours.** Completed/pending status is determined by original full-photo collector's signed terminal receipts, never the source metadata.

The fifth-wave workflow reused some old tier03 log-message prefixes `MOIS_P0016_...` in the new source job; those are stale diagnostic prefixes only. The signed source manifest actually verifies `full_source_original_tier: 4`, `selection_probability: 0.0032`, and the exact fifth-wave new photo count and source SHA. Do not infer p0016 provenance from the old print labels.

## New five-wave collector deployed and independently checked

[PR #133](https://github.com/zuizui0223/mois/pull/133) extends the existing data-only checkpoint registry with the fifth original source, pinned run 39401484641, exact 148,071-target plan and source SHA. Source/ROI terminal integrity gates are unchanged. The collector never downloads or reinfers images and does not interrupt p0008/p0016/p0032 acquisition jobs.

Production [checkpoint 39402230637](https://github.com/zuizui0223/mois/actions/runs/39402230637) passed tests and full data recovery, then published the immutable [five-source checkpoint](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39402230637-v1):

| Source | Genuine completed original-photo attempts | Frozen four-hue classified | Still pending |
| --- | ---: | ---: | ---: |
| Historical genuine FCP/MOIS/first93m photo union | 612,483 | 295,801 | n/a |
| p0004 complete | 18,119 | 7,085 | 0 |
| p0008 partial | 16,384 | 6,287 | 19,991 |
| p0016 partial | 7,168 | 2,780 | 66,704 |
| New p0032 tier04 | **0** | **0** | **148,071** |
| Completed 569-image shadow cohort | 569 | 199 | 0 |
| **Total completed unique original-photo IDs** | **654,723** | **312,152** | **234,766** |

Relative to [checkpoint 39400235870](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39400235870-v1) (653,699 attempted original photo IDs; 311,781 classifiable), this is +1,024 completed original p0016 attempts, +371 classifiable hues. The fifth-wave 0 means **no fully sealed original-image worker artifact existed at this checkpoint snapshot**; workers were actively processing and the outcome is not zero actual photograph images by assumption.

## Immutable source and reproducible accounting

Recovery ZIP artifact [12256660251](https://github.com/zuizui0223/mois/actions/runs/39402230637/artifacts/12256660251) is 28,843,351 bytes with SHA-256 `910dab4f866844a9d6f65168302baaf0bebb0319c1d79595a667a0d1ef723dc9`. Independent local checks verified the ZIP and the SHA-256 of each published data file against its receipt, and independently summed CSV source/hue rows to 654,723 terminal photo attempts and 312,152 four-hue labels. The full original-photo identity/method/source-SHA and image-terminal collector executed successfully in production; no independent local Parquet engine was available to re-read every Parquet row.

The newly measured original pixels beyond the historic 612,483-photo union had **88 download failures among 42,240 new attempted photo IDs**, with technical unknowns preserved. None are imputed coloured, native/introduced or nonflowering. The present snapshot is a **partial computational-completion cohort** and is not a valid naive global white-flower species-frequency estimator.

Next acquisition priority is the already-running original photo workers: p0008 19,991, p0016 66,704, and p0032 148,071 pending source photo attempts. Do not duplicate acquisitions. Once new worker groups publish sealed terminals, recover another source-authenticated **data-only** checkpoint without modifying the measured photo IDs or hues.
