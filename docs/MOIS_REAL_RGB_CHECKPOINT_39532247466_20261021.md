# MOIS original-photo RGB checkpoint — 21 October 2026

## Completed genuine original-photo outcomes

A [new data-only production checkpoint (run 39532247466)](https://github.com/zuizui0223/mois/actions/runs/39532247466) passed source and terminal verification and published the [immutable release](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39532247466-v1).

**751,514 distinct original-photo terminal attempts**, **348,794 automated frozen-ROI four-hue classifications**, compared with previous [run 39528207538](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39528207538-v1) (**747,418 / 347,217**): **+4,096 attempted originals and +1,577 four-hue classifications**.

| Source | Actual completed photo attempts | Frozen four-hue classified | Still unprocessed source targets |
| --- | ---: | ---: | ---: |
| Immutable FCP/MOIS and original first p0002 baseline | 612,483 | 295,801 | Not applicable |
| Second complete p0004 | 18,119 | 7,085 | 0 |
| Third complete p0008 | 36,375 | 13,874 | 0 |
| Fourth p0016, partial | 37,888 | 14,602 | 35,984 |
| Fifth p0032, partial | 30,720 | 11,847 | 117,351 |
| Geo-enriched separate tier05, partial | 15,360 | 5,386 | 6,918 |
| Completed original-photo low-resolution shadow study | 569 | 199 | 0 |
| **All authentic original photo IDs, one terminal per ID** | **751,514** | **348,794** | **160,253** |

The **4,096 newly sealed terminal rows** comprised 2,048 p0016 photo attempts (+815 four-hue classified), 1,024 p0032 (+388), and 1,024 tropical/southern geo05 (+374). Exactly **4,085** additional originals acquired successfully, **11** acquisitions failed, **2,508** successfully acquired originals were *not* four-hue classifiable. Newly labelled hues across the three waves: **759 white, 530 yellow-orange, 170 red-pink, 118 blue-purple**. These are automated photographic ROI classes, not a verified native-species white fraction or proof of botanical anthesis.

## Independently authenticated source evidence

- Immutable GitHub Actions artifact **12315773272**, size **50,907,405 bytes**, ZIP SHA256 **`c18d8ee7e5e7a88cb1fbeb99723ad9aa603a04e9a39d399391526c412ee4d626`**.
- All three output files matched the SHA256 values stored in `checkpoint_receipt.json` and the original GitHub Release asset digests. The receipt SHA256 was **`cc24bfadc578a5c824f479b4628285f2990ab15bc2785d9c3621ba24cdefccbb`**.
- Independent sum of `completed_photo_hue_counts.csv` agrees with the full reported **751,514 attempted IDs / 348,794 classifiable** and with every individual source-wave total in the receipt.
- All prior original photo worker artifact IDs, names, byte counts and SHA256 digests were unchanged between checkpoints; exactly four new complete worker artifacts were captured: **p0016 groups 35 and 36, p0032 group 29, and geo05 group 14**. No worker was reacquired and no old original hue labels were recomputed.
- The frozen original image ROI/source SHA and per-photo output terminal checks were validated in the successful production collector. Independent local checks covered ZIP and member hashes, CSV/JSON arithmetic and artifact snapshot monotonicity, **not** full independent Parquet row reading because the local environment had no Parquet engine.
- `selection_probability` is retained in the combined Parquet projection, preserving geo05's original-observation-level heterogeneous inclusion probabilities and separating it from the nested uniform observation-hash probability arms.
- Local detailed independent check: `mois_checkpoint_39532247466_independent_audit.json`, SHA256 **`b46127eb5d4094ee0aca933d48dd96c86a850691104480237b9daf306b7384bf`** (local user-facing companion).

## Continuing original-photo work

The original image jobs still run independently and were not duplicated, restarted or converted to a lossy low-resolution photo exclusion scheme. Sealed group counts at this capture: complete p0004 **18/18**, complete p0008 **36/36**, partial p0016 **37/73**, partial p0032 **30/145**, partial independent geo05 **15/22**, completed shadow cohort. Remaining **160,253 already selected photo IDs** are targets, not measured flower colours. A SHA-aware scheduled watcher requests further data-only checkpoints *only when new completed ROI worker ZIPs appear*; current checkpoint successfully passed its freshness gate.

This is a **partial computational-completion** release, not a representative global wild-flora colour census from 93,209,351 upstream metadata photo candidates. No WCVP native-only hard filter, no region-biased arm pooling and no global or evolutionary white-flower population claim was made.
