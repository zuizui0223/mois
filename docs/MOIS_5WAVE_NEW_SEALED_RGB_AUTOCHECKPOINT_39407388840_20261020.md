# MOIS real original-photo RGB: first fifth-wave pixels and source-aware checkpoint automation

**20 October 2026 — source-authenticated partial checkpoint**

## Real data confirmed, not metadata photo candidates

A data-only GitHub Actions run [39407388840](https://github.com/zuizui0223/mois/actions/runs/39407388840) finished successfully and published the immutable [mois-real-pixel-checkpoint-39407388840-v1](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39407388840-v1).

| Original-photo source wave | Completed attempted photo IDs | Frozen four-hue classified | Remaining scheduled attempts |
| --- | ---: | ---: | ---: |
| Original signed FCP/MOIS/first93m pixel baseline | 612,483 | 295,801 | n/a |
| p0004 original RGB fully complete | 18,119 | 7,085 | 0 |
| p0008 original RGB, partial | 20,480 | 7,883 | 15,895 |
| p0016 original RGB, partial | 8,192 | 3,149 | 65,680 |
| **p0032 original RGB, first completely sealed unit** | **1,024** | **400** | **147,047** |
| Completed original 569-image 320px shadow study | 569 | 199 | 0 |
| **Cumulative true original-photo terminal outcomes** | **660,867** | **314,517** | **228,622** |

Compared with [checkpoint 39405986312](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39405986312-v1), the incremental increase is **1,024 genuinely attempted original-photo IDs** and **400 more full frozen ROI-v4 four-colour classifications**. The new p0032 worker group acquired 1,021 photos, failed to acquire 3, assigned one of four hue labels to 400 and left 621 acquired originals technically unclassifiable. The source arms are mutually photo-ID disjoint and the currently pending scheduled targets are never counted as flowers or hues.

Total new terminal photo attempts beyond the immutable 612,483-photo baseline are 48,384, including 90 acquisition failures and 18,716 machine-classified four-hue photos. Historical source-photo labels and unknowns remain unchanged.

## Exact acquisition source versus model output

The original **93,209,351** photo-ID frame is a metadata and observation register. The p0032 source contains 148,983 signed candidate photo IDs in interval [0.0016,0.0032), minus 912 previously attempted originals = **148,071 new photo acquisition targets**. It is *not* a set of 148,071 measured flower-colour photographs. The original 83,835 observations remain grouped across 579 work shards; all source hashes and zero overlap with the earlier four source queues were required before original-photo acquisition.

All p0004 original photos were measured before this checkpoint. The p0008/p0016/p0032 full-resolution FCP YOLO11/EfficientSAM workers remain active. No original source photo was redownloaded or recoloured for the checkpoint. Neither native-range hard filtering nor low-resolution exclusion was introduced; the 320px pilot previously demonstrated nonzero white-flower false negatives.

This compute-completion snapshot is **not a complete probability cohort or an unbiased estimate of native global white-flower prevalence**. Floristic/taxon and geographic selection, multiple photos per observation and image missingness must be handled before macroecological inference.

## New checkpoint efficiency improvement

[PR #134](https://github.com/zuizui0223/mois/pull/134) merged the independent GitHub Actions workflow [MOIS checkpoint real RGB only when new sealed photo units arrive](https://github.com/zuizui0223/mois/actions/workflows/mois-auto-checkpoint-when-new-pixels.yml) and the pure SHA-preserving decision helper `scripts/acquisition/mois_checkpoint_fresh_artifact_gate.py`.

- A **GitHub Actions schedule at minute 17 every three UTC hours** checks the latest published checkpoint receipt against the signed worker artifact IDs and SHA-256 for all five original source arms. There is also a manual workflow trigger.
- **No new completed worker artifact → no additional checkpoint release.** A changed/missing previously sealed worker SHA or unsupported/duplicate source run is an error, not a reason to silently discard earlier results.
- A new sealed worker artifact and no other active collector → the watcher invokes the *existing* immutable, data-only collector by GitHub workflow_dispatch. It never fetches original photographs, launches ROI inference, or changes any flower labels.
- All eight source-gate unit and workflow contract tests passed in PR CI, alongside the repository source, seasonal and botanical-stage checks.
- **Production watcher run [39407356537](https://github.com/zuizui0223/mois/actions/runs/39407356537) succeeded**, observed precisely one new p0032 sealed worker artifact and no other unseen wave workers, and dispatched the above successful checkpoint 39407388840. The absence-of-new-data skip behavior is regression-tested; no scheduled 3-hour cron execution is claimed to have occurred yet.

## Immutable evidence

- Successful run: 39407388840 (data-only original photo collector)
- Immutable Release: `mois-real-pixel-checkpoint-39407388840-v1`
- Artifact ID: **12256868789**
- ZIP SHA-256: `bc0239a912ccc618780e6d0838165b9075dd70c7c8ee34012bd92df0cf26fab1`
- ZIP size: **30,246,236 bytes**
- Data files: `combined_actual_photo_checkpoint.parquet`, `new_completed_original_photo_ROIs.parquet`, `completed_photo_hue_counts.csv`, `checkpoint_receipt.json`
- Independent container-based verification recomputed the ZIP SHA-256, all three data output SHA-256 digests and independent CSV hue counts. All matched the receipt and **660,867 / 314,517** numbers. The row-by-row Parquet/source SHA reconciliation ran in the successful production collector; it was not independently repeated in the local container, where a Parquet engine was unavailable.

The next operational action is the existing p0008/p0016/p0032 worker processing of the remaining **228,622 selected** original photo IDs. New complete workers may produce new authenticated checkpoint Releases without any redundant original-image HTTP work.
