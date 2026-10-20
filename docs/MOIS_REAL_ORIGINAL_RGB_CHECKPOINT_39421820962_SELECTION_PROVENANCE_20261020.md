# MOIS: authentic original-photo colours and acquisition probabilities — 20 October 2026

## Executed result and current acquisition

**[Production checkpoint run 39421820962](https://github.com/zuizui0223/mois/actions/runs/39421820962) passed original source SHA/terminal photo identity checks and published immutable [data release 39421820962](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39421820962-v1).** It contains **676,227 unique original-photo terminal outcomes**, of which **320,352 photographs have machine-classifiable four-hue ROIs**. These are genuine original-image processing outcomes, not merely iNaturalist photo ID opportunities.

| Frozen source wave | Completed attempted photo IDs | Automatically four-hue classifiable | Original-photo attempts pending |
|---|---:|---:|---:|
| Historical source-verified FCP/MOIS/first93m RGB baseline | 612,483 | 295,801 | Not applicable |
| p0004 second source, fully sealed | 18,119 | 7,085 | 0 |
| p0008 disjoint third source, partial | 27,648 | 10,585 | 8,727 |
| p0016 disjoint fourth source, partial | 12,288 | 4,767 | 61,584 |
| p0032 disjoint fifth source, partial | 5,120 | 1,915 | 142,951 |
| geo05, independent tropical/southern enrichment, partial | 0 | 0 | 22,278 |
| Tier03 original 569-image shadow comparison | 569 | 199 | 0 |
| **Combined distinct original photo IDs** | **676,227** | **320,352** | **235,540** |

The new results since [checkpoint 39420718036](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39420718036-v1) consist of **1,024 new p0032 original-photo terminal outcomes / 360 more classified photos**, with **1,023 acquired images**, one acquisition failure, and 663 acquired but technically unclassifiable ROIs.

Across 63,744 complete post-baseline photo attempts, **111 acquisitions failed**, 24,551 were classifiable into one of four machine hue states, and 39,082 acquired images remained unclassifiable. All original hue labels and technical unknown/failure statuses were retained. Pending source photo targets are not counted as measured.

## Sampling provenance fixed in the cumulative Parquet projection

The older cumulative `combined_actual_photo_checkpoint.parquet` was a compact projection that **omitted per-photo `selection_probability`**, even though source-signed worker data contained the field. This made downstream selection-bias correction unsafe for the new enrichment cohort: it uses deliberately **heterogeneous observation-level acquisition probabilities**, unlike the global hash arm.

Merged [PR #138](https://github.com/zuizui0223/mois/pull/138) adds `selection_probability` unchanged to the cumulative metadata projection while continuing to preserve `source_arm`, `checkpoint_wave`, photo ID, observation ID, image source SHA, hue label, and acquisition failure indicator. It checks that any provided new probability is finite, strictly positive, no greater than 1, and identical for all photographs of the same original observation within the same source wave. Historical baseline rows with genuinely missing selection probabilities remain blank, not fabricated.

The regression tests were observed RED before the change, then all relevant CI checks passed after it. The production collector then passed and produced this new immutable release. It **did not** re-download original image bytes, rerun YOLO11/EfficientSAM ROI-v4 or change any machine flower colours.

**Important inferential qualification:** The retained values are the original acquisition-design source probabilities, **not automatically valid complete estimator weights** for the worldwide flora. In the independent geo05 tier [0.0032,0.0064), prespecified probabilities differ by latitude: tropical 0.0016, southern extratropical 0.0008, northern extratropical 0.00005 (unknown-latitude 0.00005). Earlier global source tiers record cumulative hash-threshold probabilities; these cannot naively be treated as each disjoint tier's incremental inclusion probability after earlier source-photo deduplication. Full baseline FCP/MOIS photographic recruitment, observation clustering, floristic/taxon sampling, photo detectability and classifier failure require separate treatment. The geographic-enrichment arm must never be pooled unweighted with the random arm to claim global or island-native white-flower frequency.

## Verified output authenticity

- GitHub release: `mois-real-pixel-checkpoint-39421820962-v1`
- Original Action artifact: `mois-actual-sealed-pixel-checkpoint-39421820962`, artifact ID **12264945380**.
- ZIP size: **33,744,565 bytes**. ZIP SHA256: `9ca1a9e371c2ee6100f01114cc7512a5203bcc23b2287452a362113957aef222`.
- File SHA256 from signed receipt all passed independent recomputation on the downloaded ZIP: `combined_actual_photo_checkpoint.parquet` `1b23393a21f7c33fb2f68b46c2572d1d00f654da626872ca07be29214d9b133a`, `new_completed_original_photo_ROIs.parquet` `d8da6995ffeeed63a7a0e5bc007440d31278e6101ccb060d48ae3882625e7f1d`, `completed_photo_hue_counts.csv` `1ff8b98ad18c22351caeaad795c389f4d7fc12ed71c790564555a36b933093b1`.
- Independent CSV aggregation reproduced **676,227 genuine terminal outcomes and 320,352 four-hue labels**, with all six waves and historical baseline accounted for.
- Full original source/terminal metadata per-row validation and Parquet materialization took place in the successful GitHub production run. Independent container re-reading of every Parquet row was not available because the local container lacks a Parquet reader.
- This release is a **partial computation-completion subset** of the selected original-photo queues; it is not a probability-complete worldwide-native-flora composition estimate.

## Mainline

Continue existing p0008/p0016/p0032 and geo05 full-resolution image workers for **235,540 pending original-photo targets**. The source-aware GitHub Action checks every three hours and publishes a new immutable verified data-only checkpoint only when one or more new fully sealed photo worker artifacts appear; no repeated photo acquisition occurs. No low-resolution image exclusion, and no native-status hard-photo filter.
