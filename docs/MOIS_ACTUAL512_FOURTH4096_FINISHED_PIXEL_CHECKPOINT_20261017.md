# Original MOIS fourth wave: source-verified 512 already processed RGB photographs, not 4,096 candidates

**Date:** 2026-10-17. **Scope:** true original photographic pixels, automated flower-like ROI-v4 hue classification, and original observation date. No new data are regarded as botanically confirmed open petals, and no FCP within-species W/D/50-km genetic-spatial inference is imported into MOIS.

The original [309,439 actually measured photographic photo-ID cohort](https://github.com/zuizui0223/mois/actions/runs/39028400690) has a signed original source-data SHA256 `546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e`, with **149,662 original automatic four-hue photo outcomes** and **149,653 authentic dated+classified ROIs**. The separate [fourth 4,096-photo source/processing workflow](https://github.com/zuizui0223/mois/actions/runs/39032258023) first selected 4,096 genuinely previously unmeasured original photo IDs from a 3,430,502-source-photo iNaturalist metadata corpus. Its frozen source selection SHA256 is `65fc139909f94448e042b6539c6749e70d1ab21bb8b13f50383bb7417180337a`; source IDs alone are NOT hue measurements. The source's observer IDs are disjoint within the fourth cohort, not necessarily across old records.

## Source pixels already genuinely processed

Two original fixed-256-photo ROI-v4 shards of the fourth cohort have already completed and produced original image-pixel classifier terminal tables. All photo IDs, each of the four original 64-image subterminal receipts, unmodified model SHA and fixed source selection metadata are verifiable by [actual Fourth4096 Actions run](https://github.com/zuizui0223/mois/actions/runs/39032258023).

| Genuinely completed photographic output | Fixed source shard 0 | Fixed source shard 1 | Actual two-shard total |
|---|---:|---:|---:|
| Real original RGB photo IDs with ROI-v4 terminal outcomes | 256 | 256 | **512** |
| Four-hue automatically classifiable original photos | 108 | 95 | **203** |
| Technically hue-unknown photos | 148 | 161 | **309** |

This is **not** an independent fourth dataset to add to the final 4,096 once it finishes; these exact 512 photos are a frozen **subset** of that still-running cohort. Upon completion, the final 4,096 result *supersedes* this interim checkpoint, with photo-ID deduplication rather than summing both.

## Exact photo-date-source join

The source-only photo IDs are joined to their originally recorded observation date, taxon and geographic location only after real RGB pixel inference. The new [reproducibility workflow](../.github/workflows/mois-checkpoint-512-actually-finished-fourth4096-photo-pixels.yml) restores both completed worker artifacts plus the full source selection and the old 309,439-photo original GZIP; verifies each source `measurement_id`, original `photo_id`, `inat_taxon_id`, `observation_id` and `observer_id`, source-worker/ROI-v4 model attestation, both full 256-photo terminal receipts, each per-shard input/output SHA256 and automatic four-hue classifiability. No original photo ID may overlap the 309,439 measured source; no old source hue/taxon/date may change.

On completion the **checkpoint** source can legitimately describe **309,951 real original photo-ID RGB/ROI attempts, 149,865 automatic four-hue labels and 149,856 date + hue labels**. The other **3,584** photo IDs of the original fourth cohort still require genuinely executed photo-pixel terminal outcomes. Any additional completed worker artifacts will be followed by a *new separately signed snapshot*, not silently counted by this specific 512-photo source-cohort receipt.

**Ecological limits:** The normal fourth-wave source is photo-ID-hashed and taxon-balanced across 1,000 taxa, not a randomized field flower-abundance census. A two-shard **execution-order** partial subset is even more susceptible to computation-completion selection; it cannot serve as an unbiased standalone estimate of global or local summer–winter pigment responses. The 309,439-photo original equal-genus climate model must not be passed off as a new 309,951-photo refit without actually reprocessing climate covariates and 999 genus bootstrap draws. The original 256-shard FCP dataset and the forth 4,096 source are distinct but can share original source-photo IDs, always deduplicated by literal photo ID. Source taxon/photo ID, visible ROI hue and local date still do NOT botanically verify open petals or distinguish flowers from bracts/fruit.

Publication implication: photograph growth counts use real source image-pixel terminal outcomes, with the 512-row checkpoint archived as a signed GitHub Actions artifact. The peer-reviewed ecological inference should continue to use the complete 309,439-photo base until the fourth wave actually completes or a bias-aware partial-sample design is reported.
