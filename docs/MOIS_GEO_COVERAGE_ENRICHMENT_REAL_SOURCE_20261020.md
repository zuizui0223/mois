# MOIS independently acquired tropical/southern flower-photo arm — 2026-10-20

## Completed source acquisition milestone and current RGB status

**[Source freeze Actions run 39417990675](https://github.com/zuizui0223/mois/actions/runs/39417990675)** passed the source registry, original photo/observation identity and five-active-wave disjointness checks. Image workers for its new source cohort are now running the exact frozen FCP YOLO11/EfficientSAM full-resolution ROI-v4 method. **The 22,278 source targets are not 22,278 measured or flower-colour-classified images.** At the subsequent six-arm checkpoint capture, no geo05 unit was yet fully sealed.

The source photo opportunity catalogue of 93,209,351 original iNaturalist photo IDs is unchanged. The new arm reads **only tier05 [0.0032, 0.0064)**, 297,967 signed original photo metadata candidates, and never rescans the full 93m original index.

### Metadata-only, fixed-before-pixels selection design

A new independent SHA-256 observation hash with frozen seed `mois-original-observation-geographic-enrichment-v1` applies the same inclusion decision to *every photograph from the same observation*. Its selection fractions within tier05 are:

| Original latitude band | Source tier05 photo candidates | Metadata hash selection probability | Admitted original photo IDs before prior-image dedup | New photo download targets |
| --- | ---: | ---: | ---: | ---: |
| Tropical (23.44°S through 23.44°N) | 18,327 | 1/2 | 9,139 | **9,040** |
| Southern extratropical | 38,465 | 1/4 | 9,658 | **9,586** |
| Northern extratropical | 241,175 | 1/64 | 3,671 | **3,652** |
| Unknown latitude | 0 | 1/64 | 0 | 0 |
| **Total** | **297,967** | *heterogeneous* | **22,468** | **22,278** |

Original photos selected by the independent geography hash but previously attempted by the existing genuine photo source: **190**, excluded from new downloading. The other **275,499 tier05 candidate originals** remain in the immutable full target registry; they were merely not included in this separate coverage-enrichment arm. The source validator found **zero image-formats unsupported by the frozen loader** among admitted candidates.

The 22,278 new image IDs belong to **11,285 distinct original observation UUIDs**. The source worker set has **88 photo shards**, no more than 256 photos per shard, split into **22 processing groups**. No per-species, per-genus, or per-taxon cap was applied.

Every RGB photo work row records its *unconditional* original observation inclusion probability: tropical `0.0032 × 0.5 = 0.0016`, southern extratropical `0.0032 × 0.25 = 0.0008`, northern extratropical `0.0032 × (1/64) = 0.00005`. The parent signed tier interval, original observation SHA seed and complete pre-/post-dedup latitude counts are frozen in `manifest.json`. **Unlike a cumulative global random photo cohort, the manifest intentionally has no single scalar `selection_probability`**, because the inclusion probability varies by region. Ignore this arm in a naive or unweighted global white-hue fraction, rather than treating enriched photo counts as a floristic prevalence estimate.

The original latitude is only a coarse acquisition variable. It is not a botanical bioregion, a taxon-wide endemicity estimate, proof of wild provenance or evidence that individual photographs are native. No native/introduced/unknown taxon source is discarded in the upstream photo register.

### Exact source signatures and collision checks

- **Source artifact**: [12262054540](https://github.com/zuizui0223/mois/actions/runs/39417990675/artifacts/12262054540), ZIP size **2,232,678 bytes**, SHA-256 `07c52958fefe5d3b9656b0953e4a2c08a431af7065f1dcbcea6ff7d587840416`.
- **Source `manifest.json` SHA-256**: `3b80efe07287c026a34eec7636038e61fa5aa4e16d49aeb139a533bb46c8feda`.
- Independent download and inspection verified all 88 source CSV SHA256 digests, 22,278 distinct original photo IDs, 11,285 unique observation UUIDs, individual region/source-arm selection probabilities and one consistent inclusion probability per original observation.
- **Exact simultaneous photo AND observation overlap**: zero versus original signed p0004 (18,119), p0008 (36,375), p0016 (73,872), p0032 (148,071), and 569-photo shadow-study source queues.

### Measurement integration and automation

- [PR #135](https://github.com/zuizui0223/mois/pull/135) added the source-only latitude selection, independent observation hash, prespecified conditional probabilities and the frozen old FCP full-resolution ROI image worker. New geographic source QA, signed original source registry QA and repository checks passed.
- [PR #136](https://github.com/zuizui0223/mois/pull/136) included source-pinned `geo05` in the **data-only** original image terminal checkpoint collector. [Checkpoint 39418633371](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39418633371-v1) published **671,107 original photo attempts / 318,460 automatic four-hue image outcomes**, with **geo05 pending 22,278** at the capture time. This checkpoint added +1,024 p0016 attempts and +398 classifiable colours; it did not count pending geo05 candidates as measured images.
- [PR #137](https://github.com/zuizui0223/mois/pull/137) added the geo05 source run and immutable artifact prefix to the previously deployed source-aware three-hour checkpoint watcher. **Production watcher run 39419053086 passed, observed no newly sealed groups, and correctly skipped an unchanged checkpoint.** The source watcher cannot reinterpret source IDs as measured flowers or dispatch duplicate image downloads.

### Next operational step

Allow the existing single geo05 measurement worker and the concurrently running p0008/p0016/p0032 workers to finish original-photo RGB units. The checkpoint watcher will only dispatch source/photo-SHA reconciled **data-only** snapshots when *new fully sealed worker artifacts* exist. Keep unknown/failed acquisitions as technical unknowns. Do not allow the geographic-enrichment arm to contaminate the original probability-arm global or species-level white-flower inference. Do not enable the loss-inducing 320px image-exclusion screen.

**Research question still unresolved**: tropical and southern enrichment corrects sample coverage opportunities, but a defensible global flower-colour environmental relationship still requires botanical flower-stage review, species/taxon/observer/photo weighting, regional nativeness annotation and sensitivity to hue-dependent ROI failure. This source file by itself is not a result about natural selection.
