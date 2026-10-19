# MOIS: full-upstream, bias-aware collection implementation

19 October 2026. Collection provenance only; no ecological regression, PCA, source photo reclassification or change to the active image executor.

## Actual source audit

The original 2026-09-27 count table was read in full and its SHA256 verified as `ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9`. It has 42,111 frame taxa, 42,095 with eligible photographs, and 93,209,351 photo candidates. The 1,000 most photographed taxa (2.375% of frame taxa) account for **45,464,071 photo candidates, 48.7763% of all candidates**; the top100 account for 12,973,635 (13.9188%). Median candidates per frame taxon: 212. These are acquisition-frame properties, not natural abundance or flower-colour results.

All eligible photos remain in scope, including multiple pictures from the same species or observation. Do not reintroduce a 32,768-total limit, new-genus-only selection, one photo/species, or one photo/observation. The upstream frame itself already restricts taxa, licences, identification and geographical eligibility; full coverage of that frame does not represent every flower in the world.

## Implemented design helpers

`scripts/acquisition/mois_sampling_metadata.py` adds two distinct operations:

1. For a FUTURE independently verified complete photo-ID index, freeze an observation-cluster priority design before looking at colours. A positive baseline probability `p` applies to every eligible observation. A second independent, metadata-stratum-specific top-up probability `q[h]` can prioritize underrepresented locations/environments. All photos within an observation inherit the same two draws. The union probability is `p + (1-p)*q[h]`, not `p+q[h]`. The helper records source and design SHA, sampling unit, both draws, the union probability and its inverse. Unselected records remain deferred in the full frame rather than excluded. Raising thresholds with the same seed, source and fixed strata gives nested milestones; setting `p=1` makes every observation/photo ready.
2. For existing deliberately species/genus-prioritized photos, create a sidecar WITHOUT overwriting original results. Global inclusion probability and inverse-probability weight are explicitly null. Neither `1` nor `8192/93209351` is a defensible retrospective probability for this legacy collection. Acquired-but-unclassifiable images and image acquisition failures remain separate. ROI failure does not prove that an image lacks a flower; botanical/human validation is `not_reviewed`, and missing flowering annotations are `not_recorded` rather than `false`.

The helper intentionally refuses a counts-only manifest. Per-taxon totals cannot supply a photo/observation/region sampling frame. A fake index-complete flag is not acceptable; the existing full-upstream scope validator must first verify the staged original photo-ID database and its per-taxon counts.

## Predeclared collection rules

- Freeze strata using source metadata: an appropriate equal-area geographic cell plus independent realm/environment information where available. Unknown environment is retained as a stratum, not rejected. One fixed stratum per observation is mandatory. Coordinate uncertainty that crosses an island boundary or fine altitude band stays uncertain; do not force it into a precise habitat.
- Set baseline/top-up probabilities from source availability and collection resources before hue outcomes. No fixed production fractions are asserted in this change because the complete geographic source index is not yet available in the inspected archive.
- Selection is at observation level; every eligible attached photo stays available. A time budget pauses collection, not the inclusion set. Do not calculate weighted results from only the fastest-finishing records of an unfinished milestone.
- These probabilities correct this designed acquisition stage only. They do not correct nonrandom iNaturalist recording, licence/Research Grade selection, annotation missingness, invisible plants, or classifier failures.
- Flower-annotation-selected refill, deliberately underrepresented-environment supplements, and legacy species/genus-first photos retain distinct source-arm labels. An independently reviewed sample must contain both classifier successes and failures; green/brown, tiny or inconspicuous flowers must not silently become leaves/background/unknown without measurement validation.
- Preserve source date, photo/observation/observer/taxon IDs, licences, coordinates and their uncertainty, classification statuses, model version, and future sampling provenance. Do not require soil or every environmental covariate before collecting pictures. Season remains supplementary.

## Actual non-destructive application

Read the original 32 completed 256-photo source CSV units for the 8,192-photo snapshot. Every source-unit CSV checksum, terminal checksum, photo ID, frozen-roster metadata field, source model signature and original hue/status/image hash was checked. The saved Parquet bytes retain SHA256 `11be4908cbd92221bf5cc3df4c1fc6d9b9fa412ccd2702b90762b9b1abfc22fe`; the source roster retains SHA256 `93585144e539ce7d4fe6e69b24bd209a292f425b7ce84329d6cd8615177f5a8d`.

All 8,192 records are retained: **2,780 classified**, **5,362 acquired but unclassifiable**, **50 acquisition failures**. A separate gzip-JSONL sampling-provenance sidecar was produced for each record, with unknown global selection probabilities left null. Sidecar SHA256: `28ba56f979fc0f31c1050e68aff0a16eeecb5438df0aeb9053325ed4c21855f7`. Original photograph results were not edited or remeasured.

Eight new local unit tests passed after a failing pre-implementation test. They cover counts-only refusal, overlapping inclusion-probability accounting, same-observation photos, retained unknown-environment support, nested progress to the full census, immutable design fingerprints and legacy missing-probability/failure distinctions.

## Status and remaining dependency

This is implemented, tested collection-design support and an executed legacy provenance audit, NOT an operational 93-million-photo extraction. The inspected authoritative upstream archive is count-only and explicitly says all photo-ID manifests are not yet frozen. The full-source photo/observation/geography index still has to be built or restored and validated before these probabilities can be applied. The active cached-photo acquisition and its prefetch/model settings are unchanged. New images acquired by this audit: **0**. No new image workflow, massive provider transfer, ecological inference or human image validation was launched by this change.
