# MOIS: full remaining tier03 acquisition, 20 October 2026

## Actual current action

PR #128 was merged as `79b59dd061a6de3625fc70dfe6b68be17eb32ad9`. Production run [39374811084](https://github.com/zuizui0223/mois/actions/runs/39374811084) passed its QA and source-freezing jobs. Workers `118470723626` and `118470723759` both finished model initialization and entered **Acquire only new photos and run original ROI without screening exclusions**. These are active acquisition/measurement jobs, NOT completed 73,872-photo outcomes.

The existing p0004 and p0008 workers were not cancelled, restarted, or changed. New tier03 processing is limited to two parallel workers. No screening model or hue classifier has been changed.

## Exact source queue verified from the saved artifact

| Quantity | Count |
|---|---:|
| Original candidate photos in disjoint tier03 | 74,883 |
| Already processed photo IDs excluded | 1,011 |
| **Genuinely new photo attempts queued** | **73,872** |
| Distinct original observation IDs | 41,745 |
| Distinct original source taxon IDs | 9,421 |
| Work units of at most 256 photos | 289 |
| Groups of at most four work units | 73 |
| New candidates unsupported by the frozen loader | 0 |

All 569 previously completed shadow-study attempts are excluded through the full published checkpoint. The 1,011 exclusions also include other previously processed source IDs; they are not a botanical or nativeness filter. The queue reads only this tier's archives from the existing complete target register, not the full 93m photo/observation index again.

Photo-ID AND observation-ID overlap checks both returned zero against each entire pinned source queue: p0004 (18,119 targets / 10,416 observations), p0008 (36,375 / 20,773), and shadow (569 / 309). The exact proof receipts are stored alongside the frozen source CSVs. Pending targets from another queue are not counted as measured data.

This tier corresponds to the fixed observation-hash interval [0.0008,0.0016). It is the remaining increment after deduplication, not the complete cumulative p<=0.0016 dataset, an all-species census, or a set of confirmed flowering pictures. The original full-resolution ROI must still determine which acquired images yield usable flower colour.

## Current measured-data boundary

The latest verified published checkpoint remains **625,340 distinct attempted photo IDs / 300,828 automatically four-hue-classified photos**: [checkpoint 39373763787](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39373763787-v1). The new p0016 target count is not added to that completed-data count. No completed new p0016 pixel release existed at this record.

The 320px pilot decision remains unchanged: no low-resolution image exclusion. No native-only hard filter is enabled. Acquisition failures and nonclassifiable ROIs retain their original technical-unknown status. The next original-photo output for this increment is `mois-full93m-observation-fourth-original-RGB-p0016-v1`, publishable only after all source IDs reconcile.

## Verification scope and receipts

- Source job `118470527511`: production prior-checkpoint Parquet SHA, actual prior-ID deduplication, full-tier count reconciliation and all three overlap proofs succeeded.
- Source artifact `12243345937`, `mois-93m-p0016-precolour-observation`: 7,317,598 bytes. ZIP SHA256 `d4cdbf049a050139546c8f369f1a336e1522df879d0525737bbdde566cc98303`.
- Source manifest SHA256: `16baf617396dd7b388debc65fd25eb8f9041f49ca5e6f78d0b3e2898f6dfeed5`.
- After downloading the artifact, all 289 worker CSV hashes, 73,872 unique photo IDs, 41,745 observation IDs and 9,421 source taxon IDs were checked independently in the local container.
- Local independent re-reading of the prior Parquet was unavailable because that container has neither pyarrow nor fastparquet. The actual prior-ID check was executed successfully by the production source job; no claim is made that this particular check was repeated locally.
- Five new wiring/syntax tests failed before the workflow was present, then passed against the actual CI-published workflow archive. The new workflow CI and existing repository CI workflows passed. This is not a claim that every repository test was run locally.

Original image bytes are not redistributed. The saved source registry is metadata, never a replacement for completed floral RGB measurements.
