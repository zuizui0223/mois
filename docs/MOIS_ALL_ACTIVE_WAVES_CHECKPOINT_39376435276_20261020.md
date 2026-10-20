# MOIS: four-wave checkpoint recovery, 20 October 2026

## Published actual results

Production recovery run [39376435276](https://github.com/zuizui0223/mois/actions/runs/39376435276) completed validation, artifact upload and immutable prerelease publication. The data release is [mois-real-pixel-checkpoint-39376435276-v1](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39376435276-v1).

The new checkpoint contains **628,412 unique original-photo terminal outcomes**, of which **302,026 are automatically classifiable into the four frozen hue categories**. Compared with checkpoint 39373763787 (625,340 / 300,828), this is **3,072 additional completed attempts and 1,198 additional classified photos**. Attempt counts include unsuccessful acquisitions; they are not counts of usable flower photographs.

| Source | Completed attempts in this checkpoint | Four-hue classified | Planned attempts not in completed artifacts |
|---|---:|---:|---:|
| Original FCP/MOIS/p0002 baseline | 612,483 | 295,801 | Not applicable |
| p0004 | 13,312 | 5,246 | 4,807 |
| p0008 | 2,048 | 780 | 34,327 |
| Completed tier03 shadow cohort | 569 | 199 | 0 |
| p0016 | 0 | 0 | 73,872 |
| **Combined measured/attempted IDs** | **628,412** | **302,026** | **113,006 pending targets across these three unfinished waves** |

Zero recovered p0016 rows means no sealed worker artifact was present in the recovery snapshot; it does not imply that its active workers have processed no images. Pending targets are never counted as measured, nonwhite or non-flowering. This is a compute-completion checkpoint, not a completed probability cohort or a global flower-colour frequency estimate.

## Exact additions since the previous checkpoint

- p0004: +1,024 attempts, +1,020 acquired images, +418 classified photos, +4 acquisition failures and +602 acquired but unclassifiable photos.
- p0008: +2,048 attempts, +2,044 acquired images, +780 classified photos, +4 acquisition failures and +1,264 acquired but unclassifiable photos.
- Total increment: **3,072 attempts = 3,064 acquired images + 8 acquisition failures**; acquired images comprise **1,198 classified + 1,866 unclassifiable**. All existing hue labels and unknowns remain unchanged.

The independently downloaded p0004 group-12 artifact contains units 48–51: 208 white, 133 yellow-orange, 31 red-pink and 46 blue-purple outcomes (418 total), plus 606 mixed/uncertain records including four acquisition failures. These are acquisition outcome counts, not ecological estimates.

## Recovery change deployed

PR [#129](https://github.com/zuizui0223/mois/pull/129) merged as `bbd17f1f2ce22a17a547980de471e086a0422ae7`. The existing checkpoint registry/downloader now includes p0016 run 39374811084, its exact 73,872-photo source plan and SHA256 `16baf617396dd7b388debc65fd25eb8f9041f49ca5e6f78d0b3e2898f6dfeed5`.

The checkpoint workflow also reacts to completion of the four exact pinned acquisition run IDs, with same-repository and main-branch checks. It recovers only independently validated completed units, even when some other units remain absent. The automatic event path is configured and covered by wiring tests; the successful run documented here was triggered by the main-branch push, not by a later acquisition-completion event.

Active acquisition workflows were not cancelled or modified. No new model, low-resolution exclusion, native-status hard filter, image downloads or colour reinference were performed by recovery.

## Verification and artifact

Artifact **12243149408**, `mois-actual-sealed-pixel-checkpoint-39376435276`, is 22,739,339 bytes. ZIP SHA256: `215d7eead4faeeee1f4ab320af0d6e72c27f6abd8d9c5f769368d962b5240da1`.

After download, the ZIP digest and all output checksums in `checkpoint_receipt.json` matched. The CSV hue counts were independently summed by source wave and checked against every receipt count; totals reconcile to 628,412 attempted IDs and 302,026 classified outcomes. The prior checkpoint ZIP was also rechecked. Full source/terminal/photo-ID/method reconciliation and Parquet row validation ran in the successful production job. This local container lacked a Parquet engine, so local verification covered file digests and CSV/JSON arithmetic, not an independent re-read of every Parquet row.

The five new recovery-wiring tests passed locally after three expected pre-change failures. The existing checkpoint test command and the new wiring test command passed in CI, and all four PR CI workflows passed. The full repository test suite was not executed locally because repository cloning was blocked by DNS; do not infer a full-suite local test claim.

Files in the checkpoint remain: `combined_actual_photo_checkpoint.parquet`, `new_completed_original_photo_ROIs.parquet`, `completed_photo_hue_counts.csv`, and `checkpoint_receipt.json`. Original raw photos are not redistributed.
