# MOIS: verified completed-pixel checkpoint request (2026-10-20)

This source-backed milestone intentionally triggers the existing **data-only** checkpoint workflow on main. It reads **previously completed and SHA-verified original-photo work units**; it neither downloads raw iNaturalist images nor reruns flower-colour inference. The ongoing p0004, p0008 and p0016 acquisition jobs are left untouched.

## Source snapshot that justifies a new immutable checkpoint

At the source check around 06:20 UTC / 15:20 JST on 2026-10-20:

| Fixed original-colour acquisition run | Completed sealed RGB worker groups | Total planned photos |
| --- | ---: | ---: |
| p0004 `39353343403` | 16 of 18 | 18,119 |
| p0008 `39371199127` | 5 completed | 36,375 |
| p0016 `39374811084` | 1 completed | 73,872 |
| shadow `39369517433` | 1 completed, full cohort | 569 |

The preceding immutable checkpoint is [39376435276](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39376435276-v1), with **628,412 unique original-photo attempts** and **302,026 machine-classified four-hue outcomes**. More worker artifacts have arrived since that checkpoint, but their new photo counts/classifiable labels MUST be recomputed by the validated production collector, not projected from the scheduled source rows.

The checkpoint source registry pins these exact four run IDs, their original photo roster SHA256 and the immutable FCP/MOIS prior baseline. New and previous photo IDs must remain disjoint, missing units must stay *pending*, acquired-but-unclassifiable images cannot become coloured flowers, and legacy unknown labels are preserved exactly.

## Operational result criterion

A successful [checkpoint workflow](https://github.com/zuizui0223/mois/actions/workflows/mois-completed-pixel-checkpoint.yml) will publish `mois-real-pixel-checkpoint-<run-id>-v1` with `checkpoint_receipt.json`, both Parquet outputs and photo/hue counts. Only after the run is successful and the new receipt matches all file digests should its totals be described as confirmed.

This milestone path is an explicit source-backed way to request an additional checkpoint during active long-running image acquisition without modifying or restarting the measurement workers. It does not schedule repeated artificial snapshots.

## Refresh request: complete p0004 source + new p0008/p0016 group terminals

This new source-verified snapshot was requested after the full p0004 terminal union was separately published as [mois-genuine630602-original-pixels-first-two-93m-waves-v1](https://github.com/zuizui0223/mois/releases/tag/mois-genuine630602-original-pixels-first-two-93m-waves-v1). The **previous partial pixel checkpoint** was [39381982167](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39381982167-v1), with 639,363 attempted photo IDs and 306,186 machine-classified four-hue images. The source-authenticated p0004 complete increment was independently verified at 18,119 attempts / 7,085 classifiable.

New archived image-worker groups are available beyond that checkpoint:
- p0004 39353343403: all 18/18 complete (unchanged since previous checkpoint)
- p0008 39371199127: 8 completed source-pinned groups; four acquisition workers active, none failed
- p0016 39374811084: 3 completed source-pinned groups; two acquisition workers active, none failed
- full shadow study 39369517433: 569 attempted and 199 classified (unchanged)

**Do not extrapolate a new classified count from group counts.** The existing collector must SHA-check and reconcile each original source photo / observation ID, byte-immutable model method and terminal outcome, then publish the exact newly classified and still-pending counts as an immutable checkpoint. Do not rerun colour inference, re-fetch images or add any low-resolution filter.

This documentation update intentionally triggers the existing one-off data-only `main` checkpoint. It does not schedule a recurring task or change image-acquisition workers.

## Further verified acquisition milestone — 2026-10-20, around 16:42 JST

After the previous immutable [checkpoint 39385699414](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39385699414-v1) published 642,435 terminal photo IDs with 307,382 frozen automatic four-hue classifications, the following additional worker files reached terminal-success status:

* p0004 39353343403: complete, all 18 of 18 archived; do NOT re-download or re-run
* p0008 39371199127: 10 of 28 completed original-photo worker groups; four currently processing, zero failures
* p0016 39374811084: 4 of 28 completed worker groups; two currently processing, zero failures
* tier03 shadow 39369517433: 569 original-photo attempts, already complete and permanently source-deduplicated

Use the pinned production collector on only the stored source SHA and actual image terminal receipts, without running image acquisition or colour inference. Publish an immutable intermediate checkpoint after verifying all newly completed units. The total new photo IDs and usable hue counts must come from exact original-photo terminal reconciliation, never assume one full group's pixels based on queued count.

This source-documented milestone is a one-off re-capture request through the existing push-triggered, data-only checkpoint workflow. It does not add an automatic recurring task.


## Source-verified checkpoint request — 18:21 JST, 2026-10-20

A fresh paginated GitHub Actions jobs read shows new completed **original pixel measurement** worker artifacts since checkpoint [39388432954](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39388432954-v1) (645,507 attempted photo IDs / 308,544 classifiable four-hue ROIs):

- p0004 original source run 39353343403: 18/18 full measurement workers successfully sealed, already in baseline checkpoint.
- p0008 new source run 39371199127: **16/36** sealed successful groups, 4 active, 16 waiting, **0 failures**; each completed group retains all original per-photo ROI outcomes.
- p0016 new source run 39374811084: **6/73** sealed successful groups, 2 active, 65 waiting, **0 failures**.
- Previous shadow pilot 39369517433: complete; its 569 original photo IDs are preserved and must not be counted twice.

The previous snapshot only contained p0008 **10** groups and p0016 **4** groups. This source-backed request triggers one more **data-only** recovery on the main branch. It must checksum-reconcile original photo/observation/source-shard IDs, immutable classifier metadata and terminal CSVs; never download the photos again, run new colour inference, or impute missing image colours.

**Outcome counts must be read from the new production recovery receipt, not inferred from these computational group counts.** Acquisition continues independently. The incremental snapshot is not a probability-complete sample or a global white-hue prevalence estimate.

## Additional sealed-pixel recovery: 2026-10-20 afternoon, fifth-wave in progress

The previous immutable [five-source checkpoint 39402230637](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39402230637-v1) was SHA/released with **654,723 attempted original photo IDs** and **312,152 four-hue classifiable photos**. It includes all 18 p0004 groups, 16 p0008 groups, 7 p0016 groups and all 569 shadow original-photo attempts; the new p0032 tier was still at zero **sealed** work units at that snapshot.

New signed worker artifacts checked against the exact original acquisition runs:
- p0004 39353343403: **18/18** fully complete, already published as original frozen true RGB and not to be downloaded again.
- p0008 39371199127: **20/36** complete worker groups; 4 actively measuring, 12 queued, **no observed job failures**.
- p0016 39374811084: **8/73** complete worker groups; 2 actively measuring, 63 queued, **no observed job failures**.
- p0032 39401484641: fifth-wave source already verified as **148,071 new candidate originals from 83,835 observations**, zero worker groups yet sealed, 2 actively measuring, with no observed failures. GitHub Jobs pagination is mandatory; counts inferred from the first page (100 jobs) are incomplete when >100 matrix jobs exist.
- tier03 569-image shadow validation already sealed; no false exclusion enabled.

Invoke the existing code- and source-authenticated **data-only** recovery exactly once via this documentation push; publish a fresh immutable checkpoint only after original source photo IDs, terminal model method, SHA and full manifest counters agree. **Never infer a new classifiable-photo count from worker-group completion alone.** Preserve pending photos and technical unknowns rather than imputing hue, and do not re-download or re-infer originals.

The already-launched p0008/p0016/p0032 workers are unchanged and remain active. This request is explicitly not another source-enumeration run, native-only photo filter, or periodic polling task.

## New sealed-image checkpoint request — 2026-10-20, about 20:13 JST

The latest published five-wave real-pixel checkpoint is [39407388840](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39407388840-v1), containing **660,867 unique terminal photo attempts / 314,517 frozen machine-classified hues**. The signed acquisition workers have newly completed source photo groups since that snapshot.

Verified with full paginated GitHub Actions Jobs and matching stored worker artifacts:

- p0004 original source run 39353343403: **18/18** complete (already in published checkpoint).
- p0008 original source run 39371199127: **23/36** complete, four workers active, nine queued, zero failed; three more complete groups compared with the preceding checkpoint.
- p0016 original source run 39374811084: **10/73** complete, two workers active, 61 queued, zero failed; two more complete groups.
- p0032 original source run 39401484641: **2/145** complete, two workers active, 141 queued, zero failed; one more complete group.
- The 569-photo independent low-resolution shadow cohort is complete and already incorporated; automatic 320px exclusion remains disabled.

This source-backed update requests exactly one **data-only** checkpoint of the newly sealed original-image terminals. Authenticate every completed source shard, original photo/observation ID and frozen ROI terminal method/checksum. Do not download or infer any new image. Do not count candidate image targets or queued unit numbers as measured hues. Original technical failures and ambiguous colours must remain unchanged.

The source queues p0008, p0016 and p0032 continue to acquire photographs independently. This update does not modify or restart those jobs. The new release's verified photo and hue counts must be read from the production collector's immutable receipt.


## New sealed fifth-wave units — 2026-10-20 evening

The latest completed source-authenticated five-wave checkpoint [39413597007](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39413597007-v1) contains **667,011 uniquely attempted original photo IDs** (including acquisition errors) and **316,897 four-hue automatically classifiable original image outcomes**.

The next paginated GitHub Actions source check finds these newly sealed original-image workers, compared with that immutable checkpoint:

- p0004 source run 39353343403: 18/18 already completed, unchanged.
- p0008 source run 39371199127: 24/36 completed groups, 4 running and 8 waiting; the preceding checkpoint used 23 groups.
- p0016 source run 39374811084: 10/73 completed groups, 2 running and 61 waiting; unchanged from the preceding checkpoint.
- p0032 source run 39401484641: 4/145 completed groups, 2 running and 139 waiting; the preceding checkpoint used 2 groups.
- Independent 569-photo shadow comparison: complete, unchanged.

**Do not assume that every group contains exactly 1024 classified original photos.** Trigger the already validated data-only five-wave recovery once. It must verify exact source-shard SHA, photo/observation IDs, original frozen ROI terminal hashes and true acquisition failures, then publish the new classified and unclassifiable counts. Neither download images nor rerun segmentation.

This is a one-time source-provenance milestone of actual sealed worker outputs. The independent three-hour watcher remains enabled for subsequent new files, and the in-flight p0008/p0016/p0032 acquisition workers are not modified or cancelled.

## Six-arm source-authenticated pixel recovery request, 20 October 2026

Last confirmed immutable original photo checkpoint: [39418633371](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39418633371-v1), with **671,107 distinct original image terminal IDs** and **318,460 classifiable four-hue automated labels**. This snapshot already had p0004 18/18 sealed groups, p0008 24/36, p0016 11/73, p0032 4/145, and the full shadow 569-photo cohort; geo05 was at 0 sealed RGB groups.

A fresh paginated GitHub Actions source check shows **p0008 27/36** sealed worker groups, **p0016 12/73**, **p0032 4/145**, **geo05 0/22**, with p0004 18/18 and shadow complete. There are no failed acquisition workers. In-progress original-image ROI workers continue independently. The new sealed p0008 and p0016 original terminal CSVs have not yet been included in the last publication.

Trigger the preexisting **source-SHA-validated, data-only** six-arm recovery once via this authenticated main-branch milestone. The collector must check every original source CSV SHA, original photo ID/observation identity, image terminal receipt/ROI model signature, preserve technical unknowns and failed acquisitions, and exclude every still-pending source target from countable measured outcomes. Neither re-download original photos nor rerun flower-colour segmentation, and do not change the ongoing image workers. **Read new photo/hue totals only from the newly sealed checkpoint receipt**, not guessed from group sizes. The geo05 enrichment probability/design labels remain distinct from the random probability-arm estimand.

The three-hour new-artifact watcher is also configured and will skip redundant checkpoints when no source-validated worker artifacts were added. This milestone merely advances the ongoing acquisition and verified integration.

## Additional six-arm real-pixel snapshot 2026-10-20-21-33-JST

Last public checkpoint [mois-real-pixel-checkpoint-39421820962-v1](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39421820962-v1) passed source/photo SHA, exact original-photo ID and full ROI terminal validation and preserves per-original-photo selection probability. The receipt confirms **676,227 unique attempted original photo IDs** and **320,352 four-hue classifiable original photos**, with p0008 27 completed groups, p0016 12, p0032 5, and geo05 0.

New photo pixel worker archives were checked directly from their *original* GitHub Actions run IDs: p0008 run 39371199127 now **28/36 completed RGB groups**, and independently sampled tropical/southern geo05 run 39417990675 now **1/22 completed RGB group**, with no failed measurement jobs. The other source snapshots are unchanged. All remaining groups are either active or queued. Do **not** infer actual classified-image counts from group totals.

This one-off `main` checkpoint request is **data-only**: reuse frozen source SHA files and sealed terminal outcomes from original worker artifact IDs, validate original observation+photo identity, preserve classifiability and acquisition-failure missingness and geo05's per-observation heterogeneous selection probabilities. Do not download original images or re-infer full-res ROI. Do not mix the geo05 coverage-enrichment arm into the original uniform probability sample as if it were representative. Only production validated terminal rows are measured. The 3-hour watcher remains enabled for subsequent newly sealed work units.


## Recovery after completed p0008 — 2026-10-21 early JST

Exact previous released source-authenticated checkpoint: [39440233813](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39440233813-v1), **694,170** distinct genuinely attempted original photo IDs with **327,144** classifiable four-hue labels. Complete p0004 (18,119 photos) and complete p0008 (36,375 photos) are both published source-SHA immutable full photo cohorts. The checkpoint contained **p0016 15/73**, **p0032 9/145**, and **geo05 2/22** sealed image worker groups, with remaining groups pending and no losses imputed.

New GitHub Actions RGB worker files have since been sealed: p0016 **17/73** (2 more), p0032 **10/145** (1 more), and separate tropical/southern geo05 **3/22** (1 more). The original workers stay active, no failed worker jobs observed. Main-branch event validation and source run IDs are frozen: p0004 39353343403, p0008 39371199127, p0016 39374811084, p0032 39401484641, geo05 39417990675, shadow 39369517433.

**Action:** run the existing data-only six-arm checkpoint collector once. Reuse the exact SHA-authenticated original source observation/photo IDs and genuine frozen ROI-v4 terminal colours, not the photo candidate counts or inferred rose/white hues; do not reacquire pixels or redo the classifier. Preserve `source_arm` and original observation-level `selection_probability` in cumulative Parquet. The geographic enrichment arm remains separate from the probability sample. Publish immutable checkpoint counts only after terminal and SHA reconciliation passes. The watcher stays enabled to skip unchanged output on later checks.


## 2026-10-21 additional completed p0016/p0032/geo05 original RGB

The last immutable verified six-arm true original-photo RGB checkpoint is [39446764842](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39446764842-v1), with **699,290 distinct actually attempted original photo IDs** and **329,072 frozen four-hue classifications**. Two probability waves (p0004, p0008) are complete; the 93,209,351 upstream photo IDs remain source metadata, not colour pixels.

New source-verified sealed worker artifacts from the live acquisition run IDs since that checkpoint:

- p0016 original tier03 run 39374811084: **18/73 completed RGB groups**, 2 active, 53 queued, 0 failed; last checkpoint had 17 groups.
- p0032 original tier04 run 39401484641: **11/145 completed RGB groups** (the first 100 GitHub Jobs response omits some matrix groups; full source plan has 145), 2 active, 132 queued, 0 failed; last checkpoint had 9 groups.
- geo05 tropical/southern enrichment run 39417990675: **4/22 completed RGB groups**, 1 active, 17 queued, 0 failed; last checkpoint had 3 groups.
- p0004 18/18 complete, p0008 36/36 complete, and the tier03 569-photo shadow are immutable and not to be reacquired.

Trigger the unchanged **six-arm data-only** recovery once. Validate original source/photo/observation identity and full-resolution ROI terminal SHA/method signatures, then add only completed measurements to the 612,483-photo immutable baseline. Preserve all failed acquisitions/unknowns, original unmodified image hue labels, original `source_arm` and **per-original-observation `selection_probability`**, including differing geo05 regional rates. Never assume all 4 newly completed groups contain classifiable flowers. Do not rerun model inference or re-download original images.

The geo05 arm is an independently selected geographic coverage-enrichment source, not a probability-complete globally representative species sample. It must remain separately labelled in any ecological analysis.


## Newly sealed original image units after automated checkpoint 39513292880 — 2026-10-21

Latest successful **automated** photo terminal release is [39513292880](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39513292880-v1), containing **737,178 unique true original-photo attempts** (including technical acquisition failures) and **343,249 genuine frozen automatic four-hue ROI classifications**. The locally inspected ZIP, `checkpoint_receipt.json`, and all output SHA256 hashes agree. Its exact signed worker snapshots included p0016 **31 groups**, p0032 **25 groups**, and geo05 **12 groups**; p0004 **18** and p0008 **36** groups are completed full source waves.

A new complete-workflow GitHub Jobs/artifacts check shows **p0016 33 sealed groups**, **p0032 26 sealed groups**, and **geo05 13 sealed groups**, with 0 failed jobs. Thus additional *previously uncollected actual source-image terminal receipts* have arrived. Existing p0016, p0032 and geo05 full-resolution ROI worker jobs continue, with no attempt to restart or redownload their pictures.

**Request exactly one independent data-only six-arm checkpoint**, using only the pinned original source and worker artifact SHA256/terminal receipts. Preserve the exact original photo and observation IDs, immutable four-hue labels, original observation-level selection probabilities, source-arm provenance and technical unknowns. The geo05 latitude-enrichment arm is distinct and cannot be pooled unweighted as if it had the global uniform probability design. Never count the 93,209,351 original metadata catalogue as 93m measured image colours.

The checkpoint's new attempted/photo-hue counts must be computed from authenticated original terminal rows, *never projected from workers or photo targets*. Do not cancel, duplicate or adjust the active acquisition waves. The hourly conditional watcher remains operational for future worker completions.

## Verified sealed-pixel refresh — 2026-10-21 10:32 JST

The immutable latest independently checked [checkpoint 39518136040](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39518136040-v1) was produced with **742,298 distinct genuine original-image terminal photo IDs** and **345,238 original full ROI-v4 four-hue classifiable photographs**. Its signed six-wave photo counts are: original baseline 612,483; full p0004 18,119; full p0008 36,375; p0016 33,792 (33 sealed groups); p0032 27,648 (27 sealed groups); geo05 13,312 (13 sealed groups); original shadow 569. The verified JSON receipt and ZIP digests match published assets. Candidate photo targets are not counted as measured images.

After this publication, the exact original workers produced further sealed RGB terminal artifacts in their pinned GitHub Actions acquisition runs, with no failed worker jobs at this check:

- p0016 source run 39374811084: **35/73** sealed successful source-colour worker groups (2 new; 2 actively running).
- p0032 source run 39401484641: **29/145** sealed groups (2 new; 2 actively running).
- separate tropical/southern geo05 source run 39417990675: **14/22** sealed groups (1 new; 1 actively running).
- Full p0004 18/18 and full p0008 36/36 are already measured and released; 569-photo lowres shadow validation is also complete.

Perform exactly one source-SHA authenticated **data-only** six-arm checkpoint on `main`, collecting only genuine sealed photo/ROI terminal files and deduplicating by original photo ID. Preserve original per-observation `selection_probability`, original source arm, acquisition failures and non-classifiable ROI outcomes, and keep the independently enriched geo05 arm distinct from the uniform observation-hash probability arm. Do **not** redownload photos, rerun frozen YOLO11/EfficientSAM colour model, reinterpret the original 93,209,351 source-photo targets as classified results, or extrapolate measured counts from completed computational groups. Publish only the exact real attempted-photo and hue counts that pass terminal/source/manifest SHA verification. Existing image-acquisition workers remain untouched. The hourly watcher is separately configured to skip checkpoints with no new signed worker files.
