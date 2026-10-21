# MOIS: 93,209,351-photo scale strategy (2026-10-21)

## Operational status

The 2026-09-27 original source contains 93,209,351 unique photo IDs representing 42,111 source taxa. The signed 128 source partitions have already been transformed once into 14 disjoint observation-hash target tiers. These are metadata opportunities, not downloaded images or flower classifications. Preserve all photo IDs and native/introduced/unknown status rather than performing a native-only hard filter.

The independently verified [RGB checkpoint 39549455436](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39549455436-v1) contains 766,874 terminal actual original-photo ID outcomes / 354,560 frozen four-hue classifications. The p0004 and p0008 original image waves are complete; p0016, p0032 and tropical/southern-enrichment geo05 remain active. Geo05's region-dependent original observation propensities must remain separate from the uniform observation-hash probability arm.

## Actual capacity and decision

GitHub Actions successful original-image group-job medians on 2026-10-21: p0016 **71.15 min** over 44 workers; p0032 **69.77 min** over 37; geo05 **68.43 min** over 18. A full work group ordinarily attempts 4 × 256 = 1,024 original images. Provisioning, model bootstrap, downloaded images, original YOLO11 + EfficientSAM ROI-v4 and results are included in these job wall times. Some short groups have fewer images; this is an approximate engineering estimate.

At 70 min / 1,024 photos per continuously occupied worker, five workers might process ~105,000 photo attempts/day, requiring **~885 days** for 93,209,351 photos even before interruptions, throttling and costs. Hypothetical 20 workers ~221 days; 100 workers ~44 days. These are capacity scenarios, not guarantees of linear scaling or permission to increase parallel provider calls.

The tested 320px confidence screen must **not** be used as an exclusion filter: actual comparison of 567 acquired photos at threshold 0.05 would lose 7 of 199 full-ROI classified photos, 5 of 93 white-labelled photographs, to save only 5.20% projected ROI-call time. The screen is not independent botanical ground truth and an observation without Flowering annotation is not proof its photo contains no open petals.

## Next full-population, staged-image architecture

1. **Full 93m target register once.** Original photo ID is the immutable deduplication key. Preserve original observation, taxon, observer, date, latitude, licence and hash tier. Add WCVP-native/introduced annotations by botanical country with uncertainty rather than deleting images. Keep separate probability sample and extra tropical/southern/island coverage arms; never compute a naive unweighted worldwide rate from a mixed selected-photo pool.
2. **Cheap image-stage screen with strong independent botanical validation.** Record candidate/negative/uncertain/error separately, model/source SHA, and original image identity. Run a fresh flower-stage classifier only after benchmarking; the failed 320px model is not approved for production exclusions.
3. **Frozen full-resolution ROI only on a controlled second stage.** Run all screen positives, unknowns and errors plus a positive-probability randomized audit of every negative-screen stratum. Original observation IDs determine reproducible random draws. Store exact second-stage probabilities and statuses, including photos not chosen for full ROI; use appropriate joint sampling weights for estimands that truly have a known probability design. In the interim SHADOW mode keeps 100% of images.
4. **Use persistent warm GPU/CPU worker services, not 90,000 short-lived cold model jobs.** Benchmark model startup, S3 pacing/request latency, decode, YOLO detection, EfficientSAM, and archival time separately. Use objectively measured end-to-end throughput and class-/region-specific missingness, not a synthetic claimed acceleration. Respect provider rate limits and photograph licenses.
5. **Partitioned append-only evidence store.** For true 93m scale, store original terminal result Parquet by immutable shard, with SHA receipt and idempotent original-photo ID index. Append newly finished workers without repeatedly downloading and rewriting every prior 93m row. Current completed-pixel checkpoints remain valid for the existing bounded waves.

## What is implemented in this change

The new script scripts/acquisition/mois_93m_stage2_probability_queue.py is a standard-library-only streaming stage1-to-stage2 allocation tool. It requires a signed stage1 CSV SHA, sorted unique original photo IDs, original observation IDs, source arm, source selection probability, audit stratum, and a prespecified JSON table of **nonzero** audit rates. It records all photographed source decisions, frozen receipt checksums, selected photo IDs, conditional ROI-selection probabilities and conditional/joint design bookkeeping. It never downloads images or assigns hues. A selected-ID file must still be joined to immutable source-photo worker shards before running a new ROI acquisition.

Default shadow mode keeps **all** source IDs. Staged mode refuses to skip screen negatives unless a version/SHA-pinned receipt of independent photo-specific open-flower review matches both the source screen and audit policy, includes at least 100 audited negative-screen photos per stratum and has a white-flower screen-recall 95% lower bound >= 0.95. This is a conservative initial safeguard, not evidence that any such production stage1 model now passes. The previous 569-photo 320px comparison does not qualify.

## Immediate next engineering gates

Continue the original p0016/p0032/geo05 acquisitions and collect newly sealed, source-SHA-verifiable original terminal files only. Separately run prospective SHADOW-only botanical flower-stage screening and full frozen ROI on a new observation/photo-disjoint audit cohort. Do not switch to staged filtering until independently botanical reviewed recall and region/hue losses justify it. For a large future batch, replace repeated full-checkpoint rewrites with immutable partitioned Parquet, and benchmark a warm worker with the *same frozen model* before any scale-up.

**Not yet accomplished:** automatic flower/nonflower truth for 93m photographs, an independently validated high-recall cheap flower-stage model, new GPU deployments, or faster live p0016/p0032/geo05 original-photo workers. No scientifically representative native-flora white-flower prevalence can be inferred from the checkpoint alone.
