# MOIS original-photo RGB: verified 753,562 terminals and batch recovery

Verified on **2026-10-21 (JST)**; source-of-truth is immutable [actual photo checkpoint 39534159275](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39534159275-v1). The independent local audit opened the original GitHub Actions ZIP, checked its SHA-256 and all three output checksums against the receipt, re-summed per-wave and four-hue CSV values, and checked exact source wave roster accounting.

## Genuine original-photo results

| Source photo wave | Terminal photo attempts (including download failures) | Frozen 4-hue machine classifications | Still pending |
| --- | ---: | ---: | ---: |
| Preexisting FCP/MOIS/p0002 original pixels (immutable) | 612,483 | 295,801 | n/a |
| p0004 (all 18 source-signed groups) | 18,119 | 7,085 | **0** |
| p0008 (all 36 groups) | 36,375 | 13,874 | **0** |
| p0016 (38/73 groups) | 38,912 | 15,003 | **34,960** |
| p0032 (31/145 groups) | 31,744 | 12,258 | **116,327** |
| Distinct regional tropical/southern geo05 (15/22 groups) | 15,360 | 5,386 | **6,918** |
| Complete low-resolution screen shadow (excluded from repeated sampling) | 569 | 199 | **0** |
| **All distinct genuinely terminal photo IDs** | **753,562** | **349,606** | **158,205** |

The preceding verified checkpoint [39532247466](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39532247466-v1) had 751,514 terminal photo IDs and 348,794 classifiable hues. This checkpoint added exactly **2,048 actual photo attempts**: 1,024 p0016 photo attempts / 401 four-hue classifications and 1,024 p0032 attempts / 411 classifications. The increment also includes two image-acquisition failures and preserves 1,234 acquired-but-unclassifiable image outcomes. The latest source-pinned workers remained active and had no failed jobs at verification.

The complete source frame of **93,209,351** photo IDs is only original photo/observation metadata and must not be described as 93m photographs whose flower colour was actually measured. The machine four-hue labels are petal-ROI candidate outcomes, not independent certified flowering plant observations. Partially completed worker groups are not complete probability samples and must not be used for global white-flower species-frequency estimates.

### Descriptive, photo-conditional hue composition

For the **349,606 classifiable photo ROI outcomes**, the four labels were **white 167,030**, **yellow-orange 110,709**, **red-pink 36,564**, **blue-purple 35,303**. White was about **47.78% of classifiable photographs**, not 47.78% of globally native flowering species.

Among the separately selected tropical/southern geo05 arm, **15,360 photos were genuinely attempted** and 5,386 had a classified ROI hue (35.07% of attempts). The automatically labelled white fraction among these classifiable regional photos was **2,587 / 5,386 = 48.03%**. Geographic selection probabilities vary by observation and are retained in both source data and the compact checkpoint Parquet; do not combine this enriched arm unweighted with the uniform probability source. The geo05 sample retains unknowns and downloaded-but-unclassifiable files as missing original image outcomes, not coloured or nonwhite plants.

### Source authenticity

- GitHub Actions checkpoint: [39534159275](https://github.com/zuizui0223/mois/actions/runs/39534159275), all validation, freshness and collector jobs **success**.
- Immutable release: [mois-real-pixel-checkpoint-39534159275-v1](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39534159275-v1).
- Source artifact **12316473051**, ZIP byte size **51,370,798** and ZIP SHA-256 **`dee78fe3e84d4862a4651d35fc855583b47be91dcca21e1e63be8480994a6c67`**.
- All three output SHA-256 hashes in `checkpoint_receipt.json` independently matched their ZIP members: combined actual original photo Parquet, newly completed original photo ROI Parquet, and photo-hue counts CSV. Independent CSV arithmetic also reproduced the 753,562/349,606 terminal/known counts and each wave's counters.
- The production collector checked every source-shard photo-ID roster SHA, original photo/observation identities and frozen original YOLO11/EfficientSAM ROI-v4 method, without image redownload or hue reinference. Local archive validation was SHA/CSV/receipt arithmetic; it was **not** an independent re-read of every Parquet row.

## Exact artifact download efficiency improvement

The previously deployed checkpoint downloaded **each of 100+ already-finished worker artifacts with a separate GitHub CLI process**, repeatedly fetching artifact metadata. [PR #142](https://github.com/zuizui0223/mois/pull/142) replaced this with **one exact-name multi-artifact `gh run download` per pinned acquisition run**, passing explicit repeated `--name` flags after freezing the original worker artifact ID/size/SHA snapshot. It never uses a broad glob that could include worker groups arriving mid-checkpoint. For multi-artifact responses the extractor validates that the produced artifact directory names match exactly the frozen inventory. Existing source/photo/ROI terminal checks are unchanged. The optional single-artifact direct-extraction behavior is accounted for.

The PR CI (checkpoint contract tests, model invariants, temporal and seasonal checks) passed. [PR #143](https://github.com/zuizui0223/mois/pull/143) subsequently confirmed the **production GitHub CLI runtime** can download p0004 exact RGB worker groups 0 and 1 in a single call, extract both with the expected directory structure, and verify **8 original shard ROI receipts / 2,048 attempts and all photos.csv SHA-256 checks**. Both pull-request and main push smoke workflows passed. **No new image classification** occurred during this optimization. The batched collector itself was merged after the latest full checkpoint, and its full production run on the next new sealed worker cohort remains to be verified; do not prematurely claim a benchmarked production speedup.

The existing three active original acquisition source runs (p0016/p0032/geo05) continue. The hourly source-aware watcher will trigger the unchanged SHA-validating data-only workflow only if additional actual completed worker archives are published, and skips unchanged snapshots. No redundant RGB inference or source-93m reread is required.
