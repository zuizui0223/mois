# MOIS 128シャード大量取得 — 実行・監視記録（2026-10-26 JST）

## Current verified execution

- [Source-only 131,072 original-photo IDs release](https://github.com/zuizui0223/mois/releases/tag/mois-full93m-bulk128-tier10-131072-original-photo-source-v1): all **128 immutable source photo partitions × 1,024** original photo attempts; 512 signed up-to-256 photo input units; 9,543,472 authenticated source-photo opportunities scanned from frozen **tier10**, not re-enumerated 93,209,351 original source photo IDs.
- Reference arm **32,768 source-photo targets**; deliberately colour-gap-targeted arm **98,304 source-photo targets**; zero P5-targeted images. Exact original photo and observation UUID duplicates = zero at source selection; source SHA archive retained. The reference is conditional on the source-selected observation tier and original opportunity frame; neither arm is a random representative census of the world's wild native flora.
- [Bulk128 source-selection success run 40338267506](https://github.com/zuizui0223/mois/actions/runs/40338267506). [Idempotent already-published source gate success run 40340510247](https://github.com/zuizui0223/mois/actions/runs/40340510247). Source manifest key mismatch on earlier repeated gate runs was corrected; no new source selection was needed.
- **Two different RGB image workers currently share the same 131,072-photo source roster:**
  - [Run 40341341397](https://github.com/zuizui0223/mois/actions/runs/40341341397): fast 500px medium RGB, source-frozen YOLO11 flower candidate boxes, sRGB-derived CIELAB box descriptors, technical missingness, **128 shards with max 12 active jobs**. Image workers are observed starting; there is **no authenticated terminal full-cohort 131,072 image pixel success receipt yet**. This produces *box-level image colour*, not anatomically validated petal colour or strict nine-palette FCP ROI-v4 results.
  - [Run 40340338433](https://github.com/zuizui0223/mois/actions/runs/40340338433): parallel 128-shard full-resolution original FCP-v4 YOLO11 + EfficientSAM ROI9, max 16 concurrent image jobs. It uses the same original image IDs but full-size image variants, and is therefore complementary in measurement type **but duplicates provider requests per photograph**. Both runs' source-gate jobs succeeded, image workers were in progress, no end-to-end completion confirmed at this snapshot.
- **Important efficiency risk:** do not dispatch another 128-shard pixel scan simply because the workflow-run API shows `queued`; inspect **job-level states** first. GitHub Actions currently exposes completed source gates plus active image jobs while the outer run can still be shown as queued. Re-running would duplicate exact photos. Avoid speculative release or worker retry. Preserve SHA-checked prior partial and complete packets before recovery.

## Recovered strictly measured source-taxon breadth

- [Original P2 source-only 4,096 unique original photos / observations / previously colour-unmeasured P2 taxa](https://github.com/zuizui0223/mois/releases/tag/mois-P2-first4096-tier07-one-taxon-one-photo-source-20261026-v1). All eight 512-photo original ROI-v4 image workers completed actual pixels in original run #40329836782, but **GitHub Release write permission** blocked publication, not colour inference.
- [DATA-ONLY recovery run #40342167711](https://github.com/zuizui0223/mois/actions/runs/40342167711) succeeded after verifying and reusing eight SHA-pinned Actions artifacts. Absolutely **no photo re-download or new original FCP inference** in recovery.
- [P2 genuine ROI9 output release](https://github.com/zuizui0223/mois/releases/tag/mois-P2-first4096-actual-originalROI9-newtaxa-20261026-v1): **4,096 actual original image attempts**, **1,467 strict nine-component original-photo ROI9 outcomes**, and **1,467 previously nine-colour-unobserved source taxa newly covered**. 16 image acquisition failures; other technical unknowns remain missing. Before: 1,045,978 photo attempts, 401,392 strict ROI9, 21,320 source taxa with strict9. After: **402,859 strict ROI9** and **22,787 of 42,111 source nominal taxa** covered, with **19,324 nominal source taxa still without a strict9 observation**. Do not combine the separate 64 image *500px Lab* pilot into strict9 counters.

## Operational monitoring contract

Hourly ChatGPT condition watch, separate from this repository, monitors **actual job statuses**, completed artifact SHA and real pixel counts, meaningful 25%/50%/75%/100% milestones, rate-limiting / 403 / 429 / 503, source identity, technical unknowns and final measured releases. Do not emit repeated notices for harmless queued jobs, CI reruns, or already known results. The watch does not itself download photos or trigger compute, and it may depend on the user's ChatGPT notification settings for proactive delivery.

Critical next actions:
1. Use [128-shard 500px run](https://github.com/zuizui0223/mois/actions/runs/40341341397) and [128-shard full-resolution run](https://github.com/zuizui0223/mois/actions/runs/40340338433) *job-level* statuses to verify incremental worker completion; beware duplication of original photo/provider requests across these two parallel measurement methods.
2. If a part fails, diagnose the specific photo source and failure first; restore SHA-sealed completed packets from the first run instead of redownloading/reinferring every photo. Current 500px workflow supports signed packet resume, but recovery must be provenance-validated and not label unknown as a flowering biological negative.
3. On terminal completion, separately verify (a) real image download/decode attempts and YOLO-box Lab missingness, (b) original strict ROI9 photo counts, (c) source taxon×realm×local season support, keeping original photo ID and observation UUID distinct and reporting the two acquisition arms separately.
4. Do not claim **131,072 measured photo colours** merely because **131,072 source targets** were selected. No 9-palette ROI success count can be inferred from 500px YOLO box colours. No measured new ecological association is established by source-selection counts alone.

This document records the observed source / Actions status at the listed snapshot, not a real-time dashboard or a claim of continuing computation within ChatGPT.


## 2026-10-27 JST — terminal **actual 131,072-photo two-pixel-method integration**

The preceding "in progress" bullet points above describe the original 26 October snapshot; they are historical, **not current job states**.

**Both photo-processing methods are genuinely finished on all 128 photo shards** and have 128 SHA-auditable GitHub Actions image output artifacts each:

1. [500px source-original image detection/Lab job 40341341397](https://github.com/zuizui0223/mois/actions/runs/40341341397): all **128 medium-photo-shards succeeded**, source gate succeeded. Original "combine" failed at script import with **`ModuleNotFoundError: No module named 'numpy'`**; no image measurement failure.
2. [Full original pixel FCP-v4 YOLO11/EfficientSAM job 40340338433](https://github.com/zuizui0223/mois/actions/runs/40340338433): all **128 image inference steps succeeded**, then all 128 workflow worker jobs were marked failed solely because **GitHub Release creation returned 403 Resource not accessible by integration**. Their numerical terminal outputs were nevertheless uploaded as all 128 real-image Actions artifacts.
3. **[Data-only unified recovery 40431588179](https://github.com/zuizui0223/mois/actions/runs/40431588179): SUCCESS**; verifies the two source runs and 128+128 authenticated worker artifacts, all 512 FCP work-unit terminal SHAs, then joins results on the exact original **photo_id, observation_id, inat_taxon_id**. No original image was downloaded, inferred, recoloured or excluded during this recovery.
4. **[Immutable merged release](https://github.com/zuizui0223/mois/releases/tag/mois-128shards-131072-fullROI9-plus-500px-sourceID-joined-20261027-v1)**: one row per original image ID, two separately preserved image-size/ROI methods, complete regional+seasonal/priority QC and explicit missingness. Full source-photo identity and image outcome provenance, audit tests and SHA sidecars are retained.

### Real picture outcomes (these are NOT planned target counts)

| Original image metric | Real result |
|---|---:|
| Signed original image + observation IDs, each represented once in the merged lake | **131,072** |
| Images actually decoded at medium500 | **130,727** |
| Medium500 high-recall flower candidate boxes | **127,600** |
| FCP-v4 actual original-image strict nine-palette ROIs | **48,784** |
| Original-size strict nine-palette outcome unknown/technically missing | **82,288** |
| Medium500 negative **but** real full original ROI9 valid | **95** |
| Medium500 flower candidate **but** strict original ROI9 missing | **78,911** |

By the *prespecified source-selection arms*:
- **Reference:** 32,768 original images; 32,699 medium RGB decoded; 31,873 flower boxes; **13,243 strict ROI9**.
- **Adaptive gap-targeted:** 98,304 original images; 98,028 medium RGB decoded; 95,727 flower boxes; **35,541 strict ROI9**.
- **P1/P2 taxon coverage, within the *frozen 26 October source-status frame*:** the targeted source arm produced authentic first ROI9 outcomes for **7,194 previously no-colour nominal taxa** (P1 184, P2 7,010). This is a verified targeted cohort taxon-coverage gain, *not* a contemporaneous global-main-branch 42,111-frame uniqueness estimate unless other active FCP/93m source arms have been reconciled too. Reference P1/P2 taxon support must not be naively added to this figure without exact taxon dedup.

### Implications for the million-to-93-million photo scan

- Fast YOLO box detection is a useful **high-recall technical pre-screen**, but **not a petal-colour measurement or perfect exclusion gate**. Even amongst box-positive images, **78,911** failed the strict FCP ROI9; this confirms that heavy ROI/petal segmentation and photographic classifiability are the present scientific/computational obstacle. **95** medium-screen negatives held genuine strict full ROI9 colours. Therefore, do not drop all medium negatives: protect a prespecified negative audit sample, or all negatives until the missingness/false-negative risk has been calibrated independently.
- The two methods were unfortunately **run separately on the same 131,072 original photo IDs**, duplicating provider requests (500px/large variants). **Integration is by photo ID, NOT 262,144 independent photos.** Future scalability should use one image acquisition per photo (image decode cached within a worker), do cheap flower-box +Lab first and selectively full ROI where biologically justified, while retaining negative technical audits and original sampling probabilities.
- The 500px box CIELAB is a colour descriptor of detected floral *boxes* and can include leaf/background; the nine-component real original FCP is an operational unvalidated floral-ROI estimate, not anatomically confirmed corolla pigment or spectral reflectance. Unknown images must not become white/nonwhite by proxy.
- The baseline [previous strict 4096 P2 original pictures: 1,467 first-colour taxa](https://github.com/zuizui0223/mois/releases/tag/mois-P2-first4096-actual-originalROI9-newtaxa-20261026-v1) is a separate **prior source arm** and must be counted once. Other active main-branch waves can have changed the latest global deduplicated count. Avoid declaring **48,784** or **7,194** to be new global net counts without cross-arm source-photo and taxon-ID join.
- A scan/worker being marked **failed** does not mean no valid photo pixels: use original SHA-verifiable Actions artifacts and final data-only receipts. Conversely, successfully frozen **source targets** do not prove actual image downloads. Monitoring should report separately the `128/128 physical original jobs finished`, `130,727 medium images decoded`, `48,784 strict original ROI9`, `7,194 distinct previously uncoloured targeted source taxa`, and `131,072 exact IDs joined`.

**Next atomic scientific pipeline integration:** update the current cross-arm global photo-ID union on main from the immutable joined parquet, checking the current accepted FCP-K1/main-wave photo and observation ledgers BEFORE promoting global cumulative counts. The 500px screen does not replace the original strict 9-palette outcome.
