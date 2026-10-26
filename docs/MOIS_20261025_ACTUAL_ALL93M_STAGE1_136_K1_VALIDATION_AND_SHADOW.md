# MOIS | true source-pixel first pass and Stage2 no-exclusion contract

**25 October 2026 — empirical results, SHA-verified original image-photo identity.**

## Decision

Full upstream frame **93,209,351** source candidate photo IDs across **128 x 14 = 1,792**
disjoint signed Stage1 scan jobs. The complete real 500px image-pixel scan has
**not** finished. First genuinely executed full-source subset: tier00 part000
**136 source photos** (no per-taxon filtering), with signed photo/gzip hash,
real image download and YOLO11 frozen weight. This paper/analysis must NOT
treat 93m source IDs as botanical petal-colour measurements.

## Source execution results, not planned counts

1. **Exact full 93m source scan plan**:
   successful [GitHub run 40163052435](https://github.com/zuizui0223/mois/actions/runs/40163052435).
   The 16 immutable source group receipts reconcile to 93,209,351
   photos across 1,792 tier×partition jobs. No actual image pixels accessed.
2. **K1 real 256-photo paired calibration**:
   [run 40162920395](https://github.com/zuizui0223/mois/actions/runs/40162920395).
   128 original classified ROI, 128 original ROI/flip-gate unknown.
   Source SHA `0d3de3374a06b1e7b9ad0281e2ab7844632ed139db871488524061191c8f34c1`.
   Genuine iNaturalist medium image bytes 256/256 acquired and decoded;
   **33,809,323 bytes**, **127.913s** end-to-end paced at <=2 requests/s.
3. **Frozen YOLO11-only 500px on the same K1 sample**:
   [run 40163518118](https://github.com/zuizui0223/mois/actions/runs/40163518118).
   Original YOLO11n model SHA `f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5`;
   actual source model class ID 0 is literally **`item`**, not `flower`.
   256/256 decoded, 251 YOLO candidate and 5 no-box; among original ROI
   classified **128/128** candidates, among original ROI/flip failed
   **123/128** candidates. Warm box inference + rough HSV proxy **16.191s**,
   image acquisition+decode **127.543s**, total **143.739s**.
   A box candidate is NOT an authenticated flower or accurate petal ROI.
4. **First actual source-tier pixels**:
   [run 40163902072](https://github.com/zuizui0223/mois/actions/runs/40163902072).
   `tier-00/part-000.csv.gz`, SHA256
   `f4fb52383147cb720367b474d8a755ef1f8181deaee57d86960d802323d82d7e`,
   contains exactly **136** globally unique source photo IDs.
   136/136 original medium images decoded, **130 candidates / 6
   negative-screen** (no assumptions of botanical flower absence), transferred
   **19,827,318 bytes**, acquired+decoded **67.774s**,
   warmed YOLO+rough box hue **5.315s**, total **73.091s**, **1.8607
   photos/sec** paced. These are genuinely image-pixel-scanned outcomes
   within the full 93m source, NOT full-frame photo-colour coverage.
5. **End-to-end Stage1 photo identity -> Stage2 safety gate**:
   successful [run 40164342354](https://github.com/zuizui0223/mois/actions/runs/40164342354).
   `mois_medium_stage1_to_stage2_shadow.py` authenticated and disk-sorted
   all 136 medium-pixel terminal original IDs and used the existing
   `mois_93m_stage2_probability_queue.py` in **SHADOW** mode.
   **136/136 kept for possible full-resolution Stage2** including **all 6**
   negative-screen images. No re-downloaded photos, no imputed 4-hue label
   and no false botanical absence. Binary detector exclusions are NOT
   deployed; the independent negative/photo-specific flowering review gate
   for production remains unsatisfied.

## Exploratory 500px HSV box-proxy vs original 1024px automatic ROI

From precisely the **128 source-classifiable K1 original photos** in the
above balanced 256-photo technical calibration, compare the argmax among
`white,yellow_orange,red_pink,blue_purple` proxy fractions inside the
500px highest-scoring YOLO box against the unchanged frozen original
four-colour ROI label. Source labels: white 61, yellow/orange 42,
red/pink 11, blue/purple 14. Proxy matches **96/128 = 75.0%**, compared to
a majority-white baseline of **61/128 = 47.7%**. Proxy white vs frozen white
binary ROC AUC **0.930** in this convenience calibration, while proxy white
pixel share averages **25.9 percentage points below** the frozen original
ROI white fraction (MAE **27.5 points**; Pearson r **0.786**).

**Inference:** simply switching to cheap YOLO bounding-box HSV loses
considerable original petal-colour information, especially the white-pixel
fraction; original YOLO candidate detection is nearly ubiquitous even among
ROI-failed photographs. Therefore full-stage efficiency requires either
inexpensive verified petal/flower segmentation with a botanically reviewed
subset, or accept per-photo uncertainty; cannot claim cheap YOLO
classifies 93m photos with adequate accuracy. No K1 unknown has been
relabelled. All measures are methodological exploratory relative to a
machine ROI baseline, *not botanical true flower-colour ground truth*.

## Actual current engineering boundary

- The first signed tier shard is fully scanned, and its original photo
  identities remain in 100% Stage2 SHADOW queue. **No real large-resolution
  Stage2 processing was launched by the shadow bridge**.
- The full 1,792-shard original target plan is verified; 1,791 original
  data shards remain pixel-unscanned under this new procedure (ignoring
  historical large-resolution samples, which have not been cross-deduplicated
  with medium outcomes).
- Extrapolating a K1 sample mean of ~132KB across 93m candidate images
  gives **~12.3TB** of 500px transfers, an illustrative capacity scenario,
  not a measured total or guarantee. At the deliberately capped **2 req/s**
  across a single worker, throughput is limited by the cap, not just YOLO.
  Real S3/provider throughput, all-part support, model misses, actual image
  extensions, photo licences and distributed-worker spend must be budgeted,
  benchmarked and centrally throttled before a 93m production crawl.
- **Next technical requirement:** validate a high-recall inexpensive
  **flower-organ** mask on independently annotated K1 failures and
  successful controls. Fast candidate BOX localization alone cannot
  resolve the principal ROI/flip failure rate. Keep screening negatives,
  error attempts and all selection/status denominators explicit.

Scripts, SHA receipts and the automated execution evidence remain on
[PR #184](https://github.com/zuizui0223/mois/pull/184). The GitHub Actions
artifacts retain derived CSV/JSON, **not redistributable raw photo bytes**.
