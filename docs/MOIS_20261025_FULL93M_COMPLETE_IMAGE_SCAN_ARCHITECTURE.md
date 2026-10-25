# MOIS full 93,209,351 image-pixel first pass (2026-10-25)

## Scope: scan ALL photographs, do not sample away the universe

Decision: MOIS will retain all 93,209,351 original iNaturalist Open Data photo IDs
(53,652,241 source observations, 42,111 nominal source taxa) for an exhaustive
**medium-image pixel scan**, not merely enumerate source metadata. This is a new
acquisition architecture, NOT evidence that 93m photographs have already been
downloaded or coloured.

Current 2026-10-25 checkpoint (historical + K1 + FCP5 completed):
1,019,272 genuine image attempts and 449,741 original four-hue ROI classifications.
The K1 ROI/flip-gate nonresponse (46,288 / 76,753) is a reason to replace
the full-resolution-all-images workflow, NOT to limit the source to five
photographs per taxon.

**One-time source work already done:** the 128 original partitions were split
into 14 disjoint original-observation-hash tiers, with source SHA receipts.
Their photo-ID total is exactly 93,209,351 and there are no per-species,
per-genus, native/introduced, geographic or monthly photo filters. The new
`mois_full93m_every_photo_scan_manifest.py` reconciles those signed receipts
into 1,792 exhaustive tier×part work items. That step opens ZERO images.
No attempt to re-enumerate the 93m original observation/photo catalogue.

## The two levels of *actual* image work

- **Stage 1, every source photo:** fetch its iNaturalist Open Data
  `medium` (max 500px) image from its validated photo ID and extension,
  decode it, record the actual image-byte SHA256 and technical failures.
  An opt-in, SHA-pinned YOLO11 flower-box detector runs **warm** in 32-image
  GPU/CPU batches, with concurrent paced image transfer and one model load per
  running shard worker. Store a screen candidate/negative/uncertain/error,
  detection confidence and rough visible HSV composition inside the detected
  bounding box. These box colours are explicitly **not petal measurements**.
  Retain original taxon, observation, day, coordinates, license and image ID.
- **Stage 2, full-resolution selectively:** fetch 1024px for candidates,
  ambiguous results and a **positive-probability validation sample of
  screen-negatives** after independent botanical reviews establish adequate
  flower and white-flower recall in every relevant geography/taxon/season
  stratum. Run stronger flower-organ segmentation / 9-palette and CIELAB.
  Reuse source-authenticated existing large-image ROI outcomes without
  redownloading. No automatic promotion of Stage-1 box proxy to source four
  original hue classes. For now, retain every screen-negative in the queue
  rather than dropping it: the older 320px 0.05 confidence shadow pilot
  missed 7/199 original classifiable photos including 5/93 white labels.

**Existing protection:** `mois_93m_stage2_probability_queue.py` defaults to
SHADOW mode (all photos go to full ROI), and will only authorize negative
subsampling with a matching source+policy SHA and independent photo-specific
flower review, >=100 audited negatives per stratum and a white-flower
screen-recall one-sided 95% lower bound >=0.95. No current 500px screen has
passed this gate. An image with no YOLO candidate is not an ecological
nonflowering observation; camera angle, petals, inflorescence, leaves and fruit
still require human/future-model assessment.

## Executable, source-frozen code now on the experimental branch

Scripts:

- `scripts/acquisition/mois_full93m_every_photo_scan_manifest.py`:
  ingest the signed full-source registry and signed 16 group receipts,
  enforce exact 128×14 complete photo coverage, write
  `all_original_photos_stage1_jobs.csv` with the exact source gzip SHA,
  expected photo count and p=1 within the 93m photographic *source*.
- `scripts/acquisition/mois_full93m_medium_pixel_scan.py`:
  one actual `tier-XX/part-YYY.csv.gz` at a time. SHA-verify the compressed
  source, read each original image from S3 (not the throttled iNaturalist API),
  decode genuine pixels, run optional frozen warm YOLO detection, keep
  all original source photo IDs including errors, write atomic
  `packet-0000000.csv.gz` + receipt files, resume only matching packet SHAs.
  It cannot complete an input shard unless its actual photo attempts equal the
  signed source count. `--max-packets` creates a **pilot only** and never a
  completed census. A full CLI scan requires a local SHA-verified flower detector.
- `tests/test_mois_full93m_medium_pixel_scan.py`: offline synthetic image
  decode, model-injection, download failure, unsupported format, source-sha
  tamper, idempotent resume, incomplete-pilot and exhaustive-manifest checks.

Example of generating the **complete work plan** after downloading existing
source receipts (no image fetching):

```bash
python -m scripts.acquisition.mois_full93m_every_photo_scan_manifest \
  --registry SOURCE/mois_all_original_photo_target_tiers_receipt.json \
  --group-receipt SOURCE/group-??-receipt.json \
  --output PLAN
```

The source release stores 224 tier-group TAR archives (14 tiers × 16 groups),
not 1,792 individually released gzip files. Extract **only the selected
signed** `tier-XX-group-YY.tar` archives to a trusted TARGETS directory
before starting those jobs; never re-enumerate the underlying full source.
The existing per-group receipt SHA verifies each extracted gzip input.

Each job row identifies an existing original tier source gzip with
`relative_source_gzip_path`, `source_gzip_sha256`, `expected_photo_rows`,
`part`, `tier`. Supply those fields to the stage-1 executor:

```bash
python -m scripts.acquisition.mois_full93m_medium_pixel_scan \
  --source TARGETS/tier-00/part-000.csv.gz \
  --source-sha256 SIGNED_SHA256_FROM_PLAN \
  --expected-count SIGNED_COUNT_FROM_PLAN \
  --part 0 --tier 0 \
  --weights /trusted/models/jrc_yolo11n_last_v4.pt \
  --weights-sha256 VERIFIED_MODEL_SHA256 \
  --flower-class-id 0 --flower-class-name item \
  --device 0 --workers 8 --batch-size 512 \
  --requests-per-second 2 \
  --output APPEND_ONLY/p000-t00
```

The SHA-pinned original FCP ROI-v4 YOLO11n weight has class ID 0 named **`item`** (verified in real GitHub Action logs). It is a flower-instance detection model, but the training label is the literal string `item`. Models with other training labels must use their actual metadata; the class name is never guessed.
No bundled weights and no automatic model download.

## Performance, source permissions and deployment gate

The prior 70-minute / 1,024-photo full 1024px ROI worker is far too slow
for 93m: five always-active workers would take approximately 885 days
*if their measured throughput stayed unchanged*. This does not imply that
medium-stage throughput is already measured or that a target speedup
has been achieved.

To scan 93,209,351 images in 30 days requires an **aggregate sustained
36.0 image attempts per second**; 60 days requires 18.0/s. If the median
medium transfer were 30 KB, the full pass would transfer ~2.8 TB; if it were
100 KB, ~9.3 TB. These are planning scenarios, **not observed media sizes,
prices or promised runtime**. CPU/GPU warm-model startup, actual average
image bytes, S3 timeouts, service throttling, per-class recall and image
decoder success must be measured on a capped **prospective pilot** before
bulk dispatch.

The image provider's website cautions that large media downloads can cause
blocks (5 GB/hour and 24 GB/day), and iNaturalist API should not be used for
bulk retrieval. The public AWS Open Data S3 dataset is the correct
bulk-source route, but **do not assume the website's or provider's bulk
traffic policies can be bypassed by concurrency or parallel IPs**.
Coordinate S3 bulk access/compute capacity and quota with the provider;
start with low concurrency and an independently measured bounded run.
Run durable warm workers near the S3 dataset (us-east-1) with centralized
rate/bandwidth budgets, not tens of thousands of cold GitHub Actions
jobs. Preserve per-photo Creative Commons/CC0 terms and attribution.

No all-93m remote job, cloud deployment, expensive model acquisition or
image fetching has been launched by these code-only commits. Existing FCP5
workers remain unaffected. The next experimental execution is a 10,000-photo
region/taxon/season-stratified **actual medium-image benchmark**, comparing
pixel decode rate, per-worker throughput, photo-specific stage-1 recall
against frozen 1024px full ROI, the stage-1 color proxy's error and 500px
YOLO boxes vs human open-flower annotations. The 93m complete scan should
then run as resumed, capped, SHA-reconciled original shards.

## Output interpretation after all shards finish

A sum of exactly 93,209,351 **source photo-ID attempts** is necessary but
not sufficient: separately report successfully decoded medium photos,
recoverable image access/format errors, model screen-negative counts,
flower-candidate counts and independent petal-validated quality. The full
93m are iNaturalist *photo opportunities*, not a global census of blooming
plants. WCVP native status and environmental/geographic covariates are
later **annotations**, never a reason to delete source photo IDs.
