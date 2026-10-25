# MOIS: 93,209,351-photo persistent scanner (25 October 2026)

## Current operational truth

The original 93,209,351 signed iNaturalist photo IDs were fully registered in
14 disjoint observation-hash tiers × 128 parts; they are PHOTO OPPORTUNITIES,
not 93 million labelled petals. Original medium 500px pixels for several small
source groups have been scanned and signed. Historical original large-image
YOLO11+EfficientSAM photo outcome releases must not be deleted or recomputed.

The previous GitHub Actions approach booted a short-lived CPU worker repeatedly,
fetched at a hard 2–3 image requests/s, ran frozen YOLO box inference on 500px
photos, then discarded the runner. The observed 2–3/s is a **configured pacing
ceiling**, not proof the S3 dataset cannot deliver faster. 93,209,351 photos
at sustained 3/s would need 360 days before outages or data retries. As the
first ~3,172 Tier0 attempts produced ~97% flower-box candidates, dropping
screen-negatives alone will not materially solve the 60% later flower-ROI
failure rate; do not pretend it does.

## Implemented structural changes

1. **Stream original signed source archives exactly once.**
   `mois_93m_prepare_persistent_source.py prepare-group` takes one original
   SHA-pinned tier×8part TAR and group receipt, verifies its S3 photo-ID
   metadata, and streams all eight source parts into bounded 4,096-photo
   chunk-gzip files. No original photograph is fetched. Call it for each of
   the 224 archives; `collect-roots` verifies all 224 metadata group results,
   totals the exact global 93,209,351 photos, and writes portable paths. These
   are compute partitions, not sampling or a native-only hard filter.
2. **Run persistent work instead of an Action per 256/1024 images.**
   `mois_93m_persistent_worker.py` reads any number of signed source chunk
   roots with one immutable YOLO detector per process. It reuses one
   `ThreadPoolExecutor` and bounded `urllib3.PoolManager` (HTTP connection
   keep-alive) across ALL selected 4,096-photo source chunks. Every completed
   512-photo packet retains a source SHA, image SHA and model identity; restart
   never fetches previously completed photos again. Result files contain no
   raw licensed image bytes.
3. **Parallelize without accidentally multiplying provider traffic.**
   Distribute original `(tier,part,chunk)` by a stable hash across
   `--worker-slots N` slots. Each worker receives an allowed
   `--approved-aggregate-rps / N` network budget. One authoritative worker
   per slot is required; cloud orchestration must not launch duplicates.
   S3 403/429/503 stops the worker. Genuine source-photo 404 is preserved as
   a technical missing result, not a nonflower or coloured label.
4. **Store once and reuse.**
   Retain signed Parquet/CSV+receipt outputs on durable user-controlled
   storage, e.g. a mounted disk or owned object store. GitHub Actions
   artifact retention is not a durable home for the full 93m census.
   Build CPU/GPU hosts in the dataset's AWS region where allowed; benchmark
   memory, throughput, retry rates and aggregate cost before allocating
   more cloud resources. Do not launch any large cloud infrastructure
   merely by merging this code.

## Scale targets, not unsupported speedup claims

| Sustained true image attempts/s (all workers) | 93,209,351 photos, uninterrupted |
| ---: | ---: |
| 3 | 360 days |
| 8 | 135 days |
| 40 | 27 days |
| 80 | 13.5 days |
| 120 | 9 days |

These are mathematical throughput requirements, **not observed worker
performance or permitted provider request rates**. A 256-photo K1 500px test
transferred ~33.81 MB (~132 KB per photo), suggesting ~12.3 TB for all photos
if the full 93m population has the same average. That assumption needs
per-region, season and photo-age validation. Costs (CPU, S3/object storage,
network, transfer and instances) are unmeasured; price actual cloud hardware
and data access before committing resources.

The iNaturalist **Open Data S3 bucket** is expressly the bulk licensed-photo
source, not the iNaturalist API/web media origin. The site's warnings about
5 GB/hour or 24 GB/day should not be mechanically interpreted as the AWS
Open Data S3 quota: the iNaturalist Help Center directs bulk photo researchers
to the Open Dataset. There is **no experimentally established S3-wide
aggregate request quota here**, so plan provider-compliant bounded benchmark
steps without falsely equating the API's request/website cap with S3 bandwidth:
https://help.inaturalist.org/en/support/solutions/articles/151000223300-how-i-download-inaturalist-photos-to-use-for-research-
https://www.inaturalist.org/pages/api+recommended+practices
https://registry.opendata.aws/inaturalist-open-data/
Amazon S3 documents connection pooling, retry and colocating compute
in the same AWS region:
https://docs.aws.amazon.com/AmazonS3/latest/userguide/optimizing-performance-guidelines.html
The Open Data source documents medium 500px images and photo-specific license
obligations: https://github.com/inaturalist/inaturalist-open-data

Provider rules and user budget still apply. No rate limit may be bypassed by
IP rotation or disguising a multi-worker job. A distributed cluster needs
operational permission/capacity and an authoritative shared traffic budget.

## Commands — examples, no hidden job launches

Create one source group's photo-ID chunks (metadata only):

```bash
python -m scripts.acquisition.mois_93m_prepare_persistent_source prepare-group \
  --registry source/mois_all_original_photo_target_tiers_receipt.json \
  --receipt source/group-00-receipt.json \
  --archive source/tier-12-group-00.tar \
  --tier 12 --group 0 --output /data/prepared/tier12-group00
```

After all 224 group outputs have been prepared, list their directories one
per line in `/data/prepared_groups.txt`, and verify full worldwide source
coverage (still zero photo requests):

```bash
python -m scripts.acquisition.mois_93m_prepare_persistent_source collect-roots \
  --group-roots-file /data/prepared_groups.txt \
  --output /data/all_signed_source_roots.txt
```

For a local plan check, again without photo downloads:

```bash
python -m scripts.acquisition.mois_93m_persistent_worker \
  --chunks-roots-file /data/all_signed_source_roots.txt \
  --output /data/pixels --worker-slot 0 --worker-slots 5
```

For an explicitly approved bounded source run, with existing signed local
model weights and mounted source/output directories:

```bash
python -m scripts.acquisition.mois_93m_persistent_worker \
  --chunks-roots-file /data/all_signed_source_roots.txt \
  --output /data/pixels \
  --weights /models/jrc_yolo11n_last_v4.pt \
  --weights-sha256 f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5 \
  --worker-slot 0 --worker-slots 5 \
  --approved-aggregate-rps 40 \
  --download-workers 16 --batch-size 512 --device cpu \
  --max-chunks 1 --execute
```

The 40/s budget in this command is *illustrative* and is not a granted
provider quota. Run multiple slots only with the same coordinator/budget,
non-overlapping slot indices, enough compute, an owned durable output store,
and provider-compliant policy. Drop `--max-chunks` only when full-run costs
and operational limits are explicitly accepted.

A durable container is defined in `infra/Dockerfile.mois93m`; neither
the image nor its build command embeds model weights or any photograph.

## What this does NOT resolve yet

- Original 500px YOLO detection boxes are not verified petal masks. The
  provisional colour composition inside a box must not become a botanically
  validated white/pigmented label.
- With ~97% screen-positive candidates in the initial source, expensive
  full-resolution SAM on every candidate would erase the image-stage cost
  savings. The next scientific estimator needs a separately validated
  one-pass flower-organ segmentation/colour output (with rigorous holdouts).
- More photographs do not remove nonrandom observability or exposure bias.
  Preserve unknown, original observation multiplicity, source sampling
  population and region/season/taxon strata in the final macroecology.
- This implementation has not scanned all 93m image pixels or provisioned
  cloud infrastructure. Extrapolate only from source-SHA-verified real-image
  benchmarks, not planning receipts.
