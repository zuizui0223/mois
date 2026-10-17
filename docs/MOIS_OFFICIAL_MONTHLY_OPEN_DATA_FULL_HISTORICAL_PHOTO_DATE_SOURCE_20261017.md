# Full historic MOIS photo-date recovery from iNaturalist's genuine monthly Open Dataset

**Research scope:** recover observation dates for *all* original historical 100,543 MOIS photo IDs (45,953 classifiable, 42,111 named taxa) **without photo-pixel reclassification**, then test whether seasonal photo-colour composition is still affected by taxon/locality/year support. This pipeline does not inherit the within-species `fcp` D/W or 50-km polymorphism estimands, and does not promote photographic automatic hue to botanist-verified open petals.

## Provider metadata source verification: real live check

The [successful MOIS official data source probe, run 39008716324](https://github.com/zuizui0223/mois/actions/runs/39008716324) independently used S3 HEAD + an exact 64-byte HTTP Range GET to verify both original gzip object identities on **17 October 2026**. **No provider photographs and no bulk metadata** were downloaded by that probe. The official files as of the probe are:

| Monthly iNaturalist Open Dataset | Original gzipped size, exact bytes | Last-Modified | Pinned ETag |
|---|---:|---|---|
| `photos.csv.gz` | **20,668,099,997 B (19.25 GiB)** | 2026-09-27 16:45:05 UTC | `"c8ded8376313d0a25357209f8cafd0db-40"` |
| `observations.csv.gz` | **13,420,521,411 B (12.50 GiB)** | 2026-09-27 16:10:27 UTC | `"441dd5c5c2ad2fd74b6dd0bdc0c6150e-26"` |

The immutable source probe receipt is available in the GitHub Actions artifact `mois-provider-official-monthly-s3-photo-observation-size-etag-source-probe` (ID `12091613604`). The two compressed objects total **31.74750 GiB**, and the complete-source stream is deliberately an *opt-in/one-shot* compute job, not a repeated timed download.

The [provider's official source schema](https://github.com/inaturalist/inaturalist-open-data/blob/main/Metadata/README.md) has **photos.csv.gz: `photo_id`, `observation_uuid`** and **observations.csv.gz: `observation_uuid`, `observed_on`, `taxon_id`, `quality_grade`**. Do not attempt the invalid direct join `historical numeric observation_id == observation_uuid`. A correct FULL-SOURCE date join is

```text
original immutable MOIS photo_id
      ↓ exact public licensed-photo ID match
public photos.csv.gz photo_id → observation_uuid
      ↓ exact provider UUID match
public observations.csv.gz observation_uuid → observed_on
      ↓ only if valid current public date exists
MOIS original photo_id + original automated 4-hue outcome + date-provenance flags
```

The dataset is only a monthly snapshot of **open-licensed photos**; historic images that have since been deleted, changed privacy/licence, moved to another observation or removed from the public dataset **must remain date-unresolved**, not artificially assigned the original observed month. The 200 independently date-verified `photo_id → numeric observation_id` API records from [run 39007591041](https://github.com/zuizui0223/mois/actions/runs/39007591041) are independent QC for current observation identity, but cannot substitute for independent numeric observation-ID verification for all 100,543 snapshot records.

## Memory/disk- and provider-conscious actual workflow

`scripts/acquisition/mois_full_historical_source_dates_official_inat_s3.py` performs an integrity-locked **single full stream per original provider metadata file**. It verifies exact `Content-Length`, `ETag` and `Last-Modified` immediately before downloading; GET requests set `If-Match` so a provider update cannot silently change the source during the scan. Compressed streams are verified by gzip CRC and rolling **full SHA256** while being decoded into rows. No 19-GiB/12-GiB file is written to disk, no image pixels are fetched, and only original historical photo-ID-matched metadata and 200-record QC objects are retained in memory. Both source records (matched or missing) are reported; a source stream failure cannot produce a successful science receipt.

The complete-source aggregation preserves **every one of the original 100,543 photo IDs** with a reason for missing provider data; leaves the actual source hue label, source taxon, original observed_on, geographic and moisture variables **immutable**; and only adds a *new column* for a verified current open-snapshot photo observation date. Comparison with 200 public numeric observation-ID pilot dates measures snapshot-to-current date disagreements without rewriting either source.

Only after the actual full stream completes does the pipeline produce:

- 100,543-row history cohort photo-level **current S3 date, original source photo ROI hue, current taxon match/mismatch and precise unmatched reason** (private-to-workflow 90-day artifact, not committed as users' original row-level data).
- Historical-only and entire **267,425-photo** source-stratified original photographic seasonal count/denominator CSVs for northern/southern temperate and tropical observation-month groups.
- Full original source-independent **same-taxon 0.25° geographic cell × local season-year** summer/winter two-observation feasibility counts after date recovery, plus original 200-source date consistency.
- Original full S3 compressed source SHA256/ETags, source photo-union SHA256, originally verified recovery numerator and classifiable four-hue/date denominators, source geographic/photo state missingness receipts.

No full-source results are labelled complete without both provider metadata streams reaching their verified source byte lengths and valid gzip EOF/CRC. GitHub Actions uploads only a labelled non-final diagnostic if provider files are unreachable or exceed its compute window.

## Ecological and publication decision boundary

Even if the full historic 100,543-photo observed-on column becomes nearly complete, automatic RGB-hue ROIs are not necessarily verified *open flowers*, and photographic absence in one season cannot imply true flower absence. Source-selected taxa, observers, taxonomic revisions and location/geographical overlap remain central inferential covariates. Changes in the raw summertime white-label ratio are not necessarily physiological pigment responses to solar radiation or water, and no actual observed-year monthly TerraClimate anomalies are created from WorldClim 1970–2000 climate normals.

The original 267,425-photo annual moisture/genus-bootstrap simultaneous-maxT slopes remain unchanged by this date-only import. The active consolidated Research Article should report **how much new historical site/year seasonal support is actually observed** before any mechanistic seasonal white-flower ecology claim.
