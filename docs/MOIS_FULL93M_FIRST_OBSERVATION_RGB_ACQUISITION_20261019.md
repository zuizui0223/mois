# MOIS: from all 93,209,351 original source photo IDs to genuine photographed hue data

This is **NOT** another capped photo candidate source. It consumes the
full, signed 2026-09-27 42,111-taxon original metadata index after the
source-index job has emitted *every one* of the 93,209,351 unique photo IDs
and corresponding original observation date/public coordinates.

The source catalogue is retained in full. Observations receive a fixed,
pre-colour hashed probability draw with p=0.0002 for the first incremental
acquisition milestone. A selected observation contributes every licensed
source photo (rather than FCP's one-per-species or 100-per-species quotas).
All photos not selected at this milestone remain in the immutable full
catalogue and are eligible at future higher thresholds, up to p=1.
The probability is an acquisition probability *within the observed
iNaturalist source*, not an unbiased probability of flowers in nature.

The first actual image runner is intentionally gated on:
1. The full 93m original photo-ID/source-observation Release, successfully
   published and source-SHA verified (not the older per-taxon count CSV).
2. Completion of the 16,384-photo cached tail already running when the
   full index job was activated. Thus a new group cannot compete with its
   original photo-loader jobs.
3. Five exact previous source-verified photo outcome Parquet releases
   (the original 317,631 core, broad 4,096 taxa source, independent wave2
   4,096, cached 8,192 and completed cached tail 16,384). Duplicate
   photo IDs are omitted from new inference, never double counted.
4. Source photo, observation and versioned per-partition SHA verification,
   including the exact original 42,111-source count digest.

Model processing then uses the original frozen location-blind
YOLO11+EfficientSAM ROI-v4 with FCP's rate-limited original photo-loader
and bounded prefetch. One runner reuses the same model bootstrap for up to
four 256-photo source units, and every actual coloured/unclassifiable/
acquisition-failed source attempt receives an original photo/observation
ID and a per-unit signed result. Completed work units are always reused
after explicit restart; no raw original licensed photo bytes are published.

GitHub Actions emits an actual-photo classified four-hue Parquet
and a provenance receipt *only if all source-signed ROI units finish*.
No unseen hues are fabricated from source metadata, no FCP
within-species polymorphism/W-D inference is imported, and no ecological
PCA/SDM or season result is calculated. Species turnover remains
legitimate for future MOIS macroecology, but all hue prevalence analyses
remain postponed per user's instructions.

**Caution:** This first p=0.0002 milestone is a nested, independently
predeclared subset of the full eligible photo universe; it is NOT a limit
on the permanent source frame or a claim that the full 93m image pixels
have already been processed. The extremely long-term p=1 full-image
census depends on separate storage/compute/resource limits and the
source's published access conditions. No recurring unlimited request
loop or uncontrolled S3 fetching is launched.

Source index:
.github/workflows/mois-build-full-photo-index-20261019.yml

Precolour source queue:
scripts/acquisition/mois_full_index_observation_queue.py

Signed original ROI executor:
scripts/acquisition/mois_full_index_real_original_photo_measure.py

First incremental genuine image workflow:
.github/workflows/mois-full-index-incremental-RGB.yml
