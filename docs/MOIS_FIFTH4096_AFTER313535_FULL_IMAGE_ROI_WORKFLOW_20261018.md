# MOIS fifth 4096 original RGB photograph pixels beyond prior313535

The primary atlas currently has **313,535 truly original RGB/photo-ROI
attempts**, **151,239 four-hue labels**, and **71,427 white-labelled**.
The 3,430,502 provider original photo IDs are opportunities, not completed
phenotypes. No number of new classified hues can be claimed until real pixels
have been scored.

This new blinded source wave excludes ALL 313,535 previously measured
original photo IDs, and also excludes the entire 4,096 original photo-ID
reservation selected for a separate supplementary local-season study
(256 of those photos had already been genuinely processed). It selects
4,096 new IDs across 1,000 original taxa by deterministic photo hash,
with original observation and observer IDs unique within the wave.
No new photo's hue, site climate, fitness, local bloom label or pigment
chemistry is read to influence selection.

The workflow performs the fixed original FCP location-blind YOLO11 +
EfficientSAM ROI-v4 model on **16 x 256 genuine original image pixels**,
four 64-photo verified worker batches per shard. Every returned colour
class or missing/uncertain state is preserved with a source and model
receipt. Original image bytes are transient and not redistributed.

Only after all 16 worker jobs are actually successful will their source
rows be joined to the previously authentic 313,535-photo original record.
The target is **317,631 distinct actual photo RGB attempts**, with the
number of four-hue classifications unknown before pixel processing.
Earlier photo_ID, taxon, hue and observation_date identities are immutable.

Next pipeline stage after source-authenticated worker success: derive
new-site WorldClim annual/19BIO/VAPR/solar/DEM values **only at 4,096
new photo IDs**, then append to existing SHA-backed 313,535 Parquet
feature stores; do not recompute all images or invent soil/terrain.
Source photo archives in GitHub Actions expire; publish an immutable
Release when actual integration has succeeded. The current source
frame is a 1,000-taxon balanced photographic sampling design, **not
an unbiased worldwide plant population census**. The main question
remains near-half WHITE versus chromatic photographed flora and their
abiotic correlates, not a maintained 50:50 adaptive equilibrium.
