# MOIS automated completed-pixel to environment to durable-release continuation

This is a **post-processing workflow** for the fifth 4,096 original
image-pixel wave in PR #67; it does not classify unmeasured provider
photo-ID metadata as hues. GitHub workflow_run triggers **only after the
named real 16-shard photographic ROI-v4 workflow fully succeeds**, from
the verified same-repository source branch. A manual workflow_dispatch
with a completed run ID is also available.

The authenticated input must contain exactly 317,631 real photo
attempts, including original 313,535 unmodified source RGB-photo/date/
hue records and another 4,096 original photos genuinely processed by
the immutable YOLO+EfficientSAM workflow. The 4,096 measurement outcomes
and hue frequencies are not asserted before they exist. It excludes
the separate season-targeted 4,096 reservation from world-white
photo frequency and the former 512/1792 partial snapshots from being
counted twice.

The follow-up reads existing immutable original photo/annual climate/
partially available static-source 313,535-Parquet snapshot and all19
WorldClim BIO/VAPR/terrain Parquet, checks source/photo IDs, outcomes
and SHA before appending. **Only the newly measured 4,096 photo sites**
are resampled from WorldClim 1970–2000 BIO1–19, 12-month mean vapour,
12-month mean solar and 2.5 arcminute DEM/terrain. Neither earlier photo
pixels nor earlier climate site values are recomputed; *source soil
remains absent* until measured. Original solar kJ m^-2 day^-1 is converted
to W m^-2 by division by 86.4.

Four durable Zstd photo-ID-keyed Parquet datasets and a SHA manifest
include the new measured visible-hue results and descriptive
photo/taxon/genus white proportions, with unknown hues preserved.
The outputs are published only after true pixel-success, in a
new versioned GitHub Release, so the only copy is not the expiring
90-day GitHub Actions artifact. Raw licensed original photo pixels
are never uploaded.

As the analysis has no independently verified open flower/petal
condition, source data cannot prove a 50:50 evolutionary balance,
pigment production costs, pollinator selection, or relative fitness.
The additional 4,096 candidate source comes from 1,000 high-depth
taxa and is not a uniform global biological census. The main-world
annual climatic PCA/conditional predictor ranking should be
recomputed only at suitable later data freezes, with missing soil
held distinct and climatic seasons in SI.
