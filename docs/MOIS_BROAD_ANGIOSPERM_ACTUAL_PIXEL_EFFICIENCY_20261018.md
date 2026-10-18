# MOIS: actual global species-first flower-photo ROI with minimal wasted compute

The primary ecological sample should cover many **angiosperm named species,
genera and environmental conditions**, not maximise the number of spatially
repeated photo observations of a few FCP high-depth species.

## Actual image acquisition after the independently frozen global source

The new official S3 monthly iNaturalist taxonomy, observations and photos
source includes a first 4,096 original species-prioritised, open-licensed,
research-grade observation photo-ID plan. The initial metadata work is
expensive and source-version pinned by ETag/complete SHA, but is performed
only once and cached for up to 32,768 subsequent photo IDs. These are
PHOTO IDs, **not** automatically classifiable flowers or actual petal hues.

This workflow waits for the same-repository source-only snapshot to
complete successfully, then processes **actual 1024px large photo RGB
pixels** with the identical location-blind historical YOLO11/EfficientSAM
ROI-v4 frozen model/detector weights and source hash.

To minimise ineffective image compute:

1. Measure **one 256-image original-photo pilot** first, with all model
   batches and scientific source SHA receipts. An arbitrary but
   prespecified **20% four-hue classification success threshold**
   stops a largely nonflower/poor-photo source from triggering another
   3,840 expensive downloads. It is a compute-quality gate, not an
   imputation of colour for an unknown photo or an adaptive choice of
   desired white versus coloured labels.
2. When pilot success meets the gate, use **four concurrent workers**
   rather than sixteen redundant model-bootstrap jobs. Each worker
   installs the SAME model once, then handles **three or four 256-image
   subshards** using the unchanged FCP photo-provider rate limiter.
   Four worker hosts have independent pacing, so transient provider
   rate-limit or access failures must abort unsealed photo units rather
   than silently become unclassified flower observations.
3. Validate all 16 authentic original 256-photo worker/source receipts,
   every photo_ID/model SHA, the observation UUID, named species,
   taxon ID, photo licence and original date. Restore the immutable
   313,535 previously photographed photo-ROI Parquet and append
   only genuine **new 4,096 source photograph attempts**, keeping each
   technical colour-unknown outcome separate. Do not reuse the
   unfinished old FCP-1,000-taxon source photo reservation or
   season-targeted SI source photo IDs.
4. Save a permanent additional photo-ID/observed-day/colour Parquet and
   a new taxon/observer/photo-licence CSV plus a signed audit of:
   new *classifiable named species*, new *classifiable genera*,
   geographic source coverage and
   **novel classifiable species/genera per 1,000 genuine RGB-photo attempts**.
   A versioned GitHub Release is published only after the actual
   16 worker results all pass. Image bytes remain ephemeral, with no
   copyrighted photo file redistributed.

The same original imagery may have leaves, buds and fruit, because
iNaturalist named-angiosperm photo metadata are not independent
floral anthesis annotations. Automatic visible hue categories are
not biochemical measurements of anthocyanin, UV flavonoids, metabolic
pigment costs or reproductive fitness. Missing classified hue is not
nonwhite; new observation/taxon source arms are not a random census
of flowering plant individuals and must be weighted accordingly.

WorldClim annual BIO1–19, VAPR, irradiance and DEM/terrain are
subsequently extracted only at new photo coordinates using the
original photo-ID incremental Parquet design; no wholesale
reprocessing of previously measured 313,535 RGB images is necessary.
Soil values are left missing unless externally/source verified.
Monthly climatic seasonality remains SI, distinct from the main
global white-versus-chromatic world flora question.

Main source:
scripts/acquisition/mois_global_angiosperm_breadth_source.py

Actual ROI worker:
scripts/analysis/mois_measure_broad_angiosperm_original_photo_roi_v4.py

Source- and outcome-verified join:
scripts/analysis/mois_integrate_broad_angiosperm_real4096.py

Workflow:
.github/workflows/mois-broad-angiosperm-actual-image-pilot-and-batched-4096.yml
