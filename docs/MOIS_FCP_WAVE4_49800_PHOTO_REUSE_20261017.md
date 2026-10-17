# Reuse FCP's additional 49,800 real original photographed flower-colour measurements in MOIS

## Main objective

Increase the **actual number of photographed visible-hue classifications** in `zuizui0223/mois`, keeping the original 100,543-photo method and particularly the global white–rainfall/atmospheric vapour-pressure question. It is scientifically acceptable for MOIS and FCP to share original **photographs**, not the FCP within-species recurrent W, D, photo pair 50-km local depletion or selection/fitness interpretation.

FCP is independently processing an **additional 49,800 original photographs from 498 species**, beyond its earlier 149,900 high-depth photos, using the original unchanged ROI-v4/EfficientSAM classifier. The authentic source is [FCP workflow run 38992543115](https://github.com/zuizui0223/fcp/actions/runs/38992543115), with 256 original source-blind terminal partitions and frozen original model/metadata receipts. Its image-pixel work is expensive; duplicating it under MOIS would waste compute and risk different sampling conditions.

## Exactly what may be admitted

`scripts/analysis/mois_import_full_FCP_wave4_completed_original_photos.py` requires the **finished 256/256 signed ROI partitions and final photo-level FCP terminal result** (not a source-only roster), matching its original model/metadata SHA256. It reads exactly 49,800 photo ID–taxon–species–site–date–hue–status rows, retains all technical unclassifiable images, and joins them **by original `photo_id`** to the MOIS photographic dataset.

The MOIS prior source is the completed 248,312-photo union plus two independently measured disjoint 4,096-photo original ROI waves. This yields **256,504 photographed source opportunities** before FCP wave4, with the actual flower-color classifier count determined from the two complete ROI source reports. Source photo-ID overlaps are intentionally allowed across FCP and MOIS but **counted only once within MOIS**. If duplicate original photo IDs yield nonidentical four-hue classifiers or taxon IDs, the original MOIS value remains unchanged and an explicit discrepancy count is included in the numeric receipt. The final total after adding wave4 is therefore **256,504 + number of truly novel wave4 photo IDs**, not automatically 306,304 without verifying overlap.

## Execution condition and compute protection

`.github/workflows/mois-import-fcp-wave4-measured-49800.yml` uses a bounded hourly source check through November 4, 2026, without opening images. Only after BOTH of these original physical operations succeed will it execute the real source/photo-level merge:

1. FCP's source/ROI run 38992543115 has **completed successfully** and has the full `fcp-wave4-149900-method-expansion-actual-measured-and-inferred` artifact.
2. MOIS's separate original second 4,096 photo wave (and its 256,504-photo source aggregation) has been completed, not merely source photo IDs selected.

Once the combined SHA-verified actual data artifact exists, later hourly checks skip it and do not re-run the expensive image or ecological models. The workflow is time-limited, not an indefinite service or assumption of successful future imaging.

The first 4,096 new real images already increased MOIS to **252,408 source photo IDs / 122,383 classifier-admitted four-hue records**, with an authentic [updated water maxT result](https://github.com/zuizui0223/mois/actions/runs/39000489558): white~BIO12 global genus-equal +1.947 pp per source SD (maxT95 +0.441,+3.453); white~atmospheric vapour +3.067 pp (maxT95 +0.865,+5.268), without support from the stricter same-genus geographic-cell results. These are source-shared observational associations, not fitness selection or independent ecological replication.

**Never add source-only 3,430,502-photo opportunities or unfinished 49,800-ROI work as classified data.** The expanded globally species/genus-weighted moisture regression must be refit on genuine new hue outcomes, using the same original 999 whole-genus bootstrap and six-anchor×four-hue maxT family, rather than promoting exploratory unadjusted p values.
