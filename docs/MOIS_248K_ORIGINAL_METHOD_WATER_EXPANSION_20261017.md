# MOIS photographed flower-colour × global moisture, same original estimator

**Research task priority, 17 October 2026:** increase measurable photo count; rerun the *original* 100,543-photo interspecific moisture method. This is **MOIS** among-species/congeneric environment inference; source photos may be shared with FCP, but FCP's within-species D/W/50-km morphological statistics are not borrowed.

## Existing measurable photographic observations

The old MOIS baseline was **100,543 original photo IDs** and **45,953 original classifier-admitted four-hue labels** across 42,111 candidate named taxa. Currently, **248,312 original photo IDs / 120,807 four-hue labels** are physically materialized using exactly original photographed hue classes. Original source files are immutable.

A separate complete date/geo census on a 1,000-taxon outcome-blind original subset contains **3,430,502 source photograph opportunities**; an ID-level comparison against the current 248,312 table on 17 October 2026 observed **9,688 literal photo-ID intersections and 3,420,814 previously unused original IDs**. The mathematically possible union is **3,669,126** unique source IDs, but the 3.42 million are **unclassified image opportunities**. They do not automatically count as measured flower hues, nor are they 3.42 million new taxa (they cover only 1,000 of the original 42,111 taxa). These are a bounded, already materialized raw candidate pool, not the entire >93M possible photos in the 42,111 taxa frame.

## Unchanged original 100k-photo moisture method

Reusing functions **directly**, not approximating their statistical structure, from `scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py` (which was copied byte-for-byte from the historical FCP code into MOIS):

1. Construct photographed four-hue *species mean* shares; retain all original photo-ID duplicates only once, and retain classifier failure in source denominators, not as a nonwhite value.
2. Normalize climate axes on the authentic original global photographic frame; remove each original named genus mean, control latitude and longitudinal sine/cosine and weight each qualifying nominal genus equally. Original per-genus minimum is three species.
3. Use the **same original six prespecified environmental anchors**: warmest-month temperature BIO5, annual rainfall BIO12, driest-month rain BIO14, precipitation seasonality BIO15, annual mean solar radiation, annual mean water vapour pressure.
4. Apply **999 multinomial whole-genus bootstraps** with the identical resampling draw shared by all 24 four-hue×six-axis contrasts; report both nominal genus-cluster 95% and the **24-contrast within-scope maxT 95%** confidence intervals. Estimate the old historical-only sample using the same raster source as a compatibility comparison, not an independent confirmation. Add genus×coarse-cell sensitivity, but never treat this as biological field-population matching.

New water-vapour data are read at the authentic photo coordinates from the official **WorldClim 2.1 10 arcminute monthly water-vapour pressure source** `wc2.1_10m_vapr.zip`. The existing 19 BIO, rainfall, elevation and long-term annual solar values were already raster-sampled for 248,312 images from original WorldClim. Solar source units are kJ m⁻² day⁻¹, correctly converted to W m⁻² before the original source-photo SD normalization. All long-term normal values reflect 1970–2000 climatology, **not the weather in the actual photo month**.

**Inference caveats:** historical MOIS strict-image-quality samples cannot be silently declared passed on the expanded photo union because original strict fluorescence/pixel-quality columns were not preserved in that current union. The new sample reuses the same source family and classifier, and a confirmatory published water-selection claim requires original quality, spatial support, MNAR and independent biological flower-organ scrutiny. The original archived 100k-paper result is retained untouched.

## Results and next measured-photo ingestion

The [independent original maxT CI pipeline](../.github/workflows/mois-original-6water-anchors-999bootstrap-248k.yml) actually downloads monthly vapor rasters, samples image sites, and calculates new original-method regression outputs for the current 248k photo census. Until that run succeeds, do **not** attribute new p-values to the 999/bootstrap method; the earlier 248k exploratory white/BIO12 +2.138 pp and BH q=0.0508 used a separate joint-adjustment estimator.

Original FCP wave4 is separately currently processing **49,800 new source photo IDs** under its own unchanged ROI-v4/EfficientSAM classifier; MOIS may import only completed original photo-level terminals when they become available, never FCP's geographic D, W, RGFCA or ecological-statistical decisions. The remaining 3.42 million original image opportunities need staged image downloads and original ROI-v4 classifiability measurements before they can be admitted to the MOIS moisture model.
