# Refit the original MOIS six environmental axes on **309,439** actually measured photographic ROIs

**Research scope:** worldwide among-taxon photographic four-hue colour versus precipitation, vapour pressure, radiation and long-term climatic normals. This is independent of FCP's within-species W/D/50-km pigment-polymorphism analyses and does not claim photo-specific open petals, pigment benefit/costs or observed-year meteorology.

## Fully source-verified genuine photo census

The [successful and merged all256 FCP photographic pixel source run](https://github.com/zuizui0223/mois/actions/runs/39028400690) has **309,439 distinct real original photo IDs**, including **149,662 frozen four-hue ROI labels** and **149,653 with a valid recorded original/public observation date**. It used all 256 genuinely finished FCP worker terminal partitions and added **7,197 original photographic pixels /3,601 real colour classifications** beyond the previously complete 302,242 measured-photo source. Previously existing photos/locations/nominal taxa/automatic hue and 99,817 historic public S3/API source dates were not modified.

The original [actually completed 295,017 original photographed-pixel genus-999 water-model run](https://github.com/zuizui0223/mois/actions/runs/39016728034) is the **environmental source baseline**, with site-level WorldClim BIO5/BIO12/BIO14/BIO15 and raw annual solar historical normals on all 295,017 original photo IDs, plus a separately sourced official vapour normal raster product. Its exact compressed climate CSV SHA256 is `84c651aeaacf0e317445593c249a85940d230a19ec8b7f7223ec0e46553013f7`. The newly enlarged 309,439-source original-photo GZIP SHA256 is `546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e`.

This gives **14,422** literal source photo IDs absent from the prior climate table, containing **7,124** original classifiable four-hue outcomes: **7,225/3,523** actual partial-FCP image additions after the earlier 295,017 cohort, followed by **7,197/3,601** from completion of every remaining source image worker partition. Only these 14,422 previously unjoined actual photo coordinates need fresh WorldClim raster point sampling. All older 295,017 original photo/taxon/hue cells are directly compared before the new model can run.

## The original statistical estimator, not a new post-hoc shortcut

`scripts/analysis/mois_refit_original_water_309439_full_FCP256.py` reuses the existing unchanged numerical model of `mois_refit_original_water_with_actual_new_original_photos.py` and historical study's `mois_refit_original_100k_water_six_anchors_248k.py`:

- six original prespecified environmental axes: WorldClim 2.1 BIO5 (warmest-month maximum temperature), BIO12 (annual precipitation), BIO14 (driest-month precipitation), BIO15 (rainfall seasonality), annual solar radiation normal, and mean annual atmospheric vapour pressure normal;
- all four original automatic photo-hue classes, including white;
- taxon-level photographic-colour contrasts standardized per **original combined source environmental SD**, genus-equal weighting, global geographic controls and an alternative genus×coarse geographic-cell comparison;
- **999 identical whole-genus bootstrap resamples**, with a 24-comparison maxT simultaneous interval family within each fixed comparison, and unchanged original source sample/seed methodology.

A physically important source warning is retained: the original WorldClim solar product stores raw **kJ m⁻² day⁻¹**, although the old source column is historically labelled W/m²; the annual standardized estimator is invariant to a positive constant conversion, but no raw solar coefficient should be interpreted as watts without division by 86.4. This run does not estimate actual observation-month anomalies.

**Independent prior benchmark:** the 295,017-photo version, with 999 shared genus resamples, had photographed white versus annual rainfall **+2.005 pp/source SD** (24-family maxT 95% **[+0.559,+3.450]**) and atmospheric vapour **+3.185 [+1.054,+5.316]**. Their corresponding genus×coarse-cell photo source intervals included zero for rainfall (**+1.010 [−2.523,+4.542]**) and vapour (**+2.217 [−2.591,+7.024]**). These are explicitly the **295,017-source benchmark, not estimates fitted to 309,439**. The full new source should be compared after its real model job passes with independent source receipts and a 309,439-photo climate-side SHA.

## Traceable source execution and inference limits

The [full data-and-water GitHub Actions workflow](../.github/workflows/mois-real309439-full256-original-genus999-maxT-water.yml) verifies the original FCP and MOIS measured pixel IDs and the full previously finished 295,017-site WorldClim source receipt, then downloads the official SHA-frozen BIO/SRAD raster ZIPs and fixed vapour archive; only completed original pixel outcomes enter a completely recomputed **genus-999 /maxT** source fit. All transformed photo-site source rows, per-comparison model numbers and source-bound receipts are archived as a GitHub Actions result artifact.

No unchanged earlier climate coefficient can automatically be claimed as a 309,439-photo coefficient. The source selection of plant taxa and photos is nonrandom, photohue classification may include nonpetal/bract/fruits, and no local within-plant pigmentation adaptation, field pollination or reproductive fitness outcome was tested. The genuine original-date **149,653-photo seasonal colour composition** results are a separate estimator and are already presented in the consolidated MOIS Research Article, with additional within-locality source-date and observer-effort caveats.
