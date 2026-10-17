# Geographic moisture gradients and seasonal taxonomic turnover in the photographed global flora

*MOIS integrated Research Article working draft — 17 October 2026. Not submitted. This draft integrates two explicitly different observational estimands and freezes the historical climate/photo sampling denominators; it is not a new pigment-selection experiment.*

## Abstract

Photo records enable global analysis of floral colour, but taxon composition and photographic effort change across space and seasons. We evaluated whether photographed white-versus-chromatic floral hue varies with annual moisture among congeneric plants and whether summer white enrichment survives same-locality/year comparison. The source contains **267,425 real original photographed ROI attempts**, including **129,795** automatic four-hue classifications. The 100,543-photo historical source received **99,732** independently date-matched photos via a SHA-verified September 2026 iNaturalist monthly S3 photo-ID→observation-UUID join and **85** additional exact public photo+numeric observation-ID API fallbacks. In total, **127,504** original source photos have both date and four-hue ROI label, representing **25,257** nominal taxa. The source-equal-genus global photographed-white response to annual precipitation was **+1.93 percentage points per environmental SD** (24-family maxT 95% **[+0.47,+3.40]**) and to atmospheric vapour pressure **+3.03 [+0.83,+5.23]**; neither direction retained simultaneous interval support under congeneric coarse-geographic-cell conditioning. Across all dated source photographs, Southern extra-tropical summer exceeded winter in photo-white fraction by **+12.49 points**. Exact photographic compositional accounting assigned **+10.60** points to changed relative photographic frequencies among taxa observed in both seasons, **+3.41** to their within-taxon photographic label-rate change and a net **−1.52** to season-exclusive taxon-photo support. Yet source photographs from the **same 0.25° cell and local season-year**, with ≥5 separate original observations and ≥2 named taxa in each season, produced a Southern summer–winter contrast of **−4.63 points** (exploratory 999-geocell-bootstrap 95% **[−16.45,+6.55]**, 18 cell-years), versus Northern **+5.04 [−8.96,+20.55]** (16 cell-years). Just **two Northern and zero Southern** taxon×cell×year groups had ≥2 original source observations in each season. The photo atlas exhibits global hue–climate gradients and seasonal taxonomic/sampling composition differences, but does **not** establish a shared local pigment response to seasonal temperature, water or sunlight, nor the flower-stage identity of each automatic ROI.


**Keywords:** floral hue, iNaturalist, white flowers, species turnover, moisture, observation bias, seasonality, macroecology, climatic normals

## 1. Introduction

The proportion of species bearing white versus chromatic flowers could depend on water limitation, solar exposure, historical diversification, photoprotection, pollinator communities and reproductive strategy. Macroecological gradients are therefore interesting, but comparing colours across places may confuse local environmental effects with **which taxa exist or get photographed there**. Time adds an analogous problem: a summer collection can contain different species from a winter collection even if no plant changes its own flower pigment.

We test two deliberately distinct questions in one source-bounded observational study:

1. **Between-taxon, geographic question:** Do photographed floral hue proportions of congeneric named taxa increase with annual precipitation or atmospheric moisture at the worldwide scale, and does a positive direction survive spatial, image-quality and alternative measurement controls?
2. **Photo-season composition question:** Does white-labelled photographic composition differ between local meteorological seasons? If so, how much of the source-photo difference is algebraically attributable to photographed-taxon frequencies, taxa lacking photo support in one season, and within-common-taxon white-label proportions, and how much local taxon×place×year counterfactual overlap actually exists?

We do not equate within-taxonomic photographic variation with segregating genetic morphs, nor source-seasonal colour with plasticity of an individual flower. FCP studies those independent within-species phenotype-space and local distribution outcomes separately; their estimates are not part of MOIS.

## 2. Materials and methods

### 2.1. Original measured photo sources and source comparability

The original historical MOIS worldwide collection contains **100,543 original photo IDs** drawn across **42,111 nominal taxa**, with **45,953** classifiable photographed hues and **54,590** unclassifiable source-photo outcomes. Photo IDs were deduplicated before any expansion. The expanded corpus combines historical MOIS with *already actually measured* original FCP photographic ROI pixel outputs and additional true original photo processing. Its current source-frozen interim size is **267,425 distinct photo IDs** with **129,795** four-hue classifiers and **137,630** remaining nonclassifiable source photos. The extra FCP wave4 subset contains only genuinely completed ROI partitions (79 of 256 at this snapshot) and is potentially biased by original compute-job completion order. Its 15,017 newly admitted IDs are not equivalent to 49,800 fully measured original photographs.

The 2026-10-17 signed original-photo union [Actions source](https://github.com/zuizui0223/mois/actions/runs/39002957961) has compressed photo-table SHA256 `fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9`. Initially **161,761** photographed source records had a dated observation, including **81,867** classifiable hues. All 100,543 historically broad photo IDs retained original photo/numeric-observation identities but lacked observed_on in the combined analysis export. A completed official September 2026 S3 monthly metadata scan verified **99,732** original historical photo dates; another **85** were independently source photo+numeric-observation-ID verified in a separate current public API cohort. The final combined corpus now has **261,578** date-bearing original photo attempts and **127,504** date+automated four-hue photo labels from **25,257 nominal taxa**. Dates remain unresolved for **726 historic** and **5,847 total** original photo records; one S3/API date disagreement retains the earlier S3 observation date as primary. See the independently signed full-source recovery receipt and the final source-ID/season ledger. These are dates of public observations, not proof of independently verified open flower petals or image EXIF timestamps.

The same underlying iNaturalist/FCP photo-source family is reused for different MOIS and FCP inferential questions, not as independent ecological replication. Photo identity and label conflicts are checked by original photo IDs and source-receipt hashes.

### 2.2. Colour outcome and photographic missingness

The automatic classifier assigns one of **white, yellow/orange, red/pink, blue/purple** when an original photographic flower ROI receives an interpretable fixed hue; otherwise its value remains missing. Our main geographic response is a nominal-taxon photographic white proportion calculated among classifiable original images. Missing hue outcomes are not implicitly designated nonwhite, and an observed-on date is not a bloom-stage confirmation. Existing small-sample phenology auditing showed taxon-dependent annotation missingness and that fruiting/flower organ identities can be mislabelled by automatic photographic ROI methods. Accordingly, statements below refer to **photographed white-labelled ROI fractions**, not white petal morphs established by botanists. We separately report the **four-hue classifiability rate among all dated original photo attempts**, the **white-labelled fraction among classifiable ROI photos**, and the yield of an observed white label among all dated attempts. Because unclassifiable ROIs are missing—not measured nonwhite/flowerless plants—differences in these rates are evaluated as separate observation-process components. Their exact symmetric product decomposition is descriptive, not a causal allocation of floral pigmentation.

### 2.3. Long-term geographic climate and taxonomic inference

We retained the original 6 anchor environmental axes (BIO5 maximum warmest-month temperature, BIO12 annual precipitation, BIO14 driest-month precipitation, BIO15 rainfall seasonality, mean annual solar radiation, mean annual atmospheric vapour pressure). Geographic and genus-centred conditioning, equal-genus weighting and 999 shared whole-genus bootstrap draws yield **24-hue×climate maxT simultaneous intervals within each fixed comparison**, not retrospective familywise protection across all decisions made during exploration. Effect sizes are *percentage points of photographed white proportion per original-source environmental SD*. The original species/support comparability rule and coarse geographic conditioning are maintained in the historical/expanded re-estimation.

WorldClim v2.1 monthly grids represent **1970–2000 climatological normals**. They are not observations of weather during a source photo's actual calendar year. The original WorldClim photo-site v1 raster join erroneously labelled raw solar radiation in **kJ m−2 day−1** as W m−2; where physical flux units are requested we divide the original numbers by **86.4** and retain the frozen source checksum. The forward WorldClim mapping now writes properly converted W m−2. Soil properties, terrain and elevation are static or slowly varying site covariates, not within-site experimental manipulations.

### 2.4. Season, photographed taxa and exact source composition identity

We assigned local meteorological seasons only beyond ±23.5° latitude and reversed them south of the equator. December is assigned to the next local winter/summer season-year. Tropical records are left in month-only summaries. Source-photo seasonal white fractions are the number of classifiable `white` labels divided by the number of classifiable four-hue labels, separately recording original attempted photos.

For summer S and winter W, let \(p_{it}\) be the white photo-label fraction of taxon \(i\) observed in season \(t\), and \(w_{it}\) its relative *photographic frequency* among shared sampled taxa. For the taxa photographed in both seasons, the exact symmetric composition identity is:

\[
\Delta_{\mathrm{common}}
=\sum_i\frac{w_{iS}+w_{iW}}{2}(p_{iS}-p_{iW})
+\sum_i\frac{p_{iS}+p_{iW}}{2}(w_{iS}-w_{iW}).
\]

The **full original-photo** difference is reconstructed without assigning unobserved taxa an opposite-season phenotype:

\[
\Delta_{\mathrm{all}}=\Delta_{\mathrm{common}}
+(1-f_S)(p_{S,\mathrm{exclusive}}-p_{S,\mathrm{common}})
-(1-f_W)(p_{W,\mathrm{exclusive}}-p_{W,\mathrm{common}}),
\]
where \(f_t\) is the fraction of season-\(t\) classifiable photos assigned to taxa photographed in both seasons. This is an **identity** for photographic composition, not a selection or causal phenology mediation analysis. We separately retain original taxonomy/site/year support and make no assumption that a taxon with zero winter photographs has zero winter flowering.

The prior common-taxon two-observation-per-season comparison also reports taxon-resampled 999-bootstrap percentile intervals; those cannot identify geographic or observation-effort independence or substitute for the annual moisture maxT family.

### 2.5. Identifiability of summer–winter locality and year comparisons

At each site width of **0.25°, 0.5°, 1°, 2°, 5°** longitude and latitude, we count nominal taxon×coarse-geographic-cell groups with both summer and winter photographic ROI colours, either pooled across years or restricted to the exact local meteorological season-year. We report thresholds of at least **one** or at least **two distinct original observation IDs** per season. Multiple photos of one observation do not satisfy the higher threshold, and distinct observations still do not establish independent photographers or identical biological populations. Coarser cells deliberately test how apparent repeated sampling changes with loss of spatial identification.

## 3. Results

### 3.1. Global moisture association persists descriptively but not geographically

In the **expanded actually measured 267,425-photo corpus**, the original taxon/genus-weighted geographic photographic-white association with BIO12 precipitation is **+1.932 pp/source SD** (24-comparison simultaneous maxT 95% **[+0.465,+3.399]**). With mean atmospheric vapour pressure it is **+3.031** (**[+0.832,+5.229]**). After genus×coarse-cell conditioning, the corresponding values are **+1.190 [−2.402,+4.782]** and **+2.561 [−2.268,+7.391]**, neither excluding zero.

The original 100,543-photo benchmark had positive global annual rainfall and vapour directions but independent strict-image-quality and congeneric local geography checks failed to confirm the rainfall sign. An exploratory 248,312-photo environment analysis found nominal white–precipitation slopes without a predeclared four-hue family q below 0.05 and a tiny within-sample prediction gain (~0.0014 R²), with no whole-genus holdout improvement. Thus the primary positive result is **global sampling-support conditional**, not a universal direct moisture response.

### 3.2. Summer photo-white enrichment in the final date-complete photographic source

Across the actual **127,504 dated and classifiable original photo ROIs**, Northern extra-tropical summer contains **18,620 white labels among 38,121** classified photographs (48.84%), versus winter **1,850/4,059** (45.58%), an unpaired contrast of **+3.27 percentage points**. In the South, summer is **2,996/5,827** (51.42%) and winter **1,589/4,082** (38.93%), an unpaired **+12.49 points**. The earlier pilot cohort of 86,172 date+colour photos gave +3.30 Northern and +17.23 Southern points because its taxon/photo-source coverage differed; those values are source-study sensitivity diagnostics rather than primary final estimates. Among original historical photos with September 2026 S3 observation dates alone, the source summer–winter photo-white contrast is **+3.20 Northern** and **+7.43 Southern** points. These are photographed four-hue labels among classifiable images, not field censuses of open petals.

### 3.3. Seasonal taxon turnover dominates the global photographic composition

Only **562 of 7,015 Northern** and **320 of 4,139 Southern** source taxa have classifiable photos in BOTH summer and winter (Jaccard **0.0801** and **0.0773**). A taxon not photographed in one season is NOT known to be reproductively absent.

The exact symmetric photographic composition identity partitions the Southern summer minus winter **+12.49-point** contrast into **+3.41 points** from white-label-rate changes within taxa photographed in both seasons, **+10.60** from shifts in their relative original photo frequencies, **+2.67** from summer-exclusive taxon-photo support and **−4.19** from winter-exclusive support. The corresponding Northern +3.27 points partitions into **+2.52**, **−1.11**, **+3.04**, and **−1.19** respectively. These are exact algebraic source photo-frequency terms, NOT biological mediators of pigmentation selection, pollinator preference or within-plant colour change.

### 3.4. Matching geography and season-year does not confirm the South-wide signed difference

At the original nominal taxon × **0.25° geographic cell × same meteorological season-year** grain, requiring at least **two distinct original observation IDs per season** retains only **two Northern comparisons (eight original photos)** and **zero Southern**. True repeated taxa cannot support a generalized local phenotype response from the current photo corpus.

A separate photo-assemblage target allows different taxa in summer and winter in the same coarse location/year. At a source threshold of **≥5 distinct original observation IDs and ≥2 source taxa each season**, the 0.25° matched sample contains **16 Northern cell-years (232 classified photographs)** and **18 Southern cell-years (367 photographs)**. The equal-cell-year summer minus winter photographed-white rate is **+5.04 points [−8.96,+20.55]** in the North and **−4.63 [−16.45,+6.55]** in the South (999 geographically clustered bootstrap percentile intervals). These intervals include zero. Enlarging cells to 1° yields **101 Northern/84 Southern** sampled cell-years but mixes widely separated habitats. The aggregate Southern +12.49-point photo-white enrichment therefore does NOT establish a common within-place/year positive signed seasonal contrast.

### 3.5. Full historical date restoration, source conflicts and taxonomic breadth

The globally broad historical photo export initially contained no observation dates for its 100,543 original photo IDs. A first photo-ID-blind 9,600-photo pilot actually verified **9,597** observation dates and **4,305** classified source hues, expanding the provisional dated/photo-colour count **81,867 → 86,172**. The subsequent complete September 2026 official metadata scan read the gzip photo and observation files sequentially using exact original photo→provider observation UUID→observed-on joining, source ETags, byte counts and full SHA256. **99,732** historical photos had a valid S3 date; an independent original photo+numeric observation ID sample supplied **85 additional** otherwise missing dates, yielding **99,817** dated historical photos and **45,637** dated/classified historic hues. **726** historic photo IDs remain undated. Among **9,512** overlapping provider dates, **one** disagreed between the September snapshot and the later public API. Both values are archived and the older source S3 date is used. The final combined 267,425-photo set has **127,504** dated/classified hues and **25,257** dated source nominal taxa without changing any original photo pixel, hue, geographic record, source taxonomy or the original annual-moisture analysis.

### 3.6. Early date-recovery controls and measurement limitations

The bounded 9,600-photo pilot increased dated/labelled nominal taxa from **2,007 to 5,315**, but over **90%** of its 3,801 newly dated historical source taxa were represented by only one photo. At that earlier stage, the strict source photo community condition above yielded merely **two Northern and one Southern** 0.25° cell-years. Full historical date recovery improves community photo support to **16 and 18**, but almost no same-taxon summer/winter source repetition is added. This separates taxonomic photo coverage from biological local temporal identification.

The earlier exploratory any-two-season match on **384 classified photographs / 61 groups / 47 taxa** used **WorldClim 1970–2000 climatological month normals**, NOT actual observation-year weather. Its white photo-label slopes were temperature (+5°C) **−0.27 [−6.93,+6.34]**, precipitation (+50 mm) **−0.93 [−10.24,+7.78]** and physical solar flux (+50 W/m²) **−0.08 [−4.81,+4.21]**, all exploratory taxon-bootstrap percentiles including zero. The full source date restoration does not convert normals to actually experienced temperature, irradiation or VPD.

The [full historical official S3 source data](MOIS_COMPLETE_100543_HISTORICAL_PHOTO_DATE_RECOVERY_REAL_RESULT_20261017.md) and [final dated photographic hue and locality analysis](MOIS_FINAL_127504_DATED_HUE_PHOTOS_GEOGRAPHIC_TEMPORAL_COMPOSITION_20261017.md) preserve source IDs, effect counts and mixed photo-anthesis uncertainty.

### 3.7. Photo-colour classifiability changes with season and geographical sampling

There are **261,578 dated original photographic attempts**, of which **127,504** have an original four-hue classification; the remaining **134,074** date-bearing photos have **unclassifiable rather than biologically nonwhite/absent** hue states. North extra-tropical four-hue classifiability is **38,121/72,857 (52.32%)** in summer versus **4,059/8,566 (47.39%)** in winter. Southern success is **5,827/12,657 (46.04%)** versus **4,082/9,421 (43.33%)**. Thus the observed photo selection process itself changes with season and sampling geography.

Let c be the fraction of dated original photo attempts receiving a four-hue ROI label and p the fraction of classifiable ROIs labelled white. The observed *white-label yield among all original photo attempts* is y = c × p, and its exact symmetric summer−winter accounting is Δy = mean(c) × Δp + mean(p) × Δc. North original photo-white label yield increases by **+3.960 percentage points**, attributable algebraically to **+1.629** conditional white-label-share points and **+2.331** four-hue classifiability points. South white-label yield increases by **+6.804** points, decomposed into **+5.580** conditional white-label composition and **+1.224** classifiability points. These terms do not measure causal selection on flower pigment; photo-specific botanical stage is unknown.

For the same **0.25° cell × original meteorological season-year** photographic-community source sample used above (≥5 classifiable original observation IDs and ≥2 classifiable named taxa each season), the mean site-year classifiability summer−winter difference actually becomes **−10.78 points** in 16 Northern and **−4.05** in 18 Southern cells. Exploratory geographically clustered 999-bootstrap 95% percentile intervals are **[−21.92,−0.73]** and **[−11.13,+4.41]**, respectively. In the Northern sample, the sign of the global image-classifiability contrast reverses after locality/year matching. This is an observational *measurement and photo sampling design* effect, not evidence that local summer flowers become harder to classify due to an identified physiological mechanism. Local cells were deliberately selected on having sufficient classified photographs, and several source supports/scales were explored: no confirmatory familywise significance is claimed. [Full original photo-denominator and locality matching protocol](MOIS_REAL127504_PHOTO_WHITE_CLASSIFIABILITY_SEASON_SELECTION_20261017.md).

## 4. Discussion

The strongest conclusion is not that moisture, temperature, or light makes an individual plant's petals turn white. Rather, broad **photographed floral-hue gradients** exist at both geographic and seasonal scales and are sensitive to the *taxonomic composition of observed images, place, and the availability of within-locality temporal comparisons*. For worldwide moisture, a positive among-taxon global direction is reproducible in the original estimator but disappears as a supported association after more restrictive geographic conditioning. For seasonality, the large Southern summer excess can be reconstructed from shifting relative photo frequencies of common taxa plus non-overlap in photographed taxa, without any identified individual-organ pigment transformation.

A further observable source limitation is that the probability of obtaining a usable four-hue photographic ROI varies across sampled seasons and reverses its apparent Northern summer direction in the narrowly paired geographic/year subset. Consequently, all ecological hue gradients should be reported with original attempt denominators and classification-success rates, rather than conditioning on apparently suitable colour images without reporting their selection process. This does not establish whether unclassifiable photos were of nonflowering plants or mislocalized petals.

Two different ecological mechanisms remain plausible but observationally unresolved: genuine seasonal changes in the flowering composition of local communities, and seasonal variation in who observes, which plants are photographed, and which flower-like organs a classifier labels. We cannot promote photo taxon nonobservation into a biological absence, or an automatic ROI label into botanically verified anthesis. A controlled flowering-species turnover paper would require independent bloom-stage labels and sampling effort, not simply a larger pool of unverified iNaturalist original photo IDs.

The linked results set an important boundary for ecological colour–climate inference: *source-level interspecific gradients are different estimands from within-lineage pigment plasticity, local adaptation, and evolutionary fitness balance*. This also rules out claiming a universal 50:50 pigment cost-benefit equilibrium from a near-half descriptive proportion among classified photographs.

### Next experiments and independent data

Site-month ecological inference will require photo-specific verified open floral organs, source-dated and accurately geolocated observations in more species×site-year strata with cross-season overlap, explicit effort/observer controls and **observation-year/month** weather (temperature, actual solar radiation, rainfall, water deficit, VPD). Spectroscopy, anthocyanin assays, common gardens, flowering surveys, pollination and reproduction would be necessary to test pigment costs and benefits, rather than their photographic correlates.

## 5. Limitations and data provenance

The historical globally broad sample, later high-depth FCP-derived taxa, and 79 completed ROI wave4 partitions are different photographic source regimes and must not be treated as independent replications; signed image/source receipts and named-taxon denominators are published in [the original MOIS refit](MOIS_267425_GENUINE_PHOTOS_ORIGINAL_MOISTURE_RESULT_20261017.md), [seasonal photo source estimator](MOIS_REAL_267425_SEASONAL_WHITE_CLIMATE_NORMAL_DIAGNOSTIC_20261017.md), and [full photographic taxon overlap audit](MOIS_REAL_PHOTOGRAPHIC_SEASON_TAXON_TURNOVER_SPATIAL_IDENTIFIABILITY_20261017.md). The original historical manuscript [dated 2026-10-14](MOIS_GLOBAL_MOISTURE_MANUSCRIPT_20261014.md) remains the frozen benchmark, not silently replaced by this exploratory 2026-10-17 integrated editorial draft.

Missing hue is nonrandom with image and floral-organ visibility; after restoration **5,847 of 267,425** source photographs still lack any verified original/public observation date, including **726** from the historical photo subset. The geographic climate analysis uses historical normals and equal-genus nominal hierarchy rather than an exhaustive time-calibrated phylogenetic covariance model. Spatial cells are source supports, not biological populations. The global association intervals adjust one six-variable-by-four-colour family inside a chosen comparison; they do not correct for all exploratory tests or selection of covariates and subregions. No experimental flower reproduction, pigment energy budget or independent systematic flowering census was observed.

## Planned figures

1. **Global white photo–moisture direction and inference scale**: annual precipitation and vapour with 999-genus maxT intervals, overall versus genus×coarse-cell, plus explicit original source/frozen expanded-photo cohorts.
2. **Seasonal photographic white fraction and observation denominators**: North/South meteorological seasons, visible colour rates separated from unclassifiable original-photo attempts, source-stratified.
3. **Exact photographic composition account**: symmetric common-taxon within-rate/frequency decomposition plus seasonal photo-only taxon support corrections, all components summing to observed whole-photo summer–winter contrast.
4. **Temporal causal-support collapse by geographical resolution**: summer/winter site×taxon matched photo groups at 0.25°, 0.5°, 1°, 2°, 5°, with and without fixed season-year, one vs two original observations per season.

## Author/editorial checklist

- Confirm contributor/coauthor scope and rights for reused iNaturalist source images; only public attribution-safe summaries should enter the eventual supplement.
- Freeze fully completed FCP wave4 partition population before updating the main annual source slopes; the present expanded source is job-completion conditional.
- Require same-organ petal verification and actual-year weather before any interpretation of **seasonal** climate plasticity.
- Distinguish descriptive photographic effect sizes, independent-source replication, null spatial identification and genuine physiological fitness outcomes in title, abstract and figures.
- Target journal, author list, final methods references, figures and supplementary tables remain **not finalized**; no submission has occurred.
