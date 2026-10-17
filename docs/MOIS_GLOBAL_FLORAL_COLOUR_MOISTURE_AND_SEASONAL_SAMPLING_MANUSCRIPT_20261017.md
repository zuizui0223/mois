# Geographic moisture gradients and seasonal taxonomic turnover in the photographed global flora

*MOIS integrated Research Article working draft — 17 October 2026. Not submitted. This draft integrates two explicitly different observational estimands and freezes the historical climate/photo sampling denominators; it is not a new pigment-selection experiment.*

## Abstract

Photographic records make it possible to investigate global floral-colour ecology, but the colour distribution among photographed species is also shaped by where and when observations are made. We evaluated whether white versus chromatic photographic flower-colour states are associated with moisture across named plant taxa, then asked how apparent seasonal white-flower enrichment partitions into changes in photographed taxonomic composition and within-taxon photographic labels. A worldwide historical archive contained 100,543 unique photo identifiers across 42,111 nominal taxa, with 45,953 successfully assigned one of four automated visible-hue states. An expanded union included 267,425 genuinely processed photographs and 129,795 classified four-hue states; 81,867 classified labels initially had an original observation date, rising to 86,172 after exact public date verification of 4,305 automatically hue-labelled original historical photographs from a separate photo-ID-blind sample. Reapplying the original genus-equal, whole-genus bootstrap, 24-family maxT estimator to the expanded union associated photographed white with annual precipitation (+1.93 percentage points per original source standard deviation; simultaneous 95% interval [+0.47,+3.40]) and atmospheric vapour pressure (+3.03 [+0.83,+5.23]); however, both intervals included zero after congeneric coarse-cell geographic conditioning. Among date-bearing extra-tropical photographs, the Southern summer–winter white-label difference was +17.23 points. An exact taxon-support decomposition apportioned +11.76 points to changing relative photo frequencies of taxa observed in both seasons, +3.27 to within-common-taxon hue-share differences, and a net +2.21 to taxa absent from photographs in one of the seasons. Despite the broad temporal pattern, only 17 Southern taxon×0.25°-cell×year groups had even one observed photo in both summer and winter, and none had two independent original observations in each season. The worldwide photographed flora therefore exhibits broad moisture and seasonal white-colour *composition* gradients, but these data cannot identify a shared local response of floral pigmentation to water availability, light or season. Photographic sampling, source taxonomic composition and spatial support are central ecological inference dimensions, not merely statistical nuisances.

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

The 2026-10-17 signed original-photo union [Actions source](https://github.com/zuizui0223/mois/actions/runs/39002957961) has compressed photo-table SHA256 `fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9`. Exactly **161,761** current photo records contain an original calendar observation date; **81,867** combine that date with a classifiable hue. In this combined artefact the entire historical 100,543-photo block and the first separate 4,096-original-pixel wave lack `observed_on`; neither is counted as a dated season record. Never conflate original photographic *opportunities*, successful image ROI classifications and independently botanically validated open flowers.

The same underlying iNaturalist/FCP photo-source family is reused for different MOIS and FCP inferential questions, not as independent ecological replication. Photo identity and label conflicts are checked by original photo IDs and source-receipt hashes.

### 2.2. Colour outcome and photographic missingness

The automatic classifier assigns one of **white, yellow/orange, red/pink, blue/purple** when an original photographic flower ROI receives an interpretable fixed hue; otherwise its value remains missing. Our main geographic response is a nominal-taxon photographic white proportion calculated among classifiable original images. Missing hue outcomes are not implicitly designated nonwhite, and an observed-on date is not a bloom-stage confirmation. Existing small-sample phenology auditing showed taxon-dependent annotation missingness and that fruiting/flower organ identities can be mislabelled by automatic photographic ROI methods. Accordingly, statements below refer to **photographed white-labelled ROI fractions**, not white petal morphs established by botanists.

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

### 3.2. A clear photographic summer-white excess crosses hemispheres

Across classifiable, date-bearing extra-tropical images, Northern summer (13,312/27,210) has **48.92%** white labels versus Northern winter (1,481/3,246) **45.63%**, a **+3.30 point** difference. Southern summer (1,627/3,030) has **53.70%** versus Southern winter (765/2,098) **36.46%**, a **+17.23 point** difference. The denominator is deliberately restricted to photographic original ROIs classified into the four fixed hue states; it is not a census of all open flowers.

### 3.3. Taxonomic turnover and photograph frequency partition the full seasonal contrast

Of 1,622 nominal source taxa photographed in Northern summer or winter, only **450** occur in both seasonal source sets (Jaccard **0.277**); of **543** in the South, only **212** occur in both (Jaccard **0.390**). In the South, **51.45%** of classified summer source photos represent taxa with no classified winter photos. Within shared taxa, the Southern summer–winter white photo proportion difference is **+15.03 points**, partitioned into **+11.76 points** from relative photographic taxon frequencies and **+3.27** from their within-taxon photographic hue proportions. Terms for taxa seen in only one season contribute **+4.00** (summer) and **−1.80** (winter) points, exactly reconstructing the total Southern **+17.23**.

In the North, within-common-taxon white rate (+1.99) almost cancels changes in shared taxon photographic frequencies (−1.95); season-exclusive photographed taxon support adds a net +3.26, reconstructing the raw **+3.30**. These identities show that no general within-taxon pigmentation shift is required to produce a large photographic seasonal pattern.

Crucially, the related high-depth `new_fcp` source has Southern summer−winter photograph-white difference **+18.80 points**, whereas the distinct incomplete `fcp_additional_ROI` source subgroup has **+0.52**. Source composition affects numerical magnitude, reinforcing the need for study-specific sample and photo-classifiability controls.

### 3.4. Geographic and observed-season-year matching largely eliminates temporal support

Matching summer and winter original photos *only by nominal taxon*, regardless of spatial distance/year, yields 450 Northern and 212 Southern common taxa. Taxon×**0.25° cell** with at least one classifiable original observation in both seasons yields **131 Northern / 70 Southern groups** pooled across all years. Requiring the **same meteorological season-year** reduces to **27 Northern / 17 Southern**; requiring **at least two distinct observations per season** leaves **one Northern group (four photos) and zero Southern groups**. Even a 5° cell×year provides only 26 Northern/13 Southern groups at the two-photo threshold, but this coarsening is geographically unsuitable for local adaptation inference.

A prior independent exploratory any-two-local-seasons fit, using 47 taxa/61 matched source photo groups/384 classifiable photos and *WorldClim month-specific normals*, found white temperature (+5°C) −0.27 pp [−6.93,+6.34], rain (+50 mm) −0.93 [−10.24,+7.78] and physical solar (+50 W/m²) −0.08 [−4.81,+4.21]; all interval families are univariate **taxon bootstrap percentiles**, not maxT confirmations and not actual-year weather effects. The stricter summer-vs-winter same-locality/year design has inadequate support for this regression altogether.

### 3.5. Recovered real observation dates expose a correctable historical source-export limitation

The historical broad-taxonomy original 100,543-photo source was exported to the combined photo-colour table with *all observed-on fields blank*, despite retaining one-to-one original photo and numeric observation identifiers. The source photo-ID-blind 9,600-photo recovery [completed on 17 October 2026](MOIS_ACTUAL_HISTORIC9600_DATE_WHITE_LABEL_RESULT_20261017.md), reusing the previous 200 exactly verified original-photo responses and retrieving current public observed-on metadata for only 9,400 additional observations. The genuine result was **9,597/9,600** exact public original photo–observation–date matches, **4,305** additional previously undated but original RGB-colour-classifiable photographic outcomes, 13 current public taxon-ID changes and 3 still missing public observations. The effective combined source therefore rises from **81,867 to 86,172** original photo IDs with both a known original/public observation date and an automatic four-hue label. The historical unmatched 90,946 remain unknown, not automatically bloom-free or undated in their original archival record.

The broader/shallowly sampled historical photographed source exhibits local summer–winter white-label fractions of **48.87% versus 43.68%** in the Northern extra-tropics, and **48.10% versus 35.09%** in the Southern extra-tropics. This shows a descriptively positive summertime photo-white pattern even within a differently taxon-selected original source frame, while continuing to confound source season with place, sampled plant taxa, open-organ validity, observer effort and classifiability. **No seasonal climatic or fitness response is established by date recovery alone.** The exact current iNaturalist public observation `observed_on` is not an independently verified camera capture timestamp; none of the 9,600 recovered observations has a current *observation-level* Flowers annotation, a separate nonrandom-missingness variable that must not be treated as biological absence.

**Taxonomic breadth and local identifiability are not interchangeable.** Source-locked reanalysis of the 4,305 newly dated historical colours found **3,801 historical photo-hue taxa**, of which **3,435 (90.37%) have only one newly dated classifiable photo**. The global dated-and-classifiable source taxa increased **2,007 → 5,315**, adding 3,308 named taxa, but the number of same-taxon × **0.25° cell × same local season-year** summer–winter photo comparisons **did not increase**: North remained 27 groups and South 17 with at least one original observation in each season; with at least two distinct observations per season North remained one four-photo group and South zero. This is an actual source-level counterexample to equating bigger worldwide flower-photo samples with a reliable local temporal pigment response. Full methods and support ledgers are in the [source date taxonomy-versus-locality audit](MOIS_HISTORICAL9600_NEW_TAXA_BUT_NO_ADDED_LOCAL_SUMMER_WINTER_COMPARISONS_20261017.md).

This pilot was intentionally bounded below the provider's bulk-data threshold. The [official 2026-09-27 monthly open-data photo/observation metadata](MOIS_OFFICIAL_MONTHLY_OPEN_DATA_FULL_HISTORICAL_PHOTO_DATE_SOURCE_20261017.md) provides a different scalable full-original-photo date recovery method, matching `photo_id → observation_uuid → observed_on` with explicit open-photo/licence coverage and source ETags. Any full-source expansion must be measured separately, not extrapolated from the 9,597/9,600 result.

### 3.6. A separate community-level temporal design is more feasible, but only when spatially coarsened

A *community of photographed taxa* is a different target from a within-named-taxon seasonal pigment change. Using the signed 86,172 date-and-colour photo data, we evaluated photographed white-hue fractions at the same geographical cell and meteorological season-year in both summer and winter while **allowing different taxa** to be photographed in each season. We required ≥5 distinct source observation IDs and ≥2 named taxa per season. At 0.25° spatial cells, just **2 Northern and 1 Southern** cell×season-year matched photo assemblages satisfy this condition. The apparent support rises to 28/12 at 0.5°, 70/31 at 1°, 108/58 at 2° and 120/74 at 5°. This expansion is driven in part by treating spatially separated floral assemblages as comparable places; a five-degree locality should never be interpreted as a population.

At a more permissive cutoff of ≥3 independent observation records and ≥2 taxa per season, 0.25° source cells support **42 Northern and 21 Southern** within-cell/year comparisons (415 and 187 genuinely classified original photos). Before the 4,305 additional historically date-recovered original hue outcomes, the corresponding source counts were 39/16 matched groups and 387/138 classified photographic labels. This is a small but genuine increase in geographic/year photographic-assemblage opportunity, although **same-taxon** summer/winter within-cell/year repeated-photo support remained virtually absent. Comparisons within selected coarse cells condition on site-level constant backgrounds by construction but do not control taxon turnover, photographic effort, unverified petal identity, or actual observation-year weather. The source method, complete geographic thresholds and exploratory clustered photo-white percentiles are reported in the [community seasonal sampling scale audit](MOIS_PHOTO_COMMUNITY_SEASON_SCALE_IDENTIFIABILITY_20261017.md).

## 4. Discussion

The strongest conclusion is not that moisture, temperature, or light makes an individual plant's petals turn white. Rather, broad **photographed floral-hue gradients** exist at both geographic and seasonal scales and are sensitive to the *taxonomic composition of observed images, place, and the availability of within-locality temporal comparisons*. For worldwide moisture, a positive among-taxon global direction is reproducible in the original estimator but disappears as a supported association after more restrictive geographic conditioning. For seasonality, the large Southern summer excess can be reconstructed from shifting relative photo frequencies of common taxa plus non-overlap in photographed taxa, without any identified individual-organ pigment transformation.

Two different ecological mechanisms remain plausible but observationally unresolved: genuine seasonal changes in the flowering composition of local communities, and seasonal variation in who observes, which plants are photographed, and which flower-like organs a classifier labels. We cannot promote photo taxon nonobservation into a biological absence, or an automatic ROI label into botanically verified anthesis. A controlled flowering-species turnover paper would require independent bloom-stage labels and sampling effort, not simply a larger pool of unverified iNaturalist original photo IDs.

The linked results set an important boundary for ecological colour–climate inference: *source-level interspecific gradients are different estimands from within-lineage pigment plasticity, local adaptation, and evolutionary fitness balance*. This also rules out claiming a universal 50:50 pigment cost-benefit equilibrium from a near-half descriptive proportion among classified photographs.

### Next experiments and independent data

Site-month ecological inference will require photo-specific verified open floral organs, source-dated and accurately geolocated observations in more species×site-year strata with cross-season overlap, explicit effort/observer controls and **observation-year/month** weather (temperature, actual solar radiation, rainfall, water deficit, VPD). Spectroscopy, anthocyanin assays, common gardens, flowering surveys, pollination and reproduction would be necessary to test pigment costs and benefits, rather than their photographic correlates.

## 5. Limitations and data provenance

The historical globally broad sample, later high-depth FCP-derived taxa, and 79 completed ROI wave4 partitions are different photographic source regimes and must not be treated as independent replications; signed image/source receipts and named-taxon denominators are published in [the original MOIS refit](MOIS_267425_GENUINE_PHOTOS_ORIGINAL_MOISTURE_RESULT_20261017.md), [seasonal photo source estimator](MOIS_REAL_267425_SEASONAL_WHITE_CLIMATE_NORMAL_DIAGNOSTIC_20261017.md), and [full photographic taxon overlap audit](MOIS_REAL_PHOTOGRAPHIC_SEASON_TAXON_TURNOVER_SPATIAL_IDENTIFIABILITY_20261017.md). The original historical manuscript [dated 2026-10-14](MOIS_GLOBAL_MOISTURE_MANUSCRIPT_20261014.md) remains the frozen benchmark, not silently replaced by this exploratory 2026-10-17 integrated editorial draft.

Missing hue is nonrandom with image and floral-organ visibility; original date is missing for 105,664 of 267,425 measured source photos. The geographic climate analysis uses historical normals and equal-genus nominal hierarchy rather than an exhaustive time-calibrated phylogenetic covariance model. Spatial cells are source supports, not biological populations. The global association intervals adjust one six-variable-by-four-colour family inside a chosen comparison; they do not correct for all exploratory tests or selection of covariates and subregions. No experimental flower reproduction, pigment energy budget or independent systematic flowering census was observed.

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
