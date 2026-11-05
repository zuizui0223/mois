# Species-evidence and same-year ALAN macro analysis

User mandate: increase real data, then restrict species, then execute global macro association analysis; inspect and use repository methods. Replaying the old 1,044-cell model is not completion.

## Approach and decisions

Use versioned Kendall & Nicholson experimental diel pollination data (Figshare 26004538 v1, CC BY 4.0), including documented flower colour, study IDs, and paired day/night reproductive outcomes. Match published species and explicit source synonyms to exact WCVP accepted IDs in MOIS. Retain ambiguity as unresolved. No species inherited from a genus-level pollination syndrome. Quantitative night-vs-day reproductive success is not moth identity or global invariant nocturnal dependence.

Alternatives considered: enlarged GloBI visits alone does not establish nocturnal effective pollination; model-based image error correction cannot manufacture missing independent classifier results. Therefore experimental species evidence is the primary addition; original ROI9 is supplementary and unknown stays unknown.

Select support without inspecting ALAN associations. Use historical_main source arm, one original observation/photo ID, 2014–2020 (intersection of observed VIIRS service and annual WorldPop service), and valid coordinates. Exact year matching is mandatory. Static climate normals/elevation are explicitly static. Preserve all known-flower-colour opportunities regardless of ROI success. Determine minimum species/location support before fetching and fitting, record an auditable selection ledger.

Global macro inference is a comparison of photographed species assemblages or species-level sampled environmental distributions. It is not species absence, extinction, abundance, a causal intervention, nor a universally representative flora census. Geographic and taxonomic replication must be explicit; do not treat many photographs as independent species. Species selection is by experimental provenance and support, not coefficient direction.

## Methodology use

Reuse original photo/observation uniqueness, SHA locks, source-arm separation, strict abstention, geographic blocks, and cloud-free VIIRS QA. Audit method/mois-hui-walter-identifiability-20261105 at 41bc2520b9313d4cfce8b3d0f9ab873cc12e6e7c. Two genuine independent white-petal classifiers do not exist; do not fit Hui–Walter or label published species flower colour as an independent same-photo classifier. Preserve experimental heterogeneity and phylogenetic dependence where supported.

## Frozen broad model before coefficient inspection

Source assembly completed: 7,722 species have precise observations; 8,667,448 precise observations in 2015–2024 yield 83,991 primary cell-years, 23,612 cells, and 556 five-degree spatial blocks. Northern/southern cell-year counts are 72,629/11,362. These are source support, not fitted environmental results.

Primary outcome is the percentage of documented single-colour species that are white in each eligible cell/year. Equal cell-year weight; log1p nonnegative VIIRS, log1p same-year WorldPop density, same-year temperature, log1p same-year precipitation, log1p original precise plant observation effort, and log documented richness are predictors. Absorb ten-degree region and year fixed effects; cluster uncertainty across five-degree blocks. Sensitivities: separate hemispheres, within-cell/year effects for cells with >=3 observed years, and genus/cell/year composition with genus/region/year effects and two-way spatial/genus clustering (each cell-year total weight one). None establishes causation or species absence. Genus controls and clustering are not a complete phylogenetic model; ecological cell association is conditional on the frozen species pool.

Same-year climate addition uses NOAA PSL NCEP/DOE Reanalysis II monthly 2 m temperature and precipitation rate, integrated over the exact observation year with twelve required months. This coarse native Gaussian grid adjusts broad climatic variation only; no claim of 0.25-degree climate resolution. High-resolution static normals, if used later, must remain separately labelled. Failed/missing environmental matches are retained in the exclusion ledger; no interpolation across years. Nine-point VIIRS is an annual sampled proxy, not official annual VNL or plant-level illumination.

## Completion contract

New published experimental and trait data downloaded with hashes; eligible species and exclusions enumerated; new same-year environmental data actually retrieved with year/quality checks; an identified macro model fit on the selected real cohort; tests and independent numerical checks; source-linked report and reproducible outputs. Unsupported models remain explicitly unestimated. No automatic merge or external communication.

## Scale amendment after explicit user correction

The 25-species/1,122-observation measured-photo subset is a technical pilot only and MUST NOT be the final global cohort. Read the entire frozen 53,652,241-observation source (93,209,351 photo IDs) through its 16 signed source archives. Extract every observation for source-matched experimental species, without per-taxon caps or ROI-success filtering. The exact source frame currently includes 90 experimental species with 1,421,838 candidate photograph IDs; these are metadata opportunities, not newly measured flowers. Retain original observer, date, coordinates and uncertainty. Extract 2014–2024 so later same-year environmental extensions do not require rescanning. Compute background observation effort across the original plant frame, never treating a non-record as an absence. Primary coordinate precision threshold is recorded <=5 km; missing precision is a separate sensitivity cohort. Geographic environmental sampling and final species support thresholds will be frozen after full-source coverage audit and before model results.

The independent public experimental corpus is Kendall & Nicholson, Pollination Across the Diel Cycle, Ecology Letters 28:e70036, DOI 10.1111/ele.70036; source data DOI 10.6084/m9.figshare.26004538.v1 (CC BY 4.0). Original raw rows, final analysis export, MD5/SHA checks, and experimental heterogeneity remain in the local source ledger. Published metadata-wide abstracts report 139 spp./135 studies; the downloaded final day_night analysis actually contains 138 spp./134 studies/400 comparisons, which is the version-specific cohort used here.

## Broad macro layer (before environmental coefficients)

User scale correction also requires a broad flora layer, with experimental diel species as a nested evidence layer rather than a claim that 90 species represent global flora. The frozen independent Island flower-colour ledger has 39,890 high/medium resolved species records, including multistate records. Only known SINGLE flower-colour states enter the binary primary analysis. Among original source taxa this yields 9,461 accepted species / 9,527 nominal taxa before the support gate. Choose one nominal source taxon per accepted species using the largest original photo-opportunity count (tie: taxon ID), and require >=100 source photo IDs before inspecting exposure. Result: 7,730 accepted species (2,991 white), 32,701,077 candidate photo IDs, not measured flowering images. Record mixed/unknown colours as excluded, never nonwhite. This threshold is an investigator decision before environmental fitting.

Primary geographic frame: 2015–2024, recorded coordinate uncertainty <=5 km, 0.25 degree cells, unique accepted species per cell/year. No per-species observation cap. Main cell-year gate: >=10 documented species and >=5 genera. Retain all-source plant observation effort for adjustment; unrecorded species are not ecological absences. Environmental variables must use the same observation year; WorldPop Global2 R2025A permits 2015–2024 consistently. Source species labels are not redistributed: only cell/year and genus/cell/year aggregates are outputs. The original full experimental-species ledger and raw metadata remain distinct from this broad aggregate layer.
