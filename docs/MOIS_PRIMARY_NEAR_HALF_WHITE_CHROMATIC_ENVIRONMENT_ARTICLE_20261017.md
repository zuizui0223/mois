# Why is white so common? Global floral photographic colour balance across environmental gradients

**MOIS primary ecological Research Article — new working draft, 17 October 2026. Not submitted.** This is a source-verified observational study and not an experimentally established pigment cost–fitness balance. Earlier integrated geographic-moisture-and-season draft remains historical context; seasonal climate results will be **Supplementary Information (SI)**, not a condition requiring within-individual seasonal flower-colour change.

## Abstract

Why is white so frequent in the world's photographed flowering flora despite diverse coloured floral signals? The question is not whether floral evolution maintains an exact 50:50 white–coloured equilibrium; rather, what ecological conditions are associated with a conspicuously large white fraction and when does its relative prevalence shift along environmental gradients? We combined **309,439 distinct genuinely image-processed original source photographs**, comprising **149,662 automatically classifiable visible-hue flower-like ROIs** and **159,777 unclassifiable photographic attempts**, with source-matched global annual environmental data. Among classified images, **70,662 (47.21%)** were white-labelled and **79,000 (52.79%)** were nonwhite across yellow/orange (46,481), red/pink (16,020) and blue/purple (16,499). The equally weighted means of source-taxon photographic white proportions and source-genus means were **46.93%** across 25,640 nominal photo-classifiable taxa and **50.01%** across 5,429 nominal genera (999-resample exploratory genus percentile 95% interval **48.97–51.02%**). This near-half figure is not fixed across sampling rules: limiting to genera with at least two or five photographed taxa gave **48.66%** and **46.55%**, and equal weighting across 207 populated 10° geographic cells gave **46.90%**. A previously completed six-variable genus-bootstrap and simultaneously adjusted long-term climate analysis found global white photographed hue positively associated with annual precipitation (**+1.99 percentage points per original environmental SD**, maxT95 **+0.55 to +3.42**) and atmospheric vapour (**+3.16**, maxT95 **+1.04 to +5.28**). Both simultaneous intervals span zero under genus×coarse geographic conditioning. These results establish **widespread photographic white colour and a geographically context-dependent environmental association**, not true adaptive equilibrium, metabolic pigment savings, selection, or flower-specific fitness. The central task is to test which geographically and taxonomically adjusted combinations of water, radiation, temperature, terrain, soil and sampling support predict white-versus-chromatic photographed composition.

**Keywords:** floral hue; white flowers; flower colour; macroecology; moisture; solar radiation; environmental filtering; observer bias; iNaturalist

## Central research question

Why is visible white floral hue *common* rather than universally replaced by chromatic hues? A globally near-half descriptive source fraction invites analysis of **different ecological environments permitting or favouring different floral appearances**, rather than a postulated global trait frequency set point.

Competing contextual mechanisms—not conclusions from the photo census—include:
- Reduced cost of producing some **visible pigments** under some resource regimes; WHITE is never assumed pigment-free, and UV-absorbing phenolics and other compounds may persist.
- Potential visual attraction or signalling gains of nonwhite floral pigments in some pollination environments; pollinator assemblages and visitation were not measured for this global source.
- Photo- and thermal-optical consequences of visible hue, and possible ultraviolet/oxidative protection from molecules contributing to floral colour or co-occurring pigments.
- Evolutionary/taxonomic history, spatial turnover and different habitat occupancy, which can create an apparent near-half global photographic average without selection maintaining that value.
- Taxon inclusion, spatial and observer effort, photography, organ identity and hue-classifier missingness, all of which can shift the observed average.

No hypothesized benefit or cost was observed as fitness or biochemical flux. The source can **reject or qualify simple universal environmental-dominance hypotheses**, but cannot by itself establish stabilizing selection on white/chromatic morphs.

## Materials and methods: immutable measured pixels and nested outcomes

**Photo source.** The immutable main sample contains 309,439 real RGB/ROI-v4 processed original photo IDs, 149,662 automatic four-hue ROIs and 149,653 original source-dated classifications. The old 100,543-photo cohort has 21,574 white among 45,953 classifiable photos, 46.95%; the much larger source is not independent replication of that earlier source. The separately site-season-targeted 256 photos are intentionally nonrepresentative and therefore excluded from the primary global fraction. Do not add partial fourth-wave source IDs until the actual RGB terminal and exact photo deduplication complete.

**Primary outcome.** Visible white among classifiable original ROIs, with yellow/orange, red/pink and blue/purple pooled as chromatic for the primary contrast. The three chromatic categories are also retained as distinct descriptive responses. Unclassifiable ROIs are not classified as nonwhite; all 309,439 attempts are the denominator for reporting classifiability and missingness. No photographed organ is assumed to be botanically validated as an open petal, no image visible hue assumed equivalent to biochemical anthocyanin expression.

**Different denominators.** Compute original photo-weighted white fraction; within nominal taxon photographed white fraction followed by equal taxon weighting; mean taxon fraction within genus then equal genus weighting; and equal 10° source-geographic cell summaries with at least 30 classifiable photographs. Cross-check genus minimum source taxa (1, 2, 5) and photo depth per taxon (1, 2, 5). The genus-replication percentile interval is descriptive uncertainty across sampled nominal genera, not a global adaptive-equilibrium test or full nonrandom source correction.

**Environment.** Reuse the frozen **original photo-ID** 309,439-annual WorldClim v2.1 six-axis model and create three independently SHA-verified Parquet tables joining original photographed visible hue, annual BIO5/BIO12/BIO14/BIO15/solar and **partially observed** 248,312-photo source soil/terrain data, with no imputation into unknown sites. Historic WorldClim solar values are kJ m^-2 day^-1 in an incorrectly named Wm2 column and must be divided by 86.4 for physical W m^-2. Long-term annual climate normals are 1970–2000, not measured weather of the photograph date. Seasonal WorldClim month normals and any separately verified observed-year monthly climate analyses belong in SI.

**Analysis hierarchy.** Test:
1. Is white consistently *common* under photo/taxon/genus/geographic weights and classifiability controls? There is no null requiring exactly 0.5.
2. How does the photographic white proportion vary across annual precipitation, atmospheric vapour, warmest-month temperature, dry-month rain, rainfall seasonality and solar radiation?
3. How does any water relationship change after controlling for elevation, latitude, longitude, coarse cell/biogeographic scale, soil and DEM landscape terrain, and complete-case site support? Never confuse a change in the measured environmental population with evidence for a causal added predictor.
4. Do four separate hues vary differently, indicating that chromatic colour is not biologically uniform and not equivalent to anthocyanins alone?
5. Do original photographic classifiability and source recruitment bias plausibly alter white-versus-colour proportions?

Use original matched-photo/genus bootstrap maxT24 water inference as **already fitted** and retain its geographic nulls. Environmental deciles are descriptive diagnostics, not new causal models. If fit/selection changes later, label them exploratory and report original source checksums.

## Source-authenticated main descriptive results

| Photographic weighting/selection rule | White fraction | Interpretation |
|---|---:|---|
| All 149,662 classifiable photo ROI outcomes | 47.21% | Photo-weighted original hue source |
| Equal nominal taxon weights, 25,640 source taxa | 46.93% | Per-taxon mean of classified photo white outcomes |
| Equal nominal genus weights, 5,429 genera | 50.01% | Mean of equal-taxon photographic genus proportions; 95% genus-resampling interval 48.97–51.02% |
| Genus-equal, genus at least 2 photo-classifiable taxa | 48.66% | More taxonomically replicated genera |
| Genus-equal, genus at least 5 photo-classifiable taxa | 46.55% | More stringent genus species depth |
| Equal source-geographic 10° cells, min 30 classifiable photos, 207 cells | 46.90% | Cell weighting, **not** area-/flora-weighted |

These findings are intentionally **not** described as a maintained 50:50 equilibrium. The genus sensitivity shows how an apparent precisely half proportion changes with taxonomic source support. The primary inference is that white is a very substantial global visible-photo hue category, with no universal chromatic dominance in these observed images.

**Existing long-term water result:** original 309,439-photo genus-equal global annual precipitation +1.99 pp per environmental SD, maxT95 [+0.55,+3.42]; mean atmospheric vapour +3.16 [+1.04,+5.28]. With coarse geographic genus conditioning, precipitation +0.90 [−2.63,+4.44], vapour +2.15 [−2.64,+6.95]. The conditional simultaneous intervals include zero; no single universal water-to-white causal response is established.

## Discussion and limits

This is a phenomenon that invites explanation, not evidence that visible pigment is invariably costly, that white flowers have superior reproduction, or that balancing selection maintains a predetermined global percentage. Many white petals have nonvisible pigments or flavonoid functions; visible classifications are not chemical trait quantification. Colour may covary with pollinator biology, climate, geography and shared plant history. The existing strong source-location and classifiability filters, photo-identification uncertainty and missingness are necessary interpretive boundaries. Distinct localities and regions can have different colours while their aggregate appears close to half.

**Figure plan:** (1) global photographic composition with photo, taxon, genus and source-cell weights and classifier-unknown denominator; (2) environmental rain/solar/temp gradients alongside fixed soil/terrain control coverage; (3) original geography-conditioned genus-equal water coefficient and uncertainty; (4) four-hue composition and classification bias; SI monthly/seasonal WorldClim normal contrasts and actual monthly weather when separately verified.

## Reproducibility

The tested image-to-photo-hue-to-annual-environment feature lake and CSV outputs are built in scripts/analysis/mois_primary_global_white_balance_feature_lake.py. Successful PR Actions verify exact photo/taxon/hue/environment identity, save 3 Zstd Parquet source tables and a SHA receipt, and a trusted main-branch Actions run publishes these derivatives as a **versioned immutable GitHub Release**. No copyrighted original image pixels are rehosted. Review rare-genus weighting, site/soil coverage and nonrandom classification before promoting ecological interpretation. This draft does not claim journal submission or external replication.
