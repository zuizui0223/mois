# MOIS — moisture × global photographed floral colour

**Active research scope (2026-11-08): worldwide variation in the composition of photographed floral colours along moisture-related environmental gradients.** The scientific focus is water availability (precipitation, atmospheric water vapour, humidity/aridity), including how water-related patterns vary across taxa, space and observation season.

## Main biological question

How does the photographed floral-colour composition of plant assemblages vary along global moisture gradients, and which part of this pattern is associated with the climate context versus the taxa and observation seasons represented in the photographic sample?

- **Primary exposure:** moisture-related climate. Models may control temperature, solar radiation, elevation, geography, soil, and source/seasonal observation composition, but these are not separate current lead programmes in MOIS.
- **Primary outcomes:** source-measured white fraction and nine-component floral-region colour composition. These are image ROI proxies, **not independently verified botanical petal colours** or unbiased worldwide species prevalence.
- **Existing empirical signal:** the executed original-photo climate refit in [MOIS draft PR #196](https://github.com/zuizui0223/mois/pull/196) found the photo-assemblage white-label rainfall coefficient **+1.314 percentage points per source-standardised SD** (nominal 95% CI +0.163 to +2.500) under technical-response weighting, using **893,253 climate-complete original observations / 350,081 strict ROI9 successes**. This is an observational association, not a within-species or causal response.
- **Sampling frame:** the original upstream index contains **93,209,351 photo candidates**, but that is *not* a set of 93 million image-measured blooming petals. A signed four-arm union reached **1,181,146 distinct attempted original photos and 451,643 strict nine-palette numeric image ROIs**; see [original source release](https://github.com/zuizui0223/mois/releases/tag/mois-original-1181146-exact-photo-crossarm-451643-strict9-taxongap-20261027-v1).

## Research programme boundary

**Island × artificial light at night (ALAN), VIIRS × flora, nocturnal flower visitors, and associated spatial/temporal analyses were moved to [zuizui0223/alan](https://github.com/zuizui0223/alan)** on 2026-11-08 from MOIS draft PR #199. All 179 source files from that PR were copied exactly, with the 54 historical Actions YAMLs made non-executing in ALAN's archive. Continue that programme in ALAN; do not append new ALAN research to MOIS.

Historical MOIS releases, earlier workflows, pre-existing geographic comparisons, and associated documentation are retained for scientific traceability, not recast as current independent research objectives. The original detailed historical README is preserved at [docs/MOIS_README_BEFORE_ALAN_SPLIT_20261108.md](docs/MOIS_README_BEFORE_ALAN_SPLIT_20261108.md).

## Where to find the moisture results

- [Global moisture manuscript](docs/MOIS_GLOBAL_MOISTURE_MANUSCRIPT_20261014.md)
- [Global floral colour, moisture and seasonal sampling draft](docs/MOIS_GLOBAL_FLORAL_COLOUR_MOISTURE_AND_SEASONAL_SAMPLING_MANUSCRIPT_20261017.md)
- [Cross-species moisture evidence and missingness-aware refit — PR #196](https://github.com/zuizui0223/mois/pull/196)
- [Original-source data provenance](docs/MOIS_20261025_ORIGINAL_9PALETTE_CLUSTER_AND_ROI_FAILURE_DECISION.md)
- [Migration and scope split record](docs/MOIS_ALAN_SPLIT_20261108.md)

## Important inference boundary

A measurable photographic colour–moisture association does **not** demonstrate pigment-level stress adaptation, evolutionary trade-offs, flower fitness, native wild-flora frequencies, or a causal humidity effect. Technical ROI nonresponse, photography-selection effects and plant identity remain central observation-process considerations.
