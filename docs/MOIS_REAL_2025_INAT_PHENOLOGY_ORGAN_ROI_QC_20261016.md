# MOIS | Real 2025 seasonal photograph × plant phenology quality audit

**Owner: zuizui0223/mois.** 16 October 2026. The source photographs are frozen from the original 27 September 2026 iNaturalist snapshot; public live observation annotations were retrieved on 16 October and may change subsequently. This is source quality, NOT a biological white-flower/temperature/moisture effect.

## Original scientific evidence

- [FCP real 2025 Asclepias photo measurement](https://github.com/zuizui0223/fcp/actions/runs/38909164068): 64 independent original observations and observers, 32 in summer and 32 in autumn. Source ROI/observed-on SHA256: 1b16cdabd2936567dee87b8344019f40f315600c6510dfc1f25a03835024b6a0.
- [MOIS measured 2025 monthly climate](https://github.com/zuizui0223/mois/actions/runs/38914456992): all 126 actual TerraClimate v1.1 CF-decoded site-variable climate time series passed, yielding 40 photo-relevant source year-month×grid units across 21 climate grid cells. The 1991–2020 same-calendar-month deviations are also available.
- [MOIS common-geocell seasonal QC](https://github.com/zuizui0223/mois/actions/runs/38917905669): only 7/21 roughly 4-km climate cells support photographs in both seasons. These contain 15 summer photo attempts (2 four-state automatic ROI labels after excluding one same-ROI white/yellow floral-organ mix) versus 24 autumn attempts (4 automatic labels).
- [Mois original-observation live phenology audit](https://github.com/zuizui0223/mois/actions/runs/38919740652): the original 64 UUID, observed calendar dates and photo IDs all match the currently public iNaturalist V2 UUID and V1 annotation records. These are *live current observations*, not historical frozen phenology.
- [Mois source-authenticated automatic ROI × live observation annotation contingency](https://github.com/zuizui0223/mois/actions/runs/38920242498): all 64 photos joined, 20 synthetic code/receipt tests passed. Automatic ROI colour must not be called independently verified open-flower colour.

## Live observation annotations, not photo-specific botanical verification

| Annotation of source observation | Summer (32) | Autumn (32) |
|---|---:|---:|
| Flowering | 2 | 0 |
| Fruiting | 0 | 4 |
| No Evidence of Flowering | 1 | 2 |
| Absent or empty phenology annotation | 29 | 26 |

The iNaturalist Plant Phenology controlled attribute is 12. Values 13, 14, 15 and 21 refer to Flowering, Fruiting, Flower Budding and No Evidence of Flowering respectively. Confirmed with [iNaturalist annotation code reference](https://forum.inaturalist.org/t/downloading-annotations/2200) and [official definitions](https://www.inaturalist.org/posts/86088-annotation-tooltips). No evidence of flowering is NOT a rigorous diagnosis that this particular photographed ROI contains no flower. Missing annotations are UNKNOWN, not biological absence. Flowering and fruiting structures can co-occur within one observed organism or observation.

| Automated image outcome (NOT validated petals) | Summer (32) | Autumn (32) |
|---|---:|---:|
| Automatic white | 7 | 2 |
| Automatic chromatic | 2 | 4 |
| Same original ROI mixes white/yellow floral organs | 3 | 0 |
| Unclassifiable ROI | 20 | 26 |

Both summer observations tagged Flowering had automatic white ROI labels, not independently confirmed open petals. One of the four autumn observations tagged Fruiting had an automatic white ROI label; the other three were unclassifiable. This creates a concrete organ-stage quality concern, not proof of a false colour: fruit and open flower can co-occur.

**Confirmed open-petaloid flower images: 0 reviewed (all 64 still UNREVIEWED).** This must never be phrased as no actual blooming flowers. The [two-botanist season-blind review form](https://github.com/zuizui0223/mois/actions/runs/38915040401) is available. Photographs with uncertain stage, fruit, bract, mixed within-flower colours or undecodable ROI remain explicit in the funnel, without colour-based replacements.

## Identifiability and next decision

The first 64-source-photo phase cannot estimate a true seasonal white-flower frequency, a universal moisture/pigment relationship, a pigment benefit/cost equilibrium, or a causal temperature/VPD response. In particular, month, climate, observer photo success and plant anthesis confound one species/site-year. The only 7 spatially shared cells have just two summer and four autumn automatic usable ROI labels, with zero independent botanical confirmation.

Next actions: independently verify actual open flowers and petal/tepal ROI in individual original photographs; increase genuinely season-overlapping photos while preserving original sampled and unclassifiable denominators; conduct observer/site/year sensitivity; fit observed-month anomalies and stable geography/soil/elevation as separate terms only when sufficient source replication is achieved. This is a MOIS seasonal programme distinct from historical MOIS genus-centred 100,543-photo BIO12 findings and from FCP W/D polymorphism inference.
