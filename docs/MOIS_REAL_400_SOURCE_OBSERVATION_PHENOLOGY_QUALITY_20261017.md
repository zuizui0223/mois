# MOIS | 400 source-photo observations: real iNaturalist phenology annotation availability

**Date of live public API read:** 17 October 2026, 01:05:48 UTC.
**Scientific scope:** ancillary observation-process bias in the MOIS seasonal global flower-colour × time-resolved temperature, light and moisture research. This is NOT a botanical photograph anthesis study outcome, and not the historical interspecific GEB paper or FCP intraspecific W/D paper.

## Authentic 400-original-photo source

The selected original image opportunities were frozen **without looking at flower hue, phenology or climate** in [the five-taxon 400-source selection](https://github.com/zuizui0223/mois/actions/runs/38929402279). Original 2026-09-27 photo selection SHA256:

b961b50ee6efe9ec08a5575a8eb1e34b1edf7b7fcb426226fb11ee901cfd3329

The selected 400 photo IDs correspond to **400 independently sourced iNaturalist observations and 400 different source observers**, with 40 photos from each of two seasons for five species at one selected original 0.25-degree geographic site/year per species. These are source photographic opportunities, not confirmed flowering photographs.

[MOIS live annotation audit #38989216871](https://github.com/zuizui0223/mois/actions/runs/38989216871) passed **20 synthetic validity tests**, then read exactly those 400 frozen observation UUIDs using the current public iNaturalist v2 UUID lookup and v1 annotation endpoint. **400/400** UUIDs, original observed_on dates and original photo IDs were independently matched in the currently public records. No photo was added, removed, or replaced, and no image pixel was acquired.

- [Authentic 10-group aggregate data](../results/mois_400_live_observation_phenology_20261017/source_siteyear_season_annotation_coverage.csv), SHA256 34f6997b85d7d98c355991ea79e9d7efe7194322da53857062ed30e7fe3f52d0.
- [Source and inference-limit receipt](../results/mois_400_live_observation_phenology_20261017/verified_source_400_annotation_coverage_receipt.json).
- [Complete saved, date-stamped read-only annotation and SHA results](https://github.com/zuizui0223/mois/actions/runs/38989216871/artifacts/12084148144), NOT committed here as 400 individual user/photo records.

The current API records were accessed in October 2026 and **must not be falsely described as 2026-09-27 immutable phenology annotations**. Later edits can change the current API response, which is why the source result and timestamp were archived.

## Major measurement result: phenology annotation missingness is taxon/recording-context dependent

Of **400 original photographed observations**, **280 (70%)** lack an interpretable Plant Phenology annotation and **120 (30%)** have at least one current annotation. A *present annotation is not proof of an open flower in that individual photograph*.

| Original sampled taxon | Local seasons represented | Original photo attempts | Current plant-phenology annotations | No interpretable annotation |
|---|---|---:|---:|---:|
| *Datura discolor* | Winter / Spring (2023 season-year) | 80 | 78 | 2 |
| *Aralia californica* | Summer / Autumn 2025 | 80 | 26 | 54 |
| *Linaria vulgaris* | Summer / Autumn 2023 | 80 | 9 | 71 |
| *Asclepias syriaca* | Summer / Autumn 2025 | 80 | 6 | 74 |
| *Achillea millefolium* | Spring / Autumn 2023 | 80 | 1 | 79 |
| **Total** | Five source site-years | **400** | **120** | **280** |

The annotation completion contrast is **78/80 (97.5%) versus 1/80 (1.25%)** in *Datura* and *Achillea* respectively. This association can reflect annotator attention, taxon recognisability, place/year composition and phenological organ visibility. It does **not** prove differential natural flowering probability.

Among 120 annotated source observations, the currently assigned **single interpreted plant-phenology values** are: 32 Flowering, 54 Fruiting, 6 Flower Budding, and 2 No Evidence of Flowering. **26 additional observations carry multiple different plant-phenology values**, represented as a separate category rather than being falsely counted as unique flowering or fruiting pictures. A plant/observation may contain flowers, buds and fruits simultaneously, so two phenology values need not be contradictory biology.

### The actual 10 frozen species×season blocks

| Species | Source local season | Photo attempts | Flowering tag | Fruiting tag | Budding tag | No-evidence tag | Multiple tags | Unknown annotation |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| *Asclepias syriaca* | Summer 2025 | 40 | 0 | 2 | 2 | 1 | 1 | 34 |
| *Asclepias syriaca* | Autumn 2025 | 40 | 0 | 0 | 0 | 0 | 0 | 40 |
| *Aralia californica* | Summer 2025 | 40 | 5 | 4 | 4 | 0 | 8 | 19 |
| *Aralia californica* | Autumn 2025 | 40 | 0 | 5 | 0 | 0 | 0 | 35 |
| *Datura discolor* | Winter 2023 season-year | 40 | 5 | 26 | 0 | 0 | 9 | 0 |
| *Datura discolor* | Spring 2023 | 40 | 13 | 17 | 0 | 0 | 8 | 2 |
| *Achillea millefolium* | Spring 2023 | 40 | 0 | 0 | 0 | 0 | 0 | 40 |
| *Achillea millefolium* | Autumn 2023 | 40 | 0 | 0 | 0 | 1 | 0 | 39 |
| *Linaria vulgaris* | Summer 2023 | 40 | 4 | 0 | 0 | 0 | 0 | 36 |
| *Linaria vulgaris* | Autumn 2023 | 40 | 5 | 0 | 0 | 0 | 0 | 35 |

**Do not pool the winter, spring, summer and autumn rows as if seasons were randomly assigned across taxa.** In this source only *Datura* contributes winter observations, and its public annotations are also unusually complete. That is taxon×site×time observation heterogeneity, not proof that winter flowers fruit more or that light/moisture controls colour.

### Month-level source opportunity asymmetry

The original pre-colour 400 photo source also demonstrates uneven months within each local season:

- *Asclepias* autumn 2025: **36/40 in September**, just 2 in October and 2 in November.
- *Aralia* autumn 2025: **34/40 in September**, 3 each in October and November.
- *Datura* winter 2023 season-year: **5/40 in December 2022, 18/40 in January 2023, 17/40 in February 2023**. The December 2022 photos belong to winter season-year 2023, not a separate winter 2022.

The seasonal photo opportunity sample has a real month-within-season imbalance. This cannot be corrected by treating a single seasonal mean air temperature or radiation as the actual climate of every photographed month.

## What the annotation evidence does and does not establish

The controlled iNaturalist Plant Phenology attribute is 12; current values 13 Flowering, 14 Fruiting, 15 Flower Budding, and 21 No Evidence of Flowering. These tags are assigned to source **observations** rather than independently verified exact photo petal/tepal ROIs. Missing values are **unknown**, not "not in bloom", and an observation showing fruit can also contain one or more flowers.

Consequently a filter keeping only observations currently tagged "Flowering" would radically overrepresent *Datura* relative to *Achillea*, change the estimand and confound photo-quality selection with botanical recognition and calendar season. It must **not** be used as the sole primary flowering-photo denominator or as a negative label for photos lacking tags.

[The canonical MOIS source-specific two-botanist independent review workflow](MOIS_400_REVIEW_BROWSER_HOWTO_20261017.md) and [private coordinator/adjudication gate](MOIS_400_REVIEW_ADJUDICATION_PROTOCOL_20261017.md) remain the appropriate next measurement step for every frozen original photograph. **Botanically reviewed open-petaloid original images: zero reviewed yet (400 awaiting review).** This is not a claim that the 400 photographs contain no actual open flowers; it simply records that no independent botanical judgements have been returned.

Only afterward may the frozen FCP ROI-v4 method classify eligible original flower organ pixels and MOIS associate verified photographic white/nonwhite states with true observed-month TerraClimate temperature, solar radiation, precipitation, water deficit and VPD; comparison requires within-species/site-year/season replication and observer/photo/classifiability sensitivity. Original source photo ID, unknown annotation, reviewer disagreements, stage failures and all denominators must remain present.

**No white-flower percentage or climate effect in the new five-taxon seasonal cohort is estimable from these metadata-only results.** The previously published-internal MOIS global between-species genus-centred BIO12 coefficients and FCP intraspecific W/D spatial results are independent observational estimates, not validations of the 400-source-photo seasonality.
