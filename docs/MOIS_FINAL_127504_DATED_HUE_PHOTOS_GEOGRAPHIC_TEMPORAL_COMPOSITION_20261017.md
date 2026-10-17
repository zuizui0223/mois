# MOIS: full 127,504 dated photographic hues — taxon turnover and the locality-scale limit

**Complete actually verified source, 17 October 2026.** The original global source contains 267,425 genuine photographic ROI attempts and 129,795 automatic four-hue labels. The original photographs, hue classifications, original taxon IDs and worldwide six-axis moisture/maxT analysis are unchanged.

## Date-source audit: two verified but different provider snapshots

The successfully completed [official public monthly S3 scan](https://github.com/zuizui0223/mois/actions/runs/39009539731) restored a valid 2026-09-27 observation date for **99,732 of the 100,543 historical original photo IDs**. The [separately completed 9,600-original-photo exact numerical observation-ID and photo-ID public API audit](https://github.com/zuizui0223/mois/actions/runs/39008584811) supplied **85 additional dates** for otherwise S3-missing original photo IDs; **46** also have original four-hue photo classifications. No photos were newly classified by either source-repair workflow.

Of 9,600 API source photos, **9,512 had a date in both source snapshots**; one original observation has discrepant current API and earlier S3 dates. We retain the older 2026-09-27 official snapshot date as canonical, preserving both versions and the disagreement flag. There are **726 historical original photo IDs still without a valid source date**. No date or flower-petal status was imputed.

| Measured original photo/date result | Count |
|---|---:|
| Historic 100,543 photographs with date from official S3 or exact API-only fallback | **99,817** |
| Historical original hue-labelled photos with date | **45,637 of 45,953** |
| All 267,425 original-photo attempts with dated original/public observation | **261,578** |
| All real original photos with BOTH observation date and one automatic four-hue ROI label | **127,504** |
| Distinct original nominal source taxa in these dated, hue-classifiable photographs | **25,257** |

Of the 25,257 dated-and-classifiable source taxa, **15,926 (63.06%) are represented by one original photographic hue**. Among the historic source subset only, 16,411 of 24,984 dated-and-classifiable source taxa (65.69%) are single-photo taxa. This dramatically improves the source-wide taxonomic scope compared with 81,867 originally dated coloured photo labels but does not create longitudinal replication of photographed petals.

## Strong *photographic* seasons but very low taxon continuity

Photographic ROIs automatically labelled white versus three nonwhite coarse colours have the following genuine original-date source composition:

| Extra-tropical source | Summer original white/photo labels | Winter original white/photo labels | Unadjusted summer−winter white-photo difference |
|---|---:|---:|---:|
| North | 18,620/38,121 = 48.84% | 1,850/4,059 = 45.58% | **+3.27 percentage points** |
| South | 2,996/5,827 = 51.42% | 1,589/4,082 = 38.93% | **+12.49 percentage points** |

Yet only **562 of 7,015 Northern** source taxa and **320 of 4,139 Southern** source taxa are photographed with classifiable colour in both local summer and winter (photographic taxon Jaccard 0.080 and 0.077 respectively). This is photographic source season/taxon turnover, not a proven field phenology or flowering absence.

An exactly additive full-photographic-census symmetric decomposition of the **Southern +12.49 points** yields: **+3.41 points** within-common-taxon photographed white-rate differences; **+10.60** from changes in their relative photographic frequencies; **+2.67** from summer-only photographed-taxon support and **−4.19** from winter-only support. The corresponding Northern +3.27 decomposes into +2.52, −1.11, +3.04, −1.19 respectively. These terms are source accounting, **not causal pigment physiology, seasonal selection or pollinator fitness effects**. Absent opposite-season photographs are never treated as evidence of an absence of flowering.

## Geographic resolution and same-season-year comparison

Strict same-named-taxon, same 0.25-degree geographical cell, same meteorological season-year and **two distinct observation IDs per summer/winter** yields only **2 Northern groups (8 original hue-labelled photos) and 0 Southern**. Hence a general within-taxon local seasonal pigment response remains unidentified.

A different target — same-cell and same-year **photographed community colour** across *possibly different taxa* — is more tractable but still source-selected. Requiring at least five different original observation IDs and two different original nominal taxa in each season gives **16 Northern and 18 Southern 0.25-degree cell×season-year comparisons**. At 1 degree the corresponding counts rise to 101 and 84, but the broader cell merges ecologically different places.

For the finest 0.25-degree matched cell-years (exploratory equal-cell-year average, 999 geographically clustered bootstrap percentile 95%), Northern summer−winter is **+5.04 points [−8.96,+20.55]** and Southern is **−4.63 points [−16.45,+6.55]**. Both intervals span zero. The apparent unpaired South-global +12.49-point summer-white photographic peak **does not recur as a robust positive signed local, same-year contrast** at the finer inspected source geographic grain. Analyses at coarser cells/alternate minimum photo counts are exploratory, not independent familywise confirmation or causal light/water effects.

These photo cell-year comparisons are not identical plants, identified breeding populations or independent observers; the grid is still tens of km wide. Flowers, floral bracts, fruits and ROI errors were not independently botanically verified as open petals in this date-only work, nor did it fit actual observed-year solar radiation, temperature, precipitation, VPD or climatic water deficit.

## Manuscript and source boundary

This provides a specific interspecific macroecology result: **strong photographed white-hue composition gradients vary by season and taxonomic sampling, and lose signed certainty when constrained to approximately co-located yearly photographic assemblages**. It is not proof of a universal flower pigment response, fitness cost–benefit balance or evolutionary 50:50 equilibrium. The annual MOIS genus-weighted water associations remain conditional on source and geographic controls and are not replicated by these temporal results. FCP's within-species W/D and 50-km photographic phenotype-space analysis remains separate.

Reproduction: the source-verified final original-photo S3/API date-merge script, GitHub Actions reanalysis and 90-day original-photo source artifact are supplied in the same PR. Original 100,543-row photo identities stay out of public Git; only the aggregate source decision and ecological manuscript report are committed.
