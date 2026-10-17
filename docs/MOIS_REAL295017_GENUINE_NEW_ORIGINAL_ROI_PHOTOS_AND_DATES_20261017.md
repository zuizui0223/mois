# MOIS: 295,017 genuinely measured original photographic ROIs, 140,247 with dates and hues

**Source scope:** 17 October 2026; a photographed global interspecific flora/observation process dataset, not FCP's separate *within-species* D/W 50-km photo polymorphism estimates.

## Combine all independently finished source pixel work, not future photo opportunities

The previous signed MOIS 267,425-photo source used **79 completed FCP wave4 ROI worker partitions** and the first 4,096-image MOIS source wave. The [real FCP 180-partition reimport](MOIS_REAL286826_FCP180_ORIGINAL_HUE_PIXELS_WATER_20261017.md) produced **286,826 original distinct photos /139,405 actually classified four-hue ROIs**, with immutable verified original hue and six-anchor genus-999 global rainfall/vapour refit. Independently, the [successful 3×4,096 genuine MOIS ROI worker waves](https://github.com/zuizui0223/mois/actions/runs/39013017991) generated **260,600 distinct original photo IDs /125,517 classifiable photo hues** from the original 248,312-photo source.

Both archives share exactly **252,409 original photo IDs**. Two independently finished new MOIS 4,096-image waves are not present in the FCP180 archive except for **one literal reused original photo ID**. That original photo has the same taxon ID, observed day and hue label in both measurements, but its observation identifier is represented as a source UUID in one workflow and a numeric iNaturalist observation ID in the other. It is **one photographed source**, counted just once, not a new plant or conflicting white/chromatic hue. No other original taxon, hue or date conflicts were found in source-comparison.

| Signed genuine source outcome | Prior 79 FCP partitions | Finished FCP180 only | All finished FCP180 + 3 MOIS waves (photo-ID union) |
|---|---:|---:|---:|
| Original photo IDs (pixels actually measured) | 267,425 | 286,826 | **295,017** |
| Original photo ROIs classifiable into 4 hues | 129,795 | 139,405 | **142,538** |
| Measured distinct photo ID increase from prior | — | +19,401 | **+27,592** |
| Measured 4-hue photo increase from prior | — | +9,610 | **+12,743** |

All 295,017 rows retain the original photographic `photo_id`, taxon identity, raw `observed_on`, frozen ROI colour and technical unclassifiable outcomes. Original colour labels are not refitted or imputed from annotation/season/climate. There are still **76 unfinished FCP worker partitions** in the FCP180 snapshot; the planned total of 49,800 new FCP photo IDs must not be confused with a fully measured dataset.

## Exact independent date provenance

The [complete original historical iNaturalist S3/API date source](https://github.com/zuizui0223/mois/actions/runs/39012110980) separately restored valid publicly documented observation dates for **99,817 of the 100,543 historical original photos**, without rewriting their original missing date column. Both newly added MOIS source waves carry their actual original `observed_on` day in the photo metadata, and the extra 19,401 finished FCP ROI source photo IDs also have source-recorded observation dates.

| Verified source photographic outcome | Signed current total |
|---|---:|
| Distinct original genuine measured photos | **295,017** |
| Real original photo labels assigned 4 hues | **142,538** |
| Original photos with a source/verified-public observation date | **289,170** |
| Original photos with **both** date and automated 4-hue ROI classification | **140,247** |
| Date-bearing original photos whose hue remains technically unclassifiable | **148,923** |
| Distinct original named taxa among date- and hue-classified source records | **25,545** |

Thus the actual photo hue + date sample increases from **127,504** to **140,247** (**+12,743 classified and dated original photograph IDs**). The original historical 726 unresolved date entries, additional earlier unknown-date photo measurements, unknown technical hues, flower-organ identification and actual observed-year weather remain explicit, rather than inventing any 2026 seasonal data.

## Current source photographic summer and winter

Photo date-only hemispheric descriptive counts (automatic white ROI labels / all classified four-hue source photo ROIs, *not botanically confirmed open flowering*) from local meteorological seasons:

| Source photo summer/winter | Summer white / classified | Winter white / classified |
|---|---:|---:|
| North extra-tropics | 18,931 /38,779 (**48.82%**) | 1,900/4,130 (**46.00%**) |
| South extra-tropics | 3,122 /6,086 (**51.30%**) | 1,659/4,266 (**38.89%**) |

These are unpaired original photo/taxon/observer/geography-seasonal descriptions. Merely adding independent global photographs does **not** establish same-named-taxon, same-locality/year pigment change. In the present source, matching taxon × **0.25° cell × one meteorological season-year**, with ≥2 separate source observations each summer and winter, remains **two Northern groups and zero Southern groups**. The expanded source is useful for global *among-named-taxon and community photographic composition* inference; source photo-organ stage, observed-year monthly climate and observer effort remain necessary before testing causal sunlight or water effects.

## Reproducibility and manuscript boundary

The [new source-locked GitHub Actions workflow](../.github/workflows/mois-genuine295017-original-photo-pixels-and-dates.yml) downloads ONLY already successful original-pixel and date-recovery artifacts, verifies three independent source photo SHA256s, photo ID/taxon/hue overlap, one numeric-vs-UUID observation ID discrepancy and exact photographic/temporal counts. It publishes original **295,017-row photographic source union**, **275k+ date-augmented photo table** and a provenance/seasonal receipt as an Actions artifact rather than committing individual photo IDs to a public repository.

The FCP180 source separately executed the original historic MOIS genus-equal **999-bootstrap maxT** white-versus-rainfall/vapour model, preserving the broad positive direction but leaving genus×geographic-cell conditioned slopes unconfirmed. This document reports the subsequent newly measured **295,017-photo union**, which has **not yet** had the full original 999-genus six-anchor model rerun. Do not substitute the preceding 286,826-photo model coefficients for unperformed 295,017-photo analysis. The original 100,543-photo historical macro benchmark and FCP's intraspecific polymorphism analyses remain independent, unmodified scientific estimands.
