# Recover original-date candidates in MOIS's previously undated historical 100,543 photographs

**Source-conserving real-data pilot, 17 October 2026.** The user requested expanding the date-bearing *actually classified photo* dataset, rather than treating millions of source photo-ID opportunities as verified flowers. The main MOIS [267,425 genuinely processed original-photo union](https://github.com/zuizui0223/mois/actions/runs/39002957961) has 129,795 automatic four-hue labels but only **81,867** with dates. Critically the full **original historical 100,543-photo MOIS block** has `observed_on` blank in the union, while the unique **original iNaturalist `observation_id` and `photo_id`** survive. All 100,543 original observation identifiers were present and unique in our source audit. Thus missing original date in an analysis export need not imply no date exists in the original public observation.

## Concrete evidence, not hypothetical upstream data

A direct current public iNaturalist observation endpoint check on **four original historical photo/observation pairs** independently recovered observed calendar dates and exact **photo-ID membership** within each returned observation, with originally named taxon IDs unchanged:

| Historical original photo ID | Original observation ID | Public observed_on | Original taxon ID still matches |
|---|---|---|---|
| `253367174` | `147261418` | 2022-10-18 | yes |
| `84410948` | `53081645` | 2020-07-14 | yes |
| `505043460` | `266179829` | 2025-03-18 | yes |
| `1860054` | `1501786` | 2015-02-05 | yes |

The fourth photo is not the first public API photo on its current observation, which proves that matching the original **photo ID across the entire photo list** is essential. Public source API: `https://api.inaturalist.org/v1/observations?id=147261418,53081645,266179829,1501786&per_page=4`. This verification concerns the **current public iNaturalist record**, not a preserved 2026-09 original API response, and can change due to edits, record deletion or taxonomic revision.

This independently establishes that the historical date column can be **enriched from verified actual observation records** rather than artificially treating historical photo dates as absent or inferring the month from climate. It is not evidence that each image contains open flower petals: one currently annotated observation with a `Flowers` tag would still require photo-specific botanical review, as public current annotations are associated with the observation, not verified independently with the target photographed organ.

## New bounded source-frozen 200-photo feasibility workflow

- Deterministically select exactly **200** real original photos **blind to colour, environment, observed date and taxon identity** by SHA256 of the original `photo_id`. Keep the original 100,543-photo sampling frame and *every* selected original ID; do not replace missing, deleted or ambiguous API observations with more favourable photos.
- Authenticate the entire **267,425-photo** original source ZIP and SHA256 `fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9`, verify the **100,543** dated-empty historical originals and one-to-one original photo/observation identity *before* contacting the provider.
- Query the official current iNaturalist public observation endpoint with **four batched requests of 50 observation IDs each**, with a mandatory >=1.25-second gap, retries with bounded backoff and strict response identity/pagination validation. This respects iNaturalist's recommended <=60 requests/minute and is far below the daily usage limit; it is only a **pilot**, not an attempted 100k API sweep.
- Require the **exact original `photo_id`** to appear within the correct original `observation_id` returned current public photo list; only then attach its current `observed_on` as a *candidate*, with original/current taxon differences and retrieval date explicitly preserved. Do not overwrite the frozen original 267,425 source photo union, original `inat_taxon_id`, historical hue classification, geo-location or original source SHA.
- Retain API nonresponses, observation-photo mismatches, invalid dates and missing current phenology as explicit known *unknowns*, not not-flowering or chromatic hues. Currently tagged `Flowers`, `Flower Buds`, `Fruits or Seeds` etc. remain **observation-level current public tags** and are never photo-specific verified flowering labels.
- The published one-shot Actions run should output a **200-photo sampled date/annotation QC CSV**, aggregated source/matched/missingness counts, source-ID SHA and current public API-response SHA for every batch. Do not fit any global white-season relation until actual successful current observation-photo date recovery is quantified.

Script: `scripts/acquisition/mois_historical_original_observed_on_inat_pilot_20261017.py`, source QA tests: `tests/test_mois_historical_original_observed_on_inat_pilot_20261017.py`, CI: `.github/workflows/mois-original100543-historical-200-real-inat-dates.yml`. The [iNaturalist API recommended practices](https://www.inaturalist.org/pages/api+recommended+practices) describe provider rate and bulk-export recommendations. For eventual full-scale historical date enrichment, prefer a compliant source snapshot or bulk export where feasible, and retain the same photo/observation-ID verification and original source IDs.

## Scientific impact and inference boundary

The existing MOIS temporal analysis currently **cannot** use its historical 100,543 photos because the analyzed export has no original `observed_on`. Reconstructing current public observed dates for even a defensible fraction would permit testing seasonal photographic composition across more of the original worldwide taxonomic breadth, not simply adding more high-depth FCP taxa. An API recovered observation date is not automatically the original-photo capture time and not verified flowering; independently matching the exact source photo ID, taxon and annotation identity remains essential.

This pilot does **not** update the 267,425-photo historical moisture estimates, the original 81,867-date-hue seasonal numerator, the 100,543-photo archival scientific conclusions, or FCP within-species W/D/local-spatial polymorphism results. It first quantifies the genuine missing date repair opportunity with a frozen photograph denominator and empirical API source QA.
