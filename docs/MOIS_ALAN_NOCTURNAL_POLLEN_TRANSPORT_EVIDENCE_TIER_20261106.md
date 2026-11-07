# MOIS ALAN — confirmed nocturnal moth pollen transport versus broad GloBI links (2026-11-06)

## Result and scientific status

A **new source-authenticated rejoin** separated three different plant-level claims in the signed 61,967-row GloBI study metadata:

1. `moth_visit` — a moth family is recorded interacting with the flower. No time of day or effective pollen transfer asserted.
2. `moth_db_pollinates` — a GloBI study asserts a `pollinates` relationship. This is not necessarily directly measured pollen deposition/plant seed output. 11,945 of 16,518 records of this relation (72.3%) derive from titles containing the Swiss trophicCH network source; source metadata must not be interpreted as 11,945 independent effective-moth-pollination experiments.
3. **`walton_night_pollen`** — same moth-plant study claim associated with Walton et al. (2020), DOI [10.1098/rsbl.2019.0877](https://doi.org/10.1098/rsbl.2019.0877). The source explicitly examined pollen borne by nocturnally collected macro-moths in agricultural England; **pollen carriage is measured, but deposition, seed set, and obligatory nocturnal dependence are NOT established**. This is a stronger direct ecological time-of-day source than moth-family identity alone.

The original paper found pollen from 47 plant species; 43 binomial plant names were found in the GloBI extract for this study, accounting for 3,650 moth-family link records. The study's original field-site conclusions must not be applied as local interactions at the unrelated locations of MOIS photographs.

## Recomputed MOIS and same-year satellite support

Exact joined signed inputs: canonical original 1,032,882 observation UUIDs, WCVP accepted botanical names, independent Island High/Medium single-state flower colour, GloBI public interaction source and the source-locked VIIRS year-specific cell tables. Historical_main strict9-success photo source, 2019–2024, 60°S–60°N. Repeated photo observations are deduplicated in the original canonical frame. Independent species white describes the species, NOT necessarily the sampled individual flower.

| Evidence tier | 2019–2024 strict9 photo rows | 2019–2024 accepted species | Species with directly documented white |
|---|---:|---:|---:|
| Recorded moth-flower interaction | 59,310 | 1,523 | 193 |
| GloBI `pollinates` moth relationship | 17,226 | 387 | 48 |
| **Walton genuinely nocturnal moth pollen-transport source** | **3,481** | **40** | **5** |

**Same year, genuine VIIRS availability** (the 2020 satellite table contains 1,044 cells; the 2024 table contains 302 cells):

| Original-photo year and VIIRS year | Pollen-transport-associated MOIS photos | Associated accepted species | Directly documented white species | Strict9 photos from these white species |
|---|---:|---:|---:|---:|
| **2020 photos + 2020 VIIRS** | 99 | 35 | **4** | **13** |
| **2024 photos + 2024 VIIRS** | 20 | 14 | **3** | **7** |

These are **different geographic sample supports**, not a panel of 1,044 cells in both years. No species loss or year-on-year changes are inferred. Use of 2020 satellite values with 2024 photos is separately flagged as TEMPORAL MISMATCH and is excluded from a same-year ecological effect.

All five known-white Walton-associated species also have GloBI butterfly flower-visit records **somewhere in the full evidence base**; only one has butterfly records from the same Walton paper. Absence of a visitor record is not biological absence; with only 4/3 independently white species under exact-year VIIRS, no ALAN × botanical white × **effective nocturnal pollination** interaction is estimable.

## Source-bias discovery

Among all 61,967 GloBI rows, `visitsFlowersOf` totals 45,449, of which 26,558 (58.4%) reference iNaturalist source metadata. Of 16,518 `pollinates` relationships, 11,945 (72.3%) have a trophicCH source title. Their number of rows is **not number of independent experimental studies** or a pollination-service frequency. The original GloBI cache is live and some source families were capped, another boundary for population-level prevalence.

The ecological implication is **not** that the proposed negative ALAN effect is absent. More general interaction records cannot substitute for plant-specific nocturnal pollen transfer/deposition and night-versus-day documented service, especially when botanical white records are sparse.

## Reproducibility and privacy/rights

Executed locally on actual previously signed source artifacts. Source SHA256 (hex):

- `GloBI_stage2_signed_extension.csv.gz`: `a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994`
- `canonical_observations_climate.csv.gz`: `f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c`
- `mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz`: `dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2`
- `mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz`: `9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3`
- `direct_species_trait_ledger.csv.gz`: `7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577`
- 2020 VIIRS: `3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3`
- 2024 VIIRS: `1c5f03e5a7c355a0f7045133140c5cf2e39a5802c548f89d1c2f66c8a235db5f`

Independent new outputs copied byte-identically to `results/mois_alan_nocturnal_pollen_tiers_20261106/`:

- `evidence_tier_MOIS_photo_coverage.csv` SHA256 **2a9e3f076477f1b408af0019084c89447c8accfe556bdff869bf4ec60e1d7528**
- `documented_white_overlap_VIIRS_support.csv` SHA256 **2e1d570a702b0a5bc8477e372daef4c5134f38c27fcaafee0f8cf60dfdf0092b**

Both original local CSV files and GitHub source blobs verified byte-identical using Git blob SHA-1. Five additional synthetic unit tests passed locally. Reproducible execution script and the test file accompany the user-facing ZIP output; they have not been committed as source files to this PR.

No source photographs fetched, inferred or recoloured. No unknown converted to white/nonwhite. No original scientific estimates overwritten, no previously forbidden MNAR sensitivity run, no private Island species-level trait records redistributed.

## Prioritized evidence acquisition roster (compiled locally, not publicly uploaded)

Of the **40** photographed accepted species from the nocturnally collected moth pollen-transport study, **5** have independently documented single-state WHITE flowers, **3** have independently documented single-state NONWHITE flowers and **32** have no eligible single-state High/Medium documented botanical flower colour. Those **32** otherwise eligible MOIS species account for **2,707 original strictROI9 successful photos** during 2019–2024, including **76 photos on original 2020 VIIRS cells** and **12 photos on actual 2024 VIIRS cells**.

An additional species-name only literature-retrieval priority queue was created for the user-facing local reproducibility bundle; it ranks those 32 species by existing same-year satellite-matched photo support and existing ROI9-success photo opportunity counts. It does NOT assign flower colour or true local nocturnal pollination from photos. The queue is held outside the public PR because the original independent individual species trait ledger has third-party redistribution constraints. Five new local unit tests passed, source SHA locks verified, and the GitHub public aggregate CSVs are byte-identical to the independently run numerical outputs.

## Next focused step

Source an *independent species-level* time-of-day / pollen-transfer candidate panel from explicitly nocturnal source studies, retaining study location, event-time verification, plant accepted ID, pollen carriage vs verified deposition/seed production, and presence of **documented day-time interactions** as separate fields. The [Diniz et al. (2025) primary pollen-metabarcoding dataset](https://datadryad.org/dataset/doi%3A10.5061/dryad.sxksn03gp) is a potential second source; some plant links are ASVs identified only at genus/family rank, so they must NOT be automatically promoted to species identities. Recheck botanical white and exact-year VIIRS support before any ALAN coefficient. A source with 500 plant-association records but zero species-level accepted names has no usable denominator for the targeted inference.
