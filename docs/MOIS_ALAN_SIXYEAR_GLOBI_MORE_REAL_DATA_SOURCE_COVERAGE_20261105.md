# MOIS × artificial light — source-complete visitor evidence and white flora, 2026-11-05

## Executed expansion of real data, not synthetic tests

After the user requested more data, the older live GloBI flower-interaction snapshot (13,913 distinct accepted study-claim combinations, sample capped at eight 250-row pages per family × relationship) was expanded with bounded, SHA-logged continuation of the **same** actual GloBI API source. No images were redownloaded or re-segmented.

- Incremental stage 1 [run 41975735211](https://github.com/zuizui0223/mois/actions/runs/41975735211) SUCCESS: **41,106** distinct validated source rows, **27,193** genuinely new distinct records. [Artifact 13392365153](https://github.com/zuizui0223/mois/actions/runs/41975735211/artifacts/13392365153), SHA256 **88ecf29dc66ad70fc198c0e5339bdb58d0da981cf0657b657b390810839b7466**. Stage 1 used 105 new provider page requests, zero API failures, but the source gives 1,001 records to a 1,000 requested limit on some pages. Its initial "0 capped" receipt was **wrong** (a bug using equality instead of >=); five queries had not finished. Stage 1 raw page SHA log proves this.
- Stage 2 [run 41976552977](https://github.com/zuizui0223/mois/actions/runs/41976552977) SUCCESS: continued **only those five** formerly capped queries from their genuine signed last offsets, **20,861** newly distinct accepted study-claim rows, for **61,967 distinct GloBI records** total. 98 additional GloBI page requests, no provider errors. Only **one** Nymphalidae visitsFlowersOf query remains capped. [Artifact 13391803415](https://github.com/zuizui0223/mois/actions/runs/41976552977/artifacts/13391803415), SHA256 **2223332cbf89ea052e82bc75487850044bcb0509a22d3c8c62fafa57f425559e**. Extended original GloBI cache SHA256 **a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994**.
- Both stages preserve botanical plant-name, genuine insect family path, interaction relation, original study-title fields, query status/offset and SHA. This is still **live/unversioned GloBI API**, NOT full pinned GloBI v0.9 2.7GB source. The full source is the future confirmatory input if the question matures; GloBI itself recommends stable releases for research [official data page](https://www.globalbioticinteractions.org/data).

No source-specific direct species-flower trait rows from Island are redistributed.

## Primary new method finding: "moth-only" fails when records accumulate

**Counts of plant names WITH at least one documented GloBI flower-visit claim**:

| Evidence category | Original limited 13,913 claims | Stage 1 41,106 claims | Stage 2 61,967 claims |
|---|---:|---:|---:|
| Moth-family recorded, butterfly-family not recorded | 563 | 719 | **519** |
| Butterfly-family recorded, moth-family not recorded | 1,611 | 2,926 | **4,350** |
| Both recorded | 831 | 1,643 | **1,843** |

A species appearing in the first category is **never proven to lack butterfly visits**. The decreasing "moth-only" count is a sampling/completeness discovery, not an ecological disappearance. A sample of records cannot supply an independently confirmed diurnal-only control group or nocturnal moth pollination. Butterflies may also visit plants at twilight; effective night pollen transfer was not recorded.

The originally proposed **"white × moth-only × artificial light"** analysis is therefore **not robust to increasing visitor evidence coverage**. Its geographic cell denominators shrink as the alleged exclusive pollinator class is reclassified to "both".

## Exact independent documented white × species-level visitor records

All source records were joined to the **frozen** MOIS original canonical one-photo-per-observation frame (1,032,882 opportunities, strict9 success 401,703, SHA f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c), exact WCVP accepted plant names and independent High/Medium direct Island species colour (SHA 7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577).

Primary source arm: historical_main, six photo years **2019–2024**, one original observation ID/photograph and photographically classifiable strict9 only, between latitudes −60° and 60°.

After stage 2, **four independent documented colour × GloBI recorded-only status categories** need all 4 source-by-colour groups present. Because "only recorded" is not true biological absence, all counts are feasibility/support diagnostics, not estimates of effective moth or butterfly pollination:

| Spatial bin | >=1 independently documented accepted species in each of four classes | >=2 in each | >=3 in each |
|---|---:|---:|---:|
| ~25 km (0.25°) | **32** | **1** | 0 |
| ~50 km (0.5°) | **66** | **3** | 0 |
| ~100 km (1°) | **105** | **7** | 0 |

Before the second extension, 25km contained 90/8/1 such cells, so source expansion reduces apparent exclusive comparative support. Coarsening grid can increase taxonomic support but conflates local ALAN at sub-cell scales; latitude/longitude degree widths are not literally identical km.

## Positive-visitor evidence without turning missing records into absence

Instead of assigning false true "night-only/day-only" statuses, **retained two independent positive flags**: source reports a moth-family flower visit? and source reports a butterfly-family flower visit? A species with both gets both flags and is never counted as two independent plant species. This is evidence metadata, not a floristic or pollinator dependency trait.

New numeric-only [run 41977526551](https://github.com/zuizui0223/mois/actions/runs/41977526551) SUCCESS. Eight unit tests, six numeric input SHA locks and result SHA verified. [Artifact 13392323286](https://github.com/zuizui0223/mois/actions/runs/41977526551/artifacts/13392323286), ZIP SHA256 **d82ef9539595c9e1cc63cde8d4edea329983ae271e1af7f9ae6035c3a7e2fdf5**.

- Original source historical_main **113,366 strict9-success photo IDs**, WCVP matched **3,601 distinct accepted plant species** with at least one of moth/butterfly documented interactions, **1,343** plant species with independently documented known single-state flower colour in Island.
- For the same six years, 0.25° cells: **4,013** have at least one documented-white species with moth evidence and at least one documented-white species with butterfly evidence. **2,221** have all four positive (white/nonwhite × moth/butterfly) counts >=1; **706** have each count >=2; **328** each >=3. **BUT** these classes overlap heavily: most cells count the SAME WHITE species as both moth- and butterfly-visited. Only **61 of the 2,221** all-four-supported cells have zero white species with both visitors. Even those can have shared NONWHITE species.
- For comparison, the 2019–2024 photo records pooled at the existing **2020-year** satellite-light 1,044 cell positions yield 388 locations with both positive white visitation groups, 270 with all four positive evidence classes >=1, 80 with all four >=2, 19 with all four >=3, but **only 4** with all four while having zero overlapping white species. **The 2019–2024 pooled photos use 2020 VIIRS brightness, a temporal mismatch**, so this cohort is **a source-support gate, not a scientific ALAN association coefficient**.
- Source-positive counts are also not proof of true nocturnality or functional dependence. White colour in independent species literature differs from original photographic ROI colour of the individual photo, and the original actual flower organ has not been verified for most automated photos.

## Ecological decision for the proposed triple difference

**Do not fit/promote a causal ALAN×white×true nocturnality 3-way coefficient from these data.** This is not a failure to obtain enough photos: biological guild identity is not measured and the comparison "no moth record" does not certify a diurnal control. Even a full GloBI database cannot turn unrecorded interactions into genuine absence without independent complete observation.

In the 2020 direct VIIRS original-photo matched analysis the white fraction against light was non-directional (full 1,044 cells: β−0.0919 pp/log-unit, CI [−2.2095,+2.0258]). The earlier 108-cell source-positive moth vs butterfly original ROI contrast was exploratory and source incompleteness made exclusive group assignment unstable. Do not overwrite either historic observation with the present source-richness findings.

The next scientifically decisive action is species-direct evidence for actual **night-opening / nocturnal EFFECTIVE pollination** and a documented diurnal reference, with study time/location and floral state; then use exact same-year source photos and VIIRS or annual pinned geospatial radiance by season. An image source larger than 1M original photos cannot substitute for this guild identity.

No original photo acquisition, no RGB/ROI re-inference, no unknown colour imputation, no historical Mois-main changes, no MNAR-tipping sensitivity tests, no species extinction claim. PR #199 stays draft.
