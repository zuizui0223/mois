# MOIS: complete p0008 authentic original-photo latitude, month and missingness audit

**Date:** 2026-10-21 (JST). **Data status:** actual frozen YOLO11/EfficientSAM ROI-v4 photo outcomes, not metadata image candidates.

## Published and independently authenticated outcome

New [PR #140](https://github.com/zuizui0223/mois/pull/140) introduced an explicit full-photo-ID, observation-aware audit of the **complete p0008 wave**. Seven synthetic regression tests and the repository CI passed; [production run 39452760985](https://github.com/zuizui0223/mois/actions/runs/39452760985) successfully re-read the original image Parquet, verified immutable source+file SHA256 and published [the complete original-photo coverage Release](https://github.com/zuizui0223/mois/releases/tag/mois-complete-p0008-original-photo-latitude-season-coverage-v1).

The source was the signed [p0008 true-original-photo full wave](https://github.com/zuizui0223/mois/releases/tag/mois-full93m-observation-third-original-RGB-p0008-v1), not the 93,209,351-photo metadata target frame. Its original source-manifest SHA256 is \`0038024ee40adfbc71aa427a112948ec5ce53cb53a7aae9879c1b32bc06eabde\`.

| Genuine observed p0008 photo metric | Count |
| --- | ---: |
| Exact unique original photo terminal ID rows | **36,375** |
| Real photographs successfully acquired | **36,331** |
| Acquisitions failed | **44** |
| Frozen four-hue ROI outcomes | **13,874** |
| Acquired but unclassifiable by frozen ROI | **22,457** |
| Distinct original observations | **20,773** |
| Observations with two or more photos | **8,191** |
| Photos belonging to multi-photo observations | **23,793** |
| Original iNaturalist taxon IDs with attempted photos | **6,549** |
| Source taxon IDs with any classifiable ROI | **4,453** |
| Taxon IDs with at least 5 ROI hue labels | **752** |
| Taxon IDs with at least 30 ROI hue labels | **18** |

The 13,874 four-hue classifications include **6,762 white**. The photographic white fraction *conditional on automatic ROI classifiability* is **48.74%**. Equal weighting of the **10,418 classifiable distinct observations**, first computing each observation's within-ROI white fraction, gives **49.00%**. This comparison is *descriptive*, not a valid correction for floristic/photo/taxon sampling selection, source-image detectability or observation-level dependence in ecological inference.

## Real-photo latitude coverage and classification missingness

Bands use original photo observation latitude only (tropics: 23.44°S through 23.44°N), not phylogeographic realms, native status or island-specific provenance.

| Latitude band | Original photo attempts | Four-hue classifiable | Classifiable per attempt | White among classifiable |
| --- | ---: | ---: | ---: | ---: |
| Northern extratropical | **29,708** | **11,678** | **39.31%** | **49.05%** |
| Southern extratropical | **4,564** | **1,479** | **32.41%** | **46.32%** |
| Tropical | **2,103** | **717** | **34.09%** | **48.68%** |
| **p0008 full genuine photo wave** | **36,375** | **13,874** | **38.14%** | **48.74%** |

Tropical source photos are **5.78% of the attempted photos**, but only **5.17% of classifiable images** because their automatic ROI classification success is below that of northern extratropics. Likewise, southern extratropical images comprise **12.55% of original-photo attempts** but **10.66% of classified outcomes**. These observations justify measuring and explicitly adjusting classifiability before interpreting geographic colour differences. They do **not** demonstrate latitude-caused classification failure: species composition, photographer behaviour, plant life stage, camera framing and local season all vary.

Real classifiable geographic coverage includes **3,342 original taxon IDs in northern extratropics**, **894 in southern extratropics**, and **435 in the tropics** (region-specific ID counts can overlap, and should not be summed as independent accepted species). Region-specific distinct classified observations: 8,835 / 1,059 / 524, respectively.

## Real-photo observation month and apparent seasonal white patterns

| Original observation month | Attempted photos | Classifiable four-hue ROIs | White share among classifiable |
| --- | ---: | ---: | ---: |
| Jan | 1,073 | 344 | 48.55% |
| Feb | 1,233 | 377 | 47.48% |
| Mar | 2,221 | 815 | 46.01% |
| Apr | 5,271 | 2,049 | 49.49% |
| May | 5,072 | 2,070 | 50.00% |
| Jun | 4,993 | 2,089 | **53.14%** |
| Jul | 4,945 | 2,038 | 50.05% |
| Aug | 3,844 | 1,473 | 45.49% |
| Sep | 3,197 | 1,182 | 45.09% |
| Oct | 2,180 | 721 | **43.55%** |
| Nov | 1,334 | 400 | 45.75% |
| Dec | 1,012 | 316 | 51.27% |

The descriptive June-to-October conditional white photo fraction differs by **9.59 percentage points**, but it is not evidence that *the same species' petals* respond to short-term temperature, moisture or solar radiation. In the northern extratropics alone, June was 1,055 white / 1,980 classifiable = **53.28%**, whereas October was 184 / 436 = **42.20%**. **Taxon and photography mixture across seasons remain confounders**, including the possibility that different species flower at different times. Observed dates span **1974-10-06 through 2026-09-27**, so this is not a single-year matched phenological experiment.

Month-specific ROI-classifiability is also variable: **February 30.58%** vs **June 41.84%**. Separately model or reweight missingness and taxon/observer/space prior to inferential claims.

## Methodological consequence for MOIS

The completed 36,375-photo p0008 sample is large enough to establish a *reproducible photograph-level seasonal pattern* and its missingness, but source biological/taxon independence is much smaller: only **18 source taxon IDs have 30+ classified photos**. Therefore the next environmental analysis should distinguish:

1. **Turnover/composition:** taxon × region × observation-month photography and the seasonally available flowering species.
2. **Within-taxon/local response:** source taxon/genus and spatial matching, conditioning on per-observation photographs; do not call cross-species seasonal shifts petal plasticity.
3. **Measurement selection:** model automatic ROIs classified/not classified by region/month/source method; preserve unknown/failed photos rather than recoding them nonwhite.
4. **Abiotic seasonality:** only after the first three distinctions, attach observation-date monthly temperature/radiation/water to photos with available location and date; leave altitude/soil/terrain as relatively static adjustment covariates.

Both p0004 and p0008 original observation-level source waves are now complete. They must not be naively summed with independently enriched geo05 photos for global white-frequency estimates, and their nominal iNaturalist taxon IDs are not necessarily unique accepted species. The ongoing p0016, p0032 and geo05 original image workers remain unchanged.

## Immutable evidence

- Source full p0008 run: \`39371199127\`; original full wave source manifest SHA: \`0038024ee40adfbc71aa427a112948ec5ce53cb53a7aae9879c1b32bc06eabde\`.
- Coverage CI and production: \`39452760985\` (success).
- Verified coverage artifact \`12277924138\`: ZIP SHA256 \`142bceb023d3429c927bc9aaf988d5bdc9beec9b5648bdfc3a44009b0291e30e\`; its four output CSVs' SHA256 all match the published \`p0008_real_photo_coverage_receipt.json\`.
- Files: \`original_photo_latitude_band.csv\`, \`original_photo_observation_month.csv\`, \`original_photo_latitude_by_month.csv\`, \`original_taxon_coverage_by_latitude.csv\` and \`p0008_real_photo_coverage_receipt.json\`.
- This full 36,375-row Parquet **was actually read row-by-row in GitHub Actions with pyarrow**, and original photo grain, source signature, classification status and observed date/coordinates verified. Local separate validation reopened the published ZIP, matched the artifact SHA, recomputed each CSV SHA and checked region/month/arithmetic against its JSON receipt. No new images were downloaded or inferred.

This is a source-backed real photo-coverage result; the new field/ROI-missingness patterns are not causal flower-colour adaptation or proof of native-flora whiteness.
