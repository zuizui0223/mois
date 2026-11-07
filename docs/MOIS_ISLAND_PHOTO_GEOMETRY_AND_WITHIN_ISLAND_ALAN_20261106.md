# MOIS: island-resolved floral colour and within-island artificial light (2026-11-06)

## Completed result, scope and provenance

The analysis now distinguishes **island versus continental photographic composition**, **differences among islands**, and **light–colour correspondence within the same island**. These are different estimands. Temporal change is not the primary question.

Original input: 1,032,882 distinct canonical photo/observation IDs, 401,703 original strict9 successes. Current spatial inventory: historical_main, years 2019–2024, latitude −60 to +60 and valid longitude, **455,089 attempts / 186,467 strict9 successes**. No new photo or original ROI measurement was made.

The old frozen island GPKG artifact was no longer downloadable. The SAME GSHHG 2.3.7 source archive was recovered and SHA verified instead: **8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41**. [Source recovery run 42175243391](https://github.com/zuizui0223/mois/actions/runs/42175243391) SUCCESS; artifact 13482747651, ZIP SHA256 **afd51f8acdecae1815b29ef1e5ab38eb841e3a28f0768ebcf8bfa592d9c9e2fc**. This source-only run did not download photographs.

GSHHG high-resolution L1 polygons are recombined by sibling_id before EPSG:6933 equal-area calculation. Use the existing island-project operational policy: 5 <= area <= 7,000,000 km² is an island, larger landmasses continental, <5 km² separate. Points within L2 inland-water polygons are separately excluded, not called continental or marine islands. Boundary/missing/outside-L1 points are not assigned to their nearest landmass. The reconstruction yields **8,260 island units plus five continental landmasses**. It is NOT claimed to be byte-identical or ID-equivalent to the older 8,265-record island GPKG; current analytical identities are GSHHG sibling IDs. Display names are assigned through stored point-in-polygon anchors, not used to determine island membership.

## 1. Actual island photographic inventory

| Spatial source class | Photo attempts | strict9 successes |
|---|---:|---:|
| Island >=5 km² | **52,737** | **19,086** |
| Continental landmass | 389,144 | 162,328 |
| Landmass <5 km² | 1,051 | 392 |
| Inside L2 / inland-island domain | 4,096 | 1,797 |
| Outside mapped L1 interior | 8,061 | 2,864 |

The island attempts occupy **1,038 islands**. Successful strict9 measurements occur on **740 islands**, with **4,166 accepted plant species** and 18,774 photographs that unambiguously resolve to accepted species. Offshore/unmatched coordinates may reflect map generalization, older coastline, reclamation, water or observation uncertainty; they are not automatically bad observations, true absence or mainland.

Largest successful island samples include Great Britain 3,368; New Zealand North Island 2,347; Taiwan 2,009; New Zealand South Island 1,996; Vancouver Island 776; Tasmania 671; Honshu 438. These are counts of measured photographs, not independent population sizes, native species richness or balanced sampling of island floras.

## 2. Matched-year and matched-landmass 2020 sample

The existing 1,044 VIIRS/WorldPop-complete 0.05° cell CENTRES comprise 72 island, 943 mainland, six inside L2 and 23 outside L1. Original photograph points were assigned independently. Retain photo–exposure joins only when the point and the cell centre lie within the SAME GSHHG landmass. Reapply the fixed source criteria: >=5 attempted photos, >=3 nominal taxa, >=3 original strict9-success photos and >=3 unambiguously accepted species in the retained cell.

**Final 2020 data: 18 islands / 66 cells versus 919 mainland cells.** Island cells have 285 strict9-success photos, of which 280 unambiguously accepted-photo records (234 species) enter the colour response. Mainland cells have 4,487 strict9 successes, of which 4,442 accepted-photo records (2,370 species) enter the response. Repeated photographs are averaged within accepted species and cell; species are equally weighted within cell, then cells equally weighted.

| Observed continuous white ROI component | Mean |
|---|---:|
| Island 66 cells | **53.14198%** |
| Mainland 919 cells | **45.81410%** |
| Unadjusted descriptive island–mainland difference | **+7.32788 percentage points** |

This is NOT the percentage of white-flowered species. The contrast is not balanced for climate, geography, habitat, species composition or photographic selection and is not an independent causal island effect. Restricting to the **85 accepted species observed in both classes**, equal-species descriptive mean difference is +4.99149 pp (median +2.01820), but those records are not matched geographically or climatically; no inferential certainty is claimed for that calculation.

Between-island source support: Taiwan 16 cells (mean white component 61.39%); NZ North 14 (58.29%); NZ South 11 (41.00%); Great Britain five (41.49%); Vancouver four (58.16%). Taiwan and the two main NZ islands contribute **41/66 cells**. Results cannot be generalized to small oceanic islands merely by treating every island record identically.

## 3. Within-island ALAN correspondence, not temporal change

Only islands with >=2 retained cells and non-identical VIIRS values can identify an island-internal association. **Eight islands / 56 cells** meet this rule. The ten single-cell islands remain in descriptive summaries but are not counted as informative clusters for the within-island slope.

Response: cell-level accepted-species-equal original white ROI mass in percentage points. Predictors: current-year log1p(VIIRS upward radiance), current-year log1p(WorldPop population density), BIO5, log1p rainfall, elevation, log1p photo attempts, and island fixed effects. All radiances are the existing Jan/Apr/Jul/Oct monthly raw values, not a new official annual composite. Nominal CR1 cluster variance by island with t(G−1) interval.

**Light coefficient = +1.96046 pp per log1p radiance unit; nominal 95% interval [−12.40851, +16.32942], p=0.7564.** This provides no clear negative within-island correspondence. Only eight informative island clusters make the interval imprecise and not confirmatory. Across the full 66 island cells the log-light/population correlation is **0.889**. Adjustment does not establish an experimentally isolated effect of illumination or eliminate land use and sampling effects.

2024 supports five islands / 20 cells under the same geographic rules; this is saved as source support only. Current 2024 population was not added, and no change model was run.

## 4. Connection to independently reviewed moth-pollen-associated species

Of the previous 32 botanical-review species and 2,707 original photos, **490 photos / 31 species / 51 islands** lie on mapped islands. Restricting to the 30 species with comparable flower-colour organs leaves **463 island photographs**. The complete same-year 2020 photo/exposure/landmass criteria leave five island photographs from this 32-species list, but two are Plantago records with noncomparable flower/anther colour. Thus the flower-colour-applicable same-year island support is only **three species / three photos in 2020**, and **zero in 2024**. UK-study moth-body pollen carriage is not proof of local effective night pollination at the NZ, Taiwan or other MOIS coordinates.

This does not overwrite the prior independently described five white Walton species, original Island trait labels, or automatic strict9 measurements. It exposes why a general island floral-composition analysis is currently much better supported than a white × truly nocturnal pollination × ALAN mechanism test.

## 5. Exact additional same-year island exposure candidates

No image acquisition is needed to prepare additional numeric-exposure joins. Existing photo support and same-island centre checks yield:

| Year | Eligible island cell-years | Islands | Already joined to current same-year VIIRS output |
|---|---:|---:|---:|
| 2019 | 53 | 19 | 0 |
| 2020 | 67 | 19 | 66 |
| 2021 | 73 | 22 | 0 |
| 2022 | 114 | 33 | 0 |
| 2023 | 146 | 36 | 0 |
| 2024 | 224 | 47 | 20 |

Total **677 island cell-years; 591 remain without the corresponding verified exposure output**. This is an acquisition manifest, not 591 new measurements. Same sites across years are not independent islands. The old WorldPop service's 2020 population must not be silently used as contemporary 2021–2024 population. The research question remains within-year spatial correspondence, not year-to-year change.

## Verification and delivery boundary

- All five canonical/WCVP/VIIRS input member hashes checked against prior source receipts.
- All extracted GSHHG L1/L2 component hashes checked against the recovered original archive receipt.
- Nine local unit tests passed: concave polygon versus bounding box, boundary handling, small islands, inland water, dateline sibling grouping, duplicate identity, cross-channel exposure mismatch, species equal weighting, unknown-colour preservation and multi-site within-island support.
- 200 actual assigned points independently rechecked by scalar point.contains and L2 intersection.
- Independent numpy least-squares and manually assembled CR1 cluster covariance reproduce the statsmodels coefficient and confidence interval to <1e−7.
- Complete local numerical reproducibility bundle: **MOIS_island_flower_ALAN_verified_20261106.zip**, SHA256 **7820fa60488771d34d8e6c5c33f6a249f24a70c2772a010845543dbeb2d097a7**. Includes scripts, tests, logs, public cell/island aggregates and acquisition candidates, not restricted individual species-colour source records or private photo-ID tables. Repository-wide analysis CI is NOT claimed for the locally executed analysis; the coastline recovery workflow itself succeeded.

No photographs downloaded, no ROI inference, no technical unknown imputation, no rejected MNAR sensitivity scenarios, no main-branch edits or PR merge. Geographic origin is not confirmed wild/native status; photo assemblages are not complete island floras; no moth dependence or local species extinction is inferred.
