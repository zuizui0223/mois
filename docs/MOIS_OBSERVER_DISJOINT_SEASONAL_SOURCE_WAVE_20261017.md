# MOIS — observer-disjoint original photo opportunities for seasonal flower-colour inference (17 October 2026)

**Canonical owner:** `zuizui0223/mois`. This is *not* an extension of the `fcp` within-species D/W ecological tests, and it does not rerun or modify the historical MOIS 100,543-photo global moisture coefficient.

## The obstacle discovered in the full original source

The [frozen dated original iNaturalist ledger](https://github.com/zuizui0223/fcp/actions/runs/38884223530) contains **3,430,502 original source photo IDs** across 1,000 outcome-blind candidate taxa, with real observed-on dates. The matching [0.25° month/5-km accuracy census](https://github.com/zuizui0223/fcp/actions/runs/38901290387) contains **835,833** site-month rows. The previous source-only two-season screen found **268** taxon×0.25°-site×season-year pairs in **39** taxa with ≥24 independent original *observations* per season.

**Distinct observations are not independent photographers.** Some high-depth taxa have more than 100 original observations per season but very few distinct observers, creating a major pseudo-replication risk. Re-reading every one of the 3,430,502 source-photo records, retaining the original site-months containing ≥2 observations and ≥2 observers, checking original positional uncertainty ≤5 km, and excluding any IDs in the previous 100,543-photo MOIS measured atlas and 64-photo Asclepias first pilot produced the following capacity audit:

| Distinct observers obtainable **in both** seasons with no cross-season reuse | Eligible original site×years | Distinct taxa |
|---:|---:|---:|
| 10 | 217 | 33 |
| 20 | 128 | 23 |
| 24 | 72 | 19 |
| 30 | 35 | 11 |
| 40 | 13 | 5 |
| 50 | 7 | 2 |
| 64 | 1 | 1 |
| 80 | 0 | 0 |

The feasible balanced independent-photographer capacity for a source site-year is `min(nA, nB, floor(|A union B|/2))`. It is a **sampling capacity**, not a number of botanically verified blooming plants and not an inference of natural flower-colour variation. This corrects the earlier tendency to treat 268 pairs as 268 high-powered biological comparisons.

## Fixed five-taxon original photo-ID sampling wave

Across the five best independent-photographer-eligible taxa, the `mois` source-only selector reserves **40 original photos in each of two seasons**: **400 distinct photo IDs and 400 distinct photographer IDs** globally, additionally one photo per original observation. An original photo and both its site-year and season are selected entirely by immutable identifier hash; photo pixels, current flower-colour classifications, ecological predictors, iNaturalist flowering annotations and image visibility do not enter this selection.

| Species | Year | Seasons sampled | Original source location (0.25° cell center) | Reserved photos |
|---|---:|---|---|---:|
| *Asclepias syriaca* | 2025 | Summer / Autumn | 45.625°N, 73.625°W | 40 + 40 |
| *Aralia californica* | 2025 | Autumn / Summer | 37.875°N, 122.625°W | 40 + 40 |
| *Datura discolor* | 2023 | Spring / Winter | 33.375°N, 116.375°W | 40 + 40 |
| *Achillea millefolium* | 2023 | Spring / Autumn | 48.375°N, 123.375°W | 40 + 40 |
| *Linaria vulgaris* | 2023 | Summer / Autumn | 42.375°N, 71.125°W | 40 + 40 |

These locations are **source metadata coarse grid cells**, not measured botanical populations. Seasons differ among taxa; the next design must model observed month, site, and taxon, rather than pooling every contrast as one universal summer-autumn comparison. The panel is concentrated in the Northern Hemisphere and is **not globally representative**. The photograph selection is conditioned on observer opportunity, which must be reported, not concealed.

### Original source identity and exclusivity

- Exact original raw photo-ID CSV SHA256 `58e037a210e23d3ff5d6f04d392665ada5ee9ac5149d0676672c239cf092df01`.
- Source 835,833-site-month file SHA256 `7691c7cf975da4e73494ad3fa0f2cf8df52302e84206017dd525051cd672b432`.
- Historical MOIS 100,543-unique-photo archive SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`.
- Previous 64 photo ID source SHA256 `d5c4bdcd632d0de8f2d20a663674e87f35e323a42b896291efefb93c46023053`.
- 268 observer-capacity audit aggregate CSV SHA256 `4dbba0d6d3b208182fcfa99471e97b2819f6602f4b05578dc7abeeb38aff74aa`.
- Deterministic 400 photo source-only CSV SHA256 `b961b50ee6efe9ec08a5575a8eb1e34b1edf7b7fcb426226fb11ee901cfd3329`; original image URL and licence attribution stay in a time-limited GitHub Actions artifact, not public Git-tracked original photo metadata.

`results/mois_original_268_observer_disjoint_400_photo_wave_20261017.json` fixes the main counts, photo source hashes and selection. `scripts/analysis/mois_observer_disjoint_photo_wave_20261017.py` and `.github/workflows/mois-observer-disjoint-season-wave-20261017.yml` reconstruct the full source using the original pinned GitHub artifacts.

## The next biological gate

**None of the 400 newly reserved IDs is a verified open flower yet.** The earlier 64-image pilot produced just 15 clean automated four-colour labels after mixed organ removal and **zero** independently botanically confirmed open-flower petal/tepal colour diagnoses. Therefore this new 400-photo wave is a pre-colour sampling product, not 400 usable floral outcomes or a demonstrated weather effect.

The next technical stage is to run blind frozen original ROI-v4 image measurement on the 400 IDs, followed by two independent botanical reviewers of photo-specific open flowers and petal/tepal ROI position, disagreement adjudication, missingness denominators, and real observed-month site-specific TerraClimate. Preserve original 400 photo IDs through all downstream missingness (no outcome-based replacement). Only after within-taxon, within-site/year verified flower colours have adequate size and biological variation should MOIS estimate white-photo fractions versus actual monthly temperature, sunlight, precipitation, water deficit or VPD. These are not a photosynthetic or pigment-fitness experiment.

The still-separate `fcp` paper is about geography and phenotype-space structure of **variation within named species**, with its own 149,900-photo inference and frozen D/W results. It does not borrow this MOIS observer opportunity as confirmation, and this MOIS study does not borrow FCP's ecological effects as a climate mechanism.
