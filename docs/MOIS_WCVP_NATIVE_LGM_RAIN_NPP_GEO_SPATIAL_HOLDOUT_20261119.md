# MOIS — WCVP native floras × LGM rainfall history and NPP (2026-11-19)

## 判定：LGM降水の全球関連はWCVP在来種・広域空間調整で頑健ではない

**完走した真正解析：** [GitHub Actions #44478550177](https://github.com/zuizui0223/mois/actions/runs/44478550177) SUCCESS。original signed GBIF accepted **37,104 photo9 species, 5,493,931 species×100km grid pairs, 396,142,825 GBIF occurrence records**, and 199 genus×occupancy-breadth hue permutations kept exactly unchanged. Source-constrained WCVP unique-exact taxa **36,762** and original fixed TDWG3 statuses. 12,295 climate/NPP/NDVI/AGB-complete cells in the signed original environmental output.

The global full-flora environment model's **all** prior 500 km FE coefficients, **12,239 fit cells and 804 blocks**, are reproduced to 1e-9 before any status-specific analysis. WCVP original 18,577 grid rows, strict-native 7,231 geographical cells (30+ native recorded species), and matching TDWG3 status counts verified. Intersection of 12,295 climate cells with >=30 strict native species: **7,201 cells**; moderate: **8,119**; same-cell strict native and introduced each >=30: **3,106 cells**.

### 1. Main seven-exposure model: current temperature/rain, LGM temperature/rain changes, NPP, NDVI, woody AGB

Each region-native floral assemblage is a subset of currently observed accepted species at the same real GBIF grid. For each subset, the **response is its own source-specific photo9 observed Shannon minus the expected Shannon under the same 199 global genus×occupancy-width label shuffles**. Negative is relative nine-hue concentration, not increased or decreased petal pigment.

Predictors and confounders: current temperature, log annual rainfall, 3-GCM median LGM-to-modern ΔT and Δlog rainfall, annual MOD17 NPP, MODIS NDVI, log ESA woody AGB; log subset species number, subset expected permuted Shannon, **log all-species GBIF record effort**, absolute latitude, longitude sine/cosine. All terms standardized in each selected cohort; 500 km or 1000 km spatial fixed effects and block-cluster SE. BH adjustment across the seven environmental predictors **within each reported model**, not globally across all exploratory model families.

| Cohort | Total status-complete cells | Fit cells | Independent block clusters | LGM rainfall beta/SD | 95% CI | p | q7 (BH) |
|---|---:|---:|---:|---:|---|---:|---:|
| **Strict WCVP native ≥30; 500 km FE** | **7,201** | 7,150 | 629 | **−0.008027** | −0.016162 to +0.000109 | .05314 | .18598 |
| Strict native ≥100; 500 km FE | 5,388 | 5,326 | 555 | −0.005909 | [see numeric output] | .07189 | .23178 |
| **Strict WCVP native ≥30; 1000 km FE** | **7,201** | 7,191 | 213 | **+0.000301** | −0.007159 to +0.007761 | .93698 | .93698 |
| Moderate native ≥30; 500 km FE | 8,119 | 8,066 | 692 | −0.006160 | [see numeric output] | .10158 | .34033 |
| Moderate native ≥30; 1000 km FE | 8,119 | 8,104 | 231 | −0.000740 | [see numeric output] | .82671 | .96491 |
| Strict introduced ≥30; 500 km FE | 3,110 | 3,020 | 374 | +0.010402 | [see numeric output] | .07700 | .17967 |
| Strict introduced ≥30; 1000 km FE | 3,110 | 3,095 | 164 | +0.004095 | [see numeric output] | .47086 | .54934 |

**Interpretation:** The previously estimated full-flora 500 km LGM rainfall association (β=−0.008781, q7=.0116) is **not robust under WCVP native restriction**. The strict-native effect is at most marginal unadjusted at 500 km and vanishes under 1000 km FE; it is unsupported at the registered seven-environment BH threshold. This does **not** prove the earlier association was caused by introduced species; the tested species pools, WCVP status, geographical support and source-colour biases differ.

### 2. Same-cell native minus introduced status comparison

For exactly the same current GBIF 100 km cells that each contain native and introduced 30+ species, the response is:
`ΔH = (H_native − E_null,native) − (H_introduced − E_null,introduced)`.

Regress on the seven exposures plus **both** native and introduced richness, their expected-null Shannon difference, all-status GBIF effort and geography. Positive ΔH means introduced assemblages have relatively stronger colour concentration under that null; this is **not** evidence of colour evolution after introduction.

| Matched subset | n source cells | n fit cells | LGM rainfall beta/SD | 95% CI | p | q7 |
|---|---:|---:|---:|---|---:|---:|
| Strict status, each ≥30; **500 km FE** | 3,106 | 3,016 | **−0.015116** | −0.027720 to −0.002512 | .01875 | .04601 |
| Strict status, each ≥30; **1000 km FE** | 3,106 | 3,091 | **−0.002408** | −0.013699 to +0.008882 | .67589 | .94247 |
| Moderate status, each ≥30; 500 km FE | 3,602 | 3,509 | −0.010592 | [see numeric output] | .06054 | .12040 |
| Strict status, each ≥100; 500 km FE | 1,299 | 1,226 | −0.007498 | [see numeric output] | .14954 | .26170 |

This paired effect is **sensitive to the geographic conditioning scale and subset richness**, and the strict/moderate cohorts overlap. The one supported q7 result is hypothesis-generating, not a stable native–introduced palaeoclimatic response difference. Same-cell comparison does not eliminate geographically structured source status and introduction histories.

### 3. Independent-region heldout test of adding palaeoclimate to current climate and vegetation

Predefined status: strict native ≥30, **7,201** environment-complete cells; **223 distinct 1000 km geographic blocks**, GroupKFold with 5 disjoint block-heldout test folds, deterministic Ridge(alpha=100, solver='cholesky'), training-fold-only standardization. Predict photo9 **native Shannon deviation** from geography/observational effort/native richness/expected null/current temperature+rain and NPP/NDVI/AGB, and add both LGM-to-modern ΔT and ΔlogP.

| Model | OOF R² | MSE |
|---|---:|---:|
| Without LGM climate history | 0.153067 | 0.00695924 |
| With LGM climate history | 0.167939 | 0.00683704 |

**ΔMSE (with-minus-without) −0.000122205**, 1000km **block bootstrap 2000 times**, 95%CI **[−0.000268661,+0.000027290]**. Direction suggests weak predictive improvement, but CI includes 0 and the improvement is **not established across heldout blocks**. The fit is a descriptive community pattern, not a past distribution model. Both original responses and covariate sources held fixed; CV deterministic.

### 4. NPP and biomass stress

- Strict native ≥30; 500 km: **NPP β=+0.009040, p=.01414, q7=.09897**. Strict native ≥100: β=+0.003438 (p=.2557); strict native 1000 km: β=+0.000860 (p=.8482).
- Strict native woody AGB 500 km: β=−0.001197 (p=.7323); 1000 km: β=−0.002515 (p=.5808). No robust independent woody biomass signal.
- Matched strict 500 km native-minus-introduced: NPP β=−0.004182 (p=.6085), 1000 km β=−0.000338 (p=.9656).
- The earlier full-flora extreme LGM displacement × NPP interaction **lost significance when the upper 25% of displacement was excluded**; that previously audited sensitivity still applies. Do not equate carbon-flux NPP or standing woody AGB with floral biomass, floral pigment cost or plant fitness.

### 5. Ecological conclusion and next decisive test

**Supported now:** The geographically structured genus×occupancy-breadth nine-hue species composition of the world's *photographed accepted plant species* contains residual colour aggregation, but this particular **LGM precipitation-history association is not general across WCVP-native assemblages and spatial conditioning scales**. Historical climate and NPP can covary with the recorded plant pool, but the available present-day photo ROI hues + GBIF cells cannot identify adaptive historical selection, true petal colours, extant native abundances or glacial refugia.

**No refugia identified.** The LGM climate values characterize reconstructed palaeoclimate at *current* GBIF sites; they are not LGM occurrences, survival areas, pollen fossil sites or species palaeoniche models. Genuine refugia require species-level projected suitability with spatial availability and 3-GCM uncertainty, and independent palaeoecological/phylogeographic records. Distinguish **current climate**, **directional LGM change**, **magnitude of LGM change**, and **historical occupied suitability**.

**Priority on next round:** test historical refugia explicitly using independent, source-verified palaeodistribution and fossil/pollen evidence and WCVP native species; block hold out whole biogeographic realms, and verify photographic floral organs before any physiological explanation. Do not turn the current historical-climate correlation into a refuge map.

### 6. Permanent, source-audited record

- Signed WCVP native geographic permutations: [Actions #44242421559](https://github.com/zuizui0223/mois/actions/runs/44242421559), ZIP SHA256 `cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd`.
- Signed climate + biomass + history: [Actions #44236695589](https://github.com/zuizui0223/mois/actions/runs/44236695589), ZIP SHA256 `4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6`.
- Signed 396,142,825 GBIF observations / accepted species grids: artifact #14166701658 ZIP SHA256 `422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2`.
- **This run** [Actions #44478550177](https://github.com/zuizui0223/mois/actions/runs/44478550177) **SUCCESS**; numeric [artifact #14507347253](https://github.com/zuizui0223/mois/actions/runs/44478550177/artifacts/14507347253); ZIP SHA256 `06eadf6e712297efe914ecce7fca5ef2ef98aa6127ffe2a6b3e366bfff618a95`. Three outputs: every model beta/CI/p/q in CSV, full model+OOF JSON, source and output SHA receipt.
- Reproduction script: `scripts/analysis/mois_WCVP_native_paleorain_NPP_20261119.py`, workflow: `.github/workflows/mois-native-paleorain-npp-20261119.yml`.
- Implemented on draft [PR #206](https://github.com/zuizui0223/mois/pull/206). **main unchanged**. This is an explicitly exploratory extension, not a separate independent dataset.

