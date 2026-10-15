# Full ecological audit: biogeography × climate × soil × elevation, with an explicit terrain-data gate

**15 October 2026 · exact original FCP source rerun · working mois/GEB supplementary research checkpoint.** The two original archive ZIPs were SHA256 verified, reopened and independently recomputed in a single analysis, with no recoding of photographed floral hue or transfer of source photo IDs, original attributed photography, taxon names or precise sites to public GitHub. Full local numerical source output passed `--write` and `--check` on exact original files; the public repository contains a safe aggregate receipt, replayer, science fixtures and a source-independent synthetic terrain test. These are **retrospective photographic source analyses, not controlled eco-evolutionary selection experiments**.

## 1. Species, genera and scope are fixed for the soil question

| Source sample | Original classifiable photographs | Nominal species with ≥3 qualifying species/genus | Qualifying genera |
|---|---:|---:|---:|
| All 26 climate features, valid photographic hue and site | **45,836** | **20,694** | **1,934** |
| All 36 (26 climate + 10 soil) features complete | **31,925** | **16,801** | **1,665** |

The soil-complete photos cover **69.65%** of the climate26-complete source photographs. The **16,801/1,665 soil36 population is kept identical in every climate–soil model and taxonomically held-out prediction comparison** below. Its own world equal-genus congeneric white/BIO12 slope was independently reconstructed as **+2.256937 percentage points per original soil36 photo-site source SD**, exactly matching the original source decision ledger. The main climate26 source uses *different* taxa and its original white/BIO12 coefficient +2.254545 has a distinct source-environment normalization; a numerical resemblance does not make the two sample sizes equivalent.

## 2. Geographic missingness in soil is an ecological limitation

Soil36 completeness was highly unequal across original **climate26-classified photograph** elevation quartiles:

| Photo-site elevation band | Source climate26 photos | Soil36-complete | Soil completion |
|---|---:|---:|---:|
| −330–96 m | 11,396 | 5,081 | **44.6%** |
| 96–289 m | 11,504 | 7,461 | **64.9%** |
| 289–941 m | 11,469 | 9,325 | **81.3%** |
| 941–6,104 m | 11,467 | 10,058 | **87.7%** |

Soil coverage was **66.8%, 67.0%, 74.3% and 68.5%** across the four 90° longitude sectors from 180–90°W source sector to eastern Asia/Oceania (90–180°E). These are source photo coverage rates, not evidence that lowland plants lack soil water or that high-altitude soils are biologically more important. The soil36 complete-case ecological population differs systematically by altitude and region: inference beyond its sampled taxa requires an explicit missingness/representation design.

## 3. Soil associations and predictive information are not interchangeable

On the fixed soil36 sample, the original geographically adjusted genus-equal rainfall model yields white/BIO12 **+2.257 pp/original soil36 SD** (new whole-genus percentile 95% interval **[+0.864,+3.732]**). Fitting precipitation, vapour, peak temperature, solar radiation and altitude together changes the rainforest/rain slope to **+1.501** with ordinary 95% **[−0.134,+3.147]**, and vapour to **+3.149 [−0.426,+6.700]**. Because source predictors are correlated, those partial coefficients are distinct conditional estimands and are not selected-familywise 24-maxT discoveries.

Adding **five soil variables** (source soil pH, SOC, N, clay and *soil available-water proxy*) to this climate/elevation model yields a source-conditional available-water-proxy coefficient **+1.907 pp/original soil36 photo-site SD**, ordinary genus-bootstrap 95% **[+0.586,+3.119]** and an exploratory *within-this-model seven-axis* maxT95 interval **[+0.199,+3.615]**. In this exact fitted model, the selected pH, SOC, N and clay individual maxT intervals contain zero. **This is a positive observed edaphic proxy association, NOT the effect of measured field capacity or an observed fitness benefit.** It does not imply better out-of-genus generalization.

The ten original soil dimensions also contain pronounced near-linear composition closure:

```text
soil_sand_raw + soil_silt_raw + 10 * soil_clay_percent ≈ 1000
```

Across the 16,801 original soil species, the left side has mean 1000.011, SD **0.287** and maximum absolute deviation from 1000 of 1.0. Unpenalized all-soil10 linear models have normal-equation condition number approximately **1.17 million**. Treating sand, silt and clay as separate, freely varying soil causes produces unstable coefficient interpretations; explicitly choose a composition-aware basis or regularized prediction model instead.

## 4. Full predictive check: five-fold whole-genus holdout

Every comparison predicts the original sampled **species-mean fraction of classifiable white photographs** among the same 16,801 soil36 taxa. Training and evaluation are split by deterministic SHA256-hashed genera (five disjoint folds, no genus shared). Each model uses training-only feature standardization, fixed `ridge penalty=0.01`, predicted fraction clipping to [0,1], and equal genus weighting. The reported loss is genus-equal **squared prediction error of photo-white fractions** (not pollen, seed or survival fitness). The 1,999 whole-genus bootstrap intervals resample already-heldout losses; they do not refit the training splits.

| Fixed soil36 evaluation model | Heldout squared error | Improvement compared with geography alone |
|---|---:|---:|
| Geography (absolute latitude + circular longitude) | 0.229810 | reference |
| Geography + elevation | 0.229306 | +0.000504 |
| Geography + **precipitation and vapour** | **0.227455** | **+0.002355 [95% +0.001125,+0.003608]** |
| Geography + moisture + peak temperature + solar + elevation | 0.227481 | +0.002329 |
| That climate model + **soil5** | 0.227916 | +0.001894 |
| That climate model + **all soil10** | 0.228101 | +0.001709 |

Relative to the **identical full climate/elevation baseline**, soil5 addition **reduces** out-of-genus prediction improvement by **−0.000435** (bootstrap 95% **[−0.000700,−0.000171]**), while soil10 addition reduces improvement by **−0.000620 [−0.000958,−0.000284]**. Thus **observing a positive within-model soil proxy coefficient is compatible with soil inputs failing to improve a heldout predictive task**. The source supports a transferable, modest water-related environmental signal more than a portable direct soil mechanism, *for these fixed predictors, fixed ridge and selected soil-complete photographs*. This is not proof that soil never affects plant colour.

## 5. Biogeography matters for transport

We performed a stricter descriptive longitude-sector transfer using the same soil36 source species: hold out an entire 90° sector and **exclude from training every genus represented in the test sector**, so evaluation holds out both geographically distinct species and source genera. The number of training species changes substantially across sectors; no generalization to spatially unseen microhabitats or Wallace biogeographic realms is claimed.

| Heldout original longitude sector | Test species | Remaining training species, excluding test genera | Added rainfall+vapour gain over geography |
|---|---:|---:|---:|
| 180–90°W | 2,031 | 8,416 | **+0.004282** |
| 90–0°W | 4,733 | 4,769 | **+0.002184** |
| 0–90°E | 5,842 | 4,618 | **−0.000233** |
| 90–180°E | 4,195 | 5,841 | **+0.002535** |

Three sectors show modest positive source-heldout water gains, but **0–90°E does not**. This is evidence that the water–flower-colour relationship is neither universally portable nor a fixed proportional advantage of white flowers on all continents. Region boundaries and original study flora affect the result.

As a separate **rough geographic composition check only**, a Natural Earth **110m country polygon** overlay was applied to each of 20,694 source species' **average photo-site coordinate**. Of these, **3,158 (15.3%) could not be matched to coarse land polygons**, including potentially small islands and coastal or cross-continent averaged coordinates; classify them as unresolved rather than assigning a continent by guess. Among the matched/source-labelled groups, photographic species-mean white fractions were approximately Africa **41.4%** (2,984 taxa), Asia **48.9%** (3,323), Europe **44.3%** (2,050), North America **47.6%** (3,286), Oceania **47.4%** (2,713), South America **43.4%** (3,180). These are **photo-source-conditioned approximate continent descriptions**; 110m continents are not a phylogenetic biogeographic realm model, true country-specific species inventories, species-range polygons, or independent colour phenotyping. The exact 110m polygon fixture shapefile SHA256 was `08e341606e8391e458c3f08deb312de664b56bfae376064c5aa0aee6681a5f55`. Re-running its descriptive overlay requires the same independently supplied Natural Earth shapefile; the main genus-heldout climate/soil analyses do not depend on that polygon file.

## 6. Elevation is measured; true terrain is still untested

The original WorldClim-linked source has **point elevation** `wc_elevation_m` and site coordinates. It has **no spatially resolved original digital elevation model (DEM) neighbourhood** and no original **slope, aspect, topographic position index (TPI), ruggedness (TRI), valley/ridge index, topographic wetness index (TWI), landform class or catchment flow accumulation**. One elevation per observation CANNOT mathematically determine these terrain properties. The altitude + geography-only response reversals after atmospheric vapour adjustment (see SI S2.9–S2.10) are already documented: **altitude is not terrain relief**.

An optional owner-executable metric-raster terrain extraction and genuine synthetic correctness gate are committed as `scripts/analysis/mois_dem_terrain_metrics_gate_20261015.py`. This tool accepts **an actual independent, metric-projected DEM GeoTIFF and private photo-site coordinates**. Its Horn 3×3 slope/downslope aspect, 3×3 Riley TRI and 3×3 TPI formulas have been tested on synthetic flat, inclined, hilltop and valley raster cells. It deliberately **rejects latitude/longitude angular cell sizes**, unprojected DEMs, grid-edge and nodata positions. It outputs only safe distributional aggregate values, not georeferenced source sites. It does **not yet** link true terrain metrics to the original classified taxon data or fit a genus-conditioned floral-trait model.

For substantive inference, obtain an independently sourced DEM (for example **MERIT DEM ~90m**; Yamazaki et al. 2017, DOI [10.1002/2017GL072874](https://doi.org/10.1002/2017GL072874), access conditions apply; or Copernicus DEM GLO-90 subject to its separate data licence). Compute terrain metrics at a prespecified *local 90m and neighbourhood kilometre scale*, keep projected metre distances, check source-elevation mismatch, and evaluate the **same photo/taxon cohorts** with terrain added to the established geography/moisture/soil CV while holding photo rights and precise sites private. Aspect must be encoded as sine/cosine, especially when interacting with radiation; TWI requires a hydrologically conditioned flow accumulation algorithm, not point altitude subtraction. Until this is done, any biological association between flower colour and slope/aspect/topographic shelter or local water accumulation is **not tested**.

## 7. Editorial decision for the GEB manuscript

**Supported:** A globally positive source-conditional *interspecific photographed-white–water-climate* gradient at broader biogeographic scales, with substantial taxonomy/region heterogeneity, and an additional positive soil available-water-proxy coefficient within an exploratory multivariate source model. The former water-related predictors provide modest transferable signal in the genus-heldout matched-soil sample; the soil panel does **not** produce a stable incremental predictive gain.

**Not supported by these inputs:** A universal geographical transfer; unique causation by atmospheric vapour versus rainfall; a genetically balanced 50:50 flower-colour system; an environmentally selected anthocyanin or pigment benefit–cost crossing; an independent causal soil fertility mechanism; or **any actual terrain slope/aspect/TPI/TRI effect**. Keep the current main GEB positive water result and inferential scope, and retain this completed factorial-style geography/climate/soil/elevation **audit in the Supporting Information**, not as an inflated additional main-paper discovery.

### File and archive provenance

- Original source SHA256 photo `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0` and environmental `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`.
- Authorised full original-source verified replayer: `scripts/analysis/mois_full_biogeography_climate_soil_audit_20261015.py`; run with `--photo-zip`, `--env-zip`, optionally `--continent-shp`, `--output` and `--write` or `--check`. The full owner-side source run and exact local-result check both **PASSED**; GitHub hosted CI has **no owner-held original photo ZIP** and therefore performs aggregate/synthetic tests only.
- Safe, non-identifying subset report `results/mois_biogeography_climate_soil_elevation_full_audit_20261015.json`; the independent full-precision internal numerical original-source report remains available as a local generated artifact in the research conversation, **not licensed as a journal's stable source-data archive**.
- Source-only synthetic and exact-aggregate regression guards `tests/test_mois_biogeography_climate_soil_terrain_20261015.py`.
- No original imagery, attributed photo rights, full source species IDs/genus names, point geolocations or derived biogeographic private mapping rows were pushed. Submission is still **Draft/unmerged/unsubmitted**; licenced DOI/reviewer access and human coauthor/author sign-off remain outstanding.
