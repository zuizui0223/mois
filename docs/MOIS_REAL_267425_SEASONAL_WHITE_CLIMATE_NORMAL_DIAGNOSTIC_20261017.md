# MOIS: actual 267,425-photo seasonal white rates and matched within-site climate-normal evidence (2026-10-17)

**Scientific owner:** `zuizui0223/mois`. This is an interspecific/worldwide photographed-flower colour study and a *seasonal photographic measurement-quality analysis*, not the `fcp` paper's within-species `D/W/50-km` polymorphism effects. No original source photograph is newly downloaded or classified by this analysis, and no genetic morph, pigment fitness or selection is inferred.

## Actual signed source, not potential photos

- [Signed successful 267,425 measured original-photo union, MOIS run 39002957961](https://github.com/zuizui0223/mois/actions/runs/39002957961): **267,425 unique actual photo IDs**, **129,795 terminal four-hue automatic photo labels**. Original photo compressed CSV SHA256: `fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9`.
- Original *observed calendar date* available for **161,761 photos**; genuine automatic four-hue label AND date available for **81,867**. An unknown hue remains unknown, not coded nonwhite.
- [Independently completed original photo-site WorldClim v2.1 source, MOIS run 38992341847](https://github.com/zuizui0223/mois/actions/runs/38992341847): includes 248,312 photo IDs, of which **74,455** have BOTH classifiable hue and observed date. All WorldClim monthly values represent *1970–2000 month-wise climatological normals*, **not** original observation-year actual weather. Immutable raster join CSV SHA256: `43dcd2608f5a25e0a5d45d35810604a62ab659fb2cb9fb1a4e4257ef756d56c3`.
- The 267,425-photo extra 19,113 photo-ID cohort does **not** magically acquire climate values from the earlier WorldClim 248,312-photo artifact. The paired weather-normal inference applies only to literal photo-ID matches with that earlier geocoded artifact. The 267,425 dataset remains valid for real-dated source-only descriptive season white labels.
- Results were independently recomputed locally with `scripts/analysis/mois_seasonal_267k_climate_normals.py` against those signed original ZIPs, with five synthetic regression tests. Complete numerical result and all monthly-seasons are produced by the GitHub workflow; durable concise receipt: [2026-10-17 result JSON](../results/mois_real_267425_season_climate_normal_decision_20261017.json).

## WorldClim solar-unit bug found and corrected

The **already frozen v1** 248,312 WorldClim joined-photo artifact mislabels *unconverted WorldClim srad* (`kJ m−2 day−1`) as `Wm2`. The actual radiation flux is `W m−2 = (kJ m−2 day−1) / 86.4`. This is an authentic **units bug**, not a new ecological observation or change of recorded image labels. We explicitly divide original v1 energy numbers by **86.4** when fitting MOIS photo-month solar effects and update the underlying WorldClim sampler so future **v2** generated `..._Wm2` columns are correctly physical. Preserve original v1 artifact SHA unmodified, do not silently relabel source bytes. Temperature (°C) and monthly rain (mm) stay in their original physical units.

## Actual northern and southern extra-tropical observed season ratios

Each numerator is *automated photographic white ROI labels*. Denominator includes *only classifiable white or three chromatic label states*, not all date-bearing photos or biologically verified flowers. Seasons are hemisphere-reversed, tropics retained as `month_only`; December is grouped in the *following* meteorological season-year for matching.

| Temperate region | Local winter | Spring | Summer | Autumn |
|---|---:|---:|---:|---:|
| Northern-hemisphere photo-white share | 1,481 / 3,246 = **45.63%** | 13,179 / 27,712 = **47.56%** | 13,312 / 27,210 = **48.92%** | 3,358 / 7,089 = **47.37%** |
| Southern-hemisphere photo-white share | 765 / 2,098 = **36.46%** | 2,213 / 5,168 = **42.82%** | 1,627 / 3,030 = **53.70%** | 758 / 1,618 = **46.85%** |

Photographic **summer** white proportions are high in both extra-tropical hemispheres. The table does **not** hold sampled species, coarse geography, original photographers, photo organ visibility, or source classifiability constant, so it is **not** evidence that individual flowers, populations or pigments causally change colour under summer sun. Global unclassified and missing-calendar-date photographs remain outside the colour-ratio denominator, without imputation.

[Complete source-photo attempt, classified, taxon counts and white numerators](../results/mois_real_267425_dated_photo_season_white_rates_20261017.csv).

## New nontrivial result: photographed taxon composition explains much of the Southern summer peak

The full dated/classifiable Southern extra-tropical photo dataset shows **+17.23 percentage points** summer minus winter, but many taxa were photographed in only one season. To distinguish a *change in which taxa get photographed/how many times* from a *within-named-taxon change in its labelled white fraction*, use only exact common original taxon IDs photographed in **both** summer and winter and decompose the pooled photographic contrast by the **symmetric Kitagawa identity**. A whole-taxon nonparametric 999-bootstrap provides *exploratory* percentile intervals, not familywise uncertainty.

| Among common photographed taxa (≥2 independent observation IDs in each season) | Southern hemisphere, **108 taxa**, 976 summer / 845 winter photos | Northern hemisphere, **238 taxa**, 2,455 summer / 1,318 winter photos |
|---|---:|---:|
| Total original-photo weighted white summer−winter difference on identical taxon support | **+14.47 pp** | **+3.69 pp** |
| **Photographic taxon-frequency mixture component** | **+11.91 pp** (taxon bootstrap 95% **[+0.19,+23.64]**) | +1.34 pp (95% [−6.07,+8.46]) |
| **Within-taxon labelled-white-rate component** | **+2.56 pp** (95% **[−1.44,+6.92]**) | +2.35 pp (95% [−0.76,+5.77]) |
| Taxon-equal mean summer−winter white share | +2.80 pp (95% [−2.29,+7.62]) | +3.69 pp (95% [+0.55,+7.12]) |

Thus **about 82% (11.91/14.47) of the Southern photographed white summer−winter difference *within the 108 common taxa*** can be algebraically assigned to changing their *relative photographed frequencies*, not to a changed individual-taxon white-photo rate. This **does not** mean 82% of the whole-atlas original +17.23-point difference has been causally explained: season-exclusive taxa lie outside the common support. With a lower one-observation threshold (212 common Southern taxa), the taxon-frequency composition term remains +11.76 pp [1.83,21.26]; with a more stringent four-observation threshold (49 taxa), its interval includes zero. This threshold sensitivity and unrestricted movement of identical taxa among sites and years forbid a universal/causal inference.

**Ecological reading:** A coarse summer photographed-white peak can be generated largely by seasonal *sampling or community/taxon-frequency composition* even if no individual taxon's flower pigment changes. It is a concrete reason to control taxon and geography before interpreting seasonal sun/water correlations. This is a photograph-sampling decomposition, not the FCP paper's within-species colour polymorphism W/D, reproductive success or proof of a taxon-season biological shift.

[Signed same-original-photo taxon-frequency result](../results/mois_real_267425_summer_winter_taxon_photo_composition_decision_20261017.json); reproducible original-source script `scripts/analysis/mois_taxon_season_composition_decomposition.py`. The CI recreates the full result from the same 267,425-photo run and preserves the complete 999-bootstrap numeric output as a run artifact.

## Geographic and species controls change the interpretation

For an observational *within-site seasonal* contrast, preserve only photographs from the **same taxon × 0.25° coarse cell × original meteorological season-year**, and require **≥2 independent observation IDs in each of at least two seasons** and all three measured authentic monthly climatology variables. This retains just **61 eligible site/taxon/season-year groups, 47 taxa, 384 actually classified photographs**, including 29 groups with *both* automated white and chromatic hues. Distinct observation IDs do not guarantee distinct observers.

Fit **three separate exploratory** within-site/year demeaned white-probability-on-climate regressions, with equal weights per site/taxon/year group; bootstrap the **whole taxon** 999 times using fixed PRNG seed `20261017`, reporting percentile 95% intervals, not familywise maxT. This is **month-specific background climatic normal**, not weather anomaly or causal adaptation.

| Within taxon/coarse site/year, background climatology | White photo-label share slope | Taxon bootstrap percentile 95% |
|---|---:|---:|
| Average *calendar-month* temperature, +5°C | **−0.27 percentage points** | [−6.93, +6.34] |
| Average monthly precipitation, +50mm | **−0.93 percentage points** | [−10.24, +7.78] |
| Average monthly incoming solar flux, +50 W/m² | **−0.08 percentage points** | [−4.81, +4.21] |

**None of these three matched local climate-normal associations has an interval excluding zero.** The preliminary global hemisphere summer peak and matched locality estimates are therefore **not a validated common seasonal temperature/rain/sun effect**. The year × month atmospheric anomalies of actual site-months have **not** been retrieved for all these photos; WorldClim 1970–2000 normals cannot resolve anomalous weather within the same calendar month across observation years.

With a looser minimum of one independent observation ID per local season, there are 852 matched groups, 359 taxa and 2,372 classified photos, but all three exploratory 95% intervals again include zero. Tightening to three IDs per season leaves only 12 groups / 8 taxa / 118 photos, much too little independently supported botanical evidence for a robust multi-testing finding. Do not select the sparse positive solar point estimate as proof from post hoc searches.

## Inference hierarchy and next paper-building work

1. Preserve the positive *historical global **between-species** annual water* trend separately (main 267,425-photo six-anchor/maxT), recognizing it fails the congeneric coarse-cell-conditional test.
2. Distinguish *time-invariant* soil, elevation, topography and coarse geography from *time-varying* observed-month climate; same-site/year comparisons absorb fixed site backgrounds automatically, so cannot separately identify soil/elevation/terrain coefficients.
3. Retain all source photo attempts, photo-classifier missingness, observation IDs and unknown original floral stage. Current iNaturalist Plant Phenology annotations are sparse and observer/taxon/season-dependent, and photographic ROI hues do not prove independently verified anthesis.
4. For a truly temporal claim, combine real **observation-year/month TerraClimate** temperature, radiation, precipitation, climatic water deficit and VPD with sufficiently sampled same-taxon/site repeated seasonal photos; use independent floral organ verification and observer/effort diagnostics. Then revisit multivariable geographic/taxonomic conditioning, positivity and ecological interpretation. Do not relabel WorldClim normals as actual weather or classify all 3,430,502 photo opportunities as blooming flowers.

**Current decision:** the global photo atlas contains a descriptive seasonal white-hue compositional pattern; the existing authentic *matched taxon and locality* photo/climate-normal evidence cannot identify a universal seasonal temperature, sunlight or moisture effect. This is **the result**, not a reason to invent nonexistent weather data or turn the FCP intraspecific W/D paper into this analysis.
