# MOIS GEB Research Article — display items, figure legends and evidence ledger

**Editorial freeze: 27 October 2026.** Use the companion machine-readable file **MOIS_GEB_20261027_EVIDENCE_NUMBERS.json** as the plotting-data contract. This is a figure production brief, not a claim that all new visual figures have already been rendered. Positive and null/disagreement results are included together. Never plot a global white-flower population mean from nonprobability-selected 1.18m photo labels.

## Main conceptual sentence for every display

**Photographed floral colour shows environmental structure, but that structure depends on the taxonomic, geographical and seasonal support of the sampled flora and the ability to measure colour from source images.**

### Table 1. Which actual original-photo sample answers which ecological question?

| Evidence source / nominal scope | Original photo ID attempts | Four-hue ROI outcomes | Strict nine-component ROI outcomes | Manuscript use |
|---|---:|---:|---:|---|
| Entire original 93,209,351 metadata/photo-opportunity source; 42,111 nominal study taxa | Not processed at full global scale | Not assessed | Not assessed | Source frame ONLY; no botanical denominator |
| Six annual climatic anchors with genus-bootstrap/maxT | **309,439** | **149,662** | Different analytical criterion, not used | Main global annual water and geography-conditional coefficients |
| Geography/habitat source with island/elevation matches, source landscape and four hue × climate cells | **313,535** | **151,239** | Not this model's outcome | Main elevation, island; four-hue supplementary gradients |
| Descriptive photo/season/taxon mixture after all tail workers | **350,399** | **164,297** | Not this model's outcome | Photo vs taxon vs genus rates; four seasons; sampling-arm comparison |
| Frozen four-arm strict original-image union (historic/K1/FCP + 2,048 + P2 4,096 + bulk128 131,072) | **1,181,146** | Do not infer from strict9 count | **451,643** | Measurement availability, 42,111-source-taxon gap, region/season nonresponse audit |
| Soil-complete 13-covariate exploratory subset of older source | **26,163 classifiable photo records** | **16,801 named taxa/1,665 genera** in its main filtered set | Not this estimand | Supplementary in-sample variance/parameter-scale sensitivity |

**Do not add rows in this table**: cohorts overlap historically, even though the 1,181,146-photo four-arm union was explicitly deduplicated internally. The later photo source does not provide refitted rainfall coefficients.

## Figure 1. The composition of photographed global floral colour depends on taxonomic weights

**(a)** Three separate photo classes plus white in the 350,399-photo source: white 78,026; yellow/orange 51,491; red/pink 17,382; blue/purple 17,398. State denominators: 164,297 machine four-hue classified / 186,102 original photos unknown. Do not claim anthocyanin from "blue/purple".

**(b)** White fractional composition among classifiable original photographs: photo-equal 47.4908%, taxon-equal 48.1891%, nominal-genus-equal 51.4631%. Use point estimates without fabricated intervals. The three weighting rules define three different observational summaries.

**(c), optional small panel rather than adding a full main chart:** 1,181,146 four-arm source ROI9 quality audit: 451,643 valid strict9 / 729,503 unknown. Mark this as a DIFFERENT eligibility criterion and later source, not a new observation of the 350,399-photo ecological outcome.

**Interpretation:** roughly half the classifiable *photo hues* are white, but representation can change its numeric share. This is not a global evolutionary colour equilibrium.

**Immutable evidence:** docs/MOIS_REAL350399_ACTUAL_RGB_FOURTERM_SEASON_AND_TAXON_20261019.md and the actual 1,181,146 strict9 release:
https://github.com/zuizui0223/mois/releases/tag/mois-original-1181146-exact-photo-crossarm-451643-strict9-taxongap-20261027-v1

## Figure 2. Positive worldwide moisture associations lose signed certainty with geographic conditioning

**Recommended: two paired forest panels, common x-axis of pp per original climate SD, vertical 0 line.**

| Variable | Model/comparison support | Estimate (pp/SD) | 24-family simultaneous 95% |
|---|---|---:|---:|
| Annual precipitation BIO12 | Genus-weighted global | **+1.9884** | **[+0.5525,+3.4244]** |
| Annual precipitation BIO12 | Genus × coarse geographical cell | +0.9034 | [−2.6350,+4.4417] |
| Annual atmospheric vapour | Genus-weighted global | **+3.1606** | **[+1.0420,+5.2792]** |
| Annual atmospheric vapour | Genus × coarse geographical cell | +2.1545 | [−2.6404,+6.9494] |

All estimates: actual **309,439 original photos /149,662 classifiable original four-hue** outcome cohort. 999 whole-genus bootstrap resamples; 24-family maxT intervals within this frozen analysis. Not a randomized climate experiment; no all-study multiplicity control. In figure, visually separate global and geography-conditioned support from mere covariate adjustment under identical support.

**Immutable evidence:** results/mois_real309439_fullFCP256_six_climate_maxT999_result_20261017.json and
docs/MOIS_REAL309439_FULL256_ORIGINAL_GLOBAL_MOISTURE_MAXT_20261017.md.

## Figure 3. Taxonomic and photographic assemblage context

**(a)** Older original identical-source **20,694 taxon/1,934 genus** comparison: photographed-white annual precipitation +8.555 pp/SD among sampled genus means vs +2.255 within nominal genera; paired difference +6.300 [3.666,8.779] (ordinary whole-genus-bootstrap 95%; NOT Fig. 2's maxT family). The "among" and "within" estimates share one original source but differ from the later 309,439-photo model; do not pool coefficients.

**(b)** Distinct 350,399-photo experimental-source expansion: old core 317,631 photo attempts (white 47.2266% of 152,827 classified) vs 32,768 deliberately breadth-targeted photos (white 51.0113% of 11,470 classified); raw difference +3.7847 pp. Genus × 10° geographically matched new-minus-old estimate +0.635 pp (exploratory genus-bootstrap 95% interval −1.847,+3.070), but only 6.05% of core versus 24.90% of new taxon–cell occurrences enter that matched support. Raw and matched contrasts are DIFFERENT populations and weights; **do not create a waterfall causal contribution** by subtracting the two.

**(c), optional SI inset:** The 13-variable soil-complete 16,801-taxon/1,665-genus sample had environment-only in-sample R²=2.000% (genus-equal) and 1.282% (species-equal), while fitting each sampled genus separately had 44.315% in-sample R² versus genus-size-preserving shuffle null 19.063%. Very different complexity: 1,665 genus parameters. The contrast is sampling-frame trait clustering, not a phylogenetic heritability estimate.

**Evidence:** older source docs/MOIS_INTEGRATED_ECOLOGICAL_FINDINGS_20261016.md on manuscript/geographic-conditional-moisture-research-20261014; docs/MOIS_REAL350399_GEOGENUS_SOURCE_SELECTION_BIAS_20261019.md; older 13var source docs/MOIS_GLOBAL_OVERALL_WHITE_13ENV_AND_GENUS_SCALES_20261016.md in older manuscript PR #1.

## Figure 4. Different environments do not share a universal white direction

Two explicitly unmatched estimands, same plotting axis in WHITE proportion difference (pp), and separate subpanel/shape notation.

**(a) Elevation**, high >=1,500m versus low below 500m one-degree assemblage cell averages, same landmass and 5° neighbourhood, 44 local neighbourhoods; high−low white = **−5.085 pp**, nominal 95% **[−8.785,−1.386]**; 31/44 same sign. 5° neighbourhoods are not mountain ranges and unequal taxon/observer support remains.

**(b) Islands**, one-to-one climate/elevation-matched 102 island–mainland one-degree source cell pairs, 65 island landmasses, 28 coarse regions; island−mainland white = **+1.396 pp**, region-cluster nominal 95% **[−1.553,+4.345]**; equal-region estimate +0.092 pp. Do not claim no island effect; say a universal island white excess is not established on this matched support.

**SI:** matched island red/pink −2.452 pp [−3.872,−1.032] in pair-equal nominal interval, but region-equal alternative includes zero; exploratory and no project-wide correction. Elevation chromatic categories jointly sum to +5.085 pp, no pigment concentration inferred.

**Evidence:** docs/MOIS_ENVIRONMENT_OVERLAP_AND_REGIONAL_ELEVATION_20261018.md, results/mois_macro_environment_matching_20261018.csv.

## Figure 5. The seasonal photographic white difference contains a strong taxonomic source-weight component

Two bar/waterfall plots, accurately reconstructing one observed white-labelled percentage-point difference with the same symmetric four-term photo composition identity. Seasonal units: summer June–August in north / December–February in south; winter reversed; extra-tropical photographs only.

| Exact photographic source component | North summer − winter, pp | South summer − winter, pp |
|---|---:|---:|
| White photo share within taxa seen in both seasons | +3.26036 | +3.57058 |
| Relative source photographic weights of common taxa | −0.85749 | **+8.03618** |
| Taxa only supported by summer photos | +2.92067 | +4.15779 |
| Taxa only supported by winter photos | −2.49426 | −5.26831 |
| **Exact photographed summer − winter total** | **+2.82928** | **+10.49624** |

For the South, 7,796 summer and 5,409 winter four-hue-classifiable original photos; original white photo rates 50.7440% and 40.2477%. The frequency term is not a causal phenological or plant-community effect; it also changes with photographer effort and image selection.

**Caption warning:** the within-taxon term pools taxa across potentially different sites and years. In same taxon×0.25° cell×same meteorological season-year source groups, just **two north + two south** have >=2 distinct original observation IDs per season. Do not infer within-plant flower-colour plasticity or actual month weather effects.

**Evidence:** docs/MOIS_REAL350399_ACTUAL_RGB_FOURTERM_SEASON_AND_TAXON_20261019.md and source release
https://github.com/zuizui0223/mois/releases/tag/mois-real350399-observed-season-taxon-photographic-mixture-v1

## Supporting Information (SI): results that strengthen validity without enlarging the central causal claim

- **S1. Cohort and source identity table.** Exact 93,209,351 photo-ID opportunities and the four-arm 1,181,146 genuine original images. Original photo = row grain; historic additional photos of the same observer observation are not independent visits. SHA receipts, licensed photo access and exclusion of the unrelated 64 photo Lab-only pilot.
- **S2. Four-hue versus strict9 classifications, missingness and botanical review boundary.** 350,399 four-hue classifiable 164,297; distinct method 1,181,146 strict9=451,643. The strict9 missingness is 729,503; geographically N temperate 305,289/726,181 (42.0%), S temperate 79,937/235,974 (33.9%), N tropical 30,579/106,008 (28.8%), S tropical 14,426/60,707 (23.8%). Lat unknown 195, season unknown 197 in georestored v2.
- **S3. Nine-mass image quality and review selection.** 512 previously strict9-invalid source photos, 495 YOLO500 candidate boxes, 14 negatives, three errors, selected for missingness/region/taxon diversity and thus NOT a global failure proportion. The botanical organ-stage labels from 128 chosen photo IDs remain blank pending human review; do not promote box candidate to verified corolla.
- **S4. Photo classifiability / explicit missing-not-at-random scenarios.** Retain all original attempted photos in denominators and disclose measurement gaps by source acquisition design, latitude and local season. On a **SEPARATE 313,535-photo original four-hue source (313,204 rainfall+solar complete; 151,107 known; 162,097 unknown)** the white-photograph wet−dry complete-case contrasts by solar level are **+3.93, +5.33 and +4.37 pp**. Under dry-side MAR, decreasing the **hypothetical** unknown wet-side white probability by **7.68, 9.74 and 7.44 pp**, respectively, suffices to erase these descriptive contrasts; solar source-attempt-standardized observed contrast +4.55 pp flips near an **8.25 pp** hypothetical wet-missing deficit. Fully signed result and detailed source MNAR assumptions: [real313535 nonresponse sensitivity](https://github.com/zuizui0223/mois/releases/tag/mois-real313535-white-solar-wetdry-MNAR-tipping-source-20261027-v1), plus [SI scientific limitations](MOIS_GEB_20261027_BIAS_ADJUSTMENT_AND_MNAR_SUPPLEMENT.md). **This is not a refit or corrected CI for the distinct genus-adjusted 309,439-photo annual moisture model.** Do not assign any unknown photo a biological white/nonwhite label.
- **S5. Genus/genera × geography positivity.** The previous taxon-broadening arms' apparent white excess shrinks under taxonomic matching, but few overlapping taxa are comparable. Distinguish genus-only, source-taxon and taxon×site requirements explicitly.
- **S6. Annual light/moisture/hot four-hue source exploratory model.** Use actual 313,535-photo true-image-feature climate mixture source: https://github.com/zuizui0223/mois/releases/tag/mois-original313535-four-hues-water-solar-genus-20261018-v1. It has four-hue wet/dry within solar-tertile composition and heldout genus/10° source-cohort Brier. Unless checked with the exact released aggregate CSV, do NOT invent a "wet red +x%" or unseen-genus Brier advantage; this pilot is supplementary to verified annual white coefficients in Figure 2.
- **S7. All 13 environment terms, within/among genus and descriptive SHAPLEY 5 blocks.** Display exact 16,801-taxon original 13var results only with explanatory/sample provenance and clear in-sample R² vs many-parameter genus model. Do not refer to soil/TPI/TRI as all observed on every global photo.
- **S8. Island matching calipers and high-elevation neighbourhood sensitivities.** Report 55 pairs under caliper 0.5; 141 under caliper 2; 86 after removing very large islands. Report region vs pair weighting and sparse high-altitude neighbourhoods rather than searching for a preferred threshold.
- **S9. Ecological colour-interpretation restrictions.** sRGB photo class labels are not UV reflectance, true petal pigment chemistry, flowers-without-bracts, genetically heritable white morphs, pollinator availability, ecological fitness, or a local reciprocal transplant experiment.
- **S10. 2026 source revisions.** Original 131,072 bulk photo lat/lon/day had been accidentally omitted from first strict9 world union; use the SHA-signed georestored v2 parquet with all 131,072 restored. Original exact photo and palette values never changed. Old geographic QC v1 must not be cited.

## Caption language to retain

- “Photo classifier white label” rather than “naturally white flower”, wherever the value is not independently botanical-anthesis validated.
- “Source SD” (from each archived analytical cohort) rather than one shared physical climate response SD across all population subsets.
- “Nominal 95%” for exploratory matching tests, “24-family maxT simultaneous” only for the corresponding 309,439-photo climate models.
- “Observed photo composition” rather than within-individual seasonal colour plasticity.
- “Genus membership” rather than molecular phylogenetic signal.
- “Source-specific association” rather than causal environmental pigment selection.
- A 500px colour BOX is not the same estimand as original full resolution FCP strict ROI9.
- No world-wide native flora estimate is available without botanically verified native status and inclusion probabilities.

## Source-to-claim index

| Main claim | Frozen source or actual outcome run |
|---|---|
| ~47.49% photographed white, four colours and taxon/genus weighting | docs/MOIS_REAL350399_ACTUAL_RGB_FOURTERM_SEASON_AND_TAXON_20261019.md |
| +1.99 annual rainfall and +3.16 vapour global, conditional intervals spanning zero | results/mois_real309439_fullFCP256_six_climate_maxT999_result_20261017.json; [original full model run 39029308911](https://github.com/zuizui0223/mois/actions/runs/39029308911) |
| Source genus-scale contrast +8.56 among vs +2.26 within | older main draft original source results/mois_genus_between_within_moisture_20261015.json on manuscript branch |
| +3.78 raw source-arm gap versus +0.64 genus/geography-matched | docs/MOIS_REAL350399_GEOGENUS_SOURCE_SELECTION_BIAS_20261019.md |
| Elevation −5.09 / island +1.40 ± nominal | docs/MOIS_ENVIRONMENT_OVERLAP_AND_REGIONAL_ELEVATION_20261018.md; results/mois_macro_environment_matching_20261018.csv |
| Southern +10.50 summer-white photographic difference with +8.04 photo-weight | docs/MOIS_REAL350399_ACTUAL_RGB_FOURTERM_SEASON_AND_TAXON_20261019.md |
| Strict9 451,643 of 1,181,146 and covered 30,048/42,111 source taxa | [exact four-arm cross-photo identity union](https://github.com/zuizui0223/mois/releases/tag/mois-original-1181146-exact-photo-crossarm-451643-strict9-taxongap-20261027-v1) |
| 131,072 source geo/day restoration; regional missingness | [real georestored v2 release](https://github.com/zuizui0223/mois/releases/tag/mois-1181146-originalROI9-128bulk-georestored-seasons-and-512QC-20261027-v2) |
| 512 old failed photos detector quality and still-unverified botany | [actual 170+342 reused technical diagnosis](https://github.com/zuizui0223/mois/releases/tag/mois-old512-realROI9unknown-technical-170cached342screen-20261027-v1) |

**Internal control:** If a plotted point or colour/source caption cannot be traced to one of these source-specific files, omit it instead of filling the panel with simulated or borrowed ecological slopes. None of these results substitutes for the rights-controlled source dataset needed for final paper archiving.
