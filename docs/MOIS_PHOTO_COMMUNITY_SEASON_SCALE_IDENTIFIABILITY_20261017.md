# MOIS: Can a local photographed floral-colour assemblage be compared between summer and winter?

**Source-identified 17 October 2026 analysis.** Original *actually classifiable* photos (267,425 original attempted pixels /129,795 four-hue labels) were combined **without changing images or colours** with **9,597 independently matched original photo-ID + numeric observation-ID public dates** from the immutable 9,600-historical-source pilot. This produces **86,172 truly dated four-hue photographic ROI labels**, of which approximately **80,083** are in extra-tropical regions. They are not independently botanist-confirmed open petal observations, and absent citizen-science photographs are not an observation of absent flowers.

The original MOIS moisture × floral hue inference remains the **among-named-taxon worldwide genus-equal historic climatic-normal** question. FCP separately owns the *within-species* phenotype-space/spatial polymorphism analyses. Here the estimand is **photographed floral-hue composition at a coarse place and date**, not whether an individual plant changes pigmentation in response to seasonal water or light.

## Compare two biologically different temporal designs

1. **Same named taxon, same place and season-year:** in the source with 86,172 date + four-colour photographic labels, only **one Northern** and **zero Southern** 0.25° geographic site × taxon × same local season-year groups possess at least two distinct original observation IDs in BOTH summer and winter. This design cannot identify a broad season-driven within-taxon pigment response.
2. **Photographed community/all taxa at the same coarse site × season-year:** this is an appropriate descriptive unit for asking whether the *mix of white-labelled photographed taxa* at a location changes with observation season. It explicitly **allows different named taxa** in summer and winter. It is therefore not an individual/persistent-population plasticity contrast. A meaningful source comparison must retain multiple distinct original observations **and** multiple named taxa in each season.

We enumerated full original photographed colours at five nominal latitude-longitude grid widths (0.25°, 0.5°, 1°, 2°, 5°), never merging across meteorological season-years, reversing seasons south of 23.5°S, and assigning December to the next winter/summer year. The strict rule is **≥5 distinct original iNaturalist observation IDs and ≥2 different named taxa in each summer/winter season**, with all source photos being four-hue classifiable.

| Geographic site-cell width (degrees) | Northern comparable cell×season-years | Southern comparable cell×season-years |
|---|---:|---:|
| **0.25°** | **2** | **1** |
| **0.5°** | 28 | 12 |
| **1.0°** | 70 | 31 |
| **2.0°** | 108 | 58 |
| **5.0°** | 120 | 74 |

The dramatic increase as cell width expands does **not** establish true local seasonality; at 5° it can combine sites several hundred kilometres apart. At 0.25°, the three usable cell×year groups comprise just **42** source photographs across the two seasons and geographically cannot establish a universal colour–weather relationship.

A slightly looser **≥3 distinct observation IDs and ≥2 taxa in each season** at 0.25° yields **42 Northern / 21 Southern** year×cell groups with **415 / 187 original source photographs** in both seasons respectively. Before the 9,600 historical photo date recovery the same method yielded **39 Northern /16 Southern** groups (the corresponding old-source photo denominators 387 /138). Thus the extra 4,305 dated original historical colour labels add a little **community-level** geographic-temporal overlap despite adding *zero* strict same-taxon local/year summer–winter repeatability. This is a quantitatively important distinction between taxon-specific and local community change.

## Sensitivity of a photographic summer–winter contrast

For a candidate site×year, estimate

\[
\Delta_{c,y}=100\left(\frac{W_{\mathrm{summer},c,y}}{N_{\mathrm{summer},c,y}}-\frac{W_{\mathrm{winter},c,y}}{N_{\mathrm{winter},c,y}}\right).
\]

Only classifiable original photos enter \(N\), and a technical unknown is not counted as coloured or white. An equally-weighted site×year contrast is a **descriptive average** over the eligible source sites, rather than an organism- or population-level photosynthetic response. An exploratory **999-resample geographic-cell cluster bootstrap** groups repeated years at the same approximate location. This does not yield a simultaneous interval correcting for the several geographic and sample-size choices, nor adjust observer effort, geographical dependence between neighbouring coarse cells, time of day, species-selection or unverified flower organ identity. Results at very small independent-cell counts are explicitly held rather than reporting misleading confidence intervals.

The full analysis and its precision/bootstrap interval output are produced, together with original-photo source hash and source denominator checks, by `scripts/analysis/mois_original_photo_community_season_geoscale.py` and `.github/workflows/mois-real-86172-community-geocell-season-year-source.yml`. They use only **already completed** official photo-colour and ID-verified metadata workflow artifacts; no additional iNaturalist photo pixels or provider API requests are made in this analysis.

## Ecological consequence for the single MOIS manuscript

The source can support **global geographic and seasonal photographic community colour composition**, with explicit sampling-effort/taxon-identity uncertainty. It cannot yet support a universal statement that *local flowering plants turn white or change pigment expression because of seasonal hotter, drier, brighter weather*. Background SoilGrids, elevation, terrain and lat/lon can be described or controlled as static site differences, but they cannot substitute for actual observation-year monthly TerraClimate radiation, temperature, precipitation, VPD or climatic-water-deficit anomalies.

**Next paper-worthy target:** after the full original 100,543-photo public month snapshot restoration, re-evaluate the **0.25° same-cell × same-year community** source photo support and the proportion of validated blooming petals, while retaining source effort and taxon composition. Enlarging date coverage may help geographic-temporal assemblage comparisons; only genuine local repeated flower-colour observations can identify *within-taxon* temporal change. FCP's independent within-species W/D and MOIS's historical genus–water maxT are not imported into this inference.
