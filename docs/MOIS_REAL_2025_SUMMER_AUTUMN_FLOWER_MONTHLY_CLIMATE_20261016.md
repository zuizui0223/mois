# MOIS | Real original summer–autumn 2025 flowers × time-resolved moisture/light

**Owner:** \`zuizui0223/mois\` (NOT FCP's within-species W/D polymorphism paper).
**Evidence date:** 2026-10-16. Method intent is prior to any seasonal climate–flower-colour association estimation. The TerraClimate source probing and first images were already inspected, so this is an internal analysis protocol, **not a preregistered result**.

## Primary ecological question

Across botanically verified open flowers, are white-petaloid versus chromatic photographed traits allocated differently among actual months at the SAME taxon and photographed locality, and do those differences covary with source-month temperature, downward solar radiation, precipitation, atmospheric vapour deficit and climate water deficit? Is there spatial/observer bias or organ-stage substitution (fruit/bracts/spikelets) masquerading as seasonal white colour?

**No evidence of adaptation, maintenance costs, genotype frequencies, or within-individual flower switching is furnished by these photographs.**

## Real FCP upstream and actual MOIS post-blind inputs

- Exact 2026-09-27 original photo source, 42,111 taxon metadata source universe; the first 1,000 outcome-blind taxa contain **3,430,502 original photo ID opportunities**. This is not a flower-colour/photo-pixel census.
- For the season-stratified first original work unit, outcome-blind FCP [source run #38908308446](https://github.com/zuizui0223/fcp/actions/runs/38908308446) frozen-selected 128 original photos of \`Asclepias syriaca\` at original 0.25° cell lat-index 529 / lon-index 435, season-year 2025, summer/autumn 64 each (first wave 32 per season, other 32 per season reserved), 64 distinct independently held-out observers in each wave and no between-season repeated observer in the first wave.
- FCP [actual first-64 blind original ROI run #38909164068](https://github.com/zuizui0223/fcp/actions/runs/38909164068) completed **64/64 genuinely newly measured original photos**, 32 summer and 32 autumn; automated four-state classes 12 summer / 6 autumn; 3 same-ROI white–yellow floral-organ mixtures in summer; after excluding the same-ROI mixtures, provisional photo labels 7 white /2 chromatic summer and 2 white /4 chromatic autumn. The true independently reviewed open-flower-petaloid count is **zero confirmed (all 64 not yet reviewed)**. These are not a biological white-frequency contrast. **Automated ROI missingness is 20/32 summer, 26/32 autumn** (62.5% vs 81.25%). This sizeable asymmetric sampling/classifiability limitation is part of the primary scientific question.
- [MOIS authentic source bridge #38914197937](https://github.com/zuizui0223/mois/actions/runs/38914197937) independently verified the SHA-locked FCP source+ROI receipts and recovered the complete 64 source observations across **21 TerraClimate ~4km cells and 40 unique source observation-month×climate-cell keys**, spanning six calendar months June–November 2025. Raw photo-colour labels are stored separately from botanical anthesis confirmation.
- [MOIS two-botanist season-blind source-review form #38915040401](https://github.com/zuizui0223/mois/actions/runs/38915040401) successfully generated a blinded review sheet for all 64 photos; **all 64 remain UNREVIEWED**. Two independent botanists should classify open bloom and petal/tepal ROI; discordance requires independent adjudication, other floral organs are a secondary stratum, and unknown/fruit-only photos remain out of the primary colour denominator. Do not substitute absent iNaturalist observation annotations for nonflowering labels.

## Truly observed 2025 monthly climate, not 2025 normals or fake units

[TerraClimate Climatology Lab](https://www.climatologylab.org/terraclimate.html) states v1.1 uses ERA5 monthly anomalies with WorldClim climatologies. The official NetCDF physical variables are \`tmin,tmax\` (monthly average daily min/max °C), \`ppt\` (mm monthly), \`srad\` (W/m²), \`def\` (mm climatic water deficit) and \`vpd\` (kPa). See [original variable definitions](https://www.climatologylab.org/terraclimate-variables.html).

**Critical verified encoding trap:** The NCSS CSV and even its single-point NetCDF response contained raw packed \`tmin=11750\` for July 2025 at a sample source locality near Boston, while declaring \`degC\`. Reading the **official original provider DAS CF metadata** showed \`scale_factor=0.01\`, \`add_offset=-99.0\` for v1.1. Therefore \`11750*0.01-99.0=18.5°C\`. Raw \`11750°C\` is not a biological climate observation. For \`srad\`, the official source separately specified \`scale_factor=0.1\` and \`add_offset=0\`. Obtain each variable's current original CF metadata and validate physical ranges rather than assuming shared scales. [Provider endpoint probe](https://github.com/zuizui0223/mois/actions/runs/38912810211); [original CF DAS verification](https://github.com/zuizui0223/mois/actions/runs/38913504194).

\`scripts/analysis/mois_fetch_terraclimate_v11_2025_by_original_month.py\` requests all 21 genuinely observed climate gridpoints and six original variables, across 1991–2025 (126 site×variable histories). It applies the real provider CF scale and offset to EVERY NCSS packed value, matches 2025 measurements to EACH photographed observation-month and calculates **1991–2020 same-grid/same-calendar-month anomalies** independently of site averages. The data are reported only if **all six variables at all 40 observed grid-month keys are available**, physically plausible and cryptographically sealed. It does not compute a colour/climate effect and must never use one year's climatology as another year's actual monthly value.

## Main white-flower ecological estimand and limitations

A photo with RGB/ROI white label and unknown flowering stage is not a white open flower. Count distinct observation UUIDs, not image/pixel counts. Preserve original sampled and excluded/unclassifiable photos and their observers by season, year and geographical cell.

Only after independently verified open petal/tepal status and a sufficient number of both confirmed chromatic and confirmed white source observations across multiple site-years/species:
- estimate species/site-year matched photographic white fraction, observer-balanced folds, missing-photo-quality/anthesis model and spatial structure;
- regress within-site/time monthly anomalies separately from between-site average climates; geography, soil, elevation and topography act as adjustment variables, not time-varying seasons;
- compare physical unit-aware solar and water estimates without mistaking predictor collinearity or the seasonality of species photographic effort for a physiological adaptation;
- report weak/negative effects when supported rather than searching for a post hoc positive coefficient.

One focal site-year (*Asclepias* 2025) has perfect confounding of seasonal calendar month with 2025 seasonal climate. Even if all photo ROI labels were true flowers, it cannot by itself identify a shared climate response. The 39 taxa /268 geocell-year outcome-blind MOIS opportunity ranking is a future PHOTO workload, not evidence of 39 taxa with confirmed flower-colour seasonality.

## Reproducibility

- Real 2025 original-photo artifact: SHA256 of FCP post-blind 64 CSV \`1b16cdabd2936567dee87b8344019f40f315600c6510dfc1f25a03835024b6a0\`.
- MOIS 2025 real source bridge: FCP run #38909164068 -> MOIS run #38914197937; signed monthly keys SHA256 \`4a181971255218aa2645313aaecd3d16739c1dd00ca7a365d6e1b4e1d974024d\`.
- MOIS original precise climate provider code + explicit no-effect CI is on this branch. Refer to a real full extraction run only once its physical-unit metadata, months, 126 provider histories, provenance SHA and availability QC all pass.
