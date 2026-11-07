# MOIS ALAN × flower colour: actual 2020 VIIRS/WorldPop result (5 November 2026)

**Status: empirically executed, exploratory, noncausal.** No claim of white-flowered species disappearance. No new flower photographs downloaded, no ROI inference, no MNAR tipping scenarios, no modifications to main.

## Main empirical comparison

- Existing 2020 original historical-main photo-observation source: 56,746 attempted photos / 23,854 strict ROI9 successes. Original source SHA256 f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c.
- Prespecified independent 0.05-degree spatial cells need ≥5 attempted observations, ≥3 attempted nominal taxa, ≥3 strict ROI9 observations, ≥3 strict ROI9 taxa. **1,059 cells / 9,966 source attempted photos / 5,087 strict9 successes**. Distribution of cells: North temperate 855, South temperate 145, North tropical 49, South tropical 10. Do not call this global flora.
- Real 2020 VIIRS *monthly* upward night-light radiance via public ArcGIS NighttimeLightsMDNB ImageServer at cell centres; actual cloud-free counts verified, require ≥2 observations in each accepted month; quarter months January/April/July/October (not true official annual composite). Units nW/cm²/sr.
- 2020 WorldPop 1km population-density service accessed independently and year-checked; units people/km². Exactly 1,044 cells had VIIRS+population measurements (9,840 photo attempts / 5,023 strict9).
- Outcome-blind same-5-degree-block 1:1 pairs require |log1p population density difference|<=0.6, BIO5 heat difference<=3°C, log rainfall<=0.5, elevation<=300m, cell separation<=3 degrees, log1p VIIRS contrast>=0.65. 350 potential pairs before exclusion/reuse, **128 independent matching cell pairs** from 36 five-degree blocks.
- Bright-group median radiance **18.19**, dark **4.33**; bright-group median population density **590**, dark **515** people/km². These are matching calipers, not experimentally randomized equal populations.
- Original automatic strict9 white ROI mass, **equal nominal-taxonomic weight per cell**, bright minus dark **+1.4649716 percentage points** (median −0.2550). Geographic block bootstrap nominal 95% interval **[−4.375,+7.018] pp**; 63 pairs bright more white and 65 bright less white. Classification success rate bright−dark **−1.2209 pp** with spatial block interval about [−4.72,+2.45] pp.
- **Conclusion:** no consistent lower white ROI fraction at high artificial night-light intensity in this selected 2020 photographed source. Not proof of no physiological/reproductive ALAN effects or lack of historic losses.

Source: [2020 actual GIS + photo run 41946071715](https://github.com/zuizui0223/mois/actions/runs/41946071715), artifact 13379405167, ZIP SHA256 49a8ee482198507d186cca45ad20aa4a8608a96ea38e49b7422ceca10c4d8f90. Independent six unit tests and all original source/output SHA checks passed.

Exact numeric artifact outputs:
- candidate_source_cells_2020.csv SHA256 d8d0b7d5327a9add2107b2c6a1fa0e5aa6358d0813ef90da6f70e478e7a4f336
- population_matched_bright_dark_pairs_2020.csv SHA256 32b53d23d4beac39846e3c618f7a4f5d8b24164a126cf63c4cc530b2e720aae6
- viirs_population_source_matched_cells_2020.csv SHA256 3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3

## Independent documented *species-level* flower colour

To avoid relying only on automatically segmented image ROI colour, species in the existing source were linked by verified exact WCVP accepted species identity to a **separate direct High/Medium Island floral trait ledger**, excluding mixed/unknown/non-direct/cultivar source records.

- Kew WCVP accession source: [run 39225659840](https://github.com/zuizui0223/mois/actions/runs/39225659840). Name mapping SHA256 dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2; accepted species names SHA256 9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3.
- Island scientific direct species trait input: [run 34191508045](https://github.com/zuizui0223/island/actions/runs/34191508045). Direct species trait ledger SHA256 7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577.
- Within the 1,059 preselected cells, **3,550 photo attempts** link to a directly documented single-state flower colour (covering 1,405 nominal original taxa). One accepted species counts only once per cell regardless of photographs/nominal synonyms. Missing and multistate direct traits are neither white nor nonwhite.
- Only **47 of the 128 identical ALAN/population matched pairs** have ≥3 independently documented accepted species each side. Across these 47, documented white species fraction bright minus dark **+6.6089466 pp** (median +6.6667). 26 pairs positive, 15 negative, 6 tied; spatial block bootstrap 95% **[+0.512,+13.180] pp**. The photographic ROI white fraction on those *same 47 matched pairs* is +4.367 pp.
- This independent measurement is still **species-level literature description**, not confirmed photographed flowers or a moth-pollination guild; its documented coverage selects a very different population. It is not an independent confirmation of ecological causal protection by artificial light. Some provenance is third-party non-redistributable: **never copy individual documented species-colour rows into this PR/artifact.**

Source: [independent 2020 result run 41947970697](https://github.com/zuizui0223/mois/actions/runs/41947970697), six new tests and all six frozen SHA pins passed; aggregate-only artifact 13378983852, ZIP SHA256 1632490802ab3083faec118d14bfaed36dbf60d259dc4850e86634a8d3fb2fc9. The only emitted independent trait output is aggregate JSON.

## Scientific scope and interpretation

Published field experiments have demonstrated reduced night pollinator visits under artificial light; **they do not imply** global loss or local extirpation of white-flowered plant species. Many white flowers have daytime pollinators, many moth-pollinated flowers are not white. This study has no species absence observations, no species temporal extinction/recolonization series, no validated nocturnal moth guild, and no natural experiment lighting onset. It tests a cross-sectional photographed assemblage association only.

Population-density matching plus temperature, rainfall, elevation, proximity, and broad region cannot remove urban land cover, habitat management, gardening/ornamentals, other pollutants, observer effort, seasonal flowering, correlated classification error or regional genus representation. VIIRS upward radiance does not directly measure ground-level ALAN exposure, spectral composition, or pollinator visual contrast. Cell-centre raster exposure is imperfect for individual flowers. If causal extinction is the target, source-backed multi-year site records, measured light installation/shutoff and moth-specialization are necessary.

All 95% bootstrap intervals describe spatial-block variation in the selected study support, not global source missingness/causal uncertainty. The work is a **retrospective exploratory** test and does not replace any previously reported positive macroecological patterns. White ROI mass and botanical species colour are different outcomes with different populations; report them separately.

## Code and reproducibility

- Photo and population/radiance analysis: scripts/analysis/mois_viirs_2020_population_match.py and associated workflow/tests.
- Independent species-colour analysis: scripts/analysis/mois_viirs_independent_white_2020.py and associated workflow/tests.
- Public GIS provenance probes: workflows mois-viirs-worldpop-source-probe.yml and mois-viirs-pop-quality-probe.yml.
- Source files and native scientific artifacts remain immutable. No photo download or ROI inference and no hypothetical MNAR tipping analysis.

Future main work: compare within-year photography composition **without interpreting absence**, improve time-series floristic sampling when independent source information permits, and add independently verified moth-pollination guild and land-cover controls. Do not fit an extinction model to missing photographs.
