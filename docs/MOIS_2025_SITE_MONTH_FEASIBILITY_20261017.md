# MOIS 2025 photo-month climate: source-verified spatial and floral eligibility (17 October 2026)

**Research owner:** `zuizui0223/mois`. This diagnostic is deliberately **not** FCP's species-conditioned spatial polymorphism inference or MOIS's historical congeneric/global rainfall slope. The data are one `Asclepias syriaca` taxon in one sampling year, and are inadequate for a flower-colour seasonality or climate-adaptation conclusion.

## What was actually executed

A deterministic re-join of source-authenticated, original observed-on dates and FCP original ROI-v4 records to TerraClimate v1.1 **real 2025 monthly climate** at the *individual photo's approximately 4-km climate grid cell*, not the earlier common site-centre climate proxy. The underlying photo source and physical climate come from:

- Genuine original 2025 FCP photos: [run 38909164068](https://github.com/zuizui0223/fcp/actions/runs/38909164068).
- MOIS source photo/date/site→month cells, [run 38914197937](https://github.com/zuizui0223/mois/actions/runs/38914197937). Original 64 source records and 40 month×cell combinations.
- MOIS original TerraClimate 2025 June–November, 126 provider variable/site requests with their complete 1991–2025 monthly history, CF packing restored into physical units, [run 38914456992](https://github.com/zuizui0223/mois/actions/runs/38914456992).
- Immutable source SHA256 values and the biological inference gate are preserved in [the concise numeric result](../results/mois_2025_exact_cell_source_feasibility_20261017.json). The [replay script](../scripts/analysis/mois_2025_exact_cell_climate_flower_stage_gate.py) uses both authentic source artifacts and refuses a mismatched checksum/date/month/climate ledger.

Every one of the **64** source photos matched its original observed-month and the correct physical climate cell. The cohort has **64 independent original observations and 64 distinct observers**, occupies **21** fine-resolution climate grid cells and **40** distinct 2025 climate month×cell combinations. This is not a dataset of 64 botanically verified flowering photos.

## Seasonal physical environment in the actual sampled photo months

Means are **original-photo/observation-weighted**, not unweighted across months or averaged at the coarse 0.25° site-centre. Precipitation and water deficit are monthly totals; solar radiation is monthly mean W/m².

| Original-month TerraClimate | 2025 summer (32 sources) | 2025 autumn (32 sources) |
|---|---:|---:|
| Mean daily air-temperature midpoint (°C) | 21.9719 | 13.4391 |
| Mean monthly precipitation (mm) | 63.6125 | 84.6063 |
| Mean solar radiation (W/m²) | 235.4156 | 132.8344 |
| Mean monthly climatic water deficit (mm) | 53.6063 | 9.6063 |
| Mean vapour-pressure deficit (kPa) | 0.9934 | 0.5559 |
| Temperature anomaly vs same-cell/same-month 1991–2020 (°C) | +0.2416 | +0.0730 |
| Precipitation anomaly vs same-cell/same-month 1991–2020 (mm) | −28.2283 | −17.8410 |

These environmental differences occur **even when nothing is known about open flowers**. They establish an authentic temporal exposure gradient among the sampled photographic opportunities, *not* a climate-linked change in flower hue.

## Design diagnosis: site, observer, flower-organ confounding

- Source observations per month: June 11, July 11, August 10, September 13, October 12, November 7.
- Of **21** approximately 4-km climate cells, **only 7** have both summer and autumn source observations (15 summer +24 autumn =39 photos). Only these seven cells can form source-level *within-cell* season comparisons. They are not verified breeding populations.
- Automated ROI colour classifications: summer **12/32**, autumn **6/32**. Excluding three same-ROI white/yellow organ mixtures leaves at most **9** summer and **6** autumn unreviewed clean source-photo colour labels.
- Independent botanical validation of an **open flower and the actual petal/tepal colour ROI**: **0 confirmed in either season**. All 64 are still marked unreviewed in the two-botanist blinded handoff.
- The previously frozen feasibility threshold for the actual flower-colour numerator/denominator is **20 botanically verified classifiable photos per season**. The current 64 photos cannot reach that threshold even if every remaining clean automated label eventually proves to be a real open flower.
- **There is no repeated observer in this 64-photo first wave.** An observer-paired seasonal hue comparison is impossible in this panel. A larger panel and independent botanical QC are required.
- Both recorded season and weather are functions of calendar time in a single year at one taxon/region. No joint seasonal weather coefficient can be identified independently of the calendar contrast using this work unit.

The immediate scientific conclusion is a **data/identification limitation**, not a positive or negative test of the seasonal biology of flowers. Unknown flowering is *not* nonflowering, and it is **not a measured white or coloured flower**. The 64-source cohort must not be pooled with the already measured 100,543-image original MOIS historical atlas or the independent 149,900-image FCP study.

## What qualifies the next ecological test

The current first wave is a functioning *source and climate* pilot but insufficient biological colour data. The next evidence-building step, owned by MOIS, is to audit the independent botanical flower-stage review for all 64 photos (two blinded reviewers with disagreement adjudication) and expand source-only sampling beyond this single site-year/taxon to the 39-taxon/268-site-season-year opportunity priority. Original source IDs, dates, observer identifiers and photo denominator must remain frozen before hue classification. In any expanded study, require verified open flowers in **both** seasons, independent observations and cross-season geographic support before attempting species/site-season white-rate or multi-taxon, within-site climate anomaly models.

The existing historical MOIS global congeneric moisture manuscript remains a **separate** scientific estimand, with its strict-photo, missing-colour and geographic nonconfirmation results intact. FCP's recurrent colour-space `W` and 50-km conspecific local-depletion tests are not used as MOIS outcomes.
