# More dated flower-colour photographs, but no more replicated summer–winter localities

**MOIS source-verified result, 17 October 2026** — original measured 267,425 photographs with immutable four-hue classifier outputs. This independently reuses the successful [200 + 9,400 recovered historical original photo records](https://github.com/zuizui0223/mois/actions/runs/39008584811) and does **not** perform another iNaturalist API request or RGB photo classification. The original 2007 source-dated taxon count is analytically separate from the much broader historical 42111 named-taxon opportunity source, which previously lacked `observed_on` in this combined export.

## Major genuinely measured coverage gain

| From signed source photo IDs, no imputed dates | Before recover | After 9,597 photo/date matches | Change |
|---|---:|---:|---:|
| Colour-classifiable + actually dated original photographs | **81,867** | **86,172** | **+4,305** |
| Distinct original named photo-hue taxa with known observation date | **2,007** | **5,315** | **+3,308** |
| All actual photographic ROI attempts with any recorded observation date | 161,761 | 171,358 | +9,597 |

Among the **4,305 newly dated historic classifiable original photos**, there are **3,801 distinct original named taxa**. Exactly **3,435 of those 3,801 taxa (90.37%)** have just **one** newly dated classifiable historic photograph, while 493 of the 3,801 taxa already had dated hues in the earlier high-depth source. Thus historical date recovery increases global taxonomic coverage dramatically but leaves most taxa without repeated within-taxon seasonal colour data.

The signed 9,600 original photo-ID cohort and unchanged original hue numerator are fully source-authenticated; [see the numerical provenance result](../results/mois_original9600_date_taxonomic_gain_vs_spatial_temporal_support_20261017.json) and the original [provider retrieval ledger](MOIS_ACTUAL_HISTORIC9600_DATE_WHITE_LABEL_RESULT_20261017.md).

## Critical falsification gate: same named taxon, same approximate locality, same season-year

For *classified original-photo* summer versus winter comparisons, use only source photos at the **same original nominal taxon × 0.25-degree geocell × original meteorological season-year**. December winter is assigned to the following season-year; southern hemisphere season labels are reversed. Counting **distinct original observation IDs per season** is stricter than counting multiple photos of one observation, but still does not prove distinct photographers or individuals.

| Summer/winter observational support (same 0.25° cell) | North before → after | South before → after |
|---|---:|---:|
| Pooled across different years, ≥1 original observation in each season | 131 → 133 | 70 → 70 |
| Pooled across different years, ≥2 original observations each season | 9 → 9 | 7 → 7 |
| **Same local season-year**, ≥1 original observation each season | **27 → 27** | **17 → 17** |
| **Same local season-year**, ≥2 independent source observations per season | **1 → 1** | **0 → 0** |

**No new taxon × coarse-locality × season-year summer/winter comparative group** was obtained from the 4,305 additional dated original source hue photos, even with a one-photo-per-season rule. The newly recovered historical photo subset alone contains **zero** same-taxon × 0.25° cell summer–winter comparisons even when years are pooled. The sole Northern site/year group satisfying two independent observations in both seasons contains four photos; none occurs in the South.

This demonstrates the key distinction between **more dated taxonomic and geographic breadth** and **genuine repeated biological seasonal contrasts at one locality**. The 0.25° cell remains coarse (tens of km) and cannot be mistaken for a known breeding population. As shown in the original MOIS seasonal turnover analysis, enlarging the cell to 5° artificially improves photographic overlap but damages biological locality identification.

## Revised publication interpretation

The descriptive Northern/Southern summer white-labelled photograph excess is useful for **seasonal turnover in the taxonomic composition of photographed flowers or flower-like organs**. The original historic-photo subset supports similar *source-level* summer ratios. Yet the stronger biological statement — **local flowers become whiter when a site experiences hotter, brighter or drier months** — still has essentially no replicable site/species/year counterfactual support.

The missingness is not just a lack of enough source photos. The newly dated original historical taxa are predominantly represented by only one classifiable photo; many exact original photos remain automatic hue-unknown; `Flowers` annotations are very sparsely available and are not photo-specific botanical petal/tepal identification. To test local causal mechanisms, source photos need replicated independent **botanically confirmed open floral organs**, original observed-year/month real weather and site/taxon season overlap. Static SoilGrids, elevation, topography and coordinates remain background geographic adjustments, not time-varying seasonal exposures.

**Separate inference boundary:** FCP's prospective within-species `W/D/50 km` polymorphism estimates and MOIS's original genus-equal annual precipitation/water vapour maxT results remain unchanged. No new summer-photo ratio is evidence for pigment fitness benefits/costs or a stable evolutionary 50:50 colour balance.
