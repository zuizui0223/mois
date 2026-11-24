# 2026-11-25 — Contemporary rainfall and native photo-white assemblages

## Scope
Exploratory same-data MOIS moisture extension. Compare rainfall slopes separately for native and introduced taxa in the same cells, then remove the introduced-richness criterion. No new photographs; no change to main or acquisition PR #200.

## Inputs and reproducibility
- Paired WCVP native/introduced artifact #14389295373, SHA256 3efaf6ee359aca4a49713d55f5bbadd0c86467ccd4d27c300f1e06936245759a.
- Environmental grid #14397855889, SHA256 4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6.
- WCVP geographical grid #14400177816, SHA256 cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd.
All outer archive and inner receipt SHA checks passed. Strict and moderate paired sample sizes 3106 and 3602. Existing difference-model ten coefficients for both modes reproduced within 1e-8 pp; native and introduced white fractions and richness cross-checked; exact linear coefficient contrast identity verified. Both scripts executed twice locally; CSV outputs byte-identical.
- scripts/analysis/mois_moisture_paired_20261125.py: blob 3014e3e75131601009183a4a76eefe6f5fd92b08; output CSV SHA256 163633e78d44facfff129e875885d09416f38e381ed6630f9012b56cfc7cb0f6.
- scripts/analysis/mois_native_grid_moisture_20261125.py: blob f00f7facadd0346be5567f276128b12482c754fb; output CSV SHA256 39f24d65f6e114b0140f68d903b40d99155a4a4646d9ace2c9d3fe6a3637de78.
Run each script with --inputs pointing at archives named github-actions-artifact-{id}.zip and a new --out directory.

## Results (pp per cohort SD of log1p annual precipitation, nominal cluster-robust 95% CI)
### Identical strict paired 3106 cells
| Spatial FE / weight | Native | Introduced | Difference |
|---|---:|---:|---:|
| 500km equal | +0.521 [0.045,0.996] | +0.404 [-0.398,1.206] | +0.117 [-0.775,1.008] |
| 500km minority richness | +0.533 [0.199,0.868] | -0.000 [-0.538,0.537] | +0.534 [-0.115,1.182] |
| 1000km equal | +0.395 [-0.157,0.948] | +0.200 [-0.394,0.794] | +0.195 [-0.549,0.939] |
| 1000km minority richness | +0.722 [0.355,1.089] | +0.001 [-0.353,0.355] | +0.721 [0.161,1.281] |
Fitted cells 3016 and 3091, clusters 374 and 164. Strict/moderate overlap, not independent.

### Broader native-only cohorts, no minimum introduced richness
| Cohort / weight | 500km FE | 1000km FE |
|---|---:|---:|
| Strict >=30, 7201 cells, equal | +0.003 [-0.549,0.555] | +0.036 [-0.698,0.769] |
| Strict >=30, native-richness weight | +0.491 [0.118,0.864] | +0.608 [0.112,1.104] |
| Strict >=100, 5388 cells, equal | +0.567 [0.133,1.001] | +0.746 [0.200,1.292] |
| Strict >=100, native-richness weight | +0.571 [0.232,0.910] | +0.769 [0.347,1.191] |
| Strict >=30, introduced <30, 4095 cells, equal | -0.675 [-1.566,0.216] | -0.193 [-1.242,0.856] |
Moderate native100 equal: +0.620 [0.249,0.992] and +0.733 [0.305,1.161]. Broader moderate native30 equal: +0.157 [-0.298,0.613] and +0.146 [-0.467,0.760].

## Interpretation and limits
Positive association occurs directionally in rich recorded native floras but is absent in the broader equal-cell sample. Distinct estimands, sampling/coverage and regional composition can explain the contrast. Only one strict matched sensitivity contrast excludes zero; no general native-introduced slope difference. No formal richness × rainfall interaction test here. Different cohorts use different SD scaling. WCVP regional status is not local organism nativeness. The GBIF parquet archive was recovered and outer SHA matched but the local execution environment lacks a parquet reader; independent cell occurrence-count effort was NOT adjusted. Models control richness only; weights are estimand sensitivities, not precision weights. Nominal CIs without experiment-wide multiplicity correction. Species photo-label composition is not measured local petals, pigment physiology, dispersal or causal moisture response.

## Next decisive test
On fixed strict-native cohort test rainfall × richness, control GBIF occurrence count, hold out geographic realms, compare genus composition and within-genus turnover. Do not infer interaction from significance differences. Keep present moisture main line, not new independent NPP/LGM narratives.
