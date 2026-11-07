# MOIS 295,017 actually photographed originals: refit the original moisture associations

**Research owner** `zuizui0223/mois` — *among-taxon photographic white vs annual temperature/moisture/solar climatology*. The `fcp` within-species W/D and 50-km phenotype-space study remains entirely independent. **2026-10-17** original source evidence and execution protocol; only actual successful water-model runs may be promoted to a numeric new result.

## Frozen source and comparability

The completed [source photo/dates result](https://github.com/zuizui0223/mois/actions/runs/39015918687) contains **295,017 genuinely measured original photo IDs / 142,538 actually classified four-hue ROI labels**, an increase of 8,191 original photographs and 3,133 classifiable ROIs relative to the independently finished [286,826 FCP180-pixel/photo water source](https://github.com/zuizui0223/mois/actions/runs/39014753277). These extra 8,191 original photos come from two source-blind MOIS real 4,096-image pixel waves, after exact original photo-ID overlap removal (one shared FCP original photo); no colour/imputed new data are introduced.

The entire previously fitted 286,826 photo-site actual WorldClim source including BIO5, BIO12, BIO14, BIO15 and monthly averaged solar is reused with SHA256 `338fbb0d50837c9cc25aece2dfc15d0bbb02e8b192728b8adf030c392c0f41f0`. Every original photographed `photo_id`, taxon ID, species, and hue is compared literally against the new 295,017 source before fitting anything. The original new union photo source SHA256 is `44e3cb6b4c6174dfb1035de208343c7c6a2e6c321e2557a509224c4d9f2ac7b3`. The extra 8,191 genuinely classified/unclassified photos receive independently sampled WorldClim monthly BIO and annual solar normal values at their actual source coordinates; mean monthly vapour-pressure climatology is sampled for all photo sites from the original official WorldClim raster source.

**The 295,017-photo genus-water model has not been reported here before a successful source-verified complete run.** The already measured 286,826-photo result alone showed +1.909 percentage points/source SD global photographed white versus annual precipitation, maxT 95% [+0.466,+3.352], and +2.992 [+0.825,+5.158] for atmospheric vapour, with both coarse-cell geographically conditioned confidence intervals spanning zero.

## Exactly the previous six-anchor × four-hue family, no proposed-outcome selection

The executable `scripts/analysis/mois_refit_original_water_295017_actual_3wave_FCP180.py` verifies strict photo-row consistency against both original photo sources, writes an authenticated compatibility receipt, then **reuses** the original source-reviewed procedure from `scripts/analysis/mois_refit_original_water_with_actual_new_original_photos.py` without altering its:

- six preexisting original annual WorldClim axes: BIO5, BIO12, BIO14, BIO15, annual solar radiation normal and mean annual atmospheric vapour pressure;
- four preexisting automatic photographic hues: white, yellow/orange, red/pink, blue/purple;
- one original source taxon/species vote in genus-centred and geographic photo-hue contrasts, including coarse-cell comparison;
- 999 resampled whole-genus bootstrap draws with existing fixed seeds and 24-family maxT simultaneous confidence intervals.

Only the input source size constants change from 286,826/139,405 to **295,017/142,538** using 8,191 authentic new original-photo pixels. The source photo-derived geography, observed date and original hue label are preserved. If any original photo ID, taxon assignment or previously classified hue differs, analysis aborts rather than silently taking the newer ROI label. If WorldClim input ZIP checksums change, the environment join aborts.

An original WorldClim v2.1 solar raster has physical units kJ/m²/day, while the earlier frozen 286,826-photo raster-join field is historically misleadingly labelled W/m². The unchanged annual source standardized regression is scale-invariant, so the new 8,191 photo sites retain the same raw kJ/day scale for estimator comparability. Do not interpret raw solar values as watts without dividing by 86.4, as already corrected for subsequent genuine physical seasonal-source analyses.

## Inference boundaries

- Original photographic four-hue outcomes remain *automatic RGB ROI image labels*, not verified spectral reflectance or independent botanical diagnoses of open petal/tepal organs. An observation date does not confirm anthesis.
- These original photo sources are partly shared across FCP and MOIS and should not be described as independent biological replications. FCP's own within-species phenotype-space and geographic polymorphism tests are not imported.
- This is an annual **long-term climate-normal association**, not actual observation-year/month rainfall, sun or atmospheric deficit experienced by an open flower.
- Original FCP wave4 image workers remain incomplete at the FCP180 source snapshot, and their compute completion order may distort geographical/taxonomic sampling. The original historical 100,543-photo benchmark and the interim 286,826/295,017 data remain separate, archived analysis populations.

**Operational source:** `.github/workflows/mois-original295017-real-pixel-genus999-maxT-water.yml` downloads existing immutable photo-source artifacts, authenticated WorldClim inputs and refits the exact historical 999-genus/maxT method. It uploads the full new 295,017-site environmental table, original water estimates and a signed, concise water-research receipt only on real execution. No numeric 295,017 global water effect is claimed from the old 286,826 model in advance.
