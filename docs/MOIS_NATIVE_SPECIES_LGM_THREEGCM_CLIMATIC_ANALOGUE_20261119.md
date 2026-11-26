# MOIS: native species-level LGM climatic analogue screening (2026-11-19)

## Current scientific scope

The preceding [native-only historical rainfall study](MOIS_WCVP_NATIVE_LGM_RAIN_NPP_GEO_SPATIAL_HOLDOUT_20261119.md) found that global photo9–LGM rainfall associations do **not** persist consistently under WCVP native restriction and 500/1000 km regional fixed effects. Therefore, **this new stage does not assume an historical refuge or an evolutionary floral-colour adaptation**. It tests whether temperature and precipitation at candidate cells within present-day recorded native regions have modern climatic analogues for recorded plant species under three independent LGM GCM outputs.

The term **climatic analogue** is essential. It is *not* an observed LGM species occurrence, palaeorefugium, historical accessibility, population persistence, or species fitness/survival.

## Real source lineage

- Original accepted photo9 species: **37,104**; historical GBIF snapshot: **5,493,931** accepted species×100 km grid presences; **396,142,825** filtered occurrences across **18,577** grid cells.
- Strict original WCVP TDWG3 native source: **2,271,852** species×100 km cells; original introduced source **436,422**. Kew accepted unique-exact name correspondence: **36,762** source species. Same source receipt and TDWG3 25-point classification revalidated: **9,288 strict** and **10,391 moderate** TDWG3-assigned cells.
- Modern WorldClim1.4 climate: 1960–1990 at **12,295** original signed environmental complete-case cells. This matches the original three LGM GCM baseline. **Do not mix it with WorldClim2.1 (1970–2000) when interpreting palaeoclimate deltas.**
- LGM: WorldClim1.4 downscaled ~22 ka CCSM4, MIROC-ESM, MPI-ESM-P and original full source SHA. The LGM grid stores **present minus LGM** change in temperature and log1p precipitation, not absolute past climate and not species occurrence. Reconstruct per GCM at the original grid without changing current baseline.
- Shared modern + all-GCM + 25-point strict native-area source domain: **7,393** 100 km cells, restricted to current recorded TDWG3 native botanical regions and valid grid climate. This is neither a complete palaeogeographic ice-age land mask nor a historically accessible area reconstruction.

Source hash constraints are embedded in `scripts/analysis/mois_native_species_LGM_3gcm_climatic_analogue_pilot_20261119.py`. Original sources are read only. New workflow `.github/workflows/mois-native-species-paleo-3gcm-pilot-20261119.yml`.

## Explicit climatic analogue estimator

For each source species with at least **50 strict WCVP-native currently occupied 100 km cells**, at least **5 distinct 500 km geographic blocks**, and available 1960–90/3GCM climate, construct a two-variable empirical rectangle:

`T_now ∈ [Q05(native T), Q95(native T)]`

`log1p(P_now) ∈ [Q05(native log1p P), Q95(native log1p P)]`.

For all candidate cells in TDWG3 regions where the species is currently recorded as native, count cells satisfying **both** conditions under current climate and under each LGM GCM. Report CCSM4, MIROC-ESM and MPI-ESM-P separately, any-GCM and **all-three-GCM consensus**, and the number of palaeo-analogous cells divided by the number of modern-analogous cells. Compute an additional 500 km block-heldout native-occurrence climatic-envelope recall and ratio to current climate-background hit rate for every candidate species. This **is not a calibrated SDM**, has no verified absences and does not model dispersal, vegetation or past suitable area outside current native TDWG3 regions.

The source species is counted once per GBIF 100 km cell, not weighted by record count. 100 km equal-area grid cells at coasts remain imperfect occupied-area units. TDWG3 native classification is an *extant regional status*, not historical glacial accessibility or specimen provenance.

## Verified pilot — 9 colour classes, deterministic source-only selection

Reconstruction, environmental projection, same-source quality and 500 km holdout successfully executed on [Actions #44513010102](https://github.com/zuizui0223/mois/actions/runs/44513010102) (**SUCCESS**). The original WCVP region status records and exact original 1960–90/3GCM climate all passed SHA checks.

- **10,882** accepted source species satisfy minimum 50 WCVP-native current cells and 5 spatial blocks.
- Deterministic outcome-independent pilot sampled up to **100 species per nine source photo9 hue classes**: **856** species (orange has only 56 eligible).
- **833/856** have at least one candidate region cell with 3GCM-consensus historical climate analogue. Zero ensemble-consensus cells for **23** species.
- Every pilot species underwent 5 spatial-block test folds. Median heldout occupied-cell recall is about **0.74–0.75** across nine groups. **580/856** have heldout climate-envelope hit rate at least **1.2 ×** the same-region current climate background proportion; climatic discrimination is uneven.
- The resulting LGM ensemble analogue count / current analogue count **differs among source nine-hue groups** but is exploratory only. Differences may be driven by genera, region pool sizes, latitude, climate-envelope width, recording process and photo-image labelling, not historical colour selection.

| Source photo hue (not verified petal colour) | pilot n | median all-3GCM LGM/current analogue-cell ratio | median current region candidate analogue cells |
|---|---:|---:|---:|
| white | 100 | 0.444 | 306 |
| yellow | 100 | 0.455 | 198 |
| orange | 56 | 0.531 | 232.5 |
| red | 100 | 0.557 | 230 |
| pink | 100 | 0.331 | 265.5 |
| magenta | 100 | 0.352 | 210.5 |
| purple | 100 | 0.248 | 200 |
| blue | 100 | 0.289 | 210 |
| bronze | 100 | 0.614 | 304.5 |

**Important quality gate:** Among 856 source species only **29** are labelled `sampling_consensus_single`, with 309 `provisional_single`, 276 `sampled_multicolour` and 242 `multicolour_candidate`. These are source ROI classes, not botanically validated floral petal colours; bronze in particular can include organ or background confounding. **Do not interpret the hue table as an ecological or phylogenetic adaptation test.**

### Next stage: all 10,882 eligible native species

The pilot demonstrates numerical feasibility. The same code has been parameterized so `--max-per-hue 0` runs **all eligible native species**, rather than the nine-class balanced pilot. All future claim testing should avoid treating the artificial pilot's 100/hue composition as worldwide plant species prevalence. The full run will still be *eligibility-truncated* (native occurrence density and geographic range). Expand species screening only after auditing this selection mechanism and the pixel ROI true-flower label gate.

## Refugia evidence cannot be claimed from this model

WorldClim LGM three-GCM correspondence does not show occupancy by any fossil or surviving plant population. A credible actual LGM refugia analysis needs independent **fossil/pollen occurrence + chronology** (for example, [Neotoma data catalogue](https://www.neotomadb.org/data)), independent phylogeography and species-level (not genus-only pollen) taxon harmonization. No Neotoma records have yet been fetched, matched or validated against these 856 species. Source geographic candidate domain remains constrained by current recorded native regions; it can miss genuinely shifted ranges.

The initial schema gate [Actions #44511325109](https://github.com/zuizui0223/mois/actions/runs/44511325109) verified 18,577 full LGM source rows and exact three-GCM delta variable definitions; the WorldClim official historical product is [WorldClim1.4 palaeoclimate](https://www.worldclim.org/data/v1.4/paleo1.4.html). Both the palaeoclimate baseline and selected caveats are source-verified, not speculative fossil evidence.

## Signed outputs

[Original genuine pilot run #44513010102](https://github.com/zuizui0223/mois/actions/runs/44513010102), artifact #14523493715, ZIP SHA256 `18c7a6c941cf0e26c1c5985e52fce1319a96aecef6f1fa3630dae0ef66486dc8`.

The source receipt records input artifact SHA256, WCVP manifest SHA, source status counts, nine-hue cohort sizes, true 3GCM candidate cell count, and individual output file SHA256. Outputs include 856 species×3GCM per-region climate analogue and heldout validation scores, per-hue summaries and a complete signed receipt. The pilot branch is [MOIS draft PR #206](https://github.com/zuizui0223/mois/pull/206), not `main`.

