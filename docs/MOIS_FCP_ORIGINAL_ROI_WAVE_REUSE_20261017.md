# MOIS measured flower-colour data, extended from 247,032 to 248,312 original photographs

The current MOIS scientific estimand is **worldwide *between-named-species* visible photographed flower hue against environment, elevation, geolocation, stable soil/terrain and real observed-time weather**. Its photo sources may overlap with FCP. **FCP within-species D/W and photograph-pair 50km local-depletion tests are not imported as MOIS outcomes.**

Previously, the exact original historical MOIS (100,543) plus high-depth FCP (149,900) merge, deduplicated at original `photo_id`, produced **247,032 photographic opportunities and 120,310 actually four-colour classifiable labels** (merged PR #15). No human botanical review gate is applied to the historical original automated photographic-hue classification; biological stage uncertainty remains an explicit limitation.

Additional source recovery found **six truly processed FCP original-photo ROI-v4/EfficientSAM technical wave results** with their own immutable GitHub Action artifacts. This is **not** 3.43 million classified pixels or a new biological sample from a different provider; each candidate image must already have a genuine frozen image-model result. Deduplicating against the 247,032 photographed source union gives:

| Original measured ROI group | Already measured original rows | New photo IDs | New four-hue labels |
|---|---:|---:|---:|
| FCP RGFCA six taxa wave 1, live API | 552 | 548 | 199 |
| FCP RGFCA six disjoint taxa wave 2, live API | 482 | 477 | 200 |
| FCP actual 2026-09-27 monthly first 64 | 64 | 64 | 27 |
| FCP actual 2026-09-27 monthly second 64 | 64 | 64 | 27 |
| FCP original exploratory site×month 63 | 63 | 63 | 26 |
| FCP original Asclepias 2025 first 64 | 64 | 64 | 18 |
| **Total extra** | **1,289** | **1,280** | **497** |

**Expanded original-photo data: 248,312 unique photos, 120,807 classifiable hues, and 127,505 unclassifiable/unknown.** The other nine sources were already used, so their second measurements are not counted as new photos. The six outcomes were matched by original ID and compared against the prior source; there were no contradictory double-classifiable four-hue labels in the preflight.

The FCP *live API* photo sources are signed to the actual observation/download workflow but are **not** automatically claimed to be IDs from the exact 2026-09-27 archived bulk snapshot; the strict monthly originals are. All sources reuse the original ROI-v4 measurement family and remain photographic phenotype measurements, not verified flower-genetic polymorphism prevalence.

### Reproduction

`scripts/analysis/mois_append_fcp_six_real_roi_waves_20261017.py` and `.github/workflows/mois-append-original-fcp-six-real-roi-waves.yml` recover **seven** original source artifacts (one already-merged 247k dataset plus six signed real FCP ROI terminal groups), check each measured terminal CSV against its original SHA256 receipt and fixed photo counts, preserve unclassifiable photos, reject conflicting classified colours, and upload the actual compressed 248,312-photo result. It does not open any new original image pixels, change FCP output, or rerun historical MOIS regressions. This is an efficient reuse of already-processed photo outcomes.

### Scientific analysis ownership

The expanded MOIS dataset is the starting point for source- and taxon-standardized global hue-versus-environment analysis. Site latitude/longitude/elevation/soil and topography must be kept conceptually static at a site. Long-term BIO baselines are neither actual observed-year monthly weather nor photo-hour of observation. Where observed-on dates exist, join original-year/month temperature, radiation, precipitation, moisture deficit and VPD separately (including year/season/geographic support); recover genuine photo time-of-day only when an authoritative timestamp/timezone exists, otherwise leave missing. More than 3.4 million *source photo IDs* are available upstream, but the unmeasured ones must still be processed through a calibrated image model before they count as new four-hue outcomes. Human review can be a technical sensitivity, not a mandatory blocking step for the entire photographic macro-analysis.
