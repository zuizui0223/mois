# MOIS authentic 334,015 photographed RGB floral-ROI cohort: 2026-10-19

## Status and provenance

**Executed, not proposed.** Full main-branch source-authenticated GitHub Action:
https://github.com/zuizui0223/mois/actions/runs/39208922929
(success 2026-10-19 UTC).
Primary released outputs:
https://github.com/zuizui0223/mois/releases/tag/mois-20261019-real-photo-source-arm-season-pre-tail-v1

Release **`cohort_white_balance.csv`** contains the values below; full
**`actual_photo_colour_and_season_support.json`** includes source photo-ID
SHA256, arm-specific hue counts and observed-local-season overlap.

Inputs (genuinely processed image pixels, not candidate photo metadata):
- Original primary 317,631 ROI photo attempts.
- Cumulative taxonomic-breadth release with original 313,535 shared ROI photos
  and 4,096 different/new actual ROI photos. The 313,535 shared photo IDs and
  original colour/taxon labels are explicitly checked and counted ONCE.
- Subsequent independently measured 4,096-photo wave2.
- Independently measured cached 8,192-photo wave.
- Additional tail 16,384 still incomplete at this checkpoint, intentionally
  excluded. Upstream 93,209,351 **candidate PHOTO IDs** have NOT yet been
  processed as image pixels; its separate full-index workflow is running.

All five input ZIP/Parquet provenance receipt SHA comparisons in the release
workflow succeeded. Strictly classified hues are only `white`,
`yellow_orange`, `red_pink`, and `blue_purple`; all alternative original
technical labels are retained as unknown. Source photo pixel outcomes are not
reclassified to reach a desired balance. **Anthesis/petal correctness remains
unvalidated by independent botanical human review.**

## Directly observed original-photo results

| Sampling arm (recruitment is NOT equal) | Actual original ROI attempts | Classified four-hue ROIs | White among classifiable | Nominal taxa classified | Nominal genera classified |
|---|---:|---:|---:|---:|---:|
| Primary core | 317,631 | 152,827 | 47.2266% | 25,651 | 5,435 |
| Separate first broad wave | 4,096 | 1,301 | 54.6503% | 1,301 | 1,264 |
| Separate second broad wave | 4,096 | 1,436 | 51.0446% | 1,436 | 1,436 |
| Cached 8,192 | 8,192 | 2,780 | 51.2950% | 2,780 | 2,308 |
| **Union: unique original photo IDs** | **334,015** | **158,344** | **47.3936%** | **31,161** | **7,518** |

At this frozen source checkpoint, classified real-image photo ROI labels:
**75,045 white** and **83,299 chromatic** (four-hue classifiable only).
**175,671 attempted pixels** are unknown/unclassified or acquisition failed,
and must not be recoded chromatic or excluded silently from quality reporting.
Observed-date coverage among classified photo ROIs is 158,335 of 158,344
(99.994%). The new 16,384 attempted photos beyond the primary core contribute
5,517 classifiable colour outcomes (33.67%). This is a different yield from the
legacy core (48.11%).

The nominal taxon-richness increase over the primary core is **+5,510 taxa**
(31,161 vs 25,651) and **+2,083 genera** (7,518 vs 5,435). These numerical
increases are for **observed named original-photo taxa/genera with classifiable
ROI labels**, NOT verified new biological species. The new waves were
deliberately prioritized for previously underrepresented taxa, so they
cannot be treated as a random independent flora sample.

White fraction sensitivity:
- Pooled source-image-weighted classifiable white: **47.3936%**.
- Taxon-equal original photographed colour: **47.8235%**.
- Equal nominal genus means: **51.2857%**.
- Primary core alone: 47.2266% photo and 50.0444% genus equal.

The near-half magnitude is intriguing and robust as a description of THESE
automatically labelled photo collections, but the difference between photo and
genus weights is real and a true worldwide flower-colour equilibrium is NOT
identified. The range of feasible white population labels if the entire unknown
photo-ROI group is adversarially white/not-white is exceptionally broad
(22.4676% to 75.0613% of all attempted images), even before geographic and
taxon observation bias.

## Season observed-day support

Computed on the REAL original observed-on calendar days, not annual climate
normals or estimated flowering times. Summer/winter are opposite by
hemisphere; December enters the subsequent local winter season-year.
Observations (not every photo) define locality×taxon×season-year replication.

North extratropical local comparison, within nominal taxon × 0.25-degree
latitude/longitude cell × local season-year:
- At least one independent observation in summer AND winter: **41 matched
  pairs, 96 observations**. Mean pair contrast summer-minus-winter white:
  **+2.439 percentage points**.
- At least two independent observations in EACH season: **2 matched pairs,
  8 observations**. Mean difference: **−25.0 percentage points**.

These are photo-label contrasts in observationally matched source groups, not
individual colour change, thermal plasticity, or a verified ecology effect;
the stricter comparison is *severely underpowered*. The exact north/south
composition and accompanying raw summer/winter rates are maintained in the
full JSON release and should be checked before interpreting global seasons.
A transparent machine-log summary is implemented separately; no reported
season coefficient is elevated to a causal result.

## Decision and next live actions

1. Use the authentic 334,015-photo checkpoint for current descriptive
   four-hue balance and taxon-coverage assessment; retain source arms in
   all tables and models, plus observed-date diagnostics.
2. The running 16,384 tail release, when completed, may be appended only by
   exact original photo ID deduplication to produce a new separately-versioned
   final result; never overwrite this pre-tail source checkpoint.
3. The full 93,209,351 source-photo-ID/observation-day index must finish
   verified metadata export before probability-sampled new RGB/photo
   collection begins. Source IDs alone are not colour observations.
4. New environmental inference must adjust for geography, nominal taxon
   selection, source arm, technical missingness, and realistic image
   flowering/petal QC. Do not claim measured fitness, pigment biosynthetic
   maintenance cost, genetics, or the seasonal weather effect from source
   photography alone.
