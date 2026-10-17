# MOIS primary hypothesis: global near-half white versus chromatic floral photography

The surprising observation is that **white visible flower-photo labels are
almost as numerous as the three chromatic categories combined**. This does
NOT mean evolution maintains an exact 50:50 evolutionary equilibrium.

In the original 100,543-photo cohort, 21,574 of 45,953 classifiable
photos (46.95%) were labelled white, and 53.05% were yellow/orange,
red/pink or blue/purple. 54,590 photographed ROI attempts were unknown.
Negative geographic compensation / stabilizing balance was not confirmed.
Source: results/source_receipts/global_classified_white_fraction_and_compensation.json

# Ecological question

Why is visible white floral colour common across the world's photographed
flora, rather than chromatic petals universally dominating, and under what
long-term environmental circumstances is WHITE proportion higher or lower?

The first environmental hypothesis is water (annual precipitation, dry-season
rain, atmospheric vapour), because the *already fitted* 309,439-source
genus-equal WorldClim BIO12 white coefficient is +1.99 percentage points
per environmental SD (within-family maxT95 [+0.55,+3.42]). The positive
water association is not confirmed after all geographic conditioning.
Other candidate gradients are solar irradiation, heat and temperature,
seasonality, broad biogeographic conditions, elevation, terrain and soil.

Pigment biosynthesis and maintenance costs, UV optics, pollinator and
thermal/oxidative benefits are **possible mechanisms** to examine later.
WHITE is not equivalent to pigment-free: UV-absorbing flavonoids and
other pigments may remain. Chromatic does not mean anthocyanin alone.
No reproductive fitness or metabolic cost was directly measured.

# Implemented photo -> colour -> environment -> cache architecture

The script scripts/analysis/mois_primary_global_white_balance_feature_lake.py
reuses the original successful YOLO/EfficientSAM pixel-ROI measurements and
checksums the exact source files before joining by original photo_id.
No new hue labels are invented, no unmeasured source ID is classified,
and unclassified ROIs never become nonwhite.

- Main immutable original source: 309,439 actual processed photo ID attempts,
  149,662 four-hue classifications, 149,653 known original days with colour,
  from GH Actions run 39028400690.
- Annual long-term WorldClim v2.1 environment per same original photo ID,
  from completed original 309,439-photo model run 39029308911.
- Prior partial 248,312-photo genuine site terrain layer and historical
  original soil variables from GH Actions run 38992341847. Missing soil or
  terrain for the additional images remains missing; no false coverage claim.

Produces compression-efficient immutable Parquet tables, separately for
(1) original photo IDs, observed days, source, hue/classifiability and
(2) annual climate and (3) partial soil/DEM terrain, plus a SHA manifest
and ecological figures for photo-weighted, named-taxon-equal and
genus-equal WHITE fractions, latitude/hemisphere breakdown and
environmental deciles among all original photo attempts. The old source
annual solar field mislabeled Wm2 is converted from kJ/m2/day by /86.4.

This workflow stores feature data without redistributing raw licensed photos.
Versioned GitHub Release assets after a trusted main-branch run are more
durable than temporary Actions artifacts alone; an off-GitHub mirror is
desirable for stronger long-term guarantees. A full fourth4096 ROI wave
will be integrated only after all image pixel workers finish and after
exact photo-ID deduplication. The separate season-gap-targeted 256
photos are NOT an unbiased addition to the global white frequency.

Analysis hierarchy: (1) near-half descriptive photo vs taxon vs genus,
(2) environmental associations under geography/taxon/source controls,
(3) four original visible hue categories, (4) SI monthly/seasonal climatic
normals and, separately, verified observed-year monthly conditions.
Neither seasonal panel requires within-individual seasonal petal changes.
