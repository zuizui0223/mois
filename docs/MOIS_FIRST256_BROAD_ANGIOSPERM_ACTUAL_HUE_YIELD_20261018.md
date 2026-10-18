# MOIS first 256 global-angiosperm RGB pilot: how many NEW flower-colour taxa?

## Authentic original-pixel pilot (18 October 2026)

The primary MOIS source contains 313,535 original photo attempts,
151,239 classified automatic visible-hue ROIs and 71,427 white.
A separate, previously unmeasured iNaturalist Angiospermae source
retains 32,768 original candidate photo IDs. Its first 4,096
source observations cover 4,096 named species, 3,788 genera and 222
coarse geographic 10-degree cells; these are metadata-only candidates,
not measured colours or botanically confirmed flowers.

The [actual first 256 photo measurement run](https://github.com/zuizui0223/mois/actions/runs/39097521005)
completed four original 64-photo YOLO11/EfficientSAM ROI-v4 work units.
It classified **69 /256 original source photo attempts (26.95%)** into
the four visible hues, with **187 unknown photo ROI outcomes**.
The frozen worker logs reported **six original image acquisition
failures** in those 187. The original photo source quality gate for
continuing the full 4,096 photos was 20% classification success,
which this pilot passed.

The independent follow-up audit reads:
1. The immutable 313,535-photo photo-ID/hue/taxon/dates Parquet source.
2. The authenticated source-only first 4,096 new species photo-ID plan.
3. The original 256 photo source SHA, actual YOLO/EfficientSAM
   photo worker terminal, four-hue classifications and missing statuses.

It revalidates model, photo ID, date, photographer, license, source taxon
and field-site identity using the same checks as the complete 4,096-photo
integration. It reports new **actually classifiable species and genera
per 1,000 attempted original image downloads**, plus detailed four-hue
counts, technical unknowns versus access failures and photographic
sampling coverage. The exact source-verified counts will be available
in the audited release only once CI has passed; the number of source
photos is NOT the number of colour-known taxa.

Automatic visible-white is not pigment free and four-hue photo ROIs
do not individually verify flowering/petal anthesis or plant fitness.
This is a deliberately taxon-breadth-focused photo sample rather than
randomly sampled global flowering individuals. A descriptive 50:50
ratio neither proves stabilising evolution nor flower pigment cost.

The continued global fourth source wave works independently in four
model-reuse groups. The accepted photo classifier results must be
source-hashed and pooled only once across the 16 completed image
subshards. All original 313,535 photo IDs, taxonomy and dates remain
immutable. Seasonality remains SI; annual climate/soil missingness
is not imputed from unknown photo locations.

Implementation: scripts/analysis/mois_broad_first256_roi_taxonomic_gain_audit.py
and .github/workflows/mois-broad-first256-true-hue-and-novel-taxa-audit.yml.
