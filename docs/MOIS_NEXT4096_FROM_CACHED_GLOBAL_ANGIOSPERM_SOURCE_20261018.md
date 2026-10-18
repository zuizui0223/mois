# MOIS second global angiosperm acquisition wave: 4096 source IDs without rescanning 500 million photos

The first official iNaturalist global Angiospermae species-first candidate
pass scanned approximately 500 million photo metadata rows once and saved
a provenance-locked **32,768 original photo ID** source pool. First 4096
source IDs cover 4096 named species and 3788 genera. That source
was archived as an immutable GitHub Release:
https://github.com/zuizui0223/mois/releases/tag/mois-angiosperm-broad-species-first-source-cache-v1-20261018

The first genuinely executed 256 original-image ROI pilot classified
**69 of 256 attempts into four visible hues**. The resulting 69
colour-classifiable original photos each belong to a previously
unrepresented source species and a previously unrepresented genus;
these are not 256 new measured colours. Fivefold and annual-climate
macroecological model estimates remain frozen to their published
sample supports.

## Source-only second wave

The implementation at
scripts/acquisition/mois_cached_angiosperm_second4096_source.py
reads only three already archived files: the source pool of 32,768,
the exact first 4096 original photo-ID roster, and its SHA-signed
monthly metadata receipt. It does NOT scan any remote iNaturalist
observation/taxonomy/photo metadata, download any original image
pixels, run any flower model, or estimate colour from a photo URL.

It SHA-verifies both source files; all photo IDs, measurement IDs,
previous source metadata and the first-wave rows must be immutable.
It excludes all 4096 first-wave PHOTO IDs *and named species* before
choosing the next 4096 original photo IDs. The original parent
source had already excluded all 313535 measured core photo IDs,
the ongoing earlier high-depth 4096 reservation, and the independent
supplementary seasonal 4096 reservation; the second wave inherits
that unchanged parent snapshot. It prefers genera absent from the
first wave, then new species in already photographed genera,
and source-photo regions previously scarce among first-wave
original photo IDs. It also prefers independent photographers,
never allows duplicate photo/observation IDs and uses a
photo-colour-blind deterministic hash for ordering.

Outputs are sixteen new 256-photo source-only worker rosters, a
4096-photo source-only roster and a SHA-verified receipt with the
real numbers of new source species, genera, observers and source
10-degree geographic cells. These are **planned photos**, not
4096 classifiable petals or 4096 successfully downloaded RGB images.

These second-wave rosters become eligible for the existing frozen
YOLO11/EfficientSAM ROI-v4 only after the first broad 4096-photo
measurement has genuinely completed and is archived. Running
the next wave after a completed first wave prevents simultaneous
uncapped requests to the original photo provider. Worker receipts
and site-climate values will be appended by literal original photo
ID with all earlier source records preserved and image/flower
classification unknowns explicitly retained.

Why this matters for MOIS rather than FCP: main ecological inference
requires more represented flowering-plant species, genera and
environmental contexts, not hundreds of photographed replicates
per single high-depth species. This still is a deliberately selected
iNaturalist photographic atlas and cannot be called an unbiased
population census; visible white is not proven metabolically
pigment-free and source frequencies do not prove a maintained
evolutionary 50:50 white-colour equilibrium or fitness costs.
Annual climate remains a historical WorldClim normal and
month-level climatic contrasts remain SI.
