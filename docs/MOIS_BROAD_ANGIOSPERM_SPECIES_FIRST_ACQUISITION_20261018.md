# MOIS: breadth-first worldwide angiosperm flower-photo source, not FCP within-species depth

MOIS is designed to describe worldwide white versus yellow/orange, red/pink
and blue/purple visible floral-photo composition across species, genera,
regions and environments. FCP's 49,800 original RGB photo cohort represents
498 high-depth species for within-species spatial polymorphism and is not
the next appropriate *primary* photo acquisition objective for MOIS.
The unfinished FCP-to-MOIS high-depth import PR #72 was therefore closed
without merge after the research-design clarification.

The original verified MOIS sample contains **313,535 real original
image-pixel attempts and 151,239 classified four-hue ROIs**. Its 47.23%
white-photograph source fraction is an observational, not evolutionary
50:50 balancing frequency. Unclassified RGB ROIs stay unknown.

## New source that actually broadens the global flora

Use iNaturalist's public monthly Open Dataset of three different
source tables: taxa (taxon_id, ancestry, rank, active, name),
observations (UUID, taxon, research-grade status, observer, calendar
date and public coordinate), and licensed photos (photo ID, observation
UUID, image extension and original file position). Official metadata:
https://github.com/inaturalist/inaturalist-open-data/blob/main/Metadata/README.md

Only active named SPECIES descending from taxon Angiospermae (47125)
are eligible for this acquisition arm, rather than repeating the same
FCP 1000 high-depth taxa. Verified research-grade records must have
a valid observer ID, original date and geographic coordinates;
positional accuracy recorded above 20 km is rejected. Coordinates
might be publicly obscured: they are not proven precise field points.
Use only photo entries with existing open CC/CC0 licenses; record license
and original observation ID; source photos remain copyrighted under
their photographer's specified terms.

The script
scripts/acquisition/mois_global_angiosperm_breadth_source.py
performs one source-locked HTTPS HEAD/ETag and gzip-range source check,
then streams three full official compressed metadata files exactly once
with byte-count/hash and gzip CRC verification. Although the monthly
metadata snapshot is large, reading it once amortizes access across
future acquisition waves. It does NOT decode a single new flower photo,
read an unseen hue or fabricate any genetic/anthocyanin measurement.

It first retains up to four outcome-blind hash-selected research-grade
observation alternatives per eligible species, chooses a valid original
photo by observation photo position and stable hash, and constructs
a reusable pool targeting up to **32,768 novel original photo IDs**.
Photo IDs already in the 313,535 genuinely processed source, the
ongoing fifth 4,096-photo FCP-like reservation, or the independent
seasonal-SI 4,096-photo reservation are excluded.

Broad acquisition priority:
1. Unrepresented ANGIOSPERM genus with a previously unmeasured species
2. New species in an already photographed genus
3. Previously photographed species from underrepresented geography
4. Independent photographer/observation repeats only after breadth

Species are intentionally cycled by genus before subsequent photo
records in the same genus. The first **4,096 PHOTO IDS** are frozen in
16 x 256 source-only shards and the full pool is retained for later
waves. All source priority stages and their real number of different
species, genera, geo-cells and original photographers are recorded
in a SHA receipt once the source has been scanned; do not claim they
are actual flower colour or 4,096 successful photo classifications.

## Actual RGB ROI measurement remains separate

Once the source-only first 4,096 and cached pool are ready, run
location-blind frozen YOLO11/EfficientSAM ROI-v4 on original 1024px
open-license photos, with genuine 64-photo worker terminals and unknown
classification outcomes. Use parallelism that respects the photo
provider's rate limits rather than blindly multiplying jobs.
Preserve source snapshot and original model/worker hashes, photo-ID
deduplication, species/genus and source arm.

The empirical **efficiency metrics** should be:
- Newly classified and previously UNSEEN named species per 1,000
  genuinely attempted original photo images
- Newly classified genera per 1,000 real photo attempts
- Newly covered 10-degree coarse geographic source cells/biomes
- Known four-hue ROI rate (not counting unknowns as coloured)
- Photo/source/taxon/genus-weighted balance stability and uncertainty

Source selection by genus/taxon is not probability-proportional sampling
of global flowering individual abundance. Inference must retain
inclusion/observer/classifier filters and use taxon/genus weights.
Not all iNaturalist angiosperm photos contain an open flower; explicit
botanical anthesis assessment still requires separate validation.
Seasonal environmental response stays supplementary, and pigment
fitness costs, ultraviolet signals and stabilising selection are
not identified by these photo source metadata.
