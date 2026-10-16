# Seasonal floral colours within the worldwide Mois source: species turnover versus repeated-species colour

**Prospective analysis specification, 16 October 2026.** Part of the ONE interspecific flower-colour × geography/climate/soil manuscript. This does not reopen the source-verified historical point estimates, does not add fruit traits and does not claim a genetic morph response to seasons.

## Main scientific question

Does the worldwide association of photographed white flowers with moisture emerge primarily from **which lineages and species flower/are photographed in different seasons**, or does the same direction remain among **the same congeneric species** across local seasons after geographic, climatic and observer/annotation support is matched?

A simple worldwide January–March versus April–June colour contrast would not answer this: a Northern spring is Southern autumn, tropical flowering can follow monsoonal/wet/dry timing, and changing photographed species composition can create a global shift even if no plant changes colour.

## Frozen local-time labels, never fabricated phenology

- Use the recorded **original observation's observed_on date**, not image-upload date, API annotation date or weather-scrape date. Missing/invalid dates remain UNKNOWN, never silently treated as spring.
- **Northern extratropics, latitude >=23.5 degrees:** spring March–May, summer June–August, autumn September–November, winter December–February.
- **Southern extratropics, latitude <= -23.5 degrees:** spring September–November, summer December–February, autumn March–May, winter June–August.
- **Tropics, absolute latitude <23.5 degrees:** retain exact month (and descriptive calendar Q1–Q4 only), **do not call Q1 dry or Q3 wet**. A real wet/dry label needs an independently byte-pinned *monthly* rainfall-climatology raster at the photo's original coordinate. Prefer within-cell seasonal rainfall anomalies or monthly ranks, with an explicitly aseasonal group where the local precipitation cycle is weak. Annual BIO12 does NOT classify a location's wet/dry month.
- Date-absent or coordinate-absent photographs are kept in UNKNOWN support groups, not counted as observed flowering in a chosen season. No direct inference about actual flowering outside the observed iNaturalist photographic opportunity frame.

## Three noninterchangeable estimands to report

**A. Observed seasonal photo/flower composition.** Within each domain × local season, report independent selected observations, taxa/genera, observers, ROI four-hue classifiability and design-weighted photographed-white conditional on classification; compare photo-opportunity weights and independent-observation weights. Separate iNaturalist Flowers-positive, nonannotated, unknown annotation, and ROI classifiability. This source composition may change through observer tourism/flower photography and seasonal species availability.

**B. Same named species in multiple seasons.** Require at least **three independently sampled four-hue-classifiable original observations in EACH compared season**, with >=5 as a sensitivity. Use the exact same eligible species in both seasons; calculate paired season differences of design-weighted within-species white probability. Paired comparison limits species turnover but NOT spatial/observer/individual-history confounding.

**C. Same species, same source geographic sector/latitude belt and original-site BIO12 stratum, multiple seasons.** Again require >=3 classified distinct source observations per season within each identical taxon × broad geography × moisture source cell, with >=5 sensitivity. This is stricter geographic support and may have *very few* qualified comparisons under a cap of 40 distinct observations/species. Report the matched taxon/cell counts **even if zero**, and avoid universal seasonal-evolution conclusions when support is thin.

**Do not add B or C to A as components of one causal variance partition**: they use different conditioning and often different species sets. Retain the original 42,111-taxon opportunity population, the separate old 100,543/45,953 measured-photo set, source photographs successfully classified, and each seasonal/matched complete-case count as DIFFERENT denominators.

## Prospectively specified environmental analysis only when the photo/model data exist

First report A–C, the cross-season species overlap, month-by-observer concentration, positive flowering annotations and classification missingness. If the selected data have adequate common genus/species/cell support, fit photo/within-species-weighted four-hue outcomes against the **same 13 original site environmental axes** with a local-season main term and a **season × annual BIO12 interaction**, clustering uncertainty by genus and explicitly checking observer/species repeat counts. For seasonal moisture mechanism, additionally fit **independently frozen photo-month climatological precipitation / atmospheric demand** (not a post hoc relabeling of annual BIO12); compare incremental held-out genus and geographic prediction to the annual 13-variable baseline, not in-sample fit alone.

The climate literature already has flower colour–moisture correlates; the novel angle here is **whether the apparent global effect travels across local flowering seasons once genus, species identity and observation bias are controlled**. A seasonal interaction is not direct evidence of pigment cost, within-individual colour plasticity or genetic selection.

Key robustness gates:
- Whole-genus bootstrap and familywise season × environmental contrast intervals, with within-species/observer structure reflected when repeat sampling permits; sparse cell intervals should be described as nonestimable.
- Observer-balanced or same-observer-in-both-season sensitivity, geographic-cell and matched-species comparison; neither a photo volume nor source observer activity is a sampled plant population abundance.
- Positive iNaturalist Flowers annotation versus the full ROI-classifiable sample, retaining NOT_VERIFIED and NO_ANNOTATION as distinct categories. Do not use flowering annotations in retrospect as if they existed in the original September source snapshot.
- Photo/observation inverse-selection probabilities, permanent acquisition failures, ROI nonclassification and missing-not-at-random sensitivity by season and water climate.
- A preregistered negative result is informative: if a global seasonal difference disappears for matched species, photo species turnover is a plausible explanation; if it persists in matched support, further observer, source site and month-climate analysis remains essential. Do not force an adaptation claim.

## Code and actual execution status

Implemented [source-attested seasonal audit](../scripts/analysis/mois_seasonal_flower_colour_opportunity_20261016.py) consumes the completed PRIVATE Mois sampled photo/ROI join and its SHA256 source receipt, optionally the approved independent annotation overlay and SHA256. It emits a private individual-source season long table, a private unsuppressed support report, and a public aggregate with small-season-cell suppression. The public evidence NEVER contains original photo IDs, observation UUIDs, site coordinates or individual taxon identities. Tropical wet/dry is explicitly NOT inferred. Its source integrity gates verify the exact original ROI-classifiable totals and do not promote partial source measurements.

Synthetic regression tests at [tests/test_mois_seasonal_flower_colour_opportunity_20261016.py](../tests/test_mois_seasonal_flower_colour_opportunity_20261016.py) cover hemispheric season inversion, tropical calendar bins, matched species versus turnover, original date/ROI/source hash integrity, annotation status and small-cell suppression. No original-date historical Mois source colour table with complete temporal identity, newly selected 42,111-taxon photo set, real phenology overlay or real seasonal all-13 environment fit has been reported. The old **16,801 taxa** global pooled all-13 colour R2=2.000% and BIO12 conditional slope=+3.770pp/source SD are **historical, NON-seasonal** estimates.

Offline example (paths are private and should NOT be committed):

    python scripts/analysis/mois_seasonal_flower_colour_opportunity_20261016.py \
      --private-joined-photo-colours /restricted/sample/private_joined_sample_flower_colours.csv \
      --sample-collection-receipt /restricted/sample/sample_collection_receipt.json \
      --out-private-season-photo-table /restricted/season/seasonal_photos.csv \
      --out-private-season-report /restricted/season/seasonal_support_private.json \
      --out-public-aggregate /restricted/season/seasonal_aggregates.json \
      --annotation-overlay /restricted/phenology/observations.csv \
      --annotation-overlay-sha256 PINNED_VERIFIED_SHA256

The annotation flags are optional; if absent, all become NOT_VERIFIED, **not** negative flowering evidence.
