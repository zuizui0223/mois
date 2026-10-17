# MOIS: recover original observation dates for 4,096 measured flower-colour photos with UUID-only observation metadata

**Why these original photographs matter:** the actual 295,017 measured-photo MOIS atlas retained **4,096** initially processed genuine photo pixels whose source observation identifier is **not numeric** and whose `observed_on` / later `qualified_date` fields were left blank. **1,576** of those same photo IDs already have a successful original four-hue ROI classification. They are genuine image measurements, but cannot yet enter a seasonal flower-colour numerator until the **correct original observation date** is matched.

A separate real-source API recovery of **1,751** missing-date original photographs with **numeric** iNaturalist observation identifiers successfully recovered 1,732 dates /706 dated four-hue photos without touching these nonnumeric 4,096 photos. This workflow complements, rather than duplicates, that exact provider API evidence.

## Existing authoritative photo/day resource, no provider API/image redownload

The [signed 2026-09-27 global FCP 1000-taxon source](https://github.com/zuizui0223/fcp/actions/runs/38884223530) contains **3,430,502 distinct original iNaturalist photo opportunity identifiers**, their source `observation_id`, source `inat_taxon_id`, and **`observed_on`** observation calendar dates. The source CSV is SHA256 locked to:

`58e037a210e23d3ff5d6f04d392665ada5ee9ac5149d0676672c239cf092df01`.

This resource is *observation/photo metadata*, NOT 3.4 million newly measured flower colours. The 4,096 genuine photographic pixels already exist in the separately [successful MOIS real295017 source](https://github.com/zuizui0223/mois/actions/runs/39015660689), whose date-and-colour CSV SHA256 is `b797d85464d2e3e39d38326196dd5ab3462d4e8e19b227dc62604b3bdf047a93`.

The reproducible `scripts/acquisition/mois_first4096_undated_uuid_photo_days_from_original3430502.py` streams the entire original 3,430,502-row source, calculating the **full raw metadata SHA256** while selecting only source lines whose *literal original photo ID* is one of the 4,096 already measured. Each accepted `observed_on` must have the same original **photo ID and taxon ID**, a nonempty actual FCP source observation identifier, a valid ISO calendar day, and the independently source-pinned 2026-09-27 metadata snapshot. The old nonnumeric observed `observation_id` is retained as original source provenance, not **falsely assumed equal to the current numeric iNaturalist observation ID**.

Every one of the 4,096 attempted photo IDs receives a matched/unmatched status. No absent source IDs, invalid dates, changed taxon IDs, colour-unknown ROIs or unverified floral-organ states are inferred. The source method never downloads a new image or makes iNaturalist API calls. It produces a 4,096-original-ID photo/date crosswalk and strict SHA+source/coverage receipt as a GitHub Actions artifact. Aggregate new photo dates and classified-photo dates are published only after the full signed source scan succeeds.

## Ecological inference

This is a **source date-enrichment** operation on 1,576 original ROIs that were already automatically classified white/yellow/orange/red/pink/blue/purple. A valid iNaturalist `observed_on` does not prove that the pixel region is an independently botanist-confirmed **open petal/tepal** (versus bract/fruit/background), nor is it camera EXIF time. The photographic sampling and original plant identifications can have region/season observer selection bias. This must not be described as a new taxon recolouring, new photograph acquisition, new field observation, or causal month-specific solar/water pigment response.

After the exact source photo/day match is complete, these dates can be merged into the separately expanded **302,242 real original photo IDs /146,061 four-hue ROI classifications** together with the independent additional 1,732 numeric-ID source-day repairs, preserving all original source photos, labels and geographic identifiers. FCP's within-species W/D/50km local original photo-colour results remain separate from MOIS's worldwide interspecific photographic colour/climate study.
