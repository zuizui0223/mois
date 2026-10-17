# Recover exact observation days for remaining original source photographs with numeric iNaturalist observation IDs

**MOIS original photograph source, 17 October 2026.** The independent earlier iNaturalist September 2026 monthly S3 photo→observation-UUID join and exact-photo/numeric-observation public API fallbacks restored the historic dates for 99,817/100,543 original photos. After the 295,017 genuine image-pixel /142,538 automatic-hue combined source was assembled, **289,170** photos had a real observation date and **140,247** had both a date and a four-hue ROI.

## Measured exact public provider date result — completed successfully

The [real GitHub Actions original-source date-recovery run #39019305535](https://github.com/zuizui0223/mois/actions/runs/39019305535) completed all **36** bounded original numerical-observation API requests. It verified the *original literal photo ID in the original observation ID's current photo array*, the current public \u0060observed_on\u0060 day, and the source taxon identity for **1,732 of 1,751** previously undated source photos. Exactly **15** old observations were absent from the current public API and **4** original photos were no longer in their recorded public observations; no dates were imputed for these 19 IDs.

Of the **715** source four-hue-classified but formerly undated numeric-observation-ID photos, **706** acquired genuinely verified public dates. The same original **295,017** photo-pixel/ROI cohort therefore rises from **289,170** to **290,902** date-bearing photos and from **140,247** to **140,953** photo IDs with both an original four-hue label and a verified calendar date. Source photographs themselves and all original 99,817 historical S3/current API date records remain unchanged.

The signed [public result JSON](../results/mois_real295017_remaining1751_numeric_photo_dates_recovered_20261017.json) identifies the artifact ID **12095777799**, its exact old and new GZIP file SHA256 and all request response hashes. The unqueried **4,096** photo rows have UUID-based original observation IDs, not numeric iNaturalist IDs. They are *not* claimed to be date recovered by this numeric API study, despite already possessing **1,576** classifiable measured ROI photo outcomes.

## Two distinct unresolved-data mechanisms

A direct inspection of the actual SHA-locked source `mois_real295017_original_photo_hues_plus_verified_observed_dates.csv.gz` (byte SHA256 `b797d85464d2e3e39d38326196dd5ab3462d4e8e19b227dc62604b3bdf047a93`, [run #39015660689](https://github.com/zuizui0223/mois/actions/runs/39015660689)) finds **5,847** previously undated photos. The reasons and targetable source denominators differ:

| Source with original photo `qualified_date` still missing | Original photo IDs | Currently four-hue classified | Source original observation ID |
|---|---:|---:|---|
| Broad historic global MOIS | **726** | **316** | iNaturalist numeric |
| Early FCP reused original photo wave1 | **548** | **199** | iNaturalist numeric |
| Early FCP reused original photo wave2 | **477** | **200** | iNaturalist numeric |
| Initial 4,096 MOIS original-photo measurement wave | **4,096** | **1,576** | source observation UUID, **not** the numeric observation ID |

Thus **1,751** original photo IDs with **715** already measured classifiable four-hue outcomes have an exact original numeric observation identifier and can be checked through the official current public iNaturalist observation API with no image/pixel redownload. The other **4,096** lack a source numeric observation ID; the numeric API cannot be applied directly to their UUID strings. Do **not** misinterpret these 1,576 classifiable but undated images as assigned photographic dates simply because pixels are measured.

## Exact technical/photographic join

The new `scripts/acquisition/mois_original295017_remaining1751_numeric_photo_dates.py` selects only all genuine missing-calendar original photos with a source numeric observation ID, independent of hue, local season, latitude, climate and organism. It checks exactly **1,751** ID-unique source photographs and calls the provider API in at most **36×50 original observation-ID groups**, with at least 1.5 seconds between groups and bounded retry handling. For every returned record, the independent iNaturalist `observation.id`, original source `photo_id` membership in the observation photo list, and current observed-on date must all be valid. If the present source taxon ID changed, the original classifiable colour and taxonomic metadata still are NOT replaced; discrepancies are reported.

A valid retrieved day fills `qualified_date` only if the original row is still missing its date; the original `observed_on`, photo ID, image-pixel colour, coordinates and source taxonomy remain byte/row coherent. The already verified 99,817 historical-photo source dates and other previously observed source calendar days stay unchanged. Unavailable/deleted/moved provider photos remain unknown, not assigned a month or flower stage. Both 1,751-source photo status rows and a newly date-augmented 295,017-photo source file are saved as an Actions artifact, while only summarized source fractions should be committed to public Git.

## Ecological meaning

This increases the number of genuinely measured photographic colour states that can be placed in seasons **only to the extent the provider confirms both original IDs and actual valid `observed_on`**. It does not create further picture-pixel colour measurements or change the globally sampled photo ROI classification denominator; it recovers **missing calendar metadata**. Current public `observed_on` might have been edited since original acquisition, and neither the original ROI nor observation-level iNaturalist `Flowers` metadata independently confirms an open petal, tepal, bract or fruit.

Any subsequent MOIS seasonal-colour vs *observed-year* weather/temperature/water/solar model must use actual site/month/year exposure sources and adjust for the selection of photos/taxa/photographers across the source sampled months. The FCP within-species W/D 50-km and earlier MOIS annual six-climate-axis genus-999 maxT study remains separate.

Reproduction: [source-verified GitHub Actions workflow](../.github/workflows/mois-original295017-recover-1751-numeric-photo-observed-date.yml). Do not report any provider success count until its finite result actually completes.
