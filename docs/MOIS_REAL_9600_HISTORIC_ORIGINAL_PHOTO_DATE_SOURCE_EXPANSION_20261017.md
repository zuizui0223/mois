# MOIS: expand verified observation dates to 9,600 historical original photographs

**Date:** 2026-10-17. **Repository:** `zuizui0223/mois`. This is an actually executable source-identified photographic calendar analysis, not a replacement for the separate FCP within-species colour space/spatial polymorphism study.

## Why this is necessary

The currently measured MOIS source has **267,425 distinct real photographic ROI attempts and 129,795 assigned four-hue labels**, but only **81,867** have both an original `observed_on` and a classifiable hue. The globally broad **100,543-photo / 42,111-named-taxon historical MOIS** block is missing observation dates *in this particular combined analysis export*, although all **100,543 original photo IDs and observation IDs are still present and unique** and 45,953 of these original image ROIs are classifiable. The original photo colours, abiotic data, and annual water/maxT inference remain unchanged.

The previous independent blind-200 source-identity [pilot](https://github.com/zuizui0223/mois/actions/runs/39007591041) actually recovered all 200 publicly recorded observation dates with the literal original photograph ID in the correct observation's full photo list, and current public taxonomic assignments unchanged. Its original photo selection salt `mois_historical_dated_source_pilot_v1:` is used unchanged, so the old **200 are exactly the immutable prefix** of the following 9,600, with no source photo re-selection.

## This step: 9,600 genuine original photo IDs, bounded provider usage

`scripts/acquisition/mois_historical9600_source_date_expansion.py` reconstructs the 267,425-photo SHA256 source union (source photo CSV SHA256 `fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9`) and selects **9,600 of its original 100,543 historical photo IDs**, deterministically by original ID SHA-256 hash **before looking at colour, date, observer, latitude, climate or species**. It then:

- Checks exact membership, taxon and observation IDs of the previous **200** recovered photographs from the old successful Actions run; uses their already retrieved dates instead of querying them again. Its source record CSV SHA256 is `5a8a4b5e5af8f4a3c2187b1e38b7dbdbafb88de035e376c28c8a0da43f946dff`.
- Retrieves metadata for only **9,400 additional original observation IDs**, in **188** bounded batches of 50 IDs, separated by at least **1.5 seconds**, with explicit retry/backoff/pagination guards. With the prior 200, the *total cohort is below the iNaturalist API recommendation of 10,000 records per bulk request*; this is not a scheme to scrape all 100,543 through repeated API batches.
- Accepts a new dated source observation *only if* the exact originally sampled `photo_id` is found in the current public `observation_id` photo list. Current taxonomic changes are reported without reassigning original taxon identity; deleted/changed/date-null records remain unknown.
- Stores all valid original image classifications, including automatic **unclassified/mixed** photo outcomes. A technical missing hue is neither a chromatic flower nor a botanically confirmed flowering state; current public observation-level `Flowers` annotation is not photo-specific proof of an open petal.
- Emits separately the **original-source 9,600 photographic monthly/season denominators** and the **whole 267,425-photo union after date-only enrichment**. The original old classifiable+dated count of 81,867 can only increase by the exact number of independently ID-verified, previously undated historical classified source photos. All original photo ID/hue rows remain unchanged; actual 2026 weather is never inferred from WorldClim monthly normals.

A real partial-batch source signature is saved if the provider returns a transient error, but an unfinished 9,600-photo run **must not** be described as having retrieved all 9,600 original dates. The result receipt includes source-ID selection SHA, original/current taxonomy agreement and date-recovery status, full source photo-CID and seasonal CSV SHA, per-batch independent API response hashes, the current public retrieval date and the old 200 source provenance.

## Why stop at 9,600 rather than 100,543 API requests?

iNaturalist explicitly asks users to prefer the [official Export tool](https://www.inaturalist.org/pages/api+recommended+practices) for large data requests above 10,000 records, and to stay at or below around 60 API requests/minute. The [iNaturalist Open Data metadata](https://github.com/inaturalist/inaturalist-open-data/blob/main/Metadata/README.md) are a separate **monthly snapshot** with `photos.csv.gz` linking `photo_id` to `observation_uuid`, and `observations.csv.gz` linking `observation_uuid` to `observed_on`. These are suitable for eventual source-wide historical photo/date recovery with **two explicit ID joins**, a dataset snapshot/source timestamp and a coverage bias audit; they should not be assumed to contain all non-openly licensed original photos. A photo's numeric `observation_id` is **not** the Open Data `observation_uuid`.

The original historical 42,111 named taxa provide breadth but often only one or two photos per species. More observation dates will enable **honest source-conditioned seasonal colour descriptives and broad taxon/continent support**, but do not automatically give enough botanically verified open flowers from the same species×fine locality×year across two seasons for environmental response or genetic flower-colour plasticity. The existing original 267,425-photo water maxT/global-vs-geocell results stay in the separate annual study section.

## Verified ecological inference gates

1. **Photo date verified:** explicit original photo and observation identity, exact current `observed_on`, time and provider source and taxonomy-change flags.
2. **Photo outcome:** original frozen ROI-v4 hue is white/chromatic or technically missing. No annotation substitutes for a photo-specific open-flower stage.
3. **Time-varying climate:** real original observed-year and month site-specific climate rather than BIO or 1970–2000 monthly normals.
4. **Sampling support:** same taxon × photo site × meteorological season-year, balanced source/observer effort, photo colour-label missingness and sensitivity to source taxonomy/season composition.
5. **Ecological conclusion:** positive photohue/climate association is not automatically a universal physiological pigmentation response, pollinator mechanism, selective cost–benefit or 50:50 equilibrium.

**Reproduction:** `.github/workflows/mois-historical-9600-observed-date-verified-source.yml` runs original photo/date authenticity and synthetic anti-pseudoreplication tests; pulls both already-successful signed upstream source artifacts, then runs all provider batches and uploads a checksum-stamped finished or partial result as a GitHub Actions artifact. The main results ledger and manuscript should be updated only after an actual completed source run.
