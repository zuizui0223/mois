# FCP × mois: independent global flower-colour mainline from related source photos

Date: 2026-10-16. **Active mois paper: global *between-species* photographic flower-colour × environment.**

The 2026-10-14 transfer was a byte-preserving migration from historical fcp macro PRs #193/#194. It was **not a new source-photograph reanalysis**, nor an independent replication of the still-active FCP New Phytologist within-species manuscript. This document and `study_scope.json` freeze the **analysis boundary**, not the historical original numeric estimates anew.

## Two related original-photo frames, two distinct estimands

| | mois (this repository) | fcp (other repository) |
|---|---|---|
| Main response | Named-species four-hue photographic composition; white proportion and its environmental directions | Within-named-species colour diversity `D`, local depletion and pre-frozen within-species colour-vector axis `W` |
| Primary frame | 100,543 unique original photographic opportunities / 42,111 named taxa; **45,953 classified, 54,590 missing hue** | 149,900 separately sampled high-depth photographs / 1,499 named species |
| Main environmental comparison | 20,694 named species / 1,934 qualified genera, climate26 | Conspecific spatial/phenotype structure in three high-depth cohorts |
| Primary analytical unit | Species hue proportion, genus-centred, geographic-support-aware | Conspecific observations/pairs and within-species phenotype variation |
| Data relation | Same iNaturalist original-photo ecosystem/source lineage; **cross-repository fixed photo-ID intersection not yet measured** | Not an independent botanical source; no cross-paper pooling |

Original source-classifiable white fraction 21,574/45,953 = 46.95% is a **conditional photographic fraction**, not white morph genetic prevalence or 50:50 fitness equilibrium.

## mois primary findings with compulsory counterevidence

- Global genus-equal photographic white–annual-rainfall (BIO12): **+2.2545 percentage points per original environmental SD**, 95% within-scope maxT **[+0.3187,+4.1903]**, among 20,694 named species/1,934 genera.
- Global vapour positive: **+3.9519 pp**, maxT **[+1.1993,+6.7044]** under the original climate26 comparison.
- The strict-image BIO12 interval **[-2.046,+4.504]** and within-genus/coarse-geographic-cell interval **[-3.556,+5.172]** include zero. The fixed East geographic-conditioning comparison 24-contrast maxT interval **[-0.258,+6.613]** also includes zero.
- Neither shared climatic measurement, nominal genus conditioning, missing-colour photos, nor photographic colour classifiers identify physiological pigment costs, universal local adaptation or a robust geographic causal moisture mechanism.
- The original same-photo *within-species* rainfall contrast in `results/source_receipts/` is a **cross-scale falsification/sensitivity audit** for the mois interspecific association. It is not an FCP prospective W result, not an independent sample source, and not a replacement for FCP's species-composition-preserving geographic null.

## Seasonality already transferred to mois as *source-only* evidence (16 October)

Since this firewall branch was created, mois main has received a separate seasonal source/QA transfer under `docs/MOIS_SEASONAL_WHITE_MOISTURE_RESEARCH_20261016.md` and `docs/MOIS_SEASONAL_FCP_SOURCE_TRANSFER_20261016.json`. The dated opportunity ledger already covers **3,430,502 original photo IDs / 1,000 taxa**; a first 64 original-pixel technical pilot yielded 27 automatic colour classifications, but **open-flower organ identity has not been independently validated**. The separately reserved within-site seasonal strata contain **no opened colour pixels**. These are source and technical measurement achievements, **not measured seasonal white-frequency/climate effects**.

The seasonal programme uses within-taxon/site contrasts to **control geographic and taxonomic turnover**, but its global target is flowering-photo hue composition and observed-month environment. It must not import FCP's photograph-pair `D`/prospective `W` as a seasonality result, and it must not pool the 3.43-million photo *opportunities* with the historical 45,953 classifiable hues. The technical source code was copied from an FCP development branch with hashes; the ecological result ownership is now mois.

## Next mois-specific analysis: seasons + more of the upstream floral photographic opportunity

Keep the **global between-species** biological target. Extend the original photo-ledger coverage upstream where permitted and identify photographed *open flowers* plus actual `observed_on` dates. Do not treat plant-present/no-open-flower photographs as classified flowering-colour observations, nor silently add unknown photo colour to the white denominator.

Build auditable `(named_species, location_cell, observed_year, observed_month)` strata:
1. Distinguish original flower-photo opportunity counts from visible open-flower/hue-classified counts, and record missing colour, image-quality and observer/study effort.
2. Estimate monthly/seasonal photographed white and four-hue composition, with species or genus-balanced reporting and explicit within-species/between-species decomposition; the **primary mois inference remains among taxa/genera**.
3. Join **year/month-specific** temperature, solar radiation, precipitation, vapour/VPD and soil moisture to actual photo dates. WorldClim/CHELSA *long-term* BIO values alone do **not** become seasonal climate by splitting photos into quarters.
4. Include spatial location/latitude/longitude, habitat and sampling opportunity in the model. Soil, elevation and terrain are relatively time-stable covariates, not seasonal exposure estimates. Audit within-genus geographic support and leave-region/genus-out performance.
5. Keep strict-photo, unclassified-hue/MNAR, taxonomic, season/source-coverage and spatial simultaneous-uncertainty comparisons. Distinguish calendar month in northern/southern hemispheres from biologically comparable season and avoid equating observed seasonality with seasonal pigment plasticity.

The underlying full source photo/environment ZIPs are **not in mois Git**. Existing CI can guard the historical transfer/receipts, but this proposed dated-photo analysis cannot honestly be reported as already refitted.

## Engineering gate

`python scripts/check_mois_migration.py` checks byte identity of original source files; `python scripts/check_study_scope.py` additionally requires the interspecific frame, 1,280-direction ledger source receipt, positive global coefficient, strict-photo and paired-geographic **nonconfirmation**, and absence of FCP W/local-depletion main outputs. Both run in mois CI. The script checks committed frozen receipts without downloading archived photos.

See the complement in `zuizui0223/fcp/docs/STUDY_SCOPE_FCP_MOIS_20261016.md`.
