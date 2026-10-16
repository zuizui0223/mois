# MOIS | Source-verified two-season sampling and flowering-stage validation

**Owner:** \`zuizui0223/mois\`. This supports the new seasonal white-flower × actual observation-month climate programme. The historical interspecific GEB moisture coefficients and the FCP polymorphism/W study remain independent.

## What can already be said

- The [FCP 2026-09-27 source ID gate](https://github.com/zuizui0223/fcp/actions/runs/38884223530) successfully reconciled **3,430,502 original licensed photo IDs** for exactly **1,000 outcome-blind selected taxa**; this is not 3.43 million flowering images. The earlier 357-taxon collision-policy mismatch was resolved before the source was admitted.
- The [complete original observation-date/cell audit](https://github.com/zuizui0223/fcp/actions/runs/38901290387) found 3,430,502 dated original photo opportunities and 2,416,225 photos with original accuracy >0 and <=5km. Source artefact: \`fcp-1000-original-monthly-observation-calendar-all-photo-IDS-NO-COLOURS\`. Its site-month source CSV SHA256 is \`7691c7cf975da4e73494ad3fa0f2cf8df52302e84206017dd525051cd672b432\`; signed source JSON SHA256 \`0442658d600c25e59289037bd311024648d1eb0185088960a0e4fd1d4dd3d658\`.
- The FCP real first-unit [ROI measurement](https://github.com/zuizui0223/fcp/actions/runs/38901816890) opened **64** genuine original source photos: 27 classifier-classifiable flower-colour regions, including **4 white–yellow same-ROI organ mixtures**; 37 classifiability failures. The remaining 8 white and 15 chromatic are *automatic image labels*, not independently verified open-flower petal/tepal phenotypes. Its two taxa are **not** a population-level or seasonally matched effect.
- The **same-taxon, same 0.25° geographic cell/year, ≥24 distinct source observations in EACH of two local seasons** rank is **268 candidate cell-years in 39 taxa**, assuming climate support through year 2025 and ≥2 observers per contributing observed calendar month. Full-reproduction target: \`season_depth_24_source_only_geocell_year_priority.csv\` SHA256 \`ef95638da60b010c8671877854163132be304f5b7245edf9c5e2f9536265102a\`.
- This first-stage rank is source **photographic opportunity**, not verified blooming, independent **between-season** observers, seasonal white fraction, weather effect, or native range. Do NOT insert the 39/268 rows as 39 colour-polymorphic taxa.

## Operational separation

**FCP only:** preserve original dated provider metadata and licensed photo provenance; choose/execute original-image blinded ROI-v4/EfficientSAM, return immutable source photo, observation and observer IDs plus SHA-sealed image-classification terminals; retain its historical 106,376-photo photo-first RGFCA sample unchanged.

**MOIS now:** regenerate the photo-opportunity/season rank from FCP's signed artifact, audit actual observed_on/site/year and independent observer overlap, choose *additional* source-only temporal photo candidates without looking at photo colours; own flowering-stage QC and monthly climate joins/uncertainty/inference. MOIS must never change FCP's historical colour method or imply reuse of old 100,543 measured photo outcomes.

The new \`scripts/analysis/fcp_source_month_season_priority.py\` is independently executable inside MOIS, with the same input schema, SHA, temperate hemisphere calendar/December season-year matching and no phenotype leakage. \`.github/workflows/mois-season-priority-real-one-shot-20261016.yml\` reads the public FCP source artifact and requires the 39/268 original result SHA; it does not re-download the 36GB provider archive or open pixels. The priority algorithm tests live in \`tests/test_mois_season_priority.py\`.

## Next evidence gate

**Independent flower-organ identity.** Automatic RGB/ROI colour alone is inadequate, particularly small-flowered inflorescences, bracts, spikelets and fruiting observations. Before a primary white-rate denominator, two independent botanical reviewers should confirm (1) an *open flower visible in the pictured photo*, and (2) colour ROI placed on a petal/tepal or explicitly distinguished floral structure. Reviewer disagreement must remain pending until adjudication; missing iNaturalist phenology annotation is **unknown**. Open flower and organ review must be outcome-independent and blind to climatic covariates, and no photo may be replaced for its observed colour/status.

**True temporal ecology.** From botanical-QC-passed seasonal photographs compute independent-observation-level white/nonwhite frequency *within the same species, site and year*, with observer-heldout checks. Join genuinely available monthly tmin/tmax, solar radiation, rainfall, water deficit and VPD; derive within-site/month climate anomalies separately from between-site averages. Soil, elevation, latitude/longitude and terrain are time-stable site adjustments, not drivers of seasonal change. Do not extrapolate from 0.25° source cells to biological populations, or use 2025 weather as a stand-in for missing 2026 months.

No seasonal ecological coefficient, pigment fitness trade-off or adaptive inference is supported by the source counts alone. Original inherited 100,543-photo MOIS interspecific results remain the separate historical main paper until this new seasonal study qualifies.
