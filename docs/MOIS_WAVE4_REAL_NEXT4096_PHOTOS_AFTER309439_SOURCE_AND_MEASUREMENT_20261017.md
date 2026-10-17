# MOIS fourth image-processing wave: source-audited next 4,096 original photographs after 309,439

**Research scope: global among-taxon photographic floral hue and seasonal sampling composition.** This is a new *actual original-photo pixel measurement* workflow, not another count of unclassified iNaturalist source candidates. It never uses the separate FCP within-species D/W/50 km photo-polymorphism results.

## Already measured source and why more is possible

All 256 FCP high-depth wave4 pixel worker partitions are actually complete. The previous exact [full original photo reassembly](https://github.com/zuizui0223/mois/actions/runs/39028400690) verified **309,439 distinct original photo-ID pixel measurement attempts**, **149,662 actual ROI-v4 four-hue classifications** and **149,653 photographic labels with observation dates**. The complete original photo/date table was SHA256-frozen as `546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e`. Only 19 records in that photograph cohort have an unverified observation day; unresolved observation/date status never makes an automatic hue a confirmed open petal.

The original independent [FCP 1,000-species source pool from 27 September 2026](https://github.com/zuizui0223/fcp/actions/runs/38884223530) contains **3,430,502 original iNaturalist photo IDs** with a checked photo provenance, originally source-identified coordinates and observed-on dates. These are **opportunities**, not 3.43 million genuinely classified flowers; photo pixels from an ID must actually be read, classified or marked technically unknown before the ID can enter MOIS's measured photo denominator.

### Fourth blinded 4,096-photo source cohort

`scripts/analysis/mois_original_fourth4096_after_309439_source_freeze.py` independently verifies both byte-level source hashes, reads all 309,439 already measured original photo IDs and scans the 3.43 million candidate IDs once. Every already measured photo ID is excluded **before** selecting a new photo. Original iNaturalist photo URL, source license, valid original observation ID/date, original taxon and genuine latitude/longitude are checked. Ranking uses **only a fixed SHA256 of the original photo ID**, stratified across 1,000 original source taxa and, subject to feasibility, 4,096 distinct source observers and observation IDs. The deterministic taxon-order hash avoids always favouring the first 96 low-numbered taxon IDs in a 4,096-of-1,000 balanced source.

No original photograph pixels, automated colour outcome, old source colour labels, observed calendar month, spatial coordinate, sunlight, climate, model result, pollinator or botanical stage is used in the hash ranking. Taxon balancing **is not field species abundance and cannot be treated as a representative new global sample without source-inclusion/effort adjustments**. Current originally measured MOIS/source taxa remain unchanged.

Only after exactly 4,096 unique original photo IDs are frozen into **16 × 256-photo source shards**, a second GitHub Actions job actually downloads each image from its license-permitted iNaturalist Open Data S3 source and runs the same original **ROI-v4/YOLO11+EfficientSAM nine palette→four hue** model used by the first three MOIS 4,096-photo expansion waves. Every original photo requires a signed terminal outcome: one of white, yellow/orange, red/pink, blue/purple, or a technical unclassifiable ROI. Image bytes are not persisted to the repository.

### Outcome replay and accurate paper denominators

The final job is **blocked until every one of the 16 image-pixel shards has a successful source-complete receipt**. It uses the original frozen model measurement attestation, verifies exact input photo/observation IDs against source SHA and each worker's 4×64-photo terminal receipts, and requires all 4,096 original photo identities be present exactly once. It refuses hue/taxon/date changes to the original 309,439 rows and refuses duplicate photo IDs. Every newly measured original source photo is attached to its **original independently source-recorded observation date**, not a day inferred from external climatological averages.

Only a **successful, source-complete** final image-inference run can report the actual new four-hue count and new date+four-hue count. The guaranteed numerical target is **4,096 additional new original photo-ID outcomes**, which should produce an enlarged 313,535 measured-original-photo union; the four-hue yield **must not be guessed** before photo RGB/palette processing. Existing 149,662 classifiable original ROIs and 149,653 date+colour photos remain unchanged.

The completed workflow emits an auditable new source GZIP plus a literal source/date/colour integrity receipt, hemisphere-corrected original photo-white season numerators, and 0.25° same-site/season-year local photographic-community support, with geographically clustered exploratory bootstrap intervals. This can directly inform the consolidated MOIS manuscript's **photo-source effect and where/when photographic taxa change** narrative without implying individual flower pigment plasticity.

## Explicit ecological limitations and open question

Earlier full 309,439-photo genus-equal maxT-999 annual rainfall and atmospheric-vapour coefficients remain a **different source cohort** after this actual new-image remeasurement. They cannot be passed off as 313,535-photo results without an actual rerun of WorldClim site sampling, the six-anchor model and all 999 genus resamples.

**Photo-specific open flower/petal/tepal identity, independent observer/plant identity across previous waves, locally measured flowering community abundance, observed-year/month meteorological variability, experimental pigmentation costs/benefits and evolutionary 50:50 equilibrium remain unidentified.** The new photo cohort is source-stratified; it may increase within-taxon photo depth but should not be assumed to be an unbiased source of true seasonal shifts without a sampling-weight/observer audit.

The finite source- and pixel-measurement workflow is `.github/workflows/mois-fourth4096-actual-fresh-original-pixel-photo-growth.yml`. The historical MOIS benchmark, all earlier photo growth ledgers and original FCP intraspecific W/D spatial results remain immutable.
