# MOIS: expand beyond 295,017 measured photographic originals using newly completed FCP wave4 ROI pixels

**Research scope:** worldwide photographed floral colour composition among taxa with observation date, not the independent FCP within-species W/D and 50-km photo-colour polymorphism hypotheses.

## Exact already-completed measured baseline

The source-verified 17 October 2026 [original photographed source run](https://github.com/zuizui0223/mois/actions/runs/39015660689) contains **295,017 distinct genuinely measured photo IDs**, **142,538 original automated four-hue ROI outcomes** and **140,247 photo IDs with both one valid source/public observation date and a four-hue label**. The latter include **99,817** original historic publicly verified photo dates (iNaturalist official 2026-09-27 S3 first, and exact numeric observation+photo ID API fallback), without assigning dates to unknown historic observations. Distinct original photo/taxon/hue records were SHA-checked against earlier MOIS original photos; FCP 180 completed authentic pixel partitions and three separately completed 4,096-photo MOIS measurement waves were deduplicated by literal photo ID.

The complete byte-level baseline hashes are:

- Original 295,017 photographed source GZIP: `44e3cb6b4c6174dfb1035de208343c7c6a2e6c321e2557a509224c4d9f2ac7b3`
- Its independently date-enriched original photo-hue GZIP: `b797d85464d2e3e39d38326196dd5ab3462d4e8e19b227dc62604b3bdf047a93`
- Actual original pixel+date GitHub Actions artifact: `mois-295017-original-photo-142538-genuine-hues-140247-verified-observation-dates` in run `39015660689`

## What this source expansion *actually* measures

The original [FCP wave4 provenance-locked pixel run](https://github.com/zuizui0223/fcp/actions/runs/38992543115) selected **49,800 original genuine photographs in 498 taxa**, completely blind to their image colour, and split photo pixels among **256 finite worker partitions**. The current MOIS 295,017 measured source already incorporates **180 genuinely completed pixel partitions** from this run. At the new source check, **212 GitHub artifacts** existed in the original FCP Actions run (one source firewall, 211 completed photo-pixel partition artifacts), i.e. approximately **31 more completed worker partitions than the previous MOIS baseline**. A workflow artifact's existence is **not** equivalent to successful four-hue classification for every image; complete pixel terminal receipts explicitly report both classified four-hue and technical-unclassifiable outcomes.

The new reproducible analysis is in `scripts/analysis/mois_append_finished_fcp_wave4_to_real295017.py`. Rather than trusting artifact names and guessing colours, the pipeline first independently reconstructs the FCP 49,800-photo source firewall and each complete worker's pixel-terminal/source-ID signatures through the existing `mois_reuse_fcp_wave4_completed_roi_partitions.py` code, generating an entirely authenticated FCP completed-partition photographic source. It then matches **each literal iNaturalist original `photo_id`** against the signed 295,017-photo measured MOIS union. It ignores every already-measured ID, refuses changed classification/taxon/original source date for overlapping IDs, retains documented alternate observation UUID-vs-numeric-ID representations without creating duplicate photographs, and admits only **new source photo IDs** with signed, completed pixels and terminal hue/unknown outputs.

For each newly admitted original FCP photo, preserve the original source `observed_on` value; add `qualified_date` only if it is an actual calendar date. Unknown technical hue stays unknown and its date never implies an independently confirmed open petal. The output source stream retains all 295,017 previously measured original photo records byte-consistently at the row-content level and never changes historical original source dates or colours.

### Publication and inference gates

The [new actual-image GitHub Actions workflow](../.github/workflows/mois-295017-next-finished-fcp-wave4-original-photo-pixels.yml) pulls the older SHA-signed source and all currently complete original FCP worker terminal receipts, runs real photo-ID/pixel verification, then produces two frozen original-photograph GZIP tables: (i) with original photo source fields, and (ii) with the existing historical-public `qualified_date` tags retained plus original date tags for *new actual FCP photographs*. It also archives a short hemispheric meteorological summer/winter photographic white-label-rate table and a numeric result receipt listing the actually added unique photo IDs, actually classified colour ROIs and additional date+colour outcomes. Counts must be recomputed from the finished-pixel source rather than inferred from 49,800 planned image IDs or from progress percentage.

The scientific main paper can then update its **separate** global annual moisture association / seasonal photographic colour-composition source denominators based on the newly verified photo pixels; previous results remain traceable source-cohort benchmarks. Every new source row has original photo/taxon/time data, but **none is photo-specifically confirmed to show an open flower**, and these photos are not randomized independent taxon/photographer/year sampling replicates. Prior annual WorldClim climatology does not represent observation-year weather. The add-on does not itself establish actual-season sun, water or pigment-fitness selection.

**Important:** FCP's full 49,800-photo run can remain not fully completed; no unmeasured remaining worker terminal, provider photo opportunity, technical missing hue, or unverified photo date enters the classified scientific denominator.
