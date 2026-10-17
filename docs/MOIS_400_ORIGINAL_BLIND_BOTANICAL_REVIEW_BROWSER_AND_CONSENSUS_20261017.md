# MOIS — original 400-photo independent botanical review: practical handoff

**Study owner:** zuizui0223/mois. **Status:** runnable reviewer interfaces and source-SHA-bound qualification logic completed. **No human reviewer responses yet.** This is not a new floral-colour effect, and is distinct from the FCP intraspecific W/D paper.

## Original photographic source and scientific boundary

The [original source-only 400-photo selection](https://github.com/zuizui0223/mois/actions/runs/38929402279) has exactly 400 distinct licensed photo IDs, 400 source observation IDs, 400 distinct photographer IDs, five taxa and ten taxon×local-season source groups of 40. Each taxon has photos from two local seasons of the same source geographic cell and source season-year. The fixed photo-ID selection is blind to visible flower colours, phenology labels and climate; it does not reuse the old MOIS 100,543 measured photos or the first 64-photo Asclepias technical ROI run.

**Photo source SHA256:** b961b50ee6efe9ec08a5575a8eb1e34b1edf7b7fcb426226fb11ee901cfd3329.

The [existing independent reviewer blank CSV handoff](https://github.com/zuizui0223/mois/actions/runs/38930414968) was source-sealed before any botanical assessments. Reviewer A and B each receive the same 400 photographs but in different hashed order. They do not receive observed calendar date, local season, fine source location, climate or automated photo-colour outcome.

## How reviewers use the browser

The [MOIS reviewer photo-browser generator](../scripts/analysis/mois_render_400_botanist_photo_webforms.py) independently rebuilt two lightweight single-page HTML files using only exact source-attested blank photo/attribution rows. [Successful proof](https://github.com/zuizui0223/mois/actions/runs/38988600457), 10 regression tests.

- [Reviewer A browser artifact](https://github.com/zuizui0223/mois/actions/runs/38988600457/artifacts/12083059923)
- [Reviewer B browser artifact](https://github.com/zuizui0223/mois/actions/runs/38988600457/artifacts/12083478803)
- [Original UI/source SHA receipt](https://github.com/zuizui0223/mois/actions/runs/38988600457/artifacts/12083369080)

Each reviewer downloads their own HTML and opens it in a desktop/mobile browser with Internet access; the page displays one public licensed iNaturalist original photograph at a time and keeps attribution visible. No photograph pixel bytes are embedded in the review artifact. Each reviewer enters a stable and distinct reviewer ID and independently fills **open_flower_stage**, **roi_organ**, and optional notes.

Allowed opening stages: open_flower, buds_only, fruit_only, no_flower_visible, uncertain. Allowed ROI tissues: petal_tepal, other_flower_organ, spikelet_glume_bract, fruit_seed, leaf_stem, background, uncertain.

The browser's **Export CSV** stores all 400 unchanged original tokens/photos, retaining empty answers as empty (never as negative flowering). **Import saved CSV** allows resuming, and refuses unknown source tokens, altered links, mixed reviewer identities, wrong answer codes and partial replies. It has no response upload server, cookies or analytics. The public image links require Internet access and could be used to look up original metadata deliberately, so date/season/colour blinding is procedural rather than cryptographically irreversible.

**Reviewer independence:** A and B should not share their answers or browse the other's work before individual CSV submission. Reviewer identifiers must be different.

## Blind response reconciliation and third adjudication

Use [the MOIS consensus/qualification script](../scripts/analysis/mois_reconcile_400_independent_botanical_reviews.py), not a manual spreadsheet join, to reattach original observation/month/site/season metadata **only after** blinded answers are received. The script checks exact original source photo manifest SHA, season/source attribution/URL/token identity, all 400 photos and 400 distinct reporters, 40 photos in each of ten source taxon-season groups, stable reviewer A/B identity and exact reviewer header schema. It rejects any review sheet containing extra date/season/colour/environment columns or altered original image links.

Both botanists must independently agree that a **visible open flower's petal or tepal** is the image region to be assessed. For genuinely discordant full reviews only, an optional separate third adjudicator file must contain these exact CSV columns in order:

review_token,adjudicator_id,adjudicated_open_flower_stage,adjudicated_roi_organ,adjudication_notes

The adjudicator must be neither reviewer A nor reviewer B. Disagreements without completed adjudication remain pending, not automatically nonflowering. Non-open buds/fruit/no-visible-flower, nonpetaloid tissues and uncertain responses are distinct QC states. A human-confirmed open petal/tepal photo is **eligible for later blind pixel-colour classification**, not yet a measured white or chromatic phenotype.

[Actual 400-source consensus test](https://github.com/zuizui0223/mois/actions/runs/38988328639) **passed 12 tests** and reassembled the original two still-unfilled forms without dropping any photo. The signed real result contains 400/400 pending human confirmation, 10 original source season groups, 0 submitted full double reviews, 0 image ROI hue measurements on the 400 new photos. [Actual original 400 botanical sample-flow receipt and CSV](https://github.com/zuizui0223/mois/actions/runs/38988328639/artifacts/12083677941).

This **does not mean zero blooming plants**. It means the selected photographs have not yet been independently botanically reviewed. The old 64-photo technical pilot also cannot establish population white-frequency seasonality.

## Next stage and reproducible statistical boundary

Once two independent reviewers have supplied photo-specific stage/tissue decisions and conflicts have third-person adjudication, Mois can report counts of **reviewed and confirmed open-petaloid original photos** by taxon×site-year×season. Then the FCP upstream frozen ROI-v4 measurement backend may be applied to those source-identified images, keeping ALL original 400 attempted photo IDs, nonflowering stages, undecodable photos and mixed original flower organs in the source funnel.

Mois, not FCP, will then join genuine observation-month TerraClimate weather, derive within-site anomalies versus same-calendar-month 1991–2020 climatology, adjust observer/site/year/phenology effort, and test genuine seasonal photographed flower-white frequencies. Soil/elevation/terrain are stable site adjustments. The current panel is five predominantly Northern-Hemisphere taxa/site-years, not a globally representative or biologically independent sample of worldwide seasonal flower colours.

No unreviewed photo, original RGB hue label, current iNaturalist observation annotation or climate-only gradient can prove flower-pigment genetic polymorphism, physiological adaptation, reproductive costs/benefits, or causal climate selection. Historical MOIS global GEB congeneric rainfall coefficients and FCP W/D intraspecific phenotypic-space results stay separate.
