# mois — provisional Global Ecology and Biogeography submission plan and cover-letter draft

**Date:** 14 October 2026. **Status:** paper draft for coauthor review, **not submitted and not ready for anonymized external review**.

## Scientific editorial decision

Working title: **The geographic and measurement dependence of flower-colour–moisture associations across a photographed global flora**

**Claim to publish:** Among 20,694 congeneric taxa, white photographic flower-colour probability has a modest positive precipitation/vapour gradient under broad spatial adjustment. Strict flower-image selection, within-genus×coarse-cell condition, matched same-species alternative spatial models, independent botanical descriptions, and geographically distant taxonomic transfer do **not** collectively establish a universal response. On exactly matched source photographs, the between-congeneric rain–white slope is positive whereas the within-species photo slope is inconclusive. This is an analysis of **conditional observational biogeography and inferential scale**, not pigment fitness or genetic morph maintenance.

**Novelty against verified precedents:** Dalrymple et al. (2020; DOI 10.1111/nph.16737) already found precipitation/other ecological correlations in 339 species; Trethowan et al. (2024; DOI 10.1111/jbi.14751) already found white flowers common in wet lowland Malesian environments. Accordingly, a result framed merely as "white in wet places" would be incremental. The potential contribution here is a much broader original species-conditional photographic comparison with **matched support, measurement nonresponse, multiple-comparison and transfer diagnostics** that show precisely where global slopes can and cannot be generalized. This is *not* a confirmed physiological mechanism or demonstration that photo quality causes the broad association.

## Provisional journal-format contract (verify again before actual submission)

Global Ecology and Biogeography (GEB) is a provisional thematic fit, not an editorial endorsement. Publicly stated guidance retrieved 14 October 2026: Research Articles are typically ~5,000 words in the main body, have roughly 50 references and six to eight total display pieces; a **structured abstract under 300 words** must include **Aim, Location, Time period, Major taxa studied, Methods, Results, Main conclusions**; references may be free-format at initial submission but should be internally consistent. Manuscripts use **double-anonymous review**, require an **anonymized manuscript** and a **separate identifying title-page**; submit supporting information separately. The publisher requests a cover letter as a separate **PDF** with a paragraph under 250 words explaining relevance. The current author-facing manuscript includes historical GitHub URLs, original SHA fingerprints and source PR descriptions, so it is **not anonymized** and must not be uploaded directly as a blinded review file.

Journal source: https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html

## Cover letter text (working draft; author review required)

To the Editors of *Global Ecology and Biogeography*,

Please consider our research article, “The geographic and measurement dependence of flower-colour–moisture associations across a photographed global flora,” for publication as a Research Article.

The study addresses a general macroecological problem: when does a trait–environment association observed across a large species set remain identifiable as sampling grain, geographic conditioning and trait measurement change? We analysed archived flower-photo hues across a global opportunity frame of 42,111 named taxa, with 20,694 species eligible for a congeneric climate comparison. Photographic white increased with annual precipitation and atmospheric water vapour under broad geographic adjustment, but the positive annual-rainfall direction was inconclusive in strict-image and local congeneric comparisons. A model comparison using exactly the same named species demonstrated dependence on spatial conditioning; a separate original-photo-matched analysis did not support the same positive slope within species. Independent botanical colour descriptions and distant geographic prediction tests further constrained generality. Our result is not the already-known association of paler flowers with wet regions in isolation; it delineates why macroecological extrapolation to a shared flower-pigmentation mechanism can fail. By reporting all negative and missingness diagnostics alongside broad positive associations, we provide a testable framework relevant beyond floral traits.

[Author to confirm that this manuscript is not submitted or under review elsewhere.] Related high-depth **within-species** floral-colour distribution research is maintained as a distinct, nonoverlapping manuscript, using partly related opportunistic observation infrastructure but targeting species-level phenotype-space and local geographic organization rather than interspecific rainfall coefficients. Both data provenance and inferential separation are disclosed.

[Corresponding author details, coauthor confirmation, conflict statement and suggested reviewers to be added.]

## Submission blockers — no substitute for human author approval

- **Authorship/title page:** complete and confirm names, affiliations, email addresses, ORCID, contributions, funding, acknowledgments and competing interests, separately from a blinded main manuscript.
- **Anonymized manuscript and supplement:** the automatically derived blinded manuscript and numbered review copy now omit explicit public repository identity, but must still be checked by human authors for contextual self-identification. **Reviewer-anonymous source data access is not yet available**. Do not upload the present copy as an anonymous review package before the data-link and metadata gate is passed.
- **Data and code:** deposit the complete, licence-compliant **original photo and environmental source** bundles or an equivalent accessible provenance/trait table with stable DOI, checksums and permission review. GitHub only contains reconstructed numerical output and original-source verifier code, not the full input archive. Full fresh source replay **has not been performed in mois**.
- **Figures:** original Figure 1 plus source-reproducible Figures 2–4 are committed. Inspect all four at journal print/readability resolution and licence/accessibility settings, and confirm the two tables and every figure caption before final review export.
- **Scientific independent measurement:** source NHM binary botanical hue with original FCP sites does not verify external geographical transfer. A second independently curated botanical/spectrometric dataset would strengthen the paper but is **not an analysis already completed**.
- **Final copyedit/reference audit:** publisher-verified citations for WorldClim, SoilGrids and EfficientSAM are inserted. Still verify the exact original iNaturalist snapshot, source environmental layer/depth specification, custom ROI-v4 implementation, third-party photo rights, all references and author metadata. The current ~4.1k-word main body, 269-word structured abstract and four figures/two tables fit typical GEB Research Article limits. Do not invent any dataset citation to close this gate.

## Production update: four visual outputs and anonymous review copy (14 October 2026)

The current independent [GEB Research Article](MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md) contains **four figures and two tables** (six display pieces); the newly source-locked Figure 4 shows all four original hue-category responses to precipitation and atmospheric vapour pressure and confirms that **only photographic white's original category-wise maxT interval excludes zero in both axes**. All four SVGs, source files and their regenerators are versioned.

Two **editable Markdown derivatives**, [blinded reviewer manuscript](MOIS_GEB_BLINDED_REVIEW_COPY_20261014.md) and [continuously numbered reviewer copy](MOIS_GEB_LINE_NUMBERED_REVIEW_COPY_20261014.txt), are derived from the exact author-facing paper using `node scripts/build_mois_blinded_review_copy.mjs --check`. The [separate identifying title-page template](MOIS_GEB_IDENTIFYING_TITLE_PAGE_TEMPLATE_20261014.md) and [cover letter for later PDF export](MOIS_GEB_COVER_LETTER_FOR_PDF_20261014.md) require human details and permission.

Seven scientific references are verified from publisher/conference pages, including WorldClim, SoilGrids and EfficientSAM ([reference and journal-guideline audit](MOIS_GEB_REFERENCE_AND_GUIDELINE_AUDIT_20261014.md)); the exact historical iNaturalist snapshot citation, photographer licences and bespoke ROI-v4 source version are **not** invented or declared complete.

**Critical GEB policy:** public GitHub tracking and successful model/output CI are *not* sufficient as a stable data/code publication archive. Anonymized reviewer data/code access is required **during review**, with suitable licensed full original source data and a stable public archive/DOI required for publication. The manuscript's private reviewer link is a clear **NOT YET PROVIDED** placeholder, not an actual URL. Our anonymous manuscript is a *working derivative*, not evidence that the required reviewer data access is ready. No empirical positive result was selected during this manuscript formatting.

## Source archive recovery route now identified

The exact historical Actions inputs were located at FCP photo artifact **11587115220** (original photo ZIP SHA256 `2d6097ba…`, approximately 34.3 MB and scheduled expiry **2027-01-06**) and all36 environmental artifact **11653183273** (original environmental ZIP SHA256 `6e287294…`). The successful original FCP run [#38502660130](https://github.com/zuizui0223/fcp/actions/runs/38502660130) independently validated both source ZIPs and source numerical calculations. [Acquisition, SHA256 and photographer licence recovery guide](MOIS_GEB_ORIGINAL_ARCHIVE_RECOVERY_AND_LICENCE_20261014.md) plus `scripts/recover_mois_original_inputs.sh` are committed. **Both exact source ZIPs have now been downloaded and verified in the working environment, but neither has been deposited in mois or in an anonymous stable reviewer archive.** Recovering the ZIPs alone would not create anonymized peer-review access or satisfy GEB's stable DOI/licence gate. Do not advertise complete reproducibility before that further work.

## Analysis stopping rule

No search for new post hoc positive rainfall hues, winner geographic sectors or image thresholds to rescue the global BIO12 sign. Existing source panels/bootstraps, 352+288 independent WLS point checks, heldout losses and unmatched-photo missingness are frozen. The next worthwhile scientific input is **source-independent standardized floral-colour/site measurement**, not further favorable p-value selection in the original iNaturalist archive.


## Verified original source point estimates and archive rights (14 October 2026)

A genuinely independent source-table reconstruction of the 20,694-species/1,934-genus world climate model now verifies **8 of 8** original four-hue×rainfall/vapour **point coefficients** (maximum discrepancy **3.51×10⁻¹⁴ pp**). Exact ZIP SHA256 and membership checks also passed, and the aggregated source licence/status census was recorded without photo IDs or coordinates. This closes the original source ZIP byte-identity and two principal environmental-axis point-estimate checks, **not** independent reruns of the original 999-bootstrap intervals, digital-image classification, or complete 20-cohort analysis. Source: `docs/MOIS_GEB_ACTUAL_SOURCE_RECOVERY_AUDIT_20261014.md` and `results/mois_original_source_independent_recovery_audit_20261014.json`.

**Current remaining hard HOLD:** license/attribution and source geoprivacy check for any released derivatives, reviewer-anonymous stable original analysis input/code access, permanent DOI, authorship and journal human approvals. Downloading the originally measured tables into a temporary working environment does not fulfil data-access requirements.
