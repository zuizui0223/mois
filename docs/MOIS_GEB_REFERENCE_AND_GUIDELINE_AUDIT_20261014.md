# GEB reference and article-format audit — 14 October 2026

**Scientific document control:** this is a journal-format and citation audit, **not** an additional biological model, independent source replication or proof of acceptance. Author-facing paper: `docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md`.

## Verified external journal submission rules

Official publisher's [Global Ecology and Biogeography author guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html), checked 14 October 2026:

| Requirement | Working Mois status | Evidence and next operation |
|---|---|---|
| Research Article type | **Appropriate provisional type**, not accepted | ~5,000 main text words, ~50 scholarly references and **6–8 display pieces** are typical, not strict absolute caps |
| Structured abstract | **PASS** | 7 required fields: Aim, Location, Time period, Major taxa studied, Methods, Results, Main conclusions; limit **300 words** |
| Running title | **Updated** | Fewer than **40 characters**; working `Floral colour and moisture patterns` is 35 characters |
| Keywords | **PASS** | Eight alphabetically arranged: citizen science, climate, floral colour, geographic sampling, macroecology, measurement bias, precipitation, trait variation |
| Main manuscript | **Draft** | Editable article text, continuous line numbering required before Word/PDF upload; running title, main text, references, Data and Code Availability, tables and figures |
| Tables + figures | **4 source-backed main figures + 2 tables** | Six display pieces, within the journal's typical range; individual SVGs require resolution/accessibility/legend review |
| Double-anonymous peer review | **Partially prepared** | A separate blinded editable manuscript and title page needed; public source-commit identifiers, FCP/this repository links and document metadata may reveal identity |
| Identifying title page | **HOLD—human author input** | Names, affiliation, correspondence email, ORCID, acknowledgements, funder/COI declarations, contributor verification, signed author consent |
| Cover letter | **Draft** | Required as a **separate PDF**; include one paragraph of under 250 words explaining relevance to globally important macroecology; coauthor sign-off and submission declarations are human tasks |
| Data and code availability for peer review | **HOLD** | Both data and code should be available for independent scrutiny during peer review, and archived in a **stable repository**. **GitHub is not treated by GEB as an adequate long-term archiving destination** |
| Anonymous reviewer link | **HOLD** | Example: reviewer-private Dryad link or another stable source-independent anonymous repository; DOI/deposit and permissions unconfirmed, never invent links |
| Raw original photo+environment input | **HOLD** | Source photo SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`; environment SHA256 `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca` verified in original FCP source replay, but ZIPs **are not present in the new mois repository** |

**Editorial positioning:** GEB evaluates general conceptual insights at broad spatial/taxonomic scales. A previously documented white–precipitation association alone is **not an adequate novelty claim**. The strongest proposed contribution is the explicit **estimand/measurement dependence** of colour–moisture relationships across genera, source localities, strict photographic conditions and between/within-species comparisons. This interpretation is bounded; it does **not** establish an ecological causal process.

## Publisher-verified scientific references currently in the text

The author-facing article's bibliographic records were verified from publisher/conference pages; these links document identity of the citation, **not evidence that the cited publication independently supports every biological inference in the current paper**.

| Article reference | Verified DOI and authority | Contextual use |
|---|---|---|
| Dalrymple et al. (2020), *New Phytologist* 228, 1972–1985 | [10.1111/nph.16737](https://doi.org/10.1111/nph.16737) | Prior flower-colour macroecology including moisture, biotic communities; do not claim novelty for a simple moisture association |
| Di Cecco et al. (2021), *BioScience* 71, 1179–1188 | [10.1093/biosci/biab093](https://doi.org/10.1093/biosci/biab093) | Citizen-science sampling/observer limitations |
| Fick & Hijmans (2017), *International Journal of Climatology* 37, 4302–4315 | [10.1002/joc.5086](https://doi.org/10.1002/joc.5086) | WorldClim climate surface methodology; actual environmental archive remains original source authority |
| Koski & Galloway (2020), *Frontiers in Plant Science* 11, 991 | [10.3389/fpls.2020.00991](https://doi.org/10.3389/fpls.2020.00991) | Within-species floral optical/thermal/geographic history precedent; no numerical FCP result substituted |
| Poggio et al. (2021), *SOIL* 7, 217–240 | [10.5194/soil-7-217-2021](https://doi.org/10.5194/soil-7-217-2021) | SoilGrids methods and spatial uncertainty |
| Trethowan et al. (2024), *Journal of Biogeography* 51, 1659–1669 | [10.1111/jbi.14751](https://doi.org/10.1111/jbi.14751) | Prior macrogeographic flower colourfulness, climate and island traits, limiting novelty |
| Xiong et al. (2024), CVPR proceedings, 16111–16121 | [10.1109/CVPR52733.2024.01525](https://doi.org/10.1109/CVPR52733.2024.01525) | EfficientSAM flower segmentation software background, **not** flower-specific accuracy validation |

## Remaining source/reference completeness gates

- **iNaturalist snapshot / image rights:** the exact observation export identifier, download citation if available, licences for original photographs, and republication permissions require source-specific verification. We do **not** fabricate an iNaturalist collection DOI or assume images are universally redistributable.
- **Custom ROI-v4 and nine-palette measurement:** the versioned source code and parameter/weight provenance need to be separately deposited; EfficientSAM's CVPR paper is not documentation for custom flower ROI training or calibration.
- **WorldClim and SoilGrids data-layer provenance:** authoritative original FCP environmental ZIP and join script must identify the exact layer versions, dates/depths, derived available-water proxy and preprocessing. General WorldClim/SoilGrids citations do not by themselves prove the specific source layers or invalidate the original checksums.
- **Date support:** original source's retained analysis-ready photo files lack a robust original observation-date column for seasonal causal claims, hence the manuscript's Time period statement is deliberately nonquantitative.
- **APA-6 and spelling:** keep names, DOI, journal volume, article-number/pages accurate; spellings may follow British English. The journal permits Free Format first submission, but internal citation formatting should remain consistent and the final review copy must have numbered lines.
- **Use of publicly available original FCP GitHub material during anonymized review:** a publicly identifiable GitHub repository can itself compromise author anonymity. During review, the blinded manuscript should use a private reviewer data-access URL not traceable to authors; the editor-facing cover letter/title page can retain source identification. **A source provenance hash isn't a substitute for accessible reviewer files.**

## Exact original ZIP archive access versus peer-review data availability

Both original Actions archives (photographic colour/site tables **11587115220** and abiotic/site tables **11653183273**) were actually recovered into the 14 October 2026 working environment, with source SHA256/ZIP integrity **PASS**. The original photo archive has a historical Actions expiry date of 2027-01-06; do not assume continued access. An independently written source join and model reconstructor agrees with the original eight principal two-axis/four-hue point coefficients (max error **3.51×10⁻¹⁴ pp**). Its aggregate-only receipt contains licence counts and source coordinate-status counts, not original images or exact source IDs. These two source ZIPs are **NOT deposited in mois or in a stable reviewer repository**. The recovery shell helper remains available but was not used for the connector download. GEB submission stays on HOLD until authorized image/site licensing, anonymous reviewer archive/DOI and human author approvals are fulfilled.

## Stop/go

**GO for internal manuscript review:** complete observed original-source result tables, frozen 1,280/576 coef ledgers, four numerical figure sources, and CI checks. **HOLD for actual GEB upload:** signed title page, reviewer-anonymous editable manuscript/metadata and reviewer data path, licence-compliant raw source archive with stable DOI, final display/layout, author-coauthor sign-off and PDF cover letter.

**No scientific p-value rescue:** additional nonlinear or split analysis must not be used to select a stronger finding from the same source after seeing outcome. This paper remains a *bounded observational inference* about generality, not demonstrated pigmentation selection or white-versus-coloured genetic morph fitness.


## Reproducibility update: source-verbatim 999-genus family rebootstrap

The original two SHA256-locked tabular source ZIPs were used for a separate independent reconstruction of the frozen world-climate26, genus-equal **24-contrast maxT family** with **999 original-seed genus resamples**. All 24 coefficient estimates, 24 percentile CIs, 24 maxT 95% intervals and support flags match original source to maximum **3.55×10⁻¹⁵ pp**; the critical value is 3.0109108050453743. See `docs/MOIS_GEB_WORLD24_GENUS_SUFFSTATS_REPLAY_20261014.md`. This is a verification of the original scoped statistical calculation, **not** a new biological sample or correction for historical post-outcome model and region selection. The internal anonymous-genus sufficient-moment intermediate has not been deposited in the public repo; source rights/site covariance must be reviewed before any external handoff or publication.


## Independent original photographic matched-geography reconstruction

The source-original taxon-matched Stage A/B/C macroenvironment model was independently recalculated with 999 original paired whole-genus draws for the **world and East region under both weighting schemes**. All **576 points, 1,152 ordinary interval endpoints and 576 paired-difference maxT interval endpoints** match the original stored source table to less than 1.60×10⁻¹⁴ pp. The original 24-maxT within-scope intervals for the world and East white/BIO12 broad−local difference include zero. This is a source-data statistical/implementation audit, not independent region-morph measurements, a study-wide multiple-testing correction or a causal result. Its public aggregate proof: `docs/MOIS_GEB_INDEPENDENT_SPATIAL576_SOURCE_REPLAY_20261014.md` / `results/mois_original_matched_spatial_full576_replay_20261014.json`. Exact original source ZIPs and anonymized reviewer dataset are still *not* published.
