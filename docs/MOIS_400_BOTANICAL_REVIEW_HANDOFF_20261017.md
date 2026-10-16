# MOIS independent botanical open-flower review handoff — 400 source photos (17 October 2026)

**Status: reviewer sheets prepared; 0 photos have been independently botanically reviewed.** This is a supporting input for MOIS global seasonal photographed flower-colour research, not the FCP New Phytologist within-species colour-diversity D/W paper.

## Source fixed *before* opening flower colours

The SHA-pinned iNaturalist source-only analysis in [MOIS PR #10](https://github.com/zuizui0223/mois/pull/10) retained 400 original licensed public photo URLs from five species × one matching 0.25° cell/site-year × two local seasons × 40 photos. There are 400 unique photo IDs, 400 unique original observation IDs and 400 distinct photographer IDs. Historic MOIS 100,543 photo IDs and the 64 previously evaluated Asclepias photos were excluded using IDs only.

The original selected photo manifest has SHA256 `b961b50ee6efe9ec08a5575a8eb1e34b1edf7b7fcb426226fb11ee901cfd3329` in the [source replay Actions run 38929402279](https://github.com/zuizui0223/mois/actions/runs/38929402279) (artifact `12058501017`).

## Actual two independent reviewer packets

`scripts/analysis/mois_prepare_400_botanical_review_handoff.py` derives two different reviewer-order shuffled CSVs and a coordinator-only identifier mapping. **No seasons, dates, coordinates, climate variables, previous automatic hue classification, or other reviewer's answers** appear in either reviewer sheet.

The 400 items have an original public photograph link, source photo licence/attribution, species name and a stable blind token. Each reviewer independently fills `reviewer_id`, `open_flower_stage`, `roi_organ` and notes.

Allowed stage codes: `open_flower`, `buds_only`, `fruit_only`, `no_flower_visible`, `uncertain`. Allowed tissue codes: `petal_tepal`, `other_flower_organ`, `spikelet_glume_bract`, `fruit_seed`, `leaf_stem`, `background`, `uncertain`.

Only a **visible open flower with verified petal/tepal ROI** can enter the eventual primary flowering-colour rate denominator, and only once its independent review agrees or has third-person adjudication. Unknown and missing labels are **unknown**, not non-flowering or white/nonwhite. Classifiability and mixed flower-organ ROI remain separate post-review technical gates.

- Reviewer A unfilled CSV SHA256: `576a55e7e613faf8c93652b4ac2d72b79f8a38bec4bfe8267c8d67bbc3fed322`.
- Reviewer B unfilled CSV SHA256: `bfbbfc92106585ca260fbd5a33573ad2598209d19907a52bbce2b941b1164203`.
- Coordinator-only original photo-to-blind-token CSV SHA256: `6119a9606b4d90eb80a8a90e7fd2576631d14063de320e788e7c274fd050a3c4`.

The independent forms are individually provided as separate time-limited GitHub Actions artifacts to reduce accidental response leakage. The coordinator mapping remains out of public Git and should not be shared with reviewers. Public source photograph URLs may permit an intentional external metadata lookup, so blinding is procedural rather than cryptographically guaranteed.

**No actual biological effects have yet been computed from these 400 URLs.** The separate photos are newly reserved source opportunities, not freshly measured RGB floral hues and certainly not genotype-specific or fitness-measured flower-colour polymorphisms. Nor do five mostly northern temperate taxa provide global generality.

The next stage is to collect two independent botanical decisions for each pictured original photo; disagreements require adjudication. Only then perform frozen ROI-v4 flower hue measurement conditioned on verified flower organs, retain all failures in denominators, join actual observed-month climate, and test ecological seasonal patterns across more geographically and phylogenetically varied replicated sites.

The original between-species global 100,543-photo MOIS BIO12 and vapour coefficients remain a separate source-specific observational result, with strict-photo and geographic negative gates unchanged. FCP's 149,900-photo species-conditioned within-species D/W and 50km geographic structuring remain independent estimands.
