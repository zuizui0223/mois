# Main-figure rendered visual check — global moisture GEB draft (updated 15 October 2026)

**Scope:** editorial readability and graphic integrity of the five existing, source-backed SVG figures; **not a new biological reanalysis**. Original numerical point estimates, normalizations, source confidence intervals and plotted marks remain unchanged.

## Actual rendering checks

The four original 2026-10-14 study SVGs were extracted from the successful source-locked Actions numerical-snapshot artifact, rasterized with CairoSVG at **1440-pixel output width**, and examined together and individually as rendered images. All four are readable at full desktop resolution, with intact headings, axes, intervals, directions and qualifiers. Source SVG generators/CI were also checked after edits. **Print proofs, colour-blindness audits, final font substitution and single-column 85-mm downsizing have not been completed**, so no journal-resolution approval is claimed.

| Main figure | Data and source | Rendered PNG check | Print-stage caveat |
|---|---|---|---|
| Figure 1 | Original broad BIO12 slope, strict-quality and matched geography, botanical descriptions; original FCP figure source | 1440 ×905 pixels; correct four original panels, zero lines and interval whiskers | Published original labels need checking at proposed print size |
| Figure 2 | Source photo classifiability/wet–dry composition and same 21,448-original-photo within-/between-species rain coefficients | 1440 ×800; three panels; long label and numerical-value text backdrops added to prevent vertical reference lines appearing through words | Rightmost panel is dense; retain full-width/two-column presentation |
| Figure 3 | Distinct outcome grains, unseen genus/spatial exclusion and nonlinear-versus-ridge source gain | 1440 ×787; three panels; masked label/value text against reference grid without moving interval points, error bars, baselines or scales | Dense interval labels; test two-column layout with actual final font embedding |
| Figure 4 | Original 1,280-row frozen CSV: four hue coefficients on BIO12/vapour and original within-family maxT 95% intervals | 1440 ×844; two panels; same source ranges/CI and clear composition/noncausality footer | Test colours and grey-white contrast under print/colour-vision simulation |
| **Figure 5** | Three independent **aggregate receipts**: photo-/species-/genus-equal white fractions, exact matched within-/between-genus moisture bootstrap summaries and positive/negative/invariant genus signed contributions | **New 1560 ×870 and 780 ×435 PNGs rasterized from the exact source-reproducible SVG; visually inspected 15 Oct.** White/other percentages, CI whiskers, 830/719/385 source count labels, 3 panel cards and limits are visible, with no clipping at full width | **Use full two-column width.** The 780-px half-scale preview has dense fine print and is unsuitable for single-column use; final actual 175–180 mm printer proof and colour-vision accessibility test remain HOLD |

### Applied non-scientific corrections

The Figure 2 source-backed renderer `scripts/render_mois_figure2_source_receipts.mjs` now draws **three header and three numerical-text background masks** behind right-panel row labels to keep the zero/reference/grid strokes from visually crossing those labels. Figure 3's `scripts/render_mois_figure3_transfer.mjs` similarly draws **nine row-label and nine coefficient-text background masks**. Their SVG outputs are regenerated from the **unchanged original JSON scientific results**. All data plot marks, coordinate/CI point values, x-axis limits, axis semantics, panel positions and source titles remain unchanged.

The CI test `tests/test_mois_geb_submission_integrity_20261014.py` guards exact mask counts. Source regen checks `node scripts/render_mois_figure2_source_receipts.mjs --check` and `node scripts/render_mois_figure3_transfer.mjs --check` verify exact source-consistent SVG output, alongside Figure 4 and the original manuscript/claim gates. Successful Actions run: https://github.com/zuizui0223/mois/actions/runs/38570933512.

**Technical PASS / journal display HOLD:** main figures have been successfully rendered and visually reviewed at full raster width, but final GEB double-anonymous Word/PDF export must still be examined at **actual journal column width** and with complete caption/supplement metadata, author sign-off and data/photograph rights. No original photos were publicly redistributed by these four numerical plots.


## Figure 5 — directly rechecked source-integrative ecology artwork (15 October 2026)

**Verified source and scope.** `node scripts/render_mois_figure5_ecological_scales.mjs --check` regenerates the exact committed `figures/mois_white_moisture_ecological_scales.svg` using only three source-verified numerical receipts:

- `results/mois_original_photo_species_genus_balance_20261014.json`: white photo-equal 46.95% (45,953 classifiable photos), species-equal 46.82% (25,064 sampled taxa), genus-equal 49.89% (5,369 classifiable genera), and ≥5-classified-species sampled-genera white 46.63% (1,185 genera)
- `results/mois_genus_between_within_moisture_20261015.json`: **the identical 20,694 species /1,934 genera** within-versus-between genus photographed-white/BIO12 and atmospheric-vapour slopes and their **ordinary paired bootstrap** 95% uncertainty; no replacement of original 24-comparison maxT limits
- `results/mois_congeneric_moisture_genus_sign_mosaic_20261015.json`: **830 positive /719 negative /385 photo-white-invariant** nominal genus contributions for BIO12; **841/708/385** for vapour

**Rendered inspection.** GitHub [run 38672349983](https://github.com/zuizui0223/mois/actions/runs/38672349983) ran all three CI jobs successfully and saved the temporary `mois-figure5-ecological-scales-visual-proof` artifact, containing **1560×870** and **780×435** PNGs rasterized from Figure 5 using CairoSVG 2.8.2. Both PNGs were directly opened and visually checked in an execution environment: full-width axis/whisker positions are clear, panels/colour legend remain separated, sign counts appear in distinct stacked segments, title/footnote are untruncated. The 780-pixel preview is **too small for comfortable unzoomed reading of subtitle/caption-sized labels**, so Figure 5 should be a **full-width (two-column) journal figure**, not reduced to a single column.

**Scope and display disclaimers.** Figure 5 contains **no original photos, named genus taxa, source sites or new biological inference**; its A-panel observation opportunities differ from the B/C climate26 source support. Percentage bars are source-classifier-conditioned photo fractions, not genetic morph frequencies; B is ordinary within-design genus bootstrap rather than original FWER control; C's per-genus signs are partial covariance contributions, not fitted genus-specific fitness/selection coefficients. Readability for the actual chosen journal printer resolution, accessible greyscale/colour-blindness simulation and final author approval **remain outstanding**. These cannot be certified from the two screen previews alone.

**Display inventory:** 5 SVG figures + 2 main tables = **7 display items** in the current author-facing GEB draft, within its provisional 6–8 figure/table target. This is editorial formatting information, not a submission.
