# Main-figure rendered visual check — global moisture GEB draft (14 October 2026)

**Scope:** editorial readability and graphic integrity of the four existing, source-backed SVG figures; **not a new biological reanalysis**. Original numerical point estimates, normalizations, source confidence intervals and plotted marks remain unchanged.

## Actual rendering checks

The 2026-10-14 study branch's four SVGs were extracted from the successful source-locked Actions numerical-snapshot artifact, rasterized with CairoSVG at **1440-pixel output width**, and examined together and individually as rendered images. All four are readable at full desktop resolution, with intact headings, axes, intervals, directions and qualifiers. Source SVG generators/CI were also checked after edits. **Print proofs, colour-blindness audits, final font substitution and single-column 85-mm downsizing have not been completed**, so no journal-resolution approval is claimed.

| Main figure | Data and source | Rendered PNG check | Print-stage caveat |
|---|---|---|---|
| Figure 1 | Original broad BIO12 slope, strict-quality and matched geography, botanical descriptions; original FCP figure source | 1440 ×905 pixels; correct four original panels, zero lines and interval whiskers | Published original labels need checking at proposed print size |
| Figure 2 | Source photo classifiability/wet–dry composition and same 21,448-original-photo within-/between-species rain coefficients | 1440 ×800; three panels; long label and numerical-value text backdrops added to prevent vertical reference lines appearing through words | Rightmost panel is dense; retain full-width/two-column presentation |
| Figure 3 | Distinct outcome grains, unseen genus/spatial exclusion and nonlinear-versus-ridge source gain | 1440 ×787; three panels; masked label/value text against reference grid without moving interval points, error bars, baselines or scales | Dense interval labels; test two-column layout with actual final font embedding |
| Figure 4 | Original 1,280-row frozen CSV: four hue coefficients on BIO12/vapour and original within-family maxT 95% intervals | 1440 ×844; two panels; same source ranges/CI and clear composition/noncausality footer | Test colours and grey-white contrast under print/colour-vision simulation |

### Applied non-scientific corrections

The Figure 2 source-backed renderer `scripts/render_mois_figure2_source_receipts.mjs` now draws **three header and three numerical-text background masks** behind right-panel row labels to keep the zero/reference/grid strokes from visually crossing those labels. Figure 3's `scripts/render_mois_figure3_transfer.mjs` similarly draws **nine row-label and nine coefficient-text background masks**. Their SVG outputs are regenerated from the **unchanged original JSON scientific results**. All data plot marks, coordinate/CI point values, x-axis limits, axis semantics, panel positions and source titles remain unchanged.

The CI test `tests/test_mois_geb_submission_integrity_20261014.py` guards exact mask counts. Source regen checks `node scripts/render_mois_figure2_source_receipts.mjs --check` and `node scripts/render_mois_figure3_transfer.mjs --check` verify exact source-consistent SVG output, alongside Figure 4 and the original manuscript/claim gates. Successful Actions run: https://github.com/zuizui0223/mois/actions/runs/38570933512.

**Technical PASS / journal display HOLD:** main figures have been successfully rendered and visually reviewed at full raster width, but final GEB double-anonymous Word/PDF export must still be examined at **actual journal column width** and with complete caption/supplement metadata, author sign-off and data/photograph rights. No original photos were publicly redistributed by these four numerical plots.
