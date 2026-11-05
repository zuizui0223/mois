# Diel same-year macro implementation plan

Goal: complete the design in ../specs/2026-11-05-diel-sameyear-macro-design.md using actual additional data.
Architecture: frozen public study inputs -> exact species bridge -> outcome-blind source support -> cached same-year exposures -> macro model and verification.
Tech stack: existing Python/pandas/statsmodels/requests; R for faithful public RData export and phylogenetic sensitivity if feasible.

## Tasks and execution ledger

- [x] Download versioned experimental data and MD5-check original files; add SHA256 receipts.
- [x] Read ALAN and Hui–Walter source code and record use boundaries.
- [x] Build tested exact-species bridge and source-support selection, with no ambiguous-name joins or ROI failure exclusions for literature-colour outcome.
- [x] Freeze estimand/support rules before examining exposure associations.
- [ ] Acquire new VIIRS and WorldPop records with explicit same-year/catalog checks, checkpointing and provenance.
- [ ] Fit macro model; audit rank, geography, independent taxa, uncertainty and source coverage.
- [ ] Test, independently review numerical results and package source-linked deliverables.

Review focus: alias conflicts; repeated species pseudo-replication; unrecorded vs absent; satellite/cloud/time mismatch; transported experimental species traits vs local pollination.

Ruling: user has explicitly instructed execution and data expansion; continue with the existing Python analysis stack and suitable current checkout. Do not stop for a further permission cycle for ordinary analysis choices. Method choices and cohort gates are investigator decisions, not a retrospective preregistration.

## Production evidence ledger (2026-11-05)

- Full source run 41982925538 succeeded: 53,652,241 observations scanned; 7,722 observed eligible species; 8,667,448 precise 2015–2024 observations. Broad primary support: 83,991 cell-years / 23,612 cells / 556 spatial blocks.
- Experimental source run 41982122333 succeeded. Broad-frame support gate retains 65 species (35 white / 30 nonwhite) and 124,157 species-cell-years; most species have one experimental study.
- Same-year coarse climate acquired and joined locally: 83,991/83,991 cell-years, twelve verified months per year. Native Gaussian grid is 94 latitudes x 192 longitudes, not 0.25-degree climate.
- Environment run 41983664064: population succeeded; VIIRS still in progress at this checkpoint. Full model run 41986239793 is waiting for verified acquisition success; no partial global model result substituted.
- Independent numerical review matched weighted absorbed estimates and one-/two-way clustered standard errors to explicit-dummy statsmodels WLS within 3e-15. Related 37 tests passed. Review fixes cover population key uniqueness, resumed raw-response caches, source/export SHA binding, and retention of successful primary results when a sensitivity is unsupported.
- Remaining: complete environmental response audit, actual fitted coefficients, sensitivity/coverage review, source-linked final scientific report and completion audit. Model confidence intervals treat literature traits and treatment ratios as given; they do not propagate all study-level uncertainty or establish a local pollination mechanism.
