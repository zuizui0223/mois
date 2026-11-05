# Diel same-year macro implementation plan

Goal: complete the design in ../specs/2026-11-05-diel-sameyear-macro-design.md using actual additional data.
Architecture: frozen public study inputs -> exact species bridge -> outcome-blind source support -> cached same-year exposures -> macro model and verification.
Tech stack: existing Python/pandas/statsmodels/requests; R for faithful public RData export and phylogenetic sensitivity if feasible.

## Tasks and execution ledger

- [x] Download versioned experimental data and MD5-check original files; add SHA256 receipts.
- [x] Read ALAN and Hui–Walter source code and record use boundaries.
- [ ] Build tested exact-species bridge and source-support selection, with no ambiguous-name joins or ROI failure exclusions for literature-colour outcome.
- [ ] Freeze estimand/support rules before examining exposure associations.
- [ ] Acquire new VIIRS and WorldPop records with explicit same-year/catalog checks, checkpointing and provenance.
- [ ] Fit macro model; audit rank, geography, independent taxa, uncertainty and source coverage.
- [ ] Test, independently review numerical results and package source-linked deliverables.

Review focus: alias conflicts; repeated species pseudo-replication; unrecorded vs absent; satellite/cloud/time mismatch; transported experimental species traits vs local pollination.

Ruling: user has explicitly instructed execution and data expansion; continue with the existing Python analysis stack and suitable current checkout. Do not stop for a further permission cycle for ordinary analysis choices. Method choices and cohort gates are investigator decisions, not a retrospective preregistration.
