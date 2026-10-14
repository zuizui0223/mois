# Internal numerical evidence bundle — non-publication staging only

**Purpose:** provide a deterministic, easily inspected **internal editorial evidence bundle** for the standalone `mois` interspecific floral colour–moisture paper. This is not the stable data archive required by *Global Ecology and Biogeography*, an anonymous reviewer-data link, or a new biological experiment.

## Exact contents and exclusions

The bundle is built by `scripts/build_mois_internal_evidence_package.py` and includes **18 previously public `mois` files** under their relative source paths: the 1,280-row world hue/environment coefficient CSV, 576-row matched spatial comparison CSV, both original decision/independent source-verification receipts, the independently rederived 8-point source audit receipt, figure source JSON, the author-facing manuscript/SI and source-verification report, original modelling/independent point-check code and the four numerical SVG figures. A machine `INTERNAL_SOURCE_MANIFEST.json` enumerates SHA256 and byte size of each packaged source file.

It excludes **both original source ZIPs**, any original iNaturalist image pixels or photo IDs, observation IDs, source taxon IDs, individual attribution, exact site latitude/longitude and per-photograph licence rows. It additionally excludes identifying title-page data. The two CSV headers are rejected if they unexpectedly include any of those source-identifying columns.

## Local and CI test

```bash
python scripts/build_mois_internal_evidence_package.py --write
python scripts/build_mois_internal_evidence_package.py --check
```

The source-content ZIP is assembled with fixed entry ordering, internal timestamps and a SHA256 file manifest, allowing an exact `--check` re-build after any source change. It is written under ignored `dist/`; the GitHub Actions `source-receipt-integrity` job produces a **14-day temporary GitHub Actions artifact** on PR checks. Those output files have the same public source contents as the repo, **not licensed raw photo/source records**. The artifact is useful for version traceability but is **NOT a permanent stable DOI or a reviewer-anonymous data deposit**.

## Still needed before journal submission

1. Acquire appropriate permission and coordinate/geoprivacy review for any *licensed row-level source derivatives* selected for anonymized scientific source reproduction. The original measurement ZIPs are recovered locally in the temporary analysis environment but have not been deposited.
2. Establish a stable, private reviewer-accessible data/code archive with non-revealing access URL, and a persistent public DOI after appropriate permission review; verify export rights per source/licence, especially the large CC BY-NC share.
3. Provide or legitimately restrict the underlying original source data if the genus-level model and bootstrap cannot be repeated from the numerical evidence bundle. **The 18-file public numerical bundle is not sufficient to refit all 20 original cohorts from scratch**, since individual source photo/site predictor rows were deliberately excluded.
4. Have coauthors approve all manuscript, source-access, acknowledgments, conflict and author-identifying records.

**Latest completed empirical check:** the independent 14 October 2026 source-table analysis reproduced the original 20,694 genus-eligible taxa, both original ZIP hashes and eight original global white/yellow-orange/red-pink/blue-purple BIO12/vapour point coefficients with maximum error 3.51×10⁻¹⁴ percentage points. The original 999-bootstrap intervals remain previously source-verified in FCP but were not newly independently recomputed by the recovery check.
