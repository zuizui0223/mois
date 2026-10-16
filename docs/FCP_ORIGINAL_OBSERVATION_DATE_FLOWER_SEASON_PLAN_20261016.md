# Original observation dates, flowering images and seasonal white-frequency analysis

Updated: 2026-10-16. Scope: within-species FCP, not the separate interspecific Mois paper.

## Question

Within a species, does the photographed frequency of white versus chromatic
flowers vary by observation month, and does within-site temporal variation
track radiation, temperature and water availability after geography and
time-stable environmental differences are controlled?

Original 2026-09-27 iNaturalist metadata: 42,111 candidate taxa,
93,209,351 licensed photo opportunities. These are NOT 93,209,351 flowering
images. The first balanced RGFCA 1,000-taxon photo selection is a screening
subset, not a biological all-photo ceiling. No old 100,543/149,900 flower
classification may replace fresh blind ROI-v4 measurement.

## Current executable, tested descriptive processing

Program: scripts/analysis/fcp_original_dated_flowering_season_qc.py.

Input A: SHA-sealed complete original photo-ID source with stable taxon,
original observation ID, observer, observed_on date and coordinates.
Input B: independently sealed newly measured original flower ROI-v4 terminal,
matched to original source only AFTER the pixel-level blind classification.
Output: original_dated_month_site_white_qc.csv plus a SHA-sealed JSON receipt.

For each species, original year-month and 5-degree source geocell report:
- Original dated candidate photos versus newly attempted real ROI photos.
- How many have classifiable flower colour, unclassifiable outcomes, and
  within-same-flower white-ray/yellow-disc colour mixture.
- Count independent original observations, not repeated photos from the
  same iNaturalist observation.
- Set a primary white indicator only when all usable pictures from one
  observation have one consistent coarse colour state and no suspect
  white+yellow mixed flower organ outcome.
- Primary descriptive fraction = white independent observations /
  independent observations with consistently white or chromatic flower ROI.
- Keep unknown-date, ambiguous and unclassifiable photos visible in missingness
  QC; never reinterpret ROI failure as proof of a nonflower.
- Southern/Northern local seasons differ; tropical sites retain month rather
  than invent a locally validated dry/wet or temperate season definition.

Six synthetic tests passed:
https://github.com/zuizui0223/fcp/actions/runs/38883506069
They verify original day (not upload day), observer-level denominator,
multi-photo contradictions, within-ROI organ mixing, hemisphere handling,
tampered source identity, changed ROI IDs and no synthetic bypass in the
production command-line interface.

**No real monthly original-source flower-colour effect has been estimated.**

## Statistical analysis after genuine full-source ROI classification

1. Retain all original candidates, record dated photo opportunities,
   measured original images, classifiable flower images and source observations
   in a transparent funnel, per taxon/site/year/month/observer.
2. First compare across months WITHIN species and geo-site-year strata where
   independent flower observations exist in more than one month. Monthly
   morph frequency is not a flower switching its colour in one individual.
3. Join time-resolved observed-month climate only after flower classification:
   air temperature, solar radiation, rainfall and moisture or VPD. Separate
   site-month deviations from site climate means. Do not assign normal-climate
   data to an observation as though it were the actual month's weather.
4. Time-stable soil, elevation, slope/aspect, geography and spatial structure
   adjust between-site contrasts. Observer, year and image/ROI classifiability
   require separate bias checks. Spatial smooths may absorb broad-scale
   climate gradients; compare constrained vs unconstrained estimates.
5. Model only when site-month replication and positivity are adequate.
   Observation-level Bernoulli outcome, species effects, cyclic day-of-year,
   within-site monthly climate anomalies, site/spatial and observer effects,
   year effects, and stable soil/topographic covariates. Partial pooling
   across species; clustered uncertainty by observation/site/observer.
6. Robustness: observer-disjoint folds, year-matched regional subsets,
   alternative geocell widths, high-confidence ROI only, white/yellow disc
   exclusions, date missingness and season-specific ROI detection rates.
   No adaptation, pigment cost, genotype-frequency or causal claim without
   additional evidence.

## Exact source blocker and reproducibility

The exact original-monthly 1,000-taxon ID job
https://github.com/zuizui0223/fcp/actions/runs/38847189688
FAILED: photo counts did not match independently authenticated full-frame
counts for 357 taxa. A new artifact-only audit, without another 36GB download,
https://github.com/zuizui0223/fcp/actions/runs/38882600217
found all 30 stored example discrepancies were negative: targeted had fewer
photo IDs than original counts, median absolute discrepancy 3.5 photos.

The failed production job checked out historical extractor commit
e2a23434ec171e3f2b7474c4972ab5cd91e31d36. Its second photo-ID
collision pass deliberately rejected duplicate photo IDs from raw rows EVEN
when the outside observation was not source-eligible. That is a real policy
difference from the full 42,111-taxon importer, which counts photo IDs
from qualified research-grade, georeferenced in-frame observations only.
This is a plausible source of the 357 taxon deficits.

The current PR branch extractor independently checks outside observations
against the same full-frame source criterion before quarantining ID
collisions. Source-identity parity regression
https://github.com/zuizui0223/fcp/actions/runs/38883913901
passed 16 synthetic tests; this does NOT prove full-real 357-taxonomy
reconciliation. The 3,430,502 original photo IDs of the 1,000-taxon group,
and eventually all 42,111 photo IDs, must be authenticated and new pixels
classified before running a real season-climate-white analysis.

## 2026-10-16 extension: season-year matching and full original observed-on opportunity coverage

**Cell/year comparison, 12 synthetic checks PASS:** `scripts/analysis/fcp_match_original_season_site_year.py`
compares the same taxon in the same original 5-degree source geocell and
meteorological season-year only when both seasons have nonzero independent,
consistent flowering ROI observation labels. Northern Dec 2025 / Jan 2026
form a single winter 2026; Southern Dec/Jan form local summer 2026;
tropical months remain month-only. Match-weighted white fractions use
`min(n_observations_season_A, n_observations_season_B)` per shared
cell-year, so surplus images from unmatched locations do not influence the
matched contrast. Publish the unmatched pooled difference separately as an
explicit observational-geography sensitivity. A 5-degree cell is FAR too
coarse to identify one population, and the aggregate table cannot establish
observer-disjoint seasonal replication. Thus results are explicitly
descriptive with no significance or climate-causality claim. [CI 12/12 pass](
https://github.com/zuizui0223/fcp/actions/runs/38891959436).

**Every source photo, not a 100-image cap, has an observed-on census ready:**
`scripts/acquisition/audit_fcp_original_rgfca_monthly_photo_opportunities.py`
takes the authentic full 1,000-taxon 3,430,502-photo-ID signed original
manifest when and ONLY when the count/identity gate passes. It partitions
EVERY original photo into exact valid observed-on day vs missing, invalid or
future-to-snapshot date. Separately provides source-date opportunities,
distinct original observations and independent observers by taxon/month;
a 0.25-degree source-cell/month diagnostic requires stated positional
accuracy >0 and <=5 km. A source-cell-year is calendar-supported only if
there are >=2 months, each with >=2 separate observations and >=2 source
observers. These are not flowering photo counts or positive-colour results.
[CI 6/6 pass](https://github.com/zuizui0223/fcp/actions/runs/38892631080).

The source-success-only actual-data runner, installed on default branch,
is `.github/workflows/fcp-rgfca-all-original-observation-month-source-after-exact-20261016.yml`.
It accepts exclusively source [run #38884223530](
https://github.com/zuizui0223/fcp/actions/runs/38884223530),
exact source-commit `b0656b6a677a8440c188939074a0fe3c8cd4ac16`,
and original full source receipt SHA. It reuses the saved 3.43M-ID
artifact; **it neither downloads the entire 36GB provider tarball again
nor opens one flower image**, and refuses partial source results.
Its counts/metadata will be a genuine result ONLY when the source scan
succeeds and source-only census workflow completes. At this document
revision, real scan #38884223530 is still scanning and no biological
season-by-climate estimate exists.

The genuine first 16–64 original ROI pixel downstream workflow on
`main` now emits BOTH the original dated descriptive month QC and this
within-geocell meteorological season-overlap table. These early outputs
are a technical first-pair unit, not prevalence estimates for all
1,000/42,111 taxa. Month and hemisphere help exclude seasonal geographic
sampling confounding, not prove that any individual flower changes
colour during the season.
