# MOIS: finish all FCP wave4 real photo pixels, preserve original source dates, extend global seasonal observations

**17 October 2026.** Scientific scope: among-named-taxon worldwide photographic floral-hue composition, not FCP's separate within-species D/W/50-km colour-polymorphism results.

## Source authenticity — the full 49,800-photo image study actually completed

The previously reconciled MOIS source contains **302,242 genuine original photo-ID pixel outcomes /146,061 automatic four-hue ROI labels /146,052 genuinely observation-dated four-hue labels**. Source: [real 302,242-photo + 4,096 original UUID-day + 1,732 numeric observation/day run](https://github.com/zuizui0223/mois/actions/runs/39020466312), with complete day-enriched photo gzip SHA256 `9c9692238f58cc61561270f7c07b1d782ce711adeade9c58550866f668fb4690`. Original historic 99,817 public source dates are retained with September 2026 S3 snapshot precedence and no new photo-specific botanist stage labels.

The independent [original FCP wave4 real pixel run](https://github.com/zuizui0223/fcp/actions/runs/38992543115) completed successfully with **all 256 planned worker partitions**; its single provider-signed final `fcp-wave4-149900-method-expansion-actual-measured-and-inferred` artifact (ID `12097805951`) contains all **49,800 genuine measured original photo pixel ROI outcomes**, of which **24,585** have an automatic four-hue class. The signed original photo/measurement source CSV SHA256 is `fc642d42b861bc79305dbd7e7e57d6b14561028e5552bc39dece93e59442f766` and the source metadata SHA256 is `8d98b14856dd2539b54e1db287aa08adea9c191e112223ff380cb49ac0f8a554`. This release was inspected and has exactly 498 distinct original FCP source taxa and all terminal validity records.

The two scientific data archives overlap in **42,603 exact original photo IDs** (not independent images). Cross-source comparison of those **42,603 duplicated photographs** found **zero disagreement** in photo hue, original taxon ID, source numeric observation ID, original latitude or longitude. For **960** shared photographs, the older MOIS raw `observed_on` field is blank, but its independently recovered `qualified_date` is **identical** to the full FCP original `observed_on` day: preserve the original raw missingness and qualified date, never overwrite the original source row or misclassify this as disagreement.

## Actually measured photo additions to the shared MOIS source

| Literal photo-ID source measure | Previously signed MOIS 302,242 photos | Finished FCP256 union (source-identical originals) | Increment |
|---|---:|---:|---:|
| **Distinct genuinely measured original photograph IDs** | 302,242 | **309,439** | **+7,197** |
| **Classifiable four-hue photo ROI outcomes** | 146,061 | **149,662** | **+3,601** |
| **Valid observation date AND four-hue ROI** | 146,052 | **149,653** | **+3,601** |
| Source photographs with observation date (including hue-unknown) | 302,223 | **309,420** | +7,197 |
| Historic/other existing original source photographs with no day | 19 | **19** | No invented date |

The 7,197 new source photo pixels are allocated as **1,786 automatic white labels, 921 yellow/orange, 493 blue/purple, 401 red/pink, and 3,596 technically unclassifiable**. All newly added 7,197 photos carry an actual source recorded date; unclassified technical hue is not interpreted as a non-white reproductive phenotype or nonflowering status. These 7,197 photographs are from taxa already in the source; **no additional original named taxon is established by appending these pixels**, though photographic depth within taxa increases.

## Spatially matched year/season source support — after source enlargement

For originally recorded dates, extra-tropical hemisphere-adjusted summer versus winter, using only automatically four-hue labelled source photos:

| Source level | Local summer white/classified | Local winter white/classified |
|---|---:|---:|
| North extra-tropics | **22,116/44,979 = 49.17%** | **2,214/4,779 = 46.33%** |
| South extra-tropics | **3,461/6,832 = 50.66%** | **1,880/4,760 = 39.50%** |

These are **unpaired** photo proportions and subject to different photographed taxa and observer-effort/classifier-selection distributions. They are not a true field census of open flowers and do not establish summer pigment changes.

More informative descriptive photo comparisons require the *same 0.25° coarse grid cell and same local meteorological season-year*, at least **five independent source observation IDs** and **two nominal taxa** with classifiable hue in summer AND winter. This admits **27 Northern /28 Southern** source cell×year groups, containing **394 /552 automatically classified photo ROIs** respectively. The equal-cell-year observed summer minus winter source-white fraction is approximately **+6.70 Northern /+4.49 Southern percentage points**. Those small and nonrandom selected groups need geographically clustered exploratory uncertainty; both labels and per-geo-cell comparisons are source photographs, not independent individual flowering transitions or actual observed-year water/light exposures. Coarsening to 1–2° admits more groups by merging distinct geographic populations, not by resolving the needed local counterfactual.

## Reproducibility and claim limits

The new original [full-FCP256 to 309,439 measured photo script](../scripts/analysis/mois_full256_fcp_original_photo_and_date_merge.py) reuses the already committed strict source-reader/ROI terminal verifier `mois_import_full_FCP_wave4_completed_original_photos.py` and the existing seasonal calendar/local community methodology. It refuses changed source SHA, inconsistent original photo/hue/taxon/observation/coordinates, a stale or incomplete 256-worker terminal receipt, photos lacking an original observed-on date, and overlapping photo IDs masquerading as new observations.

Its [new GitHub Actions source execution](../.github/workflows/mois-full256-FCP-original309439-dated-hues-source.yml) uses only already completed, signed original provider artifacts; no new images need to be downloaded or reclassified. Once successful it archives the full 309,439 original-photo/qualified-day census separately from the manuscript, together with numeric source provenance/season and locality/year checks. All former MOIS 100,543/248,312/267,425/302,242 cohorts remain reproducible as *historical source stages*, not mixed replicates.

**The primary MOIS annual six-environmental-axis / four-hue / genus-999 maxT moisture estimator has not been rerun on this exact 309,439-photo source in this merge.** Coefficients from earlier 267,425, 286,826 or 295,017 photos must not be copied as fresh 309,439-photo results. Missing photo-specific open floral petal/tepal labels, truly observed-year/month meteorology and unmeasured observer effort prohibit direct claims of seasonal pigment plasticity or sun/rain selective costs.
