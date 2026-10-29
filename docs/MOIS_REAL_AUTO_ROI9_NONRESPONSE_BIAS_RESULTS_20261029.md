# MOIS automatic bias audit: observed results, not human-validated flowers (2026-10-29)

**Status:** Executed and passed all unit/source contracts on the frozen 27 Oct 2026 original photo v2; no human review required. [True GitHub Actions run 40856639794](https://github.com/zuizui0223/mois/actions/runs/40856639794). This file is a *photographic-label and missingness* receipt for PR #196; it does not replace the climate, altitude or island model estimates in the frozen GEB main paper.

## Frozen source boundaries

- Upstream unprocessed opportunity frame: 93,209,351 photo IDs. No fresh photo downloads or original ROI model inference.
- Exact existing four acquisition-arm union: **1,181,146** genuine original photo attempts; **451,643** original strict9 success.
- One outcome-blind hash-minimized original photo per original observation UUID: **1,032,882** distinct observation UUIDs; **148,264** redundant extra photos set aside from observation-level inference (not deleted from the original photo QC ledger).
- Canonical observation outcomes: **401,703** valid original strict9 colour labels; **631,179** original strict9 technical unknowns; measured original strict9 **white hue labels = 191,823**.
- Every original strict9 unknown remains technical unknown, even if an independent unapproved 500px photo-box classifier has another prediction. **0 human annotations, 0 new photo RGB downloads and 0 new ROI inferences**.
- Exact input corrected Parquet SHA256: \`faf6951505f1002a1327f64a2ced13431b562cf13893707125fafff918160f62\`; signed georestoration receipt SHA256: \`2893308a67732678b98a20ad38c510b67db50783519bcf66535ebb0fa5a66fcb\`.

## Actual observed data and assumption-conditional correction

| Measure | Result | Identification |
|---|---:|---|
| Original strict9 technical classifiability of canonical observations | **38.8915%** | Directly observed |
| White 4-hue photo label among original strict9 successes | **47.7524%** | Directly observed; not independently verified petals |
| Source-arm × broad latitude × local-season strata | **82** | Actual source denominators |
| Strata with fewer than 30 original strict9 successes | **14** | Warning on uncertain/local extrapolation |
| Strata with zero valid original strict9 | **3** containing **4 observations** | No MAR point extrapolated to unsupported cells |
| Approximate **conditional MAR** observed-photo-label rate under stratum exchangeability | **47.5469–47.5473%** | The range is *unidentified four-record support only*, NOT a 95% CI |
| Shift relative to raw successful-photo white rate | **−0.2055 to −0.2051 percentage points** | Assumption-dependent poststandardization, NOT causal bias removed |
| No-assumption all-unknown-photo-label mathematical limits | **18.5716–79.6802%** | Not a biological flora prevalence CI; every unknown label is unobserved |

The primary MAR *point* was intentionally left \`null\`, because 3 nonempty strata contain no valid ROI9. The narrow numerical interval above is the conditional extrapolation range for four completely unsupported observations; it does **not** measure the much larger missing-not-at-random uncertainty.

### Missing-not-at-random odds scenarios

For unknown original ROI9 attempts in supported strata, hypothesize a fixed shift in log odds of WHITE versus observed strict9 labels, using Jeffreys-smoothed success odds **for the scenario only**:

| Hypothetical unknown-vs-success log odds shift | Label fraction implied by the scenario |
|---|---:|
| −2.0 | 25.242% |
| −1.0 | 33.822% |
| −0.5 | 40.196% |
| 0.0 | 47.547% |
| +0.5 | 55.086% |
| +1.0 | 61.940% |
| +2.0 | 71.674% |

Each scenario retains the 4 observations in empty-stratum identification bounds. **These are NOT measured flower colours of failures** and not estimates of hidden natural pigment state. Non-ROI photographs can contain no open corolla at all; their hypothetical binary colour is merely a mathematical sensitivity convention.

### Empirical regional/seasonal heterogeneity, conditional on photographed observation frame

| Subset | Original strict9 white fraction | Conditional MAR-supported range |
|---|---:|---:|
| North extratropics | 48.5117% | 48.48310–48.48326% |
| South tropics | 44.8197% | 44.90295% |
| Local summer | 49.4923% | 49.41133% |
| Local winter | 45.0016% | 44.85035% |

These regional and seasonal comparisons are **photographic composition**, which can differ due to sampled taxa, observers, open-flower availability, image selection and source-arm design; no native flora or temporal floral plasticity claim follows.

## Operational conclusion

**The recoverable component of technical nonresponse identifiable under coarse arm/latitude/season MAR changes the *global selected-source photographic* white fraction by only ~0.21 points, but a one-log-odds departure from MAR can move the hypothetical white fraction to roughly 33.8% or 61.9%.** Thus sample expansion and strict9 successful-photo CIs alone cannot resolve the colour dependence of technical missingness. The distinction between a *source-conditional photo-label association* and a verified natural plant-community ecological law remains essential.

No human-review requirement is included in this PR. Further work, if any, should be a separate fully automated, source-locked environmental association/blocked validation with strict source-provenance and ecological estimand labels, **not another unstructured 93m-image rescanning exercise**.

## Machine-readable outputs

- Actual output tests + source gate: [Actions run](https://github.com/zuizui0223/mois/actions/runs/40856639794).
- Immutable four-file source-analysis release: [MOIS original ROI9 fully automatic selection/MNAR result v1](https://github.com/zuizui0223/mois/releases/tag/mois-original1181146-nohuman-auto-roi9-bias-20261029-v1). All GitHub release asset SHA256 digests match the generated source receipt. [Sealing Actions run 40857151162](https://github.com/zuizui0223/mois/actions/runs/40857151162) succeeded.
- SHA-frozen generated files: \`mois_auto_roi9_bias_receipt.json\`, \`mois_auto_roi9_bias_by_scope.json\`, \`mois_auto_roi9_source_geoseason_strata.csv\`, \`mois_auto_roi9_taxon_supported_only_summary.json\`.
- Data-source script: \`scripts/analysis/mois_observation_roi9_auto_bias_20261029.py\`.
- Pure observation identity: \`scripts/analysis/mois_roi9_observation_frame_20261029.py\`.
- Tests: \`tests/test_mois_observation_roi9_auto_bias_20261029.py\` and \`tests/test_mois_roi9_observation_frame_nohuman_20261029.py\`.
