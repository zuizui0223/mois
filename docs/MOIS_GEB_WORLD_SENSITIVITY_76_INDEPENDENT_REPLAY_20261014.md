# Original-source independent rechecks of worldwide flower-colour × moisture sensitivity families

**14 October 2026 | status: PASS, reporting-only technical verification.** Three additional original FCP analysis scopes were **re-estimated in the working source-data environment** from the original exact photo-colour/site-measurement ZIP and the exact expanded 36-abiotic/site ZIP, both SHA256-verified. Their original four-colour outcomes, species identity, photo source and model definitions were unchanged. The source images' pixel data are not present. This check does **not** establish causal colour–water selection or supply independent biological observations.

## Exact original scopes and sampling denominators

The independent Python scripts used the original photo/observation/taxon one-to-one join, original frozen classification and photo technical-quality variables, original equal-area cell assignment, original species-photo means and original climate-versus-soil distinct complete-case normalization. Genus-equal weighting was by named genus, **including genus×cell source subgroups weighted by their genus's total eligible species**, not equal weighting of separate cells. The original 999 genus-resampling seeds and within-scope maxT across all six (or seven) source ecological axes × four original photo-hue outcomes were preserved.

| Original source scope (world / genus-equal) | Eligible source data | Frozen within-scope maxT critical | Verified values |
|---|---|---:|---|
| Strict photographic quality, climate26, broad coordinates | **5,848 species /727 genera /8,538 photos** | **3.0692641804**, seed 2026101455 | 24 point coefficients, 24 original percentile CIs, 24 original maxT95 CIs |
| Coarse genus×equal-area-cell condition, climate26 | **9,370 species /845 genera /14,268 photos** (1,847 eligible genus×cell subgroups) | **3.0729374107**, seed 2026102253 | 24 point coefficients, 24 original percentile CIs, 24 original maxT95 CIs |
| Soil36-complete broad geography | **16,801 species /1,665 genera /31,925 photos** | **2.9648395022**, seed 2026101397 | 28 point coefficients, 28 original percentile CIs, 28 original maxT95 CIs |

All **76 separate axis×hue coefficients**, 76 ordinary percentile intervals, 76 original scope-specific simultaneous intervals and their 76 source support indicators matched the frozen 1,280-row source result. Maximum absolute error in **strict/cell**: **4.44×10⁻¹⁶ percentage points**; soil36: **3.11×10⁻¹⁵ percentage points**. This is a numerical identity check on the same original data, not a new independent study.

### Biological focus: the two main negative sensitivity gates were reproduced

| Sensitivity population | White versus annual precipitation BIO12 pp/source SD | Source familywise maxT95 interval | Original decision |
|---|---:|---|---|
| Global climate26: strict photography | **+1.22916** | **[−2.04576,+4.50408]** | Does not establish signed association |
| Global climate26: within-genus×coarse-cell | **+0.80798** | **[−3.55638,+5.17234]** | Does not establish signed association |

Both populations also had source white/vapour-pressure point estimates with original simultaneous intervals crossing zero: strict **+1.72503 [−2.89330,+6.34335]**, genus×cell **+2.74504 [−3.56633,+9.05641]**. Across **all 24 hues×axes** in each of these two restricted populations, **none** had an original within-family maxT 95% interval excluding zero.

This does **not** mean the global signal was disproven by holding all other factors constant: image quality and spatial-conditioning restrictions change sample eligibility/support and the statistical contrast. The correctly limited claim is **lack of robust positive signed support under the stricter conditions**, not proof of environmental absence.

### Soil36 source effects are real source-model outputs, not added-soil causality

The soil36-complete photo population reidentified white/BIO12 **+2.25694 [0.13815,4.37573]** and white/mean-vapour **+3.98701 [1.03799,6.93603]** under their **separate seven-axis×four-colour within-family maxT**. Original white/soil available-water **proxy** was **+1.69770 [0.002761,3.39264]**; the simultaneous lower limit is only **+0.0028 pp**. This borderline photographic coefficient does **not** identify physical soil-water causality, because original vegetation, soil/precipitation covariates and site classification are confounded and the source sample differs from climate26. An original **same 20,786-species** joint held-out full36-versus-climate26 prediction analysis estimated soil-only added Brier gain **+0.0000757 [−0.000516,+0.000653]**; the incremental predictive evidence was not confirmed. These are different original-source estimands, **not conflicting data and not proof soil is irrelevant**.

## Reproducibility handoff (no redistributable photo-level source data)

- Machine-source original-output correspondence: `results/mois_original_world_strict_geocell_soil36_76maxT_source_audit_20261014.json`. Contains only original axis/hue aggregate coefficients, sample supports, original bootstrap seed and critical value and source numerical-replay error. **No original photographer, photo ID, taxon ID, taxon name, exact location or site measurement rows were committed**.
- Standalone independent ZIP replayer: `scripts/analysis/replay_mois_world_quality_geocell_soil36_from_original_zips_20261014.py`; requires the original **two checksum-frozen source ZIPs** and original `results/fcp_macro_direction_20261013/direction_ledger.csv`, e.g.:

```bash
python scripts/analysis/replay_mois_world_quality_geocell_soil36_from_original_zips_20261014.py \
  --photo /private/original_100543_photo_sites.zip \
  --environment /private/original_100543_photo_sites_all36_environment.zip \
  --ledger results/fcp_macro_direction_20261013/direction_ledger.csv \
  --scope strict_photo_climate26

# Repeat only for the other original frozen scope names:
# --scope within_genus_geocell_climate26
# --scope global_complete_soil36
```

These original source ZIPs are held in a separate restricted execution environment, **not GitHub public repository or stable reviewer archive**. Public CI checks exact original 1,280-source row values, the independent numeric receipts and replayer syntax, but cannot rerun original source rows without those private ZIPs. No new study-wide false-discovery correction is claimed: original maxT is **only within each separately fitted photo/population/mode/family**, not across historical exploratory scopes or regions. The two sensitivity cohorts have fewer species than the broad 20,694 panel and differing observation opportunities. No source-wide genetic polymorphic morphs, physiological pigment investment, pollinator adaptation or plant relative seed/pollen fitness were observed.

**Bottom line:** world positive white–moisture association and the strict/photo-geography nonconfirmations are **all independently numerically reproducible on original measured source tables**. Their ecological mechanisms and broad transfer remain unproven. Photo-image licensing, source coordinate privacy, stable anonymous reviewer archive/DOI, and human author sign-off remain the submission blockers.
