# Can climate predict flower-colour differences between two previously unseen congeners?

**Original-source methodological and ecological validation — 15 October 2026.** This new test addresses the macroecological question at the most relevant out-of-genus *pairwise* prediction scale. The original FCP flower-photo and environment ZIPs were independently checked with historical SHA256 and joined locally; no source photos, taxon identifiers, site points, named genus rows, photographer credits or original measured site tables were deposited into public GitHub.

## Why the target matters

There are at least **three different ecological questions** behind the same global photos:

1. **Absolute white-photo fraction**: for a species whose genus has never been seen during model training, does a general environment model predict its photographed white fraction? Earlier original-source same soil36 cohort genus-holdout showed **+0.002355** squared-error improvement from adding rainfall+vapour to broad geography, 95% whole-genus fixed-loss interval **[+0.001125,+0.003608]**.
2. **Genus-centred residual photo-white fraction**: the same taxa, centred after observation of the genus' mean white-photo fraction. Earlier source-only genus-heldout diagnostic found rainfall+vapour gain **+0.000232 [−0.000104,+0.000564]**, but its target requires test-genus observed mean phenotype and cannot prospectively forecast that genus' absolute white fraction.
3. **A prospective congeneric *pairwise difference*** (this study): when the two taxa's environmental covariates are known, can a slope from training genera predict the **difference** in their flower-colour photo-white fractions? The fitted contrast is `predicted_white_i − predicted_white_j = beta · (environment_i − environment_j)`, requiring **NO test-genus phenotype label or white-colour mean at prediction time**. Actual observed white-photo fractions of both taxa are used **only afterward for scoring**. This does not require genotypes, but the phenotype is still conditionally classified photographic colour.

These tasks must never be claimed to measure the same predictive grain, and neither indicates heritable white-morph fitness.

## Original source, holdout and exact pairwise estimator

The exact SHA256 source archives produce **16,801 soil36-complete eligible named taxa, 1,665 nominal genera, and 26,163 original classifiable flower-hue photos** in the qualified species. Each source genus contains ≥3 classifiable taxa. All candidate models have the identical taxa/white-photo outcome and genus equal representation.

Create 5 **whole-genus-disjoint** folds using original SHA256 genus hashing, then repeat with three independently salted SHA256 whole-genus fold assignments. Each test genus is disjoint from its training genera. Within each training genus, centre the original site-mean geographic/climatic/edaphic variables and classified-photo white species fraction; fit a source-genus-equal zero-intercept ridge model with fixed regularization **0.01**, standardizing environmental variables on the training genera only. Test predictions of two species' contrast require the pair's predictor difference only; neither an observed test-genus white mean nor another test-species outcome is permitted as a predictive feature.

The response is `photo-white_i−photo-white_j`; error for each pair is `[(y_i−y_j)−(ŷ_i−ŷ_j)]²`. Score **all unordered within-genus eligible species pairs**, average the squared error **within each genus** and then give each test genus **one equal weight**. There is a computational identity:

```text
mean_{i<j} {[(y_i - y_j) - (yhat_i - yhat_j)]^2}
   = [2*n_g/(n_g - 1)] * mean_i[(e_i - mean_g(e))^2]
where e_i = y_i - yhat_i .
```

Thus the full all-pair loss is efficiently calculated in O(n), and directly checked against explicit O(n²) differences in synthetic data. The mean outcome difference is **never used to generate an individual or pairwise prediction**. This is distinct from evaluating a forward absolute-white species model or centring the response as an input required at prediction time.

For uncertainty, resample the **1,665 observed genus-level fixed test-loss vectors** 1,999 times separately within each fold-assignment experiment; do not refit the models for those resamples. The four fold schemes **reuse the same original source taxa** and cannot be interpreted as four independent ecological replications.

## Actual source findings

| Whole-genus fivefold assignment | Geographic-only heldout pairwise MSE | Added rainfall+vapour: gain over geography [ordinary genus-bootstrap 95%] | Soil5 gain beyond selected full climate/elevation [ordinary genus-bootstrap 95%] |
|---|---:|---:|---:|
| Original SHA256 | **0.316285** | **+0.000665 [−0.000227,+0.001506]** | **−0.000205 [−0.000803,+0.000356]** |
| Salt 202610151 | 0.316278 | **+0.000823 [−0.000067,+0.001705]** | +0.000226 [−0.000340,+0.000773] |
| Salt 202610152 | 0.316273 | **+0.000796 [−0.000068,+0.001650]** | +0.000105 [−0.000480,+0.000686] |
| Salt 202610153 | 0.316179 | **+0.000729 [−0.000184,+0.001622]** | +0.000102 [−0.000471,+0.000653] |

The original-fold geography-only and geography+moisture pairwise MSEs are **0.316285** and **0.315620**. The positive source point gain of **+0.000665** is about **0.21% of the geography-only pairwise squared error**, and its source-conditional 95% interval **contains zero**. Adding the full selected climate/elevation block does not demonstrate stable gains over the simpler water model. Adding soil5 to the selected full climate block has an **uncertain** sign that flips with the whole-genus fold assignment; none of the four 95% genus-resampled intervals excludes zero. It is incorrect to declare the true soil coefficient zero merely because cross-genus predictive gains are not shown here.

## Interpret the macroecology together, not selectively

**Positive geographical and taxonomic-level associations are established in the source model:** white flower photo proportion covaries with larger annual rainfall and atmospheric vapour among congeneric species, with original world genus-equal +2.255 and +3.952 pp/source SD and original 24-colour×axis within-scope simultaneous positive intervals. On the same 20,694 species, genus-mean gradients are larger but per-genus contributions are heterogeneous. This is a *source-conditional* macroecological trait pattern.

**However, most congeneric differences are not explained by measured abiotic axes**: 13 selected geographic, climate, altitude and soil variables explained only **0.6202%** of genus-centred source photo-white species-fraction variation in the separately eligible 16,801-soil36-species population. The genuinely transferable, out-of-genus *pairwise phenotypic difference* gain is small and statistically inconclusive on every studied deterministic whole-genus split. Hence the observed globally positive environment–trait slope should **not** be described as a universal rule reliably forecasting which of two previously unseen congeneric species will be whiter.

**Separately**, absolute photographed-white fraction prediction across unseen genera **does** improve modestly with water climate: that reflects a different predictive target, and is compatible with poor generalization of within-genus differences. The 1999-bootstrap uncertainty resamples fixed original photo-derived predicted losses rather than refitting, so it does not solve ecological nonrandomness, geographic transfer, photo missingness or phylogenetic dependence.

This sharpens the biological question: **why does a broad geographical moisture gradient coexist with so little predictable within-genus flower-colour differentiation?** Plausible, currently unseparated mechanisms include other phylogenetically conserved traits below genus, pollinators, flower developmental pathways, species ranges/biogeographic historical structure, microhabitat, and image classifier selection. There is no direct floral pigment cost/benefit, local slope/aspect/TPI/TRI terrain measurement, genotype frequency or lifetime fitness evidence.

## Reproducibility

- Owner-only original-FCP-ZIP source replayer: `scripts/analysis/mois_pairwise_congeneric_transfer_20261015.py`. Requires original photo + environmental ZIPs, SHA256 checks, frozen source labels, and `--write` or `--check`. Synthetic all-pairs identity and privacy-safe frozen receipt assertions are tested in CI; raw photos remain outside the public checkout.
- Frozen **aggregate-only numerical result**: `results/mois_pairwise_congeneric_photo_white_transfer_20261015.json`.
- Companion original same-species all-13-environment in-sample R² and retrospectively genus-centred heldout CV: `results/mois_env_shapley_and_genus_residual_transfer_20261015.json`.
- Supplement placement: **S2.17**. The main GEB article should give no more than one scope-bounded paragraph connecting global macro pattern, low same-genus explanation and weak congeneric-pair transfer. The complete source biological interpretation remains descriptive, not evolutionarily causal.

**Scientific decision:** Broad source-photo water-climate associations are real within the original specified statistical design, but **pairwise congeneric trait differences fail to show stable, precise transfer to taxonomically unseen genera**. That distinction is more ecologically informative than further post hoc searches for a single stable climate/soil significance threshold.
