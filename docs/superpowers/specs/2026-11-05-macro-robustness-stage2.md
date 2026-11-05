# Stage 2: transportability and experimental-evidence sensitivity

User mandate: set the next goal and execute. Baseline is the completed analysis at 25a77529d1ac8556a84909a5ea7dbbaa55b90725, not a newly selected significant subset. This is a post-primary robustness investigation, not retrospective preregistration. Baseline outputs remain separate.

## Locked geographic/taxonomic checks

Reuse all 82,541 complete cell-years and the original covariates/estimand. Test the ALAN-by-southern-hemisphere interaction with ALL six predictor slopes allowed to differ and separate hemisphere/year effects; regional effects are already geographically separate. Report the interaction, not the difference between two significance decisions.

Within the southern sample run all three explicitly named longitude sectors (<-30, -30 to <60, >=60 degrees), and omit each in turn. These are coordinate sectors, not inferred country/continent membership. Omit each of the ten observation years in turn. For taxonomic influence, select the ten genera with greatest southern cell-year presence without using effect estimates, subtract their exact species/white-species counts, and reapply >=10 species / >=5 genera gates. Save every attempted estimate, failure and changed denominator. No coefficient-directed deletion or replacement of the baseline.

## Experimental evidence and phylogeny

Use the same 65 selected experimental species and the frozen Figshare v1 records and published phylogeny mapping. Audit species-to-tip uniqueness and tree support; do not fabricate missing branches. Compare the original ratio with complete-exclusion-controlled reproductive contrasts where available, preserving missing controls and negative contrasts. A contrast is not automatically proof of a specific effective nocturnal pollinator. Compare the same subset under old and new definitions before attributing a difference to adjustment.

Phylogenetic analysis is an explicitly separate sensitivity with its own sampling unit and model assumptions. Do not label the original spatial/species clustered model phylogenetically corrected. Record tree source, mapping, uncertainty and any reduced support.

Control-adjusted index: (night - complete exclusion)/(night + day), with a positive denominator and nonnegative original outcomes. Retain negative net contrasts; do not discard them to create an apparent dependence fraction. Compare all-species original, control-species original, identical control-comparison raw ratio, and identical control-comparison net index. The index is not a probability. First estimate adjusted species light profiles with species fixed effects plus the original region/year/environment controls; export the FULL spatial-cluster sampling covariance. A second-stage phylogenetic meta-regression uses that covariance plus a published-tree correlation component (metafor rma.mv); compare no-phylogeny and tree models on identical rows. This is a separate comparative sensitivity; literature-ratio measurement uncertainty remains a limitation.

## Completion

Actual sensitivity fits, audited data support and source linkage, meaningful numerical verification, all attempted-model results, and a saved stage-2 report separating robust association, conditional association, and unsupported mechanism. Failure of robustness is a valid scientific outcome; obtaining significance is not the objective.

## Completed evidence

All 28 prelisted influence models, four matched-control definition models, and six comparative/phylogenetic models estimated successfully. Southern baseline +2.279 [0.671,3.887] attenuates to +0.163 [-0.821,1.146] when the >=60E sector is omitted; single-year and top-ten-genus omissions retain positive individual intervals. This is geographically limited robustness, not proof of zero elsewhere. Formal hemisphere interaction +1.998 [0.300,3.695] remains an exploratory conditional comparison.

Control definitions use 105 comparisons/46 species, retaining all ten negative net-night contrasts. All 65 species map uniquely to the published tree. First-stage species-profile model is full rank249 over122548 rows; full spatial-cluster sampling covariance is positive definite. Only verified roundoff asymmetry6.7e-15 was symmetrized. All white-by-diel interactions remain uncertain with and without phylogeny; conditional white contrasts at index0.5 are retained rather than discarded.

Independent explicit-dummy hemisphere WLS agrees within1.11e-13. Independent GLS given the fitted variance components agrees with all six metafor fits within3.34e-15; this does not independently reoptimize REML. Two focused new tests passed, all attempted-model labels and control aggregation/source hashes audited, and two figures visually verified. Stage-2 report, exact inputs, code, support/variance tables and verification are saved under outputs/mois_global_macro_stage2_20261105. Single-tree and literature-index measurement uncertainty remain limits, not completed corrections.
