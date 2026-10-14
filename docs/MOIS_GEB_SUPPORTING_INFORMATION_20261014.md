# Supporting Information — observational global photographic flower colour and moisture

**Companion to** `docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md`. All values come from migrated original-FCP SHA-locked result files or unchanged source-verification receipts, not new image opening, p-value search or independent replication. The `mois` repository is the interspecific study; H2 high-depth within-species colour-space results remain in `fcp`.

## S1 | Data source, exact units and original analysis lineage

The source photo bundle (SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`) comprised 42,111 taxon breadth records and 85,337 taxon×cell records. The two original photographic source programmes had 26,905 shared photo IDs, so a unique-photo-ID union contained 100,543 photographs of 42,111 nominal taxa; exactly 45,953 were assigned one of four camera-visible floral colours. The original environmental bundle (SHA256 `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`) contains 36 extracted predictors at the observed photo sites. Exact matching requires taxon, observation, and photo IDs simultaneously.

| Data grain | Photo count | Named taxa | Notes |
|---|---:|---:|---|
| Original unique photo opportunities | 100,543 | 42,111 | 45,953 classified; 54,590 not classified |
| Climate26/photo-site complete | 45,836 | 25,028 | Among four-hue-classifiable original photos |
| Climate26/genus-eligible signed contrasts | Subject to 26 complete-case plus ≥3 congeneric species | 20,694 /1,934 nominal genera | Principal photographed hue–climate coefficient population |
| Soil36/photo-site complete | 31,925 | 20,786 | Ten more soil features required |
| Soil36/genus-eligible signed contrasts | Subject to soil36 complete-case plus ≥3 congeneric species | 16,801 /1,665 genera | Different taxa and genera from climate26 |
| Matched spatial B/C world | 14,268 | 9,370 /845 genera | Exactly the same B species/photo outcomes in stage C |
| Matched spatial B/C East 90–180° | 4,149 | 2,886 /280 genera | Exactly the same B species/photo outcomes in stage C |
| Same-photo congeneric versus conspecific BIO12 | 21,448 | 6,143 /787 genera | Same original image IDs, different model-centering levels |
| Independent botanical descriptions, eastern sector | 2,355 matched FCP photo sites | 670 /58 genera | Independent described *hue*, original photo-site geography |

The original macro programme evaluated 20 source photo/site model scopes, 1,280 fitted directional coefficient/sensitivity rows, six or seven ecological axes and four photographed hues; the full table is `results/fcp_macro_direction_20261013/direction_ledger.csv`. Original code is `scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py`; original source audit `independent_verified.json` independently reproduced **352 original source-row WLS coefficients**, not an independently rebootstraped confidence interval. The matched spatial source has 576 stage coefficient/difference rows and 288 independently reproduced WLS slope coefficients in `results/fcp_m1_matched_spatial_support_20261014/`. All source files are byte-verbatim clones of their historical FCP Git blobs.

The original model uses **within-nominal-genus** species-centred four-colour fractions, one source species row per original taxon, a source photo-site environmental SD, and broad geographic coordinates `abs_lat/lon_sin/lon_cos`. The strict image-quality gate is flower effective pixels ≥23,826 and `nuisance_pixel_fraction ≤0.1375717742`. The cell index is 18 longitude × nine equal-area sin(latitude) bins; it is not within a 50-km biological location. Original confounder geography, cohort membership and maxT intervals are different for the regional/strict contrasts. Original specimen organs were not expert revalidated for the complete photo source.

## S2 | Complete signed-hue result families and photo-classification missingness

### S2.1 | World genus-equal photographed hue/BIO12 coefficients

The world climate26 broad genus-equal white/BIO12 coefficient is **+2.254545 pp/source SD**, 24-comparison approx maxT95 **[+0.318741,+4.190349]**. World white/mean annual atmospheric vapour-pressure is **+3.951850 [1.199253,6.704447]**. Soil36/photo-complete genus-equal white/BIO12 under a different eligibility set is **+2.257 [0.138,4.376]**. This latter result is a contrast from a different source population; it is not a paired response to adding soil. For all four original colours and all 20 scopes, retrieve the 1,280-row `direction_ledger.csv` rather than selecting individual significant sectors.

### Four-hue compositional check and Figure 4

The original 1,280-row direction ledger retains the full **four photographed hue outcomes**, not only white. The following world-climate26 estimates have fixed genus-equal weight, broad geographic adjustment and the **same 20,694 qualified species/1,934 genera**. The intervals are simultaneous only within the original **six axes × four hues**, not across the entire historical analysis programme.

| Original classified hue | Annual precipitation BIO12 coefficient [maxT95] (pp/SD) | Annual mean vapour pressure coefficient [maxT95] (pp/SD) |
|---|---:|---:|
| White | **+2.25455 [+0.31874,+4.19035]** | **+3.95185 [+1.19925,+6.70445]** |
| Yellow/orange | −1.22082 [−2.98702,+0.54538] | −1.38543 [−3.81809,+1.04722] |
| Red/pink | −0.44593 [−1.51309,+0.62123] | −1.50897 [−3.30513,+0.28720] |
| Blue/purple | −0.58779 [−1.42880,+0.25321] | −1.05745 [−2.41338,+0.29848] |

Each column's four point slopes sums to zero within floating-point tolerance because four photo probabilities sum to one. Only white's own intervals exclude zero; **the data do not identify a unique corresponding decline in any particular coloured category** after 24-contrast adjustment. The archived original CSV, not this rounded table, controls the precise coefficients. Figure 4 regenerates from `results/fcp_macro_direction_20261013/direction_ledger.csv` with `node scripts/render_mois_four_hue_moisture.mjs --check`. This numerical figure is reporting-only, not another environmental fit.

The broad all-source climate26 genus-equal strict-photo white/BIO12 becomes **+1.229162 [−2.045756,+4.504079]**; genus×coarse-cell **+0.807976 [−3.556384,+5.172336]**. The source scientific verdict is `universal_global_abiotic_hue_direction = NOT_SUPPORTED_BY_ALL_REQUIRED_SENSITIVITIES`, while broad source-specific evidence is retained as positive.

### S2.2 | Source-original colour availability and dry/wet deciles

Overall original-photo classifiability: **45,953/100,543 = 45.705%**, unknown **54,590/100,543 = 54.295%**. Among the classifiable source photos, **21,574/45,953 = 46.94797%** are classified white. No-assumption lower/upper white fractions among *all* original photos are **21.4575% to 75.7527%**; these are extreme bounds, not confidence intervals.

| Rainfall bin | All original photo opportunities | Classified original hues | Classification success | White among classified |
|---|---:|---:|---:|---:|
| Driest source opportunity decile | 10,084 | 4,981 | 49.40% | 41.7988% |
| Wettest source opportunity decile | 9,990 | 3,904 | 39.08% | 53.2531% |

The observed wet-minus-dry classified-white fraction is **+11.454 pp**, a raw descriptive difference; it cannot be compared directly to the genus-adjusted source-BIO12 slope (+2.255 pp/global SD). The original MNAR sensitivity assumes the unclassified *dry* set is 50% white. A **41.2637% white proportion in the unclassified wet photos** would then eliminate the implied latent *all-photo* wet–dry white difference. This is a tipping condition, not a measurement of the missing photos. A classifier-success inverse-probability-weighting model depends on a missing-at-random assumption and does not establish colour-agnostic selection.

The previous species-conditioned global-latitude white-fraction evenness null was not supportive (*P*=0.214), while the hypothesised negative geographic compensation between species turnover and local colour deviations had one-sided *P*=0.997 and observed **positive** covariance +0.00018144. Thus the worldwide photo share near one half is neither unusually stable under the original null nor an evolutionary balancing equilibrium. These photographic balance tests are secondary to `mois`, not a core pigment-fitness result.

## S3 | Same taxon identity, geographic support, and photo-scale bridge

### S3.1 | Matched geographic species B/C reconditioning

At the world source scale, **845 genera /9,370 taxa /14,268 photos** qualify for matched B/C. Genus-equal white/BIO12 broad geography slope B=**+1.32586**, within-genus×cell C=**+0.80798**, B−C **+0.51788**, 24-hue–axis maxT95 **[−1.25351,+2.28927]**. East 90–180 E: **280 genera /2,886 taxa /4,149 source photos**, B=**+1.10200**, C=**−2.07568**, B−C **+3.17768**, maxT95 **[−0.25804,+6.61340]**. The eastern *nominal* paired 95% interval **[+1.07697,+5.35072]** excludes zero, but the source-defined 24-contrast interval does not. The difference in conditioning is descriptive; B and C do not fit an identical conditional geographical estimand, and a true causal spatial effect is not estimated.

### S3.2 | Same 21,448 original photo outcomes, between versus within named species

Original PR #158, `results/source_receipts/matched_congeneric_vs_withinspecies_rainfall.json`, retains **21,448 identical photographic colour outcomes from 6,143 taxa/787 genera** and BIO12 source SD **736.416 mm**. Under respective level-specific statistical fits, the original congeneric slope is **+2.72219 pp/SD** (source genus-bootstrap95 **[+0.70471,+4.71383]**); the within-species slope is **−0.51068** (**[−1.43610,+0.33512]**); difference **+3.23286** (**[+1.17542,+5.40032]**). The original retrospective one-sided paired 1,999-bootstrap *P* for a nonpositive difference was 0.0005. Positive within-species photographic white/water correspondence is not supported on these same images. The two slopes differ in statistical centering and are not comparative fitness coefficients.

A distinct originally executed PR #159 sensitivity uses *continuous segmented white-palette ROI area*, not the four-category photo label, in **14,667 source photos /3,351 taxa /489 genera**. On its own **717.253 mm** BIO12 source SD, the congeneric coefficient is **+2.229 pp/SD** (ordinary source CI +0.177 to +4.282), within-species **−0.180** (−1.011 to +0.651). A stricter source photo ROI selection gives the congeneric coefficient **+1.116** (−1.169 to +3.401), which is inconclusive. The source was not independently measured with calibrated optical equipment.

## S4 | Joint environmental prediction, scale specificity and model-form sensitivity

These estimates concern **four-hue squared-composition-error gains from adding joint environment predictors**, not directions of individual colour–BIO12 regression coefficients. The source prediction populations, folds, baselines and response levels differ. Where possible, use exact same original source taxa and same evaluation folds before comparing predictors.

| Original source prediction task | Species / genera | Genus-equal improvement (95% source bootstrap) | Decision boundary |
|---|---|---|---|
| Climate26 geographic-heldout, genus/geography baseline | 17,183 source test species | +0.002943 [+0.001749,+0.004104] | Source climate26 scope |
| Soil36 geographic-heldout, genus/geography baseline | 13,961 *different* test species | +0.001432 [+0.000216,+0.002755] | Cannot subtract first and second rows |
| Whole-genus holdout, absolute all36 photo colour | 20,786 /4,853 | +0.004518 [+0.002437,+0.006602] | Small ~0.77% source baseline gain |
| Whole-genus holdout, *within-genus species deviations* | 16,801 /1,665 | −0.000142 [−0.000748,+0.000433] | Not supported for within-genus differences |
| Whole-genus holdout, within-named-species photo differences | 5,939 /2,196 | −0.000670 [−0.001038,−0.000309] | Original predictor fails for this scoring target |
| Soil10 incremental prediction on **same 20,786 species** | 20,786 /4,853 | +0.0000757 [−0.000516,+0.000653] | Soil-only added value not established |
| All36 unseen genus only, fixed regional-fold species | 20,786 /4,853 | +0.004312 [+0.002174,+0.006357] | Source taxonomic transfer |
| Same taxa, unseen genus + longitude sector | 20,786 /4,853 | +0.000759 [−0.002011,+0.003535] | Source geographically disjoint transfer weak |
| Same taxa, unseen genus + sector + 500 km | 20,786 /4,853 | +0.000174 [−0.002629,+0.002846] | Not a robust spatially portable linear model |

Figure 3 visualizes three separate original-source contrasts: target definition (panel a), taxonomic/geographic heldout support (panel b), and nonlinear versus linear predictive baselines (panel c). These are not comparable estimates of one biological adaptation coefficient.

The original full36 whole-genus holdout **same fitted-prediction** squared-error decomposition is full absolute improvement **+0.004518 = +0.004699 genus-mean contribution −0.000181 within-genus species-deviation contribution**. This is an algebraic error partition conditional on the original fit, not a biological variance decomposition. Exact original receipts: `results/source_receipts/unseen_genus_colour_environment.json` and `unseen_genus_other_geography_500km.json`.

A pre-fixed original ExtraTrees nonlinear model under unseen genus + other sector +500km found environment gain **+0.004882** (95% [+0.002320,+0.007496]) over *its own* nonlinear geography baseline. The matched nonlinear-full versus ridge-full prediction difference was **+0.001202** (−0.001011,+0.003447), which does not establish overall nonlinear superiority. Nonlinear geography-only loss was **0.604125**, *worse* than ridge geography-only **0.600619**, so differences in environmental increments should not be interpreted as direct full-model efficacy. Original strict photo sensitivity found nonlinear environment gain **+0.000517** (−0.004122,+0.005363). Source: `results/source_receipts/nonlinear_geographic_transport.json`. All these coefficients were computed before transfer to mois and are retrospective.

The original four-hue Brier decomposition attributed **51.5%** of the source climate26 total gain to the white squared-error component and **45.4%** in soil36, but these components are algebraically constrained by four probabilities summing to one. They do not measure the proportion of pigmentation, flower-colour variance or fitness explained by environment. Source report: `docs/supporting/FOUR_HUE_ENVIRONMENTAL_PREDICTION_DECOMPOSITION.md`.

## S5 | Independent botanical hues and limits of mechanistic identification

A separate Natural History Museum text-derived floral-colour source included **59,369 named source species records**, with exactly **2,476 white or blue/purple described taxa** overlap with original FCP photographed source species sites. The original restrictive East-vs-rest support retained **670 taxa /58 shared genera**; genus-equal white-versus-blue rainfall contrast **+3.555 pp**, 95% source genus-bootstrap **[−15.340,+22.603]**, ordinary interval including zero. All original site climates were still FCP iNaturalist photographic source locations. The source classifier did not provide four comparable photo colours at the same taxon/organ grain. No independent photo-site geographic reproduction was identified.

**Inference ceiling across all files:** photographed flower hue is not guaranteed corolla reflectance, UV nectar guide, anthocyanin pathway expression, flower thermal defence, flower pigment maintenance cost or pollinator-perceived optical contrast. Nominal genus is not a dated phylogeny, and variation of historical site/climate across different congeneric taxa is not a causal within-genotype intervention. The source contains no original morph-specific seed/pollen relative fitness or plant genotypes. The separate `fcp` high-depth within-species phenotype-direction W analysis is **outside the main `mois` hypotheses**; no original new prospective H2 photos or plant traits were created for this macro source.

## S6 | Provenance, verification and bibliography boundaries

- Original photo archive SHA256: `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`
- Original environmental archive SHA256: `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`
- Original fixed worldwide direction and source: [FCP PR #193](https://github.com/zuizui0223/fcp/pull/193), original independent code+result verification (352 direct source-row WLS).
- Original identical-species geographic conditional support: [FCP PR #194](https://github.com/zuizui0223/fcp/pull/194), 288 direct source-row WLS and the 576-row paired source ledger.
- Original full photo+environment replay: [FCP Actions 38502660130](https://github.com/zuizui0223/fcp/actions/runs/38502660130).
- `mois` source file checksum inventory and historical Git blobs: `docs/MOIS_ORIGINAL_FCP_PROVENANCE_AND_TRANSFER_20261014.md`; exact byte checks in `scripts/check_mois_migration.py`.
- Full end-to-end photo- and abiotic-source replay requires the original two large ZIP inputs, which **are not shipped in mois**. Mois CI confirms unchanged tables and synthetic model behaviour; it is not independent biological replication or new original-photo bootstrap verification.

### S6.1 | Post-transfer actual ZIP recovery and independent source-table coefficient audit

Both original Actions artifacts were obtained in a separate working environment and byte-SHA256 verified as **34,282,732-byte photographed flower-colour/site tables** and **38,759,029-byte expanded abiotic/source-site tables**. There are no original camera image pixels inside either ZIP. An independent estimator, written separately from the FCP source regression implementation, recomputed 100,543 unique original photo IDs, 45,953 four-colour classifications, 45,836/25,028 climate26 photo/species cases, 31,925/20,786 soil36 cases and 20,694/1,934 genus-eligible species/genera. Eight original genus-equal point coefficients—four photographed hues on BIO12 and four on mean water vapour—matched the preserved original 1,280-result source CSV within **3.51×10⁻¹⁴ percentage points**. The confidence intervals from 999 whole-genus resamplings were **not** independently recomputed in this particular test.

The aggregate-only `results/mois_original_source_independent_recovery_audit_20261014.json` records licence counts (78,367 CC BY-NC; 13,795 CC BY; 3,229 CC0; 2,715 CC BY-SA; 2,433 CC BY-NC-SA; two CC BY-NC-ND and two unspecified) and source coordinate-status counts (100,348 public source points marked valid, 195 other status flags). No per-photo IDs, user attribution, coordinates or the original ZIPs were added to this public repository. Exact original SHA values, the optional input audit script and outstanding human source-licence/geoprivacy authorisation are described in `docs/MOIS_GEB_ACTUAL_SOURCE_RECOVERY_AUDIT_20261014.md`. Successful local byte recovery is distinct from a licensed, stable, double-anonymous reviewer-data archive, which remains **unavailable**.

### S6.2 | Fully reconstructed original world-climate26 24-colour interval family

The earlier post-transfer eight-slope audit (§S6.1) did not re-estimate bootstrap confidence intervals. **A subsequent separate original-source calculation did.** Using both original ZIPs, a new implementation averaged photographic four-hue fractions per species, applied the frozen ≥3 species/genus support rule to 20,694 species/1,934 genera, centred species values within nominal genera, computed anonymous per-genus sufficient-statistic normal equations, and executed **999 genus multinomial bootstrap resamples with source seed 2026101430**. For each of six fixed ecological axes × four mutually exclusive photographed hues, both ordinary 2.5–97.5% percentile intervals and **within-scope 24-contrast maxT simultaneous 95% intervals** were independently recomputed. All 24 point coefficients, 24 ordinary intervals, 24 maxT intervals and 24 support indicators matched the frozen original 1,280-row CSV, maximum absolute numerical difference **3.55×10⁻¹⁵ percentage points**, with original maxT critical value **3.0109108050453743**.

This numerically replays exactly the **world / climate26 / all-photo / genus_geo / genus_equal** subset; it does not rerun the other historical scopes or correct the entire source-exposed programme of geographical/model choices. A 1.06-MB anonymous-genus crossproduct/cluster-count `.npz` enables this specific 24-contrast bootstrap replay without original photo IDs, nominal genus names or exact site coordinates, but it still encodes aggregated geographical covariance and has **not passed release/licensing/geoprivacy clearance or been committed publicly**. The aggregate-only 24-result audit is `results/mois_world24_genus_moment_source_replay_audit_20261014.json`; source-only replay code is `scripts/analysis/replay_mois_world24_from_genus_moments_20261014.py`; public inference/rights contract: `docs/MOIS_GEB_WORLD24_GENUS_SUFFSTATS_REPLAY_20261014.md`. CI checks the committed result against original frozen 1,280-row coefficients but cannot rerun the original 999 draws without the restricted private moments. This improves source-data computational reproducibility, not biological replication, image-pixel validation or journal dataset availability.

### S6.3 | Independent original-source 576 spatial estimates and 999-bootstrap paired maxT intervals

A further independent source-table implementation reproduced original FCP #194 three-stage matched-spatial contrasts from the two exact original SHA256-locked ZIPs. For the **world** and **East 90–180°** data, it held each original source photo hue/site environmental value and named species fixed between stage B (broad geographic nuisance adjustment) and C (within-genus×coarse-cell conditioning), and independently recreated the source 999-replicate whole-genus bootstrap for both species-equal and genus-equal weights. Each of four region×weight cases has **3 stage estimates +3 paired stage differences**, six prespecified ecological axes and four visible-colour classes, yielding **576 fitted source rows** overall.

All **576 original points**, **1,152 ordinary 95% confidence interval endpoints**, and **576 paired-difference maxT 95% interval endpoints** matched the original frozen 576-row CSV, maximal discrepancy **1.60×10⁻¹⁴ percentage points**. For genus-equal weight, world B−C white/BIO12 **+0.517884 [−1.253507,+2.289275]** and East B−C **+3.177682 [−0.258035,+6.613398]** both have corrected intervals covering zero. The original nominal eastern comparison excluding zero is not a source-familywise finding. The portable *independently implemented* source-code replayer is `scripts/analysis/replay_mois_original_matched_spatial_576_20261014.py`; aggregate-only original-data receipt `results/mois_original_matched_spatial_full576_replay_20261014.json` and complete scope/provenance explanation `docs/MOIS_GEB_INDEPENDENT_SPATIAL576_SOURCE_REPLAY_20261014.md`. Local originals used for verification are **not released as GitHub source files**; the CI checks source receipt and original 576 rows and compiles the replayer but cannot rerun original source ZIPs without authorization. This is statistical identity on **historically exposed** data, not independent biological replication, a causal geography intervention, or programme-wide multiplicity correction.

### S6.4 | Additional global strict photo, genus-by-cell and soil36 complete source-replay families

Three further originally fitted genus-equal worldwide photo-colour families were independently reconstructed from the **same two exact source ZIPs** using the preserved photo/observation/taxon join, original photo-quality thresholds, original species means and original 999-genus-bootstrap seeds. Source denominators and the complete 24- or 28-contrast within-scope maxT family were verified as follows:

| Original source model scope | Photo support / species / genera | Within-family maxT critical value | Original familywide sign exclusions |
|---|---|---:|---|
| Strict photo quality, climate26, broad geography | 8,538 photos /5,848 species /727 genera | 3.0692641804 | **None of 24** |
| Genus×coarse-cell geography, climate26 | 14,268 /9,370 /845 (1,847 genus×cell groups) | 3.0729374107 | **None of 24** |
| Soil36 complete-case broad geography | 31,925 /16,801 /1,665 | 2.9648395022 | 3 of 28: white–BIO12, white–vapour, white–soil-water proxy |

All **76 coefficients, 76 ordinary 95% bootstrap percentile intervals and 76 original within-scope maxT 95% intervals** and their original support flags matched the frozen original 1,280-row macro coefficient file. Max absolute error was **4.44×10⁻¹⁶ percentage points** for the two climate26 sensitivity scopes and **3.11×10⁻¹⁵** for soil36. The original white/BIO12 strict-photo **+1.229 [−2.046,+4.504]** and genus×cell **+0.808 [−3.556,+5.172]** remain within-family nonconfirmations. White/vapour also has maxT intervals spanning zero in both restricted climate26 scopes.

In the *separately sampled* soil36 source, white/BIO12 **+2.257 [0.138,4.376]**, white/vapour **+3.987 [1.038,6.936]** and white/SoilGrids-derived soil available-water **proxy +1.698 [0.0028,3.393]** are source-scope discoveries under **seven×four** simultaneous intervals. The soil-water proxy finding has a lower limit barely above zero and does **not** establish an ecological mechanism or soil's *incremental predictive value* after climate: the source-verified identical-soil-population heldout comparison gave gain **+0.0000757 [−0.000516,+0.000653]**. Original source photo quality/site observation and soil-complete eligibility remain unresolved sources of biological selection.

Full aggregate-only audit with 76 original-source coefficient and CI rows: `results/mois_original_world_strict_geocell_soil36_76maxT_source_audit_20261014.json`. Standalone optional source ZIP replayer: `scripts/analysis/replay_mois_world_quality_geocell_soil36_from_original_zips_20261014.py`; complete methodology, limits and source identity: `docs/MOIS_GEB_WORLD_SENSITIVITY_76_INDEPENDENT_REPLAY_20261014.md`. Public CI checks the summary results against all relevant frozen source CSV rows and tests the replayer syntax; the original ZIPs stay outside the public repository. This is **statistical replay on the original data**, not independent species sampling, camera-image remeasurement, programme-wide multiplicity correction, trait physiology or pollinator-fitness evidence.

