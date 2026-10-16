# Original-source species and genus weighting: is white globally half of photographed floral colour?

**Independent calculation from the exact historical FCP photo measurement ZIP; 14 October 2026.** Scientific scope: **opportunistically photographed visible flower-colour labels**, not a botanical census of all species or a genetic selection equilibrium. Source archive SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0` was obtained, verified byte-for-byte and ZIP-integrity checked in the execution environment. Only aggregate numerical results and the executable replayer were published; no source photo IDs, photographer attribution or site coordinates were deposited.

## Question and weighting methods

The 42,111 breadth photo rows and 85,337 taxon-by-cell rows share exactly 26,905 photo IDs. Preserve the first breadth occurrence when dropping duplicate original photo IDs: **100,543 distinct photo opportunities from 42,111 nominal species/taxa**. Classifiable means the historical `measurement_status=classified_four_state_morph` and one of `white`, `yellow_orange`, `red_pink`, `blue_purple`. All other photos remain **unclassified**, never recoded as colourless.

For each classified photo `j` in nominal source taxon `i`, let `y_ij=1` for white and 0 for any of the three other classes; `n_i` is the number of classifiable photos for that taxon. Then define three different descriptive estimands:

- **Photo equal:** total classified white photos / total classified photos.
- **Species equal:** the arithmetic mean across observed taxa of each taxon's classified-photo white fraction `p_i = sum_j(y_ij)/n_i`.
- **Genus equal:** within each source-named genus, average those `p_i` equally across classifiable taxa, then average the resulting genus means equally. Nominal genus is the first space-delimited word of the archived source species name, as in the original FCP source analysis.

These statistics use **different weighting, not three independent data samples**. The within-genus count gate refers to the number of species with **at least one classifiable photograph**, never the true number of existing species in that genus. All source samples are non-probability observations.

We drew 1,999 whole-named-genus bootstrap resamples with seed **2026101450**, using the *same sampled genera for all three estimates* to retain their covariance. The percentile 2.5/97.5% endpoints describe the observed genus-cluster sample variability only. They are **not** phylogenetically corrected confidence intervals for the true global flora and do not remedy the many unclassifiable source photographs.

## Newly independently recalculated original-source results

| Weighting/denominator | Classified-photo source | White fraction | Whole-genus bootstrap 95% |
|---|---:|---:|---:|
| Photo equal | 45,953 photos | **46.948%** | 45.209–48.706% |
| Species equal | 25,064 classifiable taxa | **46.819%** | 45.245–48.350% |
| Genus equal | 5,369 classifiable genera | **49.891%** | 48.736–50.982% |

The **paired genus-equal minus species-equal difference is +3.072 percentage points**, source whole-genus bootstrap 95% **[+1.630, +4.602] pp**. This reflects weighting differences among the *sampled, classifiable* genera/species and cannot be attributed to natural selection, metabolic pigment costs, or ecological fitness equalization.

The original one-breadth-photo-per-nominal-taxon source alone gives **46.579%** classifiable photo white (18,457 classifiable taxa) versus **49.244%** genus equal (4,559 genera). Thus the weighting contrast is not solely produced by taxon×geographic-cell duplicate photos.

## Why near-50:50 changes with source genus richness

| Minimum classifiable species per genus | Remaining genera | Remaining species | Genus-equal white fraction |
|---:|---:|---:|---:|
| 1 | 5,369 | 25,064 | **49.891%** |
| 2 | 2,840 | 22,535 | 48.590% |
| 3 | 1,935 | 20,725 | 48.223% (bootstrap 95% 46.867–49.585%) |
| 5 | 1,185 | 18,195 | **46.626%** |
| 10 | 556 | 14,159 | 46.353% |
| 20 | 225 | 9,686 | 45.962% |

Of the 5,369 classifiable genera, **2,529 have a single classified species**, each receiving **equal genus weight** despite their minimal sampled species support. These single-classified-species genera average **51.351%** photo white. In contrast, the 225 genera with at least 20 classified species average **45.962%** white at genus-equal weight. The variation is consistent with an observed **sampling/taxonomic-composition effect**; the source does not establish that naturally species-poor lineages are more white or that biological richness causes floral colour.

The original all-hue global moisture manuscript analyses **20,694 taxa/1,934 genera** *with complete environmental features and ≥3 sampled congeners*. The source-balance population above is **25,064 classifiable taxa/5,369 genera without the same environmental filters**. Confusing these eligibility sets would create an erroneous mechanistic interpretation.

## Missingness and limits on equal-species prevalence

**17,047 of the 42,111 original nominal source taxa have no classifiable four-hue photograph.** If each of those taxa hypothetically had an unobserved scalar white-photo probability anywhere from zero to one, the mean across all **42,111 named taxa** could range **27.866%–68.347%** even before considering nonrepresentative source opportunities. For exact 50% across those original taxa, the currently unclassifiable taxa would have to average a hypothetical **54.677%** white. These are **logical missing-outcome bounds and a tipping assumption**, not empirical frequencies or confidence intervals; the classifier may have observed no biologically interpretable floral organ.

Among 25,064 classifiable taxa, **16,425 have only one classifiable photo** and **3,096 taxa have photographs assigned more than one of the four coarse colour labels**. Intra-taxon label differences may arise from organs, illumination, life stage or geographic phenotypes and must **not** be converted into genetic colour-morph polymorphism, even if the same taxon has multiple photo sites.

## Interpretation for the manuscript

**Supported:** the visually white category is close to half of classified opportunistic global flower photos, and the apparently nearly exact half arises under one particular equal-genus weighting. Restricting genus sample support reduces that percentage, while source genus-resampled variation clarifies precision for each descriptive estimand.

**Not supported:** (i) white occupies half of the world's true flowering-plant species, (ii) all plant clades or natural communities exhibit 50:50; (iii) rare white/coloured morph advantage maintains balanced genotypes; (iv) a pigment-energetic cost of coloured flowers cancels a pollination benefit; or (v) these ecological mechanisms can be determined from photo hue ratios.

A true **fitness trade-off** is the reversal in the *relative lifetime fitness* of independently confirmed morph/genotype classes across water/light/pollinator environments. A **frequency-dependent polymorphism equilibrium** additionally requires comparing relative fitness at experimentally controlled morph frequencies within the same population. The user-facing hypotheses and field experiment design remain in `docs/MOIS_WHITE_CHROMATIC_BALANCE_AND_FITNESS_HYPOTHESES_20261014.md`.

## Reproduction and rights boundary

- Aggregate-only machine result: `results/mois_original_photo_species_genus_balance_20261014.json`.
- Raw-ZIP replayer: `scripts/analysis/mois_replay_species_genus_balance_20261014.py`.
- The original photo ZIP is the FCP Actions artifact **11587115220**; it is **not in the public mois repository or any licensed reviewer archive**. Verify the SHA above and run:
  
  `python scripts/analysis/mois_replay_species_genus_balance_20261014.py --photo-zip raw_inputs/ORIGINAL_PHOTO.zip --output results/mois_original_photo_species_genus_balance_20261014.json --check`
  
  Dependencies: Python, numpy and pandas. This actual source replay and check both **PASSED** against the original historical ZIP in the working environment.
- GitHub CI without original source ZIP can check **committed aggregate values, synthetic weight invariants and source-script compilation**, but it cannot independently re-execute original source ZIP processing.

## Mechanism literature, strictly contextual

- Gigord, Macnair & Smithson (2001), **10.1073/pnas.111162598**: a manipulative test of rare-morph reproductive advantage in a rewardless orchid, not evidence of a universal global 50:50 ratio.
- Tenhumberg, Dellinger & Smith (2023), **10.1111/1365-2745.14057**: a mechanistic model of pollinator-mediated selection versus vegetative survival and pigment traits; not a parameter estimate from the present dataset.
- Narbona et al. (2025), **10.1038/s41598-025-94709-4**: plant flower pigment measurements show visibly white flowers need not lack UV-absorbing phenylpropanoids.

**Decision:** Keep taxonomic weighting and near-parity as a descriptive source-backed Supporting Information result, and keep balancing selection/pigment fitness crossing as explicit **unverified follow-up hypotheses**, not a new causal claim in the GEB Abstract.
