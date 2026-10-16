# Two taxonomic scales of the global photographed white–moisture association

**15 October 2026 · independent numerical analysis of the already used original FCP source tables.** Not new flower measurements, a new biological cohort or evidence of pigment-fitness compensation. Original photo and environment ZIP SHA256 were verified, original source photo IDs joined one-to-one, and no source identifiers or precise coordinates were released. The GEB main article's prior coefficient, geographic/source negatives, and licensing gates remain unchanged.

## The ecological question

Prior `mois` analysis found two patterns:

1. Globally, photographic white and coloured flowers are relatively close to half when assigned equal taxonomic genus weight; however, most sampled genera are strongly white-biased or nonwhite-biased, rather than internally half-and-half. On the 20,725 eligible taxa in 1,935 genus-sampled groups the earlier genus-clustering R² was 34.38%, versus an expected 9.33% under genus-size-preserving randomisation.
2. White photographed flowers are more likely among nominal congeneric species in wetter environments: original broad genus-equal BIO12 coefficient +2.255 percentage points per source SD, original within-scope maxT95 [+0.319,+4.190]; atmospheric vapour +3.952 [+1.199,+6.704].

**Can the second pattern be ascribed exclusively to genus turnover?** The original genus-centred effect already implies the answer is no for that source and geographical model. The new ecological comparison quantifies the *genus-mean* gradient **on exactly the same eligible species** and contrasts it to the archived congeneric response, with paired genus resampling.

## Matched source and models

Main climate26 four-hue source: 45,836 photo-site climate26-complete observations across 25,028 named taxa. Restrict each genus to ≥3 climate26-qualified taxa to keep **20,694 species, 1,934 genera, 38,498 classified images** in both estimates. For each taxon, white is the proportion of its classifiable original photos; all analysis weights are taxon- and genus-based, not raw photographic effort. BIO12 source SD is 736.415681 mm and annual water-vapour source SD 0.6755593 in original extracted units. Both models use the same three broad site-mean geographical covariates (absolute latitude and circular longitude sin/cos), genus-equal weighting and one original-source axis SD; neither models causal selection or a phylogenetic covariance tree.

- Within-genus: centre source covariates and white response by sampled genus and fit a broadly geographically adjusted white–moisture slope.
- Between-genus: regress equally weighted genus means in photographic white, moisture and broad geography with an intercept, again one source environmental SD.
- Paired uncertainty: same original 1,934 genus units drawn 1,999 times in shared multinomial bootstrap samples, **ordinary** percentile 95% CIs for each slope and between−within difference; original 24-contrast simultaneous maxT within-scope interval for the legacy within-genus result remains the principal source finding.

## Results

| Source matched target | White/BIO12 pp per source SD (95% bootstrap) | White/annual vapour pp per source SD (95% bootstrap) |
|---|---:|---:|
| Within genus among congeneric species | **+2.255 [+1.030,+3.563]** | **+3.952 [+2.134,+5.815]** |
| Between equally weighted genus means | **+8.555 [+6.220,+10.786]** | **+9.370 [+5.958,+12.893]** |
| Between minus within, paired | **+6.300 [+3.666,+8.779]** | **+5.419 [+1.724,+9.380]** |

The within-genus original full source coefficients were **exactly reconstructed** (numerical absolute discrepancy below 1×10⁻⁹ pp) independently of the historical regression code. For contrast, a pooled genus-equal species regression without genus centring gave +5.303 pp for rain and +6.163 for vapour on these same taxa; these are not arithmetic combinations of the scale-specific coefficients.

Retaining only genera with ≥5, ≥10 or ≥20 climate26-eligible taxa preserves the **positive point direction at both scales** for both water axes. For example, with ≥5 taxa/genus (18,164 taxa; 1,184 genera), within BIO12 = +1.756 and between BIO12 = +10.154 pp/source SD. Under ≥10, within = +2.822 and between = +12.870 (14,130 taxa, 555 genera). These are **descriptive subset coefficients**, not newly corrected significance families; changing support can also change spatial/environmental variation.

## Interpretation: what advances, what remains unknown

**Source-backed macroecological answer:** Wetter environments host a higher mean photographed white fraction **among sampled genera** *and*, within a sampled genus, species with wetter photo-site environments also tend to have higher white fraction. Therefore **genus membership replacement alone cannot completely describe** the positive world genus-centred moisture signal in this fitted observational source. A genus-mean association stronger than the congeneric association is compatible with geographic species-pool turnover, clade/habitat history and classification bias, but not direct proof that lineage turnover is the **causal mechanism** behind the pooled trend.

**Not demonstrated:** An adaptation of petals to rainfall, conserved pigment synthesis cost, UV protective trade-off, equal lifetime reproductive fitness, genetic frequency-dependent selection maintaining a white/coloured equilibrium, actual local community-level colour balance, or biologically distinct white morph alleles. Local geographic and stringent photo-quality tests remain imprecise, as already disclosed in the manuscript. The global near-half is still **a statistic of classifiable photo-labelled taxa under source-specific weighting**, not a fixed 50:50 evolutionary steady state.

**Sharper biological follow-up:** Verify independent reflectance/UV-pigment labels and source-independent occurrence sites on a phylogenetically resolved species sample. Test whether the genus-mean association persists with lineage history and independent soil/light/pollinator variables, and separately whether the within-genus moisture gradient is consistent across clades, floristic regions and image standards. Only in confirmed polymorphic populations should one manipulate colour-morph frequencies and environmental treatments, measuring pollen export, pollen receipt, seed/offspring fitness and survival across years. Species turnover and within-genus sorting are ecological patterns; **relative fitness crossing requires direct reproductive or survival data**.

## Reproducibility and publication position

- Original source ZIPs: photo SHA256 `2d6097ba…` and environment SHA256 `6e287294…`, recovered and kept outside the public `mois` repository pending photographer licensing and source-location/geoprivacy clearance.
- Full original-source code: `scripts/analysis/replay_mois_genus_between_within_moisture_20261015.py`. Original ZIPs are required for `--photo-zip ... --env-zip ... --output ... --check`; exact original source replay and original 1,999-draw percentile result passed locally.
- Public **aggregate-only** numeric receipt: `results/mois_genus_between_within_moisture_20261015.json`. No source photo IDs, taxon names, per-genus rows, attribution, coordinates, or individual source geographic moments included.
- GEB main working manuscript: [Results 3.2 and Discussion 4.1](MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md). Detailed inference and support gates [Supporting Information S2.6](MOIS_GEB_SUPPORTING_INFORMATION_20261014.md).
- **Draft PR not submitted.** Publisher-accessible licensed reviewer archive/DOI, coauthor approval, anonymized source access and final full-resolution figure/cover-letter proof still unresolved.

**Scientific decision:** Keep the interspecific positive moisture–white association as the main GEB-scale result. Report the new identical-species genus-mean versus congeneric comparison as a test of **taxonomic-scale ecological heterogeneity**, not as evidence that the two slopes mathematically partition a climate-selection process.
