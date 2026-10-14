# Global half-white photographs mask marked nominal-genus clustering

**Original-source retrospective ecological analysis — 15 October 2026.** This is a supplementary, exploratory analysis of previously classified iNaturalist flower photographs, using original FCP source photo ZIP SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`. Two independently defined source-photo opportunity sets were audited. No flower pixels were reclassified, no genotypes or plant fitness were measured, and **no fitness compensation or balancing selection is established**.

## Research question

Does a worldwide white–nonwhite photograph fraction close to one half describe many genera that themselves exhibit intermediate colour proportions, or does it aggregate genera with very different photographed colour compositions? A near-half pooled mean cannot distinguish those ecological structures.

## Source-derived statistical estimands

Exactly 100,543 unique photo opportunities remain after removing 26,905 photo-ID overlaps between 42,111 breadth-photo and 85,337 taxon×cell records. The original four-hue label gate retains 45,953 classifiable photos, grouped into 25,064 nominal named taxa and 5,369 genera. Per taxon `i`, the outcome `p_i` is the classified-photo fraction labelled visible white, on [0,1]. The other three original visible hues are grouped as nonwhite for this secondary question. The genus is the first token of the original species name, not an externally validated species-tree node.

For genera with **at least three classifiable taxa**, compute the species-count-weighted **between-genus fraction of total taxon-level `p_i` variance**:

```text
R2_genus = [sum_g n_g * (mean(p_i in g) - mean(p_i))^2] / [sum_i (p_i - mean(p_i))^2]
```

This is a retrospective ANOVA-like descriptive fraction, *not* heritability, phylogenetic Pagel lambda or the fraction of colour variation causally explained by evolutionary history. We also count nominal genera with averaged classified-photo white fraction **<=10% or >=90%**; the latter are `near-pure photograph-hue genera`, not genetically monomorphic lineages.

The **global permutation null** reassigns each observed species-level `p_i` to a new genus while keeping *exactly the observed genus sample size structure* and the full distribution of source species colour fractions. We performed 1,999 frozen-seed permutations. An alternative **coarse source-conditioning null** reassigns species colour fractions only within source BIO12 rainfall-quintile × absolute-latitude (0–30°, 30–60°, 60–90°) × classified-photo-depth (1,2,>=3) blocks, including explicit missing-value strata. It holds those distributions fixed but **does not condition on true biome, complete plant phylogeny, local habitat, pollinator guild, observer traits or measurement bias**. The complete source was also independently checked using the **original breadth-only one-photo-per-nominal-taxon programme**, with its own permutation families and smaller eligible species/genera. All nulls are exploratory and were run after seeing source-level near-parity.

## Original-source results

| Source and randomisation | Classifiable taxa / genera (>=3 taxa) | Observed genus R2 | Null mean [2.5–97.5%] | Exploratory upper-tail P |
|---|---:|---:|---:|---:|
| Unique-photo union; unrestricted shuffling | 20,725 / 1,935 | **34.38%** | 9.33% [8.84–9.88%] | 0.0005 |
| Unique-photo union; rainfall/latitude/depth strata | 20,725 / 1,935 | **34.38%** | 9.62% [9.08–10.15%] | 0.0005 |
| Original one-breadth-photo programme; unrestricted | 14,639 / 1,529 | **33.86%** | 10.45% [9.82–11.11%] | 0.0005 |
| Original one-breadth-photo programme; conditioned | 14,639 / 1,529 | **33.86%** | 10.64% [9.97–11.32%] | 0.0005 |

There are **571/1,935 (29.51%)** photographed near-pure nominal genera among the source union's genera with >=3 classifiable taxa, compared with **7.70% expected** from global reassignment. Under a stricter >=5 classifiable taxon threshold, **277/1,185 (23.38%)** are near-pure, compared with a global null mean of **2.19%** and conditional mean **2.21%**; the two classes contribute **101 near-white** and **176 near-nonwhite** genera. Original breadth-only >=5 produces **21.79%** near-pure, compared with **2.27%** under global rearrangement. Every source-defined primary and sensitivity permutation upper-tail result is `P=0.0005` (minimum with 1,999 randomizations; not pre-registered programme-wide error control).

The source's descriptive **49.891%** white among 5,369 equally weighted sampled genera remains true, but it **does not mean that a typical genus is evolutionarily maintained at a half-white/half-coloured equilibrium**. Classifiability and taxonomic sampling are strongly unequal: 2,529 genera have only one classified taxon, which cannot inform within-genus variance. The >=3 and >=5 clustered cohorts change taxon support, and their mean photographed whiteness need not equal the unfiltered genus-equal fraction.

## Ecological alternatives and a discriminating next test

1. **Taxonomic/geographic species-pool sorting (supported as a pattern, mechanism unresolved).** The substantial nominal-genus clustering is consistent with lineage retention, habitat occupancy, and/or shared photographic classification conditions. Genuine phylogenetic signal requires a resolved species tree, colour measurement cross-validation and a phylogenetic comparative null; present genus ANOVA is **not** that test.
2. **Local community ecological sorting (untested).** At explicitly matched floral sites/plant assemblages, evaluate coexistence and floral trait dispersion conditional on habitat and phylogenetic distances. The present global genus memberships and coarse source coordinates do not establish flower-colour competition.
3. **Within-population frequency-dependent maintenance (untested).** Find independently morph-verified segregating populations, manipulate known white/coloured *individual* frequencies and measure pollen removal, pollen deposition, fruit/seed success and survival; estimate relative fitness as a function of morph frequency and environment. A worldwide tally of white versus nonwhite species is *not* a test of this.
4. **Stress/pollination benefit–cost trade-off (untested).** Manipulate water and radiation across genotypes/families, quantify VIS/UV floral chemistry, growth and lifetime reproductive fitness, then ask whether the relative fitness difference genuinely crosses zero between treatments. Nonwhite visible colour is not itself a directly measured energetic pigment cost.

Existing study context is mixed rather than a universal claim of conserved colour: McEwen & Vamosi (2010; [DOI 10.1098/rspb.2010.0501](https://doi.org/10.1098/rspb.2010.0501)) observed weak phylogenetic signal in spectrometrically measured flower colours in a limited subalpine flora. Gigord, Macnair & Smithson (2001; [DOI 10.1073/pnas.111162598](https://doi.org/10.1073/pnas.111162598)) experimentally demonstrated rare-floral-morph advantage within a deceptive orchid. Neither setting estimates the present worldwide photo-label genus R2 or validates global 50:50 balancing selection.

## Reproducibility, identifiers and current manuscript placement

- Source-only raw-photo ZIP replayer: `scripts/analysis/mois_replay_genus_clustering_20261015.py`, with `--photo-zip raw_inputs/ORIGINAL_PHOTO.zip --output results/mois_original_genus_clustering_20261015.json --check`. The reader checks original source SHA256 and ZIP integrity, deduplicates original photo IDs, preserves frozen colour labels, and redoes each of four 1,999-permutation families.
- Public **aggregate-only receipt**: `results/mois_original_genus_clustering_20261015.json`, including frozen seeds, permutation null means/quantiles and sample denominators. No named taxa, photo/observation IDs, precise geographic positions or photographer attribution are written. The original ZIP remains outside public GitHub and requires a license/position/privacy decision before any reviewer archive.
- Journal manuscript: succinct 3.1/Discussion mention; full measurement and null assumptions in **Supporting Information S2.5**; this memo and executable test fixture. No replacement of the original primary global BIO12/water-vapour signed effects, no biological evidence increment at pigment/fitness level.

**Scientific decision:** The new result is **strong source-conditional genus-level colour clustering**, not a stable white/coloured evolutionary balance and not an estimated pigment synthesis–pollination benefit trade-off. Treat pooled near-half and clustered genus composition as different levels of description; only independent morph/genotype and relative-fitness data could establish balancing selection.
