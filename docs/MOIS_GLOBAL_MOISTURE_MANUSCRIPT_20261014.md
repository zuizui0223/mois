# Moisture-associated floral colours across the worldwide photographed flora: taxonomic and spatial limits

**Global interspecific study — transferred editorial draft, 14 October 2026.** This is a research manuscript *draft*, not a submission or independent replicated biological experiment. All source measurements/bootstraps were performed in the historical FCP branches and source-locked; we do not rename retrospective evidence prospective. This manuscript excludes the FCP high-depth intraspecific polymorphism / H2 white-axis study.

## Summary

Global flower-colour datasets often assign one colour per species, but a photographic opportunity atlas contains substantial within- and among-species sampling heterogeneity. We ask whether white and chromatic photographed flower states correspond to precipitation and atmospheric water vapour among congeneric taxa worldwide, whether this relationship generalizes geographically, and whether its direction survives measurement and model limitations. An original-source union comprised 100,543 unique images of 42,111 named species; 45,953 could be assigned one of four visible flower hue states. In 20,694 species from 1,934 qualified genera, white photographic colour increased with BIO12 annual precipitation by +2.255 percentage points per source SD (within-family maxT 95% interval [+0.319,+4.190]) and mean atmospheric water vapour by +3.952 ([+1.199,+6.704]). Strict image quality and within-genus×coarse-geographic-cell comparisons did not confirm the signed precipitation effect; an identical-species paired geographic conditioning comparison also had a multiple-testing interval spanning zero. Original photographic colour classifiability, genus history and site sampling remained important alternatives. Thus the worldwide photographic flora shows a **broad interspecific moisture–colour association under its observed source distribution**, not an identified common pigment-response, natural selection or morph-maintenance mechanism.

## Introduction

Floral colours can vary across plant lineages, climates and local communities. A broad photographic association between floral colour and annual precipitation could reflect many distinct processes: historical species turnover, physiological differences, pollinator environments, photo-organ availability and environmental or photographic selection. We therefore treat **named-species photographed hue proportions**, not genetic polymorphic morph frequency, as the observable response. The objective is to estimate and aggressively test interspecific white–chromatic moisture relationships rather than infer adaptation from a map.

A separate FCP project investigates **within-species** photographic flower-colour distributions and the recurrent white-versus-nonwhite direction W; it is deliberately excluded from this manuscript as a primary result. The high-depth W result does not establish the environmental mechanism for the global association.

## Materials and Methods

### Source, sample support and colour measurement

Source opportunity collection: iNaturalist Research Grade flowering photographs with original photo/species/observation/site identifiers; the deduplicated archive contains 100,543 unique photo IDs for 42,111 nominal named taxa. An original fixed visible-light flower detection/segmentation pipeline assigned coarse white, yellow/orange, red/pink or blue/purple to 45,953 images; 54,590 nonclassifiable outcomes remain **missing**. The colour label is an RGB-photo descriptor, **not a reflectance-spectrum or biochemical pigment assay**. The original photo archive SHA256 is `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`.

A checksum-verified original-photo identifier comparison with the separately sampled FCP high-depth within-species colour study (149,900 photo IDs) found **3,411 shared original photo IDs**: 3.39% of this global MOIS photo archive and 2.28% of FCP's high-depth observations. The overlaps comprised 1,248 FCP discovery, 1,326 validation and 837 prospective photographs; none of the three FCP cohorts shared photo IDs with each other. Four-colour labels agreed for all 1,693 shared photos classifiable in both measurements. This describes a *shared photograph and measurement system*, not independent external colour validation. The archived ID-level matching report and exact source SHA256 values are preserved in `results/fcp_mois_original_photo_id_overlap_20261016.json`; individual photo identifiers are not included in the published audit summary. These two papers test distinct *between-species* and *within-species* estimands, and do not treat one another as independent biological replicates.

Climate26 nonsoil predictors are complete for 45,836 classifiable photos / 25,028 species, of which 20,694 /1,934 genera satisfy the fixed ≥3 named species per nominal genus signed-comparison rule. Soil36 is complete for 31,925 photos /20,786 species, with 16,801 /1,665 genera qualified for analogous hue regressions; these are **different complete-case supports**. The environmental archive SHA256 is `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`.

### Inter-specific hue and environment

Responses are species photo-colour proportions (four categories summing to one); covariates include BIO5, BIO12, BIO14, BIO15, annual solar radiation, mean atmospheric vapour pressure (plus an available-water soil proxy where complete). These six/seven ecologically interpretable axes were selected retrospectively. Within-nominal-genus centring and adjustment for absolute latitude and sine/cosine longitude isolate source-conditional gradients; slopes are in photographed colour probability percentage points per original-global environmental SD. Fit species-equal and genus-equal weighting in distinct, verified scopes. Within each of 20 geographic/photo/scope combinations, 999 whole-genus bootstrap draws and approximately simultaneous maxT intervals adjust that scope's 6/7×4 hue contrasts, but not the history of model and region selection.

To challenge the broad geography, a separate same-genus×same-original coarse spatial cell comparison requires ≥3 congeneric named species within each relevant group. A matched A/B/C stage further distinguishes changes in species support from changes in spatial conditioning; within its B/C comparison *the original named species and their measured photo-colour/environment values are identical*. This is model conditioning, not a causal geography intervention. Independent NHM botanical colour descriptions were considered a distinct phenotype measurement, but their original geographical site source still overlaps FCP.

### Reproducibility

Original code, paired spatial analysis, synthetic tests, independent original-row WLS checks and source-result CSV/JSON were **copied intact** from original FCP #193/#194. The original source-verified Actions run and file SHA fingerprints are linked in the provenance manifest. No original photographs were reclassified by this repository transfer.

## Results

### World-scale association, with interval support

For the genus-equal global climate26 panel, photographed-white versus annual precipitation BIO12 yielded **+2.255** percentage points per one original source SD with maxT 95% **[+0.319,+4.190]**; mean atmospheric vapour pressure yielded **+3.952** with maxT 95% **[+1.199,+6.704]**. Global soil36 under its separately complete panel yielded BIO12 white **+2.257** ([+0.138,+4.376]). The global photo source was 45,953 four-colour classifiable observations among 100,543 unique opportunities, so all results condition on measured hues.

Strict photographic quality yielded BIO12 white **+1.229** (maxT95 **[−2.046,+4.504]**) and same-genus×same-coarse-cell yielded **+0.808** (**[−3.556,+5.172]**): neither signed direction was confirmed within those source populations.

### Same-species comparison exposes spatial-support dependence

The paired spatial support audit fixes 9,370 original species/845 genera in its world B/C comparison and 2,886 species/280 genera in the East 90–180° E. White–BIO12 under broad geographic conditioning is +1.326 (world), +1.102 (East); under genus-by-cell conditioning +0.808 (world), −2.076 (East). East B−C +3.178, paired within-family 24-contrast maxT95 **[−0.258,+6.613]**, includes zero. This rules out treating the nominal post hoc East point estimate as a robust independent climatic cause.

### Small predictive gains and independent-hue nonconfirmation

The climate26 genus-equal geographical heldout Brier gain was +0.002943 on 17,183 test species and soil36 +0.001432 on a **different** 13,961 test species. They cannot be subtracted to infer soil importance. An original, separately executed same-soil-cohort unseen-genus comparison estimated incremental soil gain +0.0000757 (95% [−0.000516,+0.000653]) on identical 20,786 source species: an inconclusive added soil effect in that design. An independent NHM floral-description comparison of 670 photo-site-matched East named taxa/58 genera estimated +3.56 pp (95% [−15.34,+22.60]), without a convergent directional conclusion.

## Discussion

The broad congeneric photo-white association with annual rain and vapour is a reproducible **historical source-data/model output** with bounded generalizability. Stronger photo-quality, more local geographic conditioning, alternative botanical hue and geographic transfer leave a shared species-independent water–colour response unestablished. The source photographs are selected for colour classifiability and regional observation opportunity; climate variables are highly correlated, and neither photometric pigments nor genotype/relative fitness were observed.

Across the original-classifiable photos, 21,574/45,953 (46.95%) were called white. That descriptive source proportion is not an evolutionary 50:50 fitness equilibrium; 54,590 photo-colour outcomes were missing and source species-conditioned latitude/compensation nulls did not establish adaptive stabilization. Same-original-photo interspecific-versus-intraspecific rainfall contrasts, archived in the source receipts, further bound extrapolation from average species photo-colour to individual or population environmental change.

The fcp high-depth intraspecific D, W and spatial tests answer a **different biological question**. Their recurring white–chromatic phenotype-space direction is not independent confirmation of a moisture effect and should not be pooled with genus-centred macro associations. Experimental floral spectroscopy, pigment chemistry, common-garden water treatment and reproductive fitness would be needed for the proposed benefits/costs of pigmentation.

## Data, code and limitations

Canonical transferred code/results/figure/data files are listed in `docs/MOIS_ORIGINAL_FCP_PROVENANCE_AND_TRANSFER_20261014.md`. Historic input zipped archives are **not stored here**; reconstructing them requires the preserved original FCP source workflow/release identity and original SHA256 checks. This paper remains **draft** and requires graphics, external source availability validation, authorship and editorial approval before submission.

**Limitations that may not be dropped:** missing-classifier MNAR, photo-white near-clipping, coarse flower/background misclassification, complete-case changes between 26/36, conditional genus rather than dated-tip phylogenetic control, spatial nontransport, and no pigment/metabolic/reproductive experimental measurement.
