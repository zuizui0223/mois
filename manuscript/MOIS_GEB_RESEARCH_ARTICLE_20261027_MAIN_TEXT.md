# Global floral-colour associations with climate depend on taxonomic composition and geographical support

**Submission-oriented Research Article draft | 27 October 2026 | Provisional target: Global Ecology and Biogeography**

**Status:** Complete first-pass narrative and source-identified main results, **not submitted**. The data stem from different, partially overlapping SHA-frozen, original-photo cohorts and must not be represented as one giant re-fitted climate model. Double-anonymous manuscript version, coauthor approval, final independent source rights, reference check and finalized figures remain outstanding.

**Running title:** Context dependence of global floral colour

## Structured abstract

**Aim:** Associations between flower colour and climate may reflect environmental filtering, geographical differences among plant assemblages, or the taxa and photographs that observers sample. We asked whether moisture-related variation in photographed floral colour is globally consistent, and whether apparent colour differences among elevations, islands and seasons remain interpretable across their distinct sampling supports.

**Location:** Worldwide photographed flowering-plant assemblages.

**Time period:** Historically heterogeneous community-science observations; long-term WorldClim 1970–2000 climatic normals for the annual environmental models.

**Major taxa studied:** Named flowering-plant taxa represented in source-audited iNaturalist photographs; 37,107 nominal taxa with classifiable four-hue photographs in the 350,399-photo descriptive cohort.

**Methods:** We synthesized separately authenticated original-image studies: a 309,439-photo six-climate-axis analysis with whole-genus resampling; a 350,399-photo four-hue and seasonal composition analysis; geographically matched elevation and island–mainland photographic assemblage comparisons; and a later 1,181,146-photo, four-acquisition-arm audit of nine-component colour measurability. We explicitly separated photo-, taxon- and genus-weighted results and retained technically unclassifiable images as missing outcomes.

**Results:** White comprised 47.49% of classifiable four-hue photographs, versus 48.19% and 51.46% when sampled taxa and genera were weighted equally. Photographed-white frequency was positively associated with annual precipitation (+1.99 percentage points per source SD, simultaneous 95% interval +0.55 to +3.42) and atmospheric vapour (+3.16, +1.04 to +5.28), but geographically conditioned intervals included zero. In 44 nearby elevation comparisons, white was 5.09 points lower at high elevation, whereas an environmentally matched island–mainland difference (+1.40 points) was uncertain. Southern summer photographs were 10.50 points more often white-labelled than winter photographs, with large contributions from changing photographic weights of taxa.

**Main conclusions:** Broad climatic and topographic gradients coexist with taxonomic and photographic assemblage structure. They do not imply one universal locally transferable floral-pigmentation response, or an inferred pigment fitness trade-off. Separately audited original-image nonresponse places important limits on generalization.

**Keywords:** floral colour; flower colour; macroecology; climatic moisture; taxonomic composition; community science; geographic matching; seasonality; photographic nonresponse.

---

## 1 | Introduction

Flower colour is an ecologically conspicuous reproductive trait. Its geographical pattern can arise from pigmentation, petal optics, inherited or developmentally constrained trait differences, pollinator interactions, and the climatic and biogeographic distributions of the plants that bear it. Visible floral colours have been associated with both abiotic environments and biotic communities in previous macroecological studies (Dalrymple et al., 2020). Nevertheless, a worldwide relationship between photographically visible white flowers and annual rainfall does not establish that exposure to more rainfall induces or favours white flowers within the same population. Broad-scale colour patterns can emerge when communities, lineages and the locations at which they are observed differ systematically across environments.

Three forms of scale dependence require particular attention. First, the world-wide proportion of white photographs need not equal the proportion of white-flowered taxa or botanical individuals. Frequently photographed species, heavily represented genera and easily segmented floral objects contribute disproportionately to opportunistic imagery. Second, moisture, elevation, temperature and geographical position covary. Positive global associations may weaken when comparisons are restricted to taxa photographed within geographically comparable environments. Third, the flora visible in photographs changes over time: summer and winter source photographs can represent different taxa and markedly different photographic weights among taxa shared across seasons, even if no individual plant changes flower colour. Community-science observation patterns are known to be spatially, temporally and taxonomically uneven (Di Cecco et al., 2021).

We address these issues using genuine original photographic image measurements and multiple source-frozen ecological analyses. Rather than treating every exploratory contrast as independent evidence of selection, we ask three questions. **(1)** Does a broad positive white–moisture association persist in an enlarged photographed flora, and how dependent is its direction on geographic conditioning? **(2)** Are apparent floral colour shifts shared across distinct settings—higher elevations, islands and local seasons—or do their effects and sampling support differ? **(3)** How much can differences in observed plant taxa and image-classification success alter the meaning of a worldwide colour composition?

We use *photographed floral colour* throughout to designate operationally classified visible flower-like regions, rather than independently confirmed corolla reflectance. White and the three chromatic categories refer to visible-sRGB classifier labels, not absolute absence or presence of anthocyanins or other pigments. This deliberately limited phenotype definition makes the ecological estimands explicit and avoids silently extrapolating to experimentally measured floral colour, reproductive fitness or evolutionary costs.

## 2 | Materials and methods

### 2.1 | Original image sources and units of analysis

Our analyses use source-audited iNaturalist/FCP original photo identifiers, recorded observation identities, source taxa and georeferences. A frozen upstream source registry contains 93,209,351 photographic opportunity records mapped to 42,111 nominal target taxa; these records constitute an opportunity frame, **not 93.2 million measured flowers**. The ecological models and descriptive comparisons were completed on several actual-image subsets under shared but not identical eligibility rules.

The original four-hue global moisture model was re-fitted after measuring **309,439 independent original photo IDs**, of which **149,662** had a classifiable four-hue label. We summarize colour frequencies and actual observation-date composition in a separate, subsequently enlarged **350,399-photo** source, with **164,297** four-hue classifications. Regional elevation and island analyses use the independently signed **313,535-photo** feature-lake source and aggregated, species-equal geographical cells. The 13-covariate species-and-genus comparison is a separate soil-complete subset of **16,801 named taxa and 1,665 nominal genera**, drawn from older original photographic analysis resources; its numerical in-sample explanatory fractions are treated as secondary context, not refitted on the 350,399 or later 1.18-million-photo frame.

To audit phenotype availability independently of the fitted environmental models, a later four-acquisition-arm union reconciled **1,181,146 distinct genuine original photo IDs** with **451,643 valid strict nine-component original ROIs**. It preserves the heterogeneous 1,043,930-photo historical/FCP/K1 source and the separate 2,048, 4,096 and 131,072-photo replenishment arms. The large-source geo- and observation-date fields were corrected from 512 SHA-signed source work units after detection of missing source metadata, without changing any previously measured image colour. For analyses of geographic nonresponse the corrected geography release, not the preceding geography-incomplete export, is used. This expanded union supplies a data-quality assessment, **not updated environmental regression coefficients**. A distinct 64-photo 500-pixel flower-box pilot is omitted from the strict original-ROI denominator.

The principal identity unit in all image analyses is a unique original photo ID. Photos sharing an original observation UUID are not separate ecological visits. Models and summary comparisons average by source taxon or by geographical cell where described; a photo-level denominator must not be confused with the number of independent botanical individuals or taxa. Sampling arms selected specifically to enrich unmeasured taxa remain distinct from source-reference cohorts.

### 2.2 | Image-derived floral colour and technical missingness

The frozen FCP original-image ROI-v4 workflow uses a YOLO11 flower detector and EfficientSAM segmentation, followed by a visible nine-category colour composition and admissibility gates. Operational four-hue labels are **white**, **yellow/orange**, **red/pink** and **blue/purple**; an image whose ROI or colour is unsuitable remains **unknown**. The nine-category fractions (white, yellow, orange, red, pink, magenta, purple, blue and bronze) are distinct from four mutually exclusive photo labels. The four-hue classifications and the stricter nine-component admissibility criterion have different definitions and denominators, even when they derive from the same source photo.

For descriptive floristic contrasts we use each classifiable photo's four-hue label, and compare photograph-equal, taxon-equal and genus-equal white frequencies. For the 1.18-million-photo quality audit, every source image remains in the original attempted-photo denominator and a strict nine-component ROI is either valid or technically missing. Technically missing colour is **never** interpreted as a biologically white or chromatic flower or as absence of bloom.

A separately developed 500-pixel YOLO flower-box/CIE-L*a*b* measurement is an auxiliary detector/colour-box proxy, potentially contaminated by foliage and background. It is not numerically interchangeable with the full-resolution segmented original image's nine-component flower ROI. In particular, 500-pixel detector negatives are not used to remove original-image records from the phenotype denominator.

### 2.3 | World-wide moisture association and geographic conditioning

The primary annual-climate analysis uses six source-prespecified long-term variables: maximum temperature of the warmest month (BIO5), annual precipitation (BIO12), precipitation of the driest month (BIO14), precipitation seasonality (BIO15), mean annual solar radiation and atmospheric vapour pressure. The original source-specific observation-site standard deviations are the effect units. These gridded climatic normals do not represent weather measured during a photo's own month or year. The original solar grid encodes kJ m−2 day−1; physical conversions to W m−2 require dividing by 86.4.

We retain the existing genus-weighted original source estimator, the geographically conditioned genus × broad geographical-cell comparison, and the identical pre-frozen 999 whole-genus resamples. The reported intervals are simultaneous 95% maxT bounds **within a 24 hue × environment comparison family for each model**, not a project-wide correction for all exploratory hypotheses. Positive coefficients indicate increased percentage points of the white photographic label per one source environmental SD among the taxa meeting the comparison's support requirements. We report globally pooled and geographically conditioned slopes as distinct estimands, not as a causal treatment-control pair.

The original smaller climate-complete source also permitted within- versus between-genus slopes estimated on **the same 20,694 nominal species and 1,934 genera**, after the original broad geographic controls. We use these as supporting source-frozen scale diagnostics; their environmental standardization and sample differ from the later 309,439-photo fitted model. The historical soil-complete 13-predictor weighted least-squares analysis is included as an in-sample secondary check of overall environmental association and observed taxonomic clustering, not a cross-validated causal partition.

### 2.4 | Elevation and island geographical contrasts

An independent 313,535-photo geographic feature-lake analysis summarized **species-equal photographic hue proportions within one-degree cells**, requiring at least five sampled taxa and ten classifiable photo labels per cell. Its matched island comparison initially included 177 island and 1,217 mainland cells, and paired islands and mainland cells within the same 10-degree latitude × 30-degree longitude coarse stratum on standardized warmest-month temperature, log precipitation, annual solar radiation and elevation. The primary one-to-one environmental caliper of 1 yielded **102 pairs across 65 island landmasses and 28 source regions**. Pair-equal white differences were evaluated with nominal region-cluster intervals; additional calipers, equal-region weighting and island area restrictions were sensitivity comparisons. The matching design does not require the same species on an island and mainland, and source island assignment does not establish a sampled plant's native island status.

For the elevation diagnostic, low-elevation cells (mean below 500 m) and high-elevation cells (mean at least 1,500 m) were compared within the same landmass and 5° × 5° geographical neighbourhood. **Forty-four neighbourhoods** met the dual-elevation support criteria. Contrasts compare photographed cell assemblages, not the same species ascending a mountain or independently validated alpine habitats. Their nominal neighbourhood-based intervals do not jointly account for every climate or phylogenetic confounder.

### 2.5 | Photographic local seasons and taxonomic composition

The separate **350,399-photo** descriptive source retains original observation dates. Summer is June–August in the north and December–February in the south, with opposite definitions for winter. Source records within ±23.5° latitude are excluded from the extratropical summer–winter comparison. We compare white-labelled proportions among classifiable original photographs and decompose an unpaired summer-minus-winter difference exactly into: **(i)** a photographed-white frequency change for taxa represented in both seasons, **(ii)** a change in the relative source-photographic frequencies of those shared taxa, **(iii)** taxa contributing classifiable photo support only in summer, and **(iv)** taxa contributing only in winter. This is a symmetric algebraic accounting identity, not a causal phenological mediation model. Taxa with no photos in a season are not deemed biologically absent or colourless.

To assess whether the unpaired contrast identifies a response within a taxon and place, we audited exact taxon × 0.25° source cell × local meteorological season-year groups with independently identified source observations in both seasons. Multiple photos from a single observation do not meet replication thresholds. We report the available matched support rather than filling missing taxa, locations or seasons with presumed colour states.

### 2.6 | Observation and selection-bias checks

To evaluate whether colour patterns depend on which taxa were photographed, the 350,399-photo source compares its original 317,631-photo core with 32,768 independently source-selected photographs enriched toward previously underrepresented taxa. We distinguish a raw photo-weighted arm difference from matched genus × geographical-cell estimates and report how many source taxon × cell records enter each matching procedure. Source selection is outcome-conditioned in some supplementary arms; therefore, coefficients and white fractions from these cohorts are not freely pooled to estimate world-wide biological abundance.

The newer 1,181,146-photo strict9 availability audit additionally summarizes classification rates by latitude band and local season, separately retaining original source design. A selected 512-image diagnostic subset of already-unclassifiable strict9 records, spanning 512 nominal source taxa, was evaluated for 500-pixel flower-box availability. The subsample was purposefully balanced by source arm, locality, season and unmeasured taxon, so detector yield from these images is **not** a prevalence estimate of all unclassifiable photographs. An independently prepared 128-photo human flower-organ review has blank annotations at the time of this draft and is **not** reported as validated botanist ground truth.

## 3 | Results

### 3.1 | A near-half photographic white fraction changes with taxonomic weighting

In the **350,399-source-photo** descriptive cohort, **164,297** photographs were successfully assigned one of four visible hue categories, while **186,102** attempts remained unknown or unclassifiable. Among the classifiable images, the counts were **78,026 white**, **51,491 yellow/orange**, **17,382 red/pink** and **17,398 blue/purple**, making white **47.49%** of eligible photographed colour labels. When equally averaging named sampled taxa rather than their photo counts, the white fraction was **48.19%**; equal weighting of nominal genera yielded **51.46%**. The approximately four-point difference between photo and genus weighting shows that the photographic hue composition depends on taxonomic sampling depth. These proportions are conditional on the four-hue classification process and should not be described as a world-wide 50:50 selective equilibrium.

A further check of source recruitment supports this sensitivity. Relative to the classifiable core photographs, a deliberately taxon-broadening supplementary arm had a white-labelled photo fraction **3.78 percentage points higher**. Restricting to common sampled nominal genus × 10° geographical cells and equal-weighting genera reduced the *source-arm contrast* to **+0.64 points** (exploratory genus bootstrap interval **−1.85 to +3.07**). The raw and matched contrasts refer to different sampled support and different weights, and hence the reduction is **not** an identified causal taxonomic-turnover effect. Nevertheless, their divergence makes the composition of sampled floras a practical determinant of estimated white-photo frequencies.

### 3.2 | Moisture covaries with white photographs globally, but does not establish the same local direction

The completed **309,439-photo** original-image six-axis model retained a positive global photographed-white association with **annual precipitation: +1.99 percentage points per original source SD**, within-model 24-comparison simultaneous 95% interval **+0.55 to +3.42**. Mean annual atmospheric vapour showed **+3.16 points (+1.04 to +5.28)**. These source-specific estimates are consistent in direction with the preceding 295,017-photo version, but do not demonstrate within-population selection.

On the corresponding **genus × coarse-geographical-cell support**, the precipitation point estimate was **+0.90 points (−2.63 to +4.44)** and vapour **+2.15 (−2.64 to +6.95)**. Both intervals included zero. A worldwide genus-weighted direction therefore coexists with uncertainty about its generality when taxa are compared within more restricted geographical support.

Earlier, source-frozen and independently fitted between/within-genus comparisons offered a compatible perspective on scale. On the **same 20,694 nominal taxa across 1,934 genera**, the broad geographically adjusted photographed-white rainfall slopes were **+8.56 points between sampled genus means** and **+2.26 within nominal genera** per that analysis's source SD; the paired between-minus-within difference was **+6.30 points** (ordinary whole-genus-bootstrap 95% interval **+3.67 to +8.78**). These are neither direct genetic phylogenetic effects nor observations of one plant species changing its pigment; they demonstrate that colour-associated climatic sorting occurs on more than one sampled taxonomic scale.

An additional complete-case, **16,801-taxon** exploratory 13-environment regression had in-sample white-photo R² of only **2.00%** under equal-genus weights (**1.28%** equal-species weights). Fitted genus means on the same source had much larger in-sample grouping R² (**44.32%**), but a 1,665-genus high-parameter grouping model is intrinsically flexible (size-preserving permutation mean **19.06%**). These values describe observational structure under unequal model complexity, **not** the percentage of flower colour heritability explained by lineages versus climate.

### 3.3 | Higher elevation differs from matched islands rather than sharing one white-colour syndrome

The geographically paired elevation comparison encompassed **44 nearby landmass neighbourhoods**. The high-elevation side had **5.09 percentage points less white photographic composition** than the low-elevation side (nominal 95% interval **−8.79 to −1.39**), and **31 of 44** local comparisons had the same direction. Relative chromatic photographic categories increased jointly on the high-elevation side, but this is not evidence of greater measured pigment concentration or a specific protective anthocyanin mechanism.

By contrast, the same broader source's **102 environmentally matched island–mainland pairs** yielded an island-minus-mainland white difference of only **+1.40 points**, nominal region-cluster interval **−1.55 to +4.35**, with equal-region weighting close to zero. A simple island–mainland white excess observed without matching therefore should not be elevated to a universal island syndrome. No island-specific pollinator shortage, endemic or native origin was measured in this comparison.

Both setting analyses are assembled-photograph comparisons, not independent experimental tests of identical populations; the contrast shows that “more extreme environmental conditions” cannot be assigned one universal predicted direction of white floral-colour composition without specifying the ecological and geographic context.

### 3.4 | Seasonal white-photo differences are strongly dependent on photographed taxa

Among classifiable original photographs outside the tropics, Southern summer was **10.50 percentage points** whiter-labelled than Southern winter; the equivalent Northern difference was **+2.83 points**. In the South, **+8.04 points** of the exact decomposition were assigned algebraically to changing relative photographic weights of taxa observed in both seasons. The within-shared-taxon white-photo fraction term was **+3.57**, while the summer-only and winter-only taxon-photo support terms had opposite directions (**+4.16** and **−5.27**). These values precisely reconstruct the unpaired source difference but do not isolate an underlying biological flowering-community shift from seasonal observer behaviour.

Stringent within-locality temporal replication was scarce. For the same named taxon, 0.25° cell and local meteorological season-year, only **two Northern and two Southern matched taxon-cell-years** had at least **two independent original observation IDs in each season**. Consequently, neither the large Southern photographic summer-minus-winter contrast nor its component terms are evidence that individual flowers plastically change colour with seasonal temperature or light. It remains a strong empirical illustration of how phenology-associated **photographic composition** can change broad-scale floral-colour estimates.

### 3.5 | Technical colour missingness has large and geographically structured support

The separate quality-assessment union contained **1,181,146 genuinely attempted original photo IDs**, of which **451,643** were admissible strict nine-component ROIs (38.24%) and **729,503** remained strict9 unknown. Of the frozen 42,111 nominal target taxa, **30,048** had at least one strict9 image; **12,063** remained unsupported at this stricter colour threshold. These are operational source taxa and photo measurements, not the true number of flowering species globally.

After recovering original locations and observed dates for the 131,072 previously geometadata-incomplete bulk-shard photos from SHA-pinned source work units, true strict9 classification success remained heterogeneous: **42.0% for Northern temperate**, **33.9% Southern temperate**, **28.8% Northern tropical** and **23.8% Southern tropical** original photo attempts. The restored world-wide audit had only **195 coordinate-unknown** and **197 local-season-unknown** original photos. Because success varies across study strata, the classifiable subset's composition cannot simply be assumed to represent the unclassified photographs.

In a deliberately stratified existing-photo diagnostic sample of **512 strict9-unknown original images**, a 500-pixel YOLO flower-box detector returned **495 candidate boxes**, 14 negatives and three technical errors. This does **not** establish that 495 truly depicted open petals or could all be rescued by a changed segmentation gate. It does establish that the detector-negative category alone cannot plausibly be used as a sufficient explanation for the ROI9 failure of the selected photographs. The botanical truth annotations for the separately frozen 128-photo review remain pending. No inferred white or coloured label was written to any technical unknown.

## 4 | Discussion

### 4.1 | Broad ecological associations are real photographic patterns, not universal local laws

The strongest consistent environmental result is a positive global association between moisture and photographed-white composition. This supports the ecological relevance of water-related gradients in broad-scale source assemblages, while the simultaneous intervals that cross zero after geographical restriction explicitly limit the claim to the support on which it is measured. The association can reflect physiological selection, community-level environment filtering, divergence among sampled taxa or observational visibility; the present evidence does not identify their relative causal contributions.

The high-elevation white deficit supplies a separately structured contrast, but does not indicate that altitude acts independently of temperature, ultraviolet radiation, moisture or biotic assemblages. The matched island result, conversely, weakens a simple rule that geographical isolation universally increases white flowers. Taken together, these patterns call for hypotheses formulated at explicit taxonomic, regional and source-support scales rather than a globally homogeneous colour-response coefficient.

### 4.2 | Floristic composition matters across environmental and temporal comparisons

Changing the unit of representation from photographs to taxa or nominal genera shifts the estimated overall fraction of white photographic ROIs. Broad between-genus moisture gradients are stronger than within-genus gradients in an earlier identical-taxon comparison. Seasonal photo-white contrasts also contain prominent relative-weight and season-specific taxonomic support components. These are different endpoints, but they share one macroecological implication: photographed trait distributions are jointly structured by environmental context and the identities/weights of organisms included in the observation.

The seasonal decomposition is an exact accounting of the photographed source rather than direct evidence of phenological trait plasticity. Likewise, source-arm genus/geography matching demonstrates sensitivity of the white-label denominator, but the matched strata are a small and selected overlap population. Taxonomic composition is thus important, but it has not been isolated as the only cause of the white–climate associations. Genus membership here is not independently calibrated molecular-phylogenetic relatedness.

### 4.3 | Detection and classification are part of a trait-macroecological observation process

Large image numbers improve taxonomic and geographic opportunity, but cannot automatically make classifiable photographs representative. Our modern 1.18-million-photo quality audit found more than three-fifths of genuine original images without admissible strict nine-component ROIs, with lower success in the sampled tropics. A future paper claiming an absolute global frequency of naturally white-flowered taxa would need phenotype validation, defined taxonomic sampling probabilities and explicit sensitivity to missing classification. This article instead consistently states the narrower, supported target: observed colour composition **conditional on photographed source inclusion and colour measurability**.

The contrast between strict9 ROI success and a broad flower-box candidate is particularly informative. A box can detect a potential flower-like region and still fail to locate evaluable corolla pixels, while a detector negative does not establish that a photo contains no flower. Therefore, neither image absence nor pigment class may be imputed from detector status. The purposefully failure-enriched 512-photo diagnostic is an audit of likely technical failure modes, not an independent probability sample of source photographs.

### 4.4 | Biological mechanisms and predictive scope

The observed results are consistent with several, potentially simultaneous explanations: water and irradiance regimes may affect optical/pigment functions; pollinator assemblages may change geographically; taxa differ in their ancestral colour tendencies and geographic ranges; and colour-dependent image acquisition and segmentation may modify observed class proportions. These mechanisms have different units of intervention and cannot be distinguished without dedicated flower-organ validation, phylogeny, pollen/pollinator/fitness data or source-aware repeated within-population observations.

The modest explanatory R² of the 13-variable ecological comparison and geographically uncertain moisture contrasts caution against using the global slope as a local species-colour predictor. The contribution of this study is not to identify an ultimate cost–benefit equilibrium maintaining white and coloured flowers, but to **define what macroecological floral-colour signals are reproducible at which observation, taxonomic and geographic scales**. This distinction matters whenever dense opportunistic occurrence records are used to infer broad trait biogeography.

### 4.5 | Implications for future comparative tests

Comparisons of islands, mountains and seasonal plant assemblages should track three noninterchangeable objects: source-photo sampling opportunities, photo-specific floral-colour measurability, and the sampled taxa's actual visible colour distributions. Matching environments while ignoring which taxa are available tests a source assemblage contrast; matching taxa while discarding most sites tests a distinct, narrower support. More imaging alone is not a replacement for preserving the selection probabilities and auditing sources with no valid flower masks. Mechanism-focused work can build on these global associations by measuring independently validated petals, biotic interaction contexts and within-taxon replication.

## 5 | Conclusions

Photograph-derived worldwide floral-colour composition covaries with annual moisture, and geographically comparable higher-elevation photographed floras show a lower share of white-labelled flower-like regions. Yet the apparent worldwide white frequency changes with taxonomic weighting, the moisture direction weakens under more geographically constrained comparison, an island-white excess is uncertain when environmental matches are required, and season-associated photographic colour differences are strongly shaped by how taxa are represented. These findings jointly support a **scale-dependent, assemblage-conditional view** of floral colour biogeography rather than a single universally transferable pigmentation rule. Independent original-photo outcome audits define how far these inferences can currently be generalized.

---

## Data accessibility and research transparency — working submission text

Source-identifiable, SHA-locked numerical summaries and code are available in the MOIS study repository and its immutable releases. A submission-anonymous source/code package and clear licensing/access provisions for the original photo-level climatic tables are being prepared; the paper must not claim that privately held original ecological analysis inputs or third-party licensed images have already been publicly deposited. Original photos have not been redistributed by this manuscript. Directly reconstructible primary ecological coefficients and exploratory follow-up summaries are enumerated in the accompanying evidence ledger.

The principal environment association was fitted on 309,439 original photo attempts and is **not** a regression refitted on all 1,181,146 signed original photo attempts. The larger union supports image measurement, nominal-taxon coverage and technical nonresponse statements only. The 350,399-source seasonal descriptive analysis and 313,535-source environmental geographical matches use their own frozen source populations. These independently source-audited panels are **not independent full-world replications**.

## Author contributions, funding and declarations — author completion required

Author names, contribution statement, funding/permit confirmation, licence acknowledgements, competing interests, and identification-bearing title page are pending author/coauthor review. The anonymous submission copy must omit GitHub owner-identifying URLs and metadata; do not insert fictional affiliations or contributor declarations.

## Selected references — verified and to be completed before submission

Dalrymple, R. L., Kemp, D. J., Flores-Moreno, H., Laffan, S. W., White, T. E., Hemmings, F. A., & Moles, A. T. (2020). Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist*, 228, 1972–1985. https://doi.org/10.1111/nph.16737

Di Cecco, G. J., Barve, V., Belitz, M. W., Stucky, B. J., Guralnick, R. P., & Hurlbert, A. H. (2021). Observing the observers: How participants contribute data to iNaturalist and implications for biodiversity science. *BioScience*, 71, 1179–1188. https://doi.org/10.1093/biosci/biab093

Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: New 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37, 4302–4315. https://doi.org/10.1002/joc.5086

Koski, M. H., & Galloway, L. F. (2020). Geographic variation in floral color and reflectance correlates with temperature and colonization history. *Frontiers in Plant Science*, 11, 991. https://doi.org/10.3389/fpls.2020.00991

Xiong, Y., Varadarajan, B., Wu, L., Xiang, X., Xiao, F., Zhu, C., ... Chandra, V. (2024). EfficientSAM: Leveraged masked image pretraining for efficient segment anything. *Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition*, 16111–16121. https://doi.org/10.1109/CVPR52733.2024.01525

**Bibliography completion gate:** verify all source-year/volume/DOI entries, add remaining relevant macroecology/plant colour and citizen-science bias studies, and provide full methodological primary source attribution before external submission. No uncertain citation should be silently invented.

## Figure plan

Figures 1–5, their numerical provenance, caption-level inferential restrictions, and the supplementary source-to-estimand mapping are in the companion manuscript figure and evidence ledger. Do **not** draw a photo-weighted global ecological white-rate curve using the new adaptively enriched 1.18m union.

**Working draft only: not reviewed by coauthors or submitted.**
