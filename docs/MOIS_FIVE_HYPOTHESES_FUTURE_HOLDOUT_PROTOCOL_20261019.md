# MOIS protocol freeze: five global flower-colour hypotheses, photographic source and future unseen outcomes

**Version:** 2026-10-19, `v1`. **Status:** an analysis protocol for a **newly unseen, independent photo-outcome holdout**, **NOT a retrospectively legitimate preregistration of already examined 350,399 photos**. The machine-readable immutable choice of directions, source exclusions, outcomes and gates is `analysis_specs/mois_five_global_photo_hypotheses_future_holdout_v1.json`. All modifications after the frozen commit must create a new version, with a visible changelog and no replacement of the original.

## 1. Scientific question

**Across worldwide photographed plant assemblages, when does the balance of apparent WHITE versus chromatic flowers vary with abiotic stress, native island assembly, latitude, and the independently measured availability of pollinators?**

The primary response is a **taxon-equal, geographically balanced probability of a white *visible open-flower photograph*** among actually observable, classifiable flower photographs, not the unconditioned percentage of all vascular-plant individuals in nature. An iNaturalist photograph's original automated ROI-v4 hue is an *image annotation*, not yet verified open petals. Where independent botanical review is absent, analyses must be labelled **photographic ROI-label analyses**, not literal botanical flowering. Absent photographs do not demonstrate absence of flowering.

The `fcp` within-species diversity/polymorphism project is distinct: MOIS examines **species and spatial assemblage composition**, using an external FCP polymorphism `D` estimate only after independent taxon ID/source verification. Do not generalize an FCP sample's reported median dominance to all sampled MOIS species or deterministically force multicoloured species to a single class.

## 2. Existing results already seen: cannot be preregistered retrospectively

At protocol creation, [350,399 authentic attempted original photo pixels](https://github.com/zuizui0223/mois/releases/tag/mois-20261019-real-photo-source-arm-season-including-tail16384-v1) and [164,297 four-hue classifiable ROIs](https://github.com/zuizui0223/mois/blob/main/docs/MOIS_REAL350399_ACTUAL_RGB_FOURTERM_SEASON_AND_TAXON_20261019.md) have **already been examined**. Photo-white 47.49%, nominal taxon-equal 48.19%, genus-equal 51.46%. The global species-first addition changes the photographed-white fraction and [matching genus×geography attenuates the between-source difference](https://github.com/zuizui0223/mois/releases/tag/mois-real350399-geogenus-matched-source-selection-v1). These are *discovery* results.

| Hypothesis | Direction for a future untouched holdout | Already exposed evidence (discovery, NOT a new test) |
|---|---|---|
| H1 High elevation | Higher elevation → lower WHITE | Yes: [previous actual-image elevation contrasts](MOIS_ISLAND_MAINLAND_ALTITUDE_FIRST_RESULTS_20261018.md), non-monotonic alternatives remain |
| H2 Dry-season stress | More negative log-dry-month rainfall → lower WHITE | Yes: inverse of exploratory moisture directions, not an independent replication of H5 |
| H3 Native islands | Native island assemblages → higher WHITE than matched native continental assemblages | Yes: [earlier raw and geography-conditioned island contrasts](MOIS_ISLAND_MAINLAND_ALTITUDE_FIRST_RESULTS_20261018.md) and [support-limited sensitivity](MOIS_ENVIRONMENT_OVERLAP_AND_REGIONAL_ELEVATION_20261018.md) |
| H4 Arctic | Northern Arctic native photo assemblages → higher WHITE than matched adjacent non-Arctic assemblages | Yes: prior Arctic-circle subset is observed and *does not show a universal white excess* |
| H5 Atmospheric humidity | More atmospheric water vapour → higher WHITE | Yes: previous positive global model directions lose sign support under geographic matching |

All five directions are **hypotheses informed by prior observation**. A statistically significant fit to the same previously studied original pixels is **not independent confirmation**, regardless of this file's timestamp. The ongoing 93,209,351 source-photo-ID index is *metadata*, not 93 million measured flower colours.

## 3. Independent observation-level, pre-colour future outcome split

Only a new, genuinely unmeasured **observation** cohort from a complete, signed 42,111-taxon / 93,209,351-photo source index can be called prospective for the **new photographs' colour outcomes**. The existing `mois-full-index-observation-queue` defines a fixed pre-colour p=0.0002 hashed selection, retains *all attached licensed photos* of selected observations and does not impose a species/genus cap. This source frame has iNaturalist source-filter biases even if the selection is internally random.

Before accessing its hue results, freeze and SHA-sign:

1. Every source photo/observation index shard; full 93,209,351-photo count and original observation IDs; original sampled observation hash threshold, seed, source inclusion probability, source-arm labels and the complete selection manifest.
2. The 350,399 previously observed **photo IDs** and all their **observation IDs**. Remove from the unseen-outcome evaluation **every observation already represented among previously seen photo outcomes**, not merely duplicate photo IDs. Otherwise additional images of a known observation contaminate the untouched holdout. This conditions the estimand on previously *unseen* observations and the exclusion itself must be reported.
3. Taxonomic accepted-name, native-status and geographical joins and source SHA; fixed plan plus model/feature engineering, predicted directions, quality calibration sample, holdout region/genus blocks and outcome code before viewing any held-out colour labels.
4. End-to-end completion of the actual ROI-v4 original-pixel jobs, including failures, unknown hues and botanical-stage audit. A partial batch of only fastest-classified photos is *not* a released unbiased cohort.

The [protocol-gate script](../scripts/analysis/mois_five_hypotheses_protocol_gate.py) checks that the declared split is non-overlapping at both photo and observation levels and *refuses confirmation* without full-index signatures and provenance. A holdout generated using another photography platform is truly independent of observer ecology only when taxon/region matching and platform permissions permit that claim.

## 4. Estimand, equal-area geography, colour labels, effort

Aggregate all genuine source photos to **one source observation**, then accepted taxon × **100-km equal-area (EPSG:6933) grid cell** × an annual observation window. Average colour from multiple validated original open-flower photographs in one observation; estimate the taxon's probability of WHITE, **not a binary forced majority**. One-pixel success is uncertain; report denominator and classification/flower-stage validation, sensitivity at ≥2 independently annotated observations per taxon×cell. Primary weights: equal taxa within each cell and equal supported cells within biogeographic realm. Descriptive photo-, observer- and genus-equal estimates must remain separate, not conflated.

For seasonally assembled/photo-detectable floras, separate hemisphere-specific meteorological seasons and year (December assigned to following local season year) and decomposed `source-photo` taxon-label frequency:

- Apparent WHITE rate change among taxa **photographed in both seasons**.
- Changes in relative photographic frequency of shared taxa.
- Support corrections for taxa photographed **only in summer** and **only in winter**, without imputing a missing opposite-season hue.

The already observed Southern unpaired summer minus winter **+10.50 pp**, including **+8.04 pp common-taxon relative photographic frequency**, is a *photo mixture*, **not a measured 8.04-point phenological process**. To approach flowering-phenology inference, assemble counts of **all identifiable plant observations**, per-cell/month **distinct observers**, photo source effort and imaging/classifiability outcome denominators. Non-flower observations can act as an imperfect observation-effort denominator if their ascertainment and anatomical labelling are validated; they are **not verified absences of blooming**. Independent field checklists would strengthen biological interpretation. Observer-disjoint analyses, local cell/year overlap and a negative-control non-flower snapshot sensitivity are required before suggesting seasonal biological availability rather than shooting behaviour.

## 5. Signed external environmental, range and pollinator joins

- **Native/introduced**: Kew [WCVP](https://powo.science.kew.org/about-wcvp) records accepted names and native/introduced ranges **at TDWG level 3**. It cannot automatically determine if one particular small island population is native. `introduced`, `native_at_TDWG3_but_island_uncertain`, `native_on_specific_island_documented`, `not_resolved` are *separate states*. Never treat missing WCVP distribution as native. Island-specific native/endemism main tests require an additional authoritative island flora/range check with a provenance tag; otherwise H3 remains exploratory/unresolved.
- **Islands**: use original GSHHG exact landmass polygons and area, preserve coastline/coordinate-uncertainty unassigned, distinguish small vs huge islands, geographical realms and latitude. Greenland/newly emerged oceanic islands are not assumed to have uniformly few pollinators. Flora native and introduced strata must be shown separately.
- **UV**: use independent [glUV/UVB1](https://doi.org/10.1111/2041-210X.12168) **15-arc-minute 2004–2013 UV-B** normals, explicitly not the WorldClim *short-wave solar W/m²* proxy. Spatial scale and period differ from the existing **WorldClim 2.1 1970–2000** atmospheric moisture/temperature/elevation fields.
- **Bee richness**: [Orr et al. 2021](https://doi.org/10.1016/j.cub.2020.10.053) combines checklist/occurrence sources and documents extreme geographic sampling bias (particularly for Asia/Africa). This is a candidate external ecological *context*, **not a validated grid of local bee visits, abundance, or pollination effectiveness**. A real source-backed, appropriately resolved and effort-audited global surface must be acquired before an H3/H4/abiotic “pollinator pathway” comparison is run. If absent, report **pollinator component not estimable** rather than plugging in a fabricated raster or a plant-photography-co-sampled pseudo-proxy.
- **Photographic brightness**: blind exposure/petal ROI visual-quality metrics must be computed on actual source imagery in a calibration subset (classifiable and unknown outcomes); until available, light-dependent white misclassification is an *unmeasured sensitivity*, not an automatically fitted adjustment or guaranteed conservative-bias direction. Unverified OR=1.44 is not accepted as an MOIS source result.
- **Other factors**: include annual temperature/BIO05, BIO12, BIO14, atmospheric VAPR, solar, roughness/elevation, realm and season/year, source cohort and photographic effort; soil only if source verified (no imputation). Independent bee and UV sources are held at *unknown* until checksum/coverage/resolution validation.

## 6. Joint test and no-causality gate

Estimate a *conditional association*, not a fitness/anthocyanin result. Use hierarchical binomial/beta-binomial photo/observation counts with **accepted genus** (family where reliable) and spatial covariance/random effects, accepted taxon×cell composition, standardization fitted **within training folds**, and train/test blocked jointly by geography and genus where identifiable. Distinguish source recruitment arms. A non-identifiable/empty native island, Arctic, UV or bee join is **not** silently substituted with an easier sample.

The principal specification contains altitude, dryness (negative log1p BIO14), VAPR, BIO05, BIO12, UVB1, solar, island/native, Arctic, observation year/effort, realm/latitude/longitude, verified independent bee richness (if supported) and image-quality calibration. Because drought and VAPR often co-vary, if the competing H2/H5 coefficients are unstable, report **shared water-related predictor block** without claiming separately identified stress and bee-mediated mechanisms. Distinguish the full fitted model from geography-only, taxonomy-only, source/effort-only, abiotic-only and bee-added models on **exactly the same holdout support**.

The five named directional contrasts form a prespecified Holm-adjusted family (α=0.05) **only in a genuinely untouched future observation/photo outcome cohort**. One-sided signs do not grant licence to ignore opposite results; signed effects and blocked uncertainty are mandatory. H2/H5 represent related moisture gradients and cannot count as two mechanistically independent replications. All historical source results and source-dependent model explorations remain explicitly labelled exploratory.

**Negative outcomes must be reported**: complete-case taxon/cell/observation count, review failure, unmatched biogeographic cells, native status unavailable, unknown hue fractions, taxon-genus-holdout transfer failure, bee coverage and source-photo sampling rates. A very small overlapping native island subset is a limited-support test, not evidence of no island effect.

## 7. What this protocol does not establish

The study does **not** measure biochemical anthocyanin absence/presence from RGB WHITE, UV floral contrast, actual pollen transfer, per-plant fitness, evolved pigment maintenance costs, or a single shared causal pathway. A multivariable association in non-experimental macroecology cannot choose among direct plant stress, community composition, ecological interactions and sampling artifacts by statistical adjustment alone.

### Existing real-data discovery basis

- [Current 350,399 original-pixel data and four-term season composition](MOIS_REAL350399_ACTUAL_RGB_FOURTERM_SEASON_AND_TAXON_20261019.md)
- [Actual genus×geography source-selection matching](MOIS_REAL350399_GEOGENUS_SOURCE_SELECTION_BIAS_20261019.md)
- [Existing island/elevation/Arctic discovery](MOIS_ISLAND_MAINLAND_ALTITUDE_FIRST_RESULTS_20261018.md)
- [MOIS independent plant-botany review QC](MOIS_400_BOTANICAL_REVIEW_HANDOFF_20261017.md)
- [Upcoming observation sampling data acquisition protocol](MOIS_FULL93M_FIRST_OBSERVATION_RGB_ACQUISITION_20261019.md)

**Next data dependency:** the complete source-index Release of photo/observation IDs and source-authenticated future observation cohort; the WCVP/UV/bee joins cannot be claimed complete before retrieval, name resolution, licence checks, scale alignment and SHA receipts.
