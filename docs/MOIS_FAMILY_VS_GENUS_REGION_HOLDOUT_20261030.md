# MOIS: independent Kew families vs accepted genera in six geographic holdouts (2026-10-30)

## Authoritative frozen input and result

- Existing original photographic strict9 cohort, **not a new photo sample**: previous numerical stage [41018155739](https://github.com/zuizui0223/mois/actions/runs/41018155739).
- Accepted WCVP species/genus previous data-only step: [41042569706](https://github.com/zuizui0223/mois/actions/runs/41042569706), 30,054 region×accepted species entries after resolving 30,158/30,634 nominal taxon×region entries; 476 unresolved and 104 synonym collapses. Genus source [document](MOIS_GENUS_REGION_HOLDOUT_REAL20261030.md).
- **New true independent official family lookup:** Original Kew WCVP source `wcvp.zip` downloaded once by this data-only workflow and SHA-pinned to `d32ea2b3a85e489b14e83bcc9eae7274532e1d113753f7be290d4b2dfde573fa`. Official June-2026 WCVP PIPE-separated `wcvp_names.csv` contains `family`, accepted `plant_name_id`, `taxon_name`, `taxon_status`, `taxon_rank` and `genus`. The original complete Kew archive is NOT redistributed.
- [Family run 41043422845](https://github.com/zuizui0223/mois/actions/runs/41043422845) **SUCCESS**; [numeric artifact 12970763367](https://github.com/zuizui0223/mois/actions/runs/41043422845/artifacts/12970763367), ZIP SHA256 `3eb62f10607c13893e410dc1b7abdf9f90d4f1e730b08edecb245669823d1c63`.

**Source audit:** 1,448,984 original Kew name rows streamed, exactly **22,712/22,712 required accepted species IDs** resolved to accepted species name and original `family` without genus conflicts. All **30,054 sector×accepted species rows**, representing **417,111 measured successful source-photo observations across overlapping sector definitions**, carry a family; **312 distinct accepted families**. The original source's 476 unresolved nominal taxon×region entries remain excluded, never assigned families. Summing regional photo supports is **not** the union of independent photos, because island/continent sector labels and latitudinal sectors are overlapping classifications.

## Different equal-weight definitions give different regional signals

These are original automatically segmented strict9 ROI white-mass average differences (percentage points), **not fractions of botanical white-flowered individuals**.

| Region A minus B | Equal accepted species from previous study | Equal accepted genus | Equal accepted family | Common family only |
|---|---:|---:|---:|---:|
| Arctic minus north subpolar | −0.620 | −2.664 | **+0.650** | −2.763 |
| Southern temperate minus northern temperate | −1.436 | −1.644 | **−1.267** | −0.878 |
| Island minus mainland (historic exact landmass cache) | +1.465 | +0.934 | **+4.178** | +3.192 |

Changing from equal-genus to equal-family weighting changes the estimand; the +4.178 island family contrast is **not stronger evidence for island selection**. Species turnover across families and within families, imbalance in species per family and source geographic support all matter. Requiring the same family again changes the target, not merely removes confounding. The independent original photo sampling remains preferential, and familial means are not phylogenetic contrast tests.

For each comparison, 5 metrics (white/blue/purple/bronze ROI components and nine-palette entropy) were decomposed into common family trait differences, common-family relative weights, focal-only families and reference-only families; **15 exact numerical closures passed**. Original reference region taxa do not necessarily occur in the focal region.

## Heldout accepted species: how much can phylogenetic rank alone transfer?

Only accepted species **absent from the training region** are scored. Family model requires at least **5 different training-region accepted species** with a WCVP family; genus model requires at least **3 accepted species** in the same training-region genus. The 6 comparisons evaluate the same broad sector definitions with opposite training/target direction, but are not independent random biome holdouts.

| Target ← training region | Novel target species | Family-supported | Genus-supported (previous) | Identical species scored by both | Family white MSE | Genus white MSE |
|---|---:|---:|---:|---:|---:|---:|
| Arctic ← north subpolar | 132 | **131** | 99 | 99 | 0.1217 | **0.0915** |
| North subpolar ← Arctic | 1,544 | **981** | 306 | 286 | 0.1301 | **0.1219** |
| South temperate ← North temperate | 8,492 | **7,488** | 2,030 | 2,013 | **0.1112** | 0.1136 |
| North temperate ← South temperate | 12,452 | **11,486** | 3,555 | 3,493 | 0.1039 | **0.1009** |
| Island ← mainland | 162 | **150** | 72 | 71 | 0.1128 | **0.1051** |
| Mainland ← island | 1,538 | **1,346** | 446 | 435 | 0.0990 | **0.0905** |

**The two MSE columns are computed only on the same intersection of previously unseen accepted target species in each row**, not on all family-only or genus-only supported taxa. In these observations, genus prediction has smaller numerical squared error in five of six directional comparisons. No inference/uncertainty interval is computed for their difference; this is a reproducible descriptive cross-region test, not evidence of statistically established superiority. Family provides broader coverage because clades span more species and regions, but the average within-family trait may be less specific. For islands, family can score 150/162 novel species, but only 150/1,457 = **10.3% of all focal accepted species** are newly testable under this exclusion of shared species, so avoid interpreting this as 92.6% of island-flora predictability.

## Safeguards

- Measured traits: only existing strict ROI9 **visible colour components** and normalized nine-palette entropy; no direct botanical confirmation of flower petals, chlorophyll investment, dark petals, ultraviolet reflectance, calibrated chroma, display area or background contrast.
- No green independent component; bronze cannot be called verified brown flower or inconspicuousness.
- WCVP family = taxonomic rank, **not a pollination guild or flowering syndrome**. No ecological guild imputed from family.
- Original Kew TDWG3 native-range flags were **not** interpreted as native on a particular island. Geographic results apply to already cached island/mainland photo assignment only.
- No new iNaturalist images, inference, source photo counts, unknown-colour imputation, MNAR scenario, or change to main manuscript.
- This family reclassification is independent botanical metadata added to the previously measured source, not new independent flower observations.

## Verification and next step

Nine synthetic provenance/heldout unit tests and the original Kew SHA256 gate passed on GitHub. Every result file's SHA256 matches the released receipt; an independent local reimplementation verified all six family-vs-genus same-species MSE values to max 1.25e−16 and all 15 family partitions to numerical closure.

The remaining functional-syndrome hypothesis requires independent species-level measurement of **corolla/display size, perianth reduction and background-relative reflectance**, plus scientifically sourced pollination-system records, not inferred from family or the 500px YOLO box. Until then, the valid result concerns *taxonomic structure of photographed floral visible-colour mixtures* rather than an estimated world distribution of 'dull flowers'.
