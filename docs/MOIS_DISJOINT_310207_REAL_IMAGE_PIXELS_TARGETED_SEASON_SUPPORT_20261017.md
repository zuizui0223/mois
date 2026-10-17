# MOIS | Disjoint additional image pixels and targeted season support (17 October 2026)

## Scientific and data decision

The **primary observational macroecological article denominator is frozen** at
**309,439 actually image-processed original photographs**, **149,662 ROI-v4
four-hue labels**, and **149,653 original photo labels linked to verified observation dates**.
The six-axis, equal-genus precipitation/vapour maxT999 WorldClim 1970–2000
association and the unpaired regional summer-winter photo-colour rates continue
to refer to that original source. They must not be naively recomputed after
targeted enrichment.

Two further completed original RGB/model-photo arms exist:

* **Fourth source wave, partial two-worker checkpoint**: 512 unique original
  photographed pixel outcomes and 203 date-bearing four-hue ROIs;
  [real-photo source and model attestation](https://github.com/zuizui0223/mois/actions/runs/39034832999).
  These are the FIRST TWO SHARDS OF the still-running 4,096-photo cohort,
  **not** an additional independent 512 to be summed with its final 4,096.
* **Observer/observation-distinct, taxon/cell/season-year TARGETED arm**:
  256 actually processed new photographic RGB pixel outcomes; **71** have a
  date-bearing four-hue label, 185 remain technically unclassifiable;
  [true pixel-work, model, source-shard and site/year receipts](https://github.com/zuizui0223/mois/actions/runs/39033985014).
  Candidates were selected on previous classifiability support, named taxon,
  authentic observed-on and 0.25-degree cell only, **not new hue**. All 256
  photo IDs were excluded from the separate fourth-4096 reservation.

The [reconciliation workflow](../.github/workflows/mois-reconcile-disjoint-fourth512-targeted256-real-pixels.yml)
compares both complete cohort files against the identical immutable
309,439-photo GZIP SHA, checks previous prefix row-for-row, compares every
new photographic photo ID across all three arms, checks completed signed
worker receipts, and saves each arm's actual classification yield separately.
It never imputes pixel outcomes from the remaining source-only candidate IDs.

### Counts once photo-ID-disjointness is verified

| Source arm | Pixel-processed photo IDs | Four-hue ROIs | Dated four-hue ROIs |
|---|---:|---:|---:|
| Original macroecology/season source | 309,439 | 149,662 | 149,653 |
| Ordinary fourth4096 first-two-worker checkpoint (subset!) | 512 | 203 | 203 |
| Purposefully selected seasonal sampling arm | 256 | 71 | 71 |
| **Combined photo-ID audit only** | **310,207** | **149,936** | **149,927** |

The combined 310,207 is a legitimate *unique image-pixel measurement
inventory*, **not** a new globally representative flowering census. Check
source signatures and exact numerical total on successful action before
claiming the combined photo-archive artifact complete.

### Empirical spatial/seasonal support

Using only named-taxon × 0.25-degree cell × identical meteorological
season-year original photographic observations that have a real classified
four-hue ROI, the completed targeted 256 run reported:

| | Northern extratropics | Southern extratropics | Total |
|---|---:|---:|---:|
| Baseline >=1 photo observation in BOTH summer/winter | 40 | 34 | 74 |
| With 256 enriched real image outcomes >=1 in BOTH | 87 | 54 | 141 |
| Added >=1 groups | +47 | +20 | **+67** |
| Baseline >=2 distinct original observations in BOTH | 2 | 0 | 2 |
| With 256 enriched real image outcomes >=2 in BOTH | 2 | 1 | 3 |

The ordinary first 512-photo shard checkpoint added **zero** exact
same-named-taxon/site/season-year summer-winter groups; it is not the
replication-targeted selection scheme. The 256 targeted image pixel work
has improved **source support**, but 141 minimally matched photo groups are
**not 141 adequate within-plant seasonal effect estimates**. The genuinely
replicated >=2 subset contains only 3 groups.

### Limits for the manuscript

- ROI-v4 visible hue identifies a flower-like photographic target, **not
  photo-specifically botanically confirmed open petals**. Fruits, bracts,
  stage and colour/missingness remain relevant audit targets.
- A 0.25-degree cell spans kilometres; a named taxon and source observation
  ID do not establish repeated individuals, flowering population abundance,
  independent photographers across arms, or local weather exposure.
- The targeted high-opposite-season source selection is intentionally
  **unequal-probability ascertainment**. Never add its white images as if an
  unbiased global atlas, nor infer seasonal pigmentation, flower availability,
  pigment cost/benefit or fitness from the +67 sampling-support increase.
- Actual observed-year monthly temperature, solar radiation, precipitation
  and water balance have not been mapped to all new original photographs.
  The article's current annual WorldClim predictors are historic climate
  normals, and soil/elevation/terrain remain slow/static geographic controls.
- No FCP within-species whiteness/darkness/50-km phenotype estimator is
  transplanted into this global among-taxon MOIS study.

**Next data-acquisition criterion:** After the full fourth4096 real RGB job
completes, replace the 512 subset with all 4096 measured photo outcomes
using literal photo ID, preserving this independent targeted256 arm only
once. To advance actual seasonal ecology, target >=2 independent field or
source observations per season/taxon/cell/year, and separately verify
open-petal stage; metadata-only candidates do not close that gap.
