# MOIS — the seasonal selection of classifiable original flower-colour photographs

**Date:** 17 October 2026. **Source:** 267,425 real original photo IDs, 129,795 original ROI four-hue labels, and 127,504 classifiable labels linked to exact iNaturalist public observation dates. This is an observation-process quality analysis on the same original source photos, **not** a new FCP within-species W/D study or a claim about flower pigment plasticity.

## The original photo denominator matters

The full date-enriched original source has **261,578 date-bearing photographic attempts**. Of those, **127,504** carry an automated label of white, yellow/orange, red/pink or blue/purple, while **134,074** have an unknown technical hue. We cannot classify the latter as chromatic, nonflowering, white, or a missing biological morph. The independent botanical certainty of an **open floral organ** remains unverified in the enlarged source.

To avoid conflating two source processes, define for each sampled local season:

- \(c_t\) = fraction of all *date-bearing actual original photo attempts* with a classifiable original four-hue ROI;
- \(p_t\) = fraction of classifiable ROIs automatically labelled white;
- \(y_t = c_t p_t\) = fraction of all attempted original date-bearing photos with an observed *white label* (not true white-flower prevalence).

The symmetric **exact photographic detection algebra** is

\[
100(y_S-y_W)=100\bar c(p_S-p_W)+100\bar p(c_S-c_W),
\]
with seasonal means \(\bar c=(c_S+c_W)/2\) and \(\bar p=(p_S+p_W)/2\). This is NOT a causal decomposition of true open flower colour, biology or observation effort; missingness may depend on flower organ identity and source/region/taxon.

## Whole-source extratropical summer−winter differences

| Original source photo measure | North summer | North winter | South summer | South winter |
|---|---:|---:|---:|---:|
| Actually dated original photo attempts | 72,857 | 8,566 | 12,657 | 9,421 |
| Four-hue classifiable photo ROIs | 38,121 | 4,059 | 5,827 | 4,082 |
| Automated white photo labels | 18,620 | 1,850 | 2,996 | 1,589 |
| **Four-hue classification success** \(c\) | **52.32%** | **47.39%** | **46.04%** | **43.33%** |
| **White share among classifiable** \(p\) | **48.84%** | **45.58%** | **51.42%** | **38.93%** |
| White-labelled outcome per original photo attempt \(y\) | **25.56%** | **21.60%** | **23.67%** | **16.87%** |

The photo-classifiability rate itself changes across photographed seasons: **North +4.94 percentage points summer−winter**, **South +2.71 points**. The original automatic-white fraction *conditional on classifiability* changes by **+3.27** and **+12.49**, respectively.

The exactly reconstructed summer−winter difference in the chance of *recording a white label per original photographed attempt* is:

| Photographic detection identity (percentage points) | North | South |
|---|---:|---:|
| Actual summer − winter white-label yield among **all** source photo attempts | **+3.960** | **+6.804** |
| Conditional automated-white-photo-share component | +1.629 | +5.580 |
| **Four-hue classifiability / technical selection component** | **+2.331** | **+1.224** |
| Identity check | 1.629 + 2.331 = 3.960 | 5.580 + 1.224 = 6.804 |

Thus in the Northern source, *more than half of the observed-white-label **per-attempt yield** contrast* is algebraically linked to different rates of technically classifiable photographs. This does **not** say that more than half of the original conditional white-hue-share difference is a causal artifact, nor that the unavailable colours are nonwhite. In the South, the much larger conditional white-label share component is dominant, but source selection remains a concern.

## Geographic sampling changes even the sign of the classification process

The main MOIS seasonal locality comparison only admits original **0.25° site cells and the same meteorological season-year** where **at least five distinct original observations and two original nominal taxa** have classifiable colours in *each* summer and winter. It keeps the original set of **16 Northern and 18 Southern matched cell-years** used for the earlier white-label composition comparison. Recounting *all* originally dated photo attempts, including technical-hue-unknown photos, at those exact source cells reveals:

| Summer − winter per-site/year (equal-cell-year mean) | North (16 cells) | South (18 cells) |
|---|---:|---:|
| **Four-hue classifiability percentage-point difference** | **−10.78** | **−4.05** |
| Exploratory 999-geocell bootstrap percentile 95% | [−21.92, −0.73] | [−11.13, +4.41] |
| Conditional white-label share difference on that same source support | +5.04 | −4.63 |

The North-wide summer increase in classification success **reverses direction** within these limited, matched original sites and observation years. This is strong descriptive evidence that the source's geographical/time and photo-eligibility composition confounds a simple global summer classification gradient. However, the within-cell/year subset is deliberately *selected because photos with classifiable ROIs exist in both seasons*, contains only 16 and 18 unique original coarse cells, and is not a randomized observer- or flower-stage-balanced sample. The Northern exploratory CI excludes zero in this one chosen comparison, but no maxT correction was applied to the multiple sample/site scales explored; do not turn this post hoc test into a universal seasonal observation mechanism.

## Method and reproducibility

The original source [run #39002957961](https://github.com/zuizui0223/mois/actions/runs/39002957961) contains frozen photographed original pixels and four-hue attempts; [full-date source reassembly run #39012110980](https://github.com/zuizui0223/mois/actions/runs/39012110980) supplies source-verifiable final 99,817 historical dates (S3 September 2026 first; 85 independently numeric observation/photo-ID verified API fallbacks) without changing original photo IDs/ROIs. The new script reads these source archives and refuses a changed photo/source-URI/date/ROI/taxon identity, then exports a full-source original photo classification/white-yield decomposition and locality-year matched classification efficiency sensitivity to a source-signed aggregate JSON/CSV artifact. Four-hue classification status is always based only on the original measured ROI-v4 `morph` state.

The original historical **annual** white-photo vs rainfall and atmospheric moisture genus-equal 999-bootstrap maxT inference remains separate and unchanged. The global summer/winter ratios are unpaired photo composition estimates, and locality comparisons are coarse cell-averaged photo-source outcomes. None of these original photographs were newly botanically verified as *open petals or tepals*, and no observed-year weather, genetic morph, pollinator selection, pigment metabolic cost or reproductive fitness was measured here.

**Paper implication:** main article should distinguish (i) photographed white among classifiable images, (ii) whether original photos receive a valid colour classification, and (iii) actual biological flowering. This makes a global seasonal colour association less likely to be mistaken for a universal within-locality colour adaptation.
