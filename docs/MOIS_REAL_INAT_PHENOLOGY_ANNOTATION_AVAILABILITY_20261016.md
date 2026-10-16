# iNaturalist photographed-plant observation bias: real flowering-annotation audit

**2026-10-16. Actually executed public API count study — NOT the original dated FCP photo-ID population, not measured flowers and not the new 13-variable model.**

## Why this check matters

The worldwide Mois paper uses photographic flower colours. iNaturalist supports Flowers and Fruits observations (term ID 12) with Flowers (value 13), No Flowers or Fruits (value 21), buds and fruit values. It is tempting to use these as automatic flowering labels. But iNaturalist annotations are *per observation*, not per photo. An observation with multiple images may have a flower in some photos but not in the representative photograph chosen for measurement. Conversely, a photographed flower may exist in an observation without any user-added Flowers annotation.

The official [iNaturalist annotation search documentation](https://help.inaturalist.org/en/support/solutions/articles/151000170738) and [API recommended practices](https://www.inaturalist.org/pages/api+recommended+practices) permit the small count-only diagnostic below. The platform advises approximately one request per second and pre-generated export rather than high-volume observation scraping.

## Source and fixed design

[Actual successful GitHub Actions collection](https://github.com/zuizui0223/mois/actions/runs/38843210801), [verified independent artifact receipt recovery](https://github.com/zuizui0223/mois/actions/runs/38843566842), complete [35-query JSON count/audit receipt](../results/mois_real_inat_photo_phenology_annotation_count_pilot_20261016.json), and its [reproducible code](../scripts/analysis/mois_real_inat_phenology_count_bias_pilot_20261016.py).

The filter is taxon_id=47126 (Plantae), photos=true, per_page=1 with only **total_results** persisted; no image pixels, observer names, individual record IDs or raw observation JSON were retained. The actual API snapshot was collected in October 2026 for **observed_on date windows in 2024–2025**. Regions/date bins were fixed before reviewing counts: the world for all of 2025; Northern East Asia (31–45°N, 129–145°E) spring March–May versus summer June–August 2025; southern Australia (22–38°S, 113–153°E) local spring September–November 2025 versus summer December 2024–February 2025; and tropical South America (10°S–10°N, 75–45°W) calendar Q1 and Q3 of 2025. The two tropical quarters are deliberately **NOT** equated with local wet and dry months.

Each cohort was queried for five category counts: all photographic plant observations, any Flowers and Fruits annotation, specifically Flowers, specifically No Flowers or Fruits and without any Flowers and Fruits annotation. Every request URL and SHA256 of its original response JSON are stored in the complete result. All 35 requests ran sequentially, with >=1.35 seconds between requests, with no individual observation archive or provider image download.

## Actually counted photographic plant observations

| Fixed publicly searchable observation group | All photographed-plant observations | Any Flowers/Fruits annotation | Specifically Flowers annotation | Flowers share of all |
|---|---:|---:|---:|---:|
| World: observations in 2025 | 28,565,884 | 2,550,229 (8.93%) | 1,615,314 | **5.65%** |
| N. East Asia: local spring 2025 | 52,403 | 6,568 (12.53%) | 5,351 | **10.21%** |
| N. East Asia: local summer 2025 | 38,805 | 4,481 (11.55%) | 3,245 | **8.36%** |
| S. Australia: local spring 2025 | 395,531 | 84,266 (21.30%) | 67,646 | **17.10%** |
| S. Australia: local summer 2024–25 | 143,217 | 25,636 (17.90%) | 16,502 | **11.52%** |
| Tropical S. America: calendar Q1 2025 | 26,576 | 2,497 (9.40%) | 1,519 | **5.72%** |
| Tropical S. America: calendar Q3 2025 | 28,687 | 1,648 (5.74%) | 1,105 | **3.85%** |

For the worldwide 2025 query, explicitly **No Flowers or Fruits** was annotated on 424,670 photographic observations, about **1.49%** of the total. This is a different, explicit statement from **without** an annotation. The 'without Flowers/Fruits annotation' filter returned **26,015,038** observations (~91.07% of all). A roughly **91% unannotated photographic-opportunity share does NOT mean 91% of plants were not flowering**, nor does it prove that a particular photographed image lacks coloured petals.

**Non-atomic API uncertainty:** the five searches are distinct provider queries, not a transactionally frozen API database export. For the global cohort, 'any annotation' + 'without annotation' totalled **28,565,267**, **617 less than** the overall search count. Five of the six regional filters have smaller discrepancies (1–14 records); only Northern East Asian spring had an exact complement. These are recorded as actual search inconsistencies, rather than silently reassigning observations or assuming perfect partitioning. The denominator/missingness ratios remain descriptive individual query snapshots, not independently linked observation-level records.

## Scientific interpretation — what is and is not identified

1. **Observation annotation availability is heterogeneous in space and time.** The share with explicit Flowers differs across the fixed pilot strata and dates (e.g. 10.21% vs 8.36% in Northern East Asia, 17.10% vs 11.52% in Southern Australia, 5.72% vs 3.85% in the tropical Q1/Q3 comparison). These differences incorporate **which species flower and are photographed, where observers spend time, which observations get annotated, and when**, not just annotation behaviour. Thus the results identify heterogeneous ascertainment, **not its causal components**.
2. **The user-proposed Flowers filter is informative but not a complete measurement frame.** Flowers annotation is high-confidence positive **observation-level** evidence where present; the much larger unannotated denominator can contain both flowering and nonflowering photo opportunities. Filtering only positively annotated observations would select on user annotations and season/geography/plant composition. It would produce a differently ascertained population than the original FCP source.
3. **iNaturalist API counts are not the Mois source frame.** The worldwide observed-year 2025 count of 28,565,884 denotes current platform photographic plant observations across ALL Plantae taxa and includes observations outside the 42,111 predefined named species, potential source licence differences and changed/today's annotations. FCP's independently frozen 2026-09-27 bulk metadata **93,209,351 source photo opportunities** are *photographs*, often multiple per observation, and are NOT interchangeable with the API's year-2025 observation counts.
4. **No white-flower seasonal effect is estimated here.** There are zero newly downloaded plant image pixels, zero new photo-colour classifications, no ascertainment probability estimated for the fixed 42,111 taxa, no observer-level covariance available and no climate value joined to these API counts. No selection fitness, pigment cost or population incidence can be inferred.

## Mandatory follow-up in the single worldwide flower-colour paper

Keep the prospective region × precipitation stratum photo selection outcome-blind (or separately frozen local-season-balanced V2), retaining original inclusion probabilities and missing photographs. **Never exclude an unannotated observation from the primary photographic opportunity frame solely due to absent Flowers annotation.** After the complete original photo-ID freeze and colour ROI measurement, perform the linked observation-level **positive-Flowers versus unannotated/unverified diagnostic** in the independently pinned sample, and estimate photo-colour classifier response by observation source season/geography/genus and observer. Hold out entire genera/geographic cells in the 13-variable extension; do not explain a white–BIO12 coefficient using independently sourced annotation counts with no joined photograph hue.

In the Results/Discussion, this actual API finding is a **secondary platform-data ascertainment pilot**, not an independent confirmation that wetter real sites have more white flowers. The main paper's already source-verified 16,801-taxon/1,665-genus same-13-environment soil36 source and positive conditional BIO12 photograph-white slope are untouched.


## Taxonomically appropriate primary sensitivity: Angiospermae (flowering plants)

Because the whole-Plantae API denominator also includes **nonflowering plant lineages** (e.g. ferns and mosses), the previous 5.65% worldwide Flowers-annotation number is a **broad-plant diagnostic**, not the best denominator for a flower-colour question. A predeclared follow-up queried the identical 35 count filters and seven original spatial/date bins for Angiospermae, iNaturalist taxon ID **47125**, without any original source plant media retrieval.

[Real successful Angiospermae query run](https://github.com/zuizui0223/mois/actions/runs/38844098039) and [full 35-query immutable response-URL/SHA record](../results/mois_real_inat_angiosperm_phenology_annotation_count_pilot_20261016.json) use the same public provider limits, with original output artifact bytes SHA256 **67c22183e6919be4ed5f7d466f6047c178ca64d18a34479cee71c51b95ccc03e**.

| Flowering-plant observation cohort | All photographed Angiospermae observations | Any Flowers/Fruits annotation | Explicit positive Flowers | Positive Flowers / all |
|---|---:|---:|---:|---:|
| World 2025 | 26,069,247 | 2,550,231 (9.78%) | 1,615,314 | **6.20%** |
| Northern East Asia spring | 47,063 | 6,568 (13.96%) | 5,351 | **11.37%** |
| Northern East Asia summer | 35,072 | 4,481 (12.78%) | 3,245 | **9.25%** |
| Southern Australia spring | 379,937 | 84,266 (22.18%) | 67,646 | **17.80%** |
| Southern Australia summer | 134,793 | 25,636 (19.02%) | 16,502 | **12.24%** |
| Tropical South America calendar Q1 | 24,649 | 2,497 (10.13%) | 1,519 | **6.16%** |
| Tropical South America calendar Q3 | 26,648 | 1,648 (6.18%) | 1,105 | **4.15%** |

The global no-phenology-annotation filter returned **23,518,398 (90.21%)** of the observed photo-backed flowering-plant records, while **424,672 (1.63%)** were explicitly annotated **No Flowers or Fruits**. The former should **never** be interpreted as proof of no flowers. The sequential API queries were not one atomic snapshot; any-stage plus without-stage counts fell short of the separate whole-group total by **618 records**. The fixed regional/season cohorts also span different numbers and compositions of species and observers, so neither within-region differences nor between-region contrasts isolate an annotation-behaviour causal effect.

**Single-paper decision:** The Angiospermae-only count is the taxonomically appropriate *platform-context* sensitivity for the main manuscript's discussion of photographed flowering-annotation ascertainment. The earlier Plantae count remains in this document as the broad-plant baseline. Neither platform cohort is a probability sample of the original 42,111 FCP named source taxa, original licensed photo IDs, images containing a classifiable bloom or wild flowering-plant individuals. The existing Mois historical white–BIO12 association cannot be reweighted from aggregate API counts. After actual source-ID admission and ROI measurement, use the *independently joined original observation UUID annotation overlay* and report Flowers-positive/ROI-classifiable versus unannotated/ROI-classifiable cohorts by genus, site, season and observer.
