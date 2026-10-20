# MOIS: completed p0004 image-level coverage and missingness audit (2026-10-20)

## Empirical status and source

This audit independently re-read **all 18 finished GitHub Actions worker ZIPs** and the immutable full p0004 pre-colour source artifact, not merely a release-summary claim.

* Frozen source: [p0004 Action 39353343403](https://github.com/zuizui0223/mois/actions/runs/39353343403), source artifact 12234737795 (ZIP SHA256 \`90db5e1335ee4c9cd5c6c33377b6f112da47867e1ee9519e003a7b56cb371750\`).
* Frozen original source manifest SHA256 \`c9e2cb8b80e696bf34c0e5f1cdba2bd3579e64c5126781c9aa75989222a437eb\`.
* The 18 worker ZIPs contain exactly the **71 original source-photo shards**. All 71 source CSV SHA256 signatures, every completed \`photos.csv\` SHA256, every terminal CSV SHA256, every photo and observation ID relationship, and every per-shard classification/failure counter were checked. The original per-photo source columns match the signed pre-colour source.
* The authentic complete RGB output and the **630,602-source-photo historical+second full wave union** are released at [mois-genuine630602-original-pixels-first-two-93m-waves-v1](https://github.com/zuizui0223/mois/releases/tag/mois-genuine630602-original-pixels-first-two-93m-waves-v1).
* No original image downloaded or recoloured for this diagnostic. Raw licensed photos were not distributed.

## Full p0004 incremental image outcomes

| Measure | Photos / IDs |
| --- | ---: |
| Newly sampled and genuinely attempted original image IDs | **18,119** |
| Acquired original images | **18,065** |
| Acquisitions failed | **54** |
| Automatically classifiable into four frozen hues | **7,085** |
| Acquired but unclassifiable ROIs | **10,980** |
| Distinct original observation IDs | **10,416** |
| Distinct source taxon IDs with at least one attempted photo | **4,390** |
| Distinct source taxon IDs with at least one four-hue classification | **2,867** |
| Source taxon IDs with at least five classified photos | **326** |
| Source taxon IDs with at least thirty classified photos | **2** |

All 18,119 records had parseable observation dates and latitude/longitude coordinates. **4,037 observations** contributed two or more photographs (**11,740 photos**), so photo rows are not independent biological observations. The 7,085 classified photos came from **5,277 distinct observations**.

The automatic photo hues were white 3,504, yellow-orange 2,316, red-pink 637, blue-purple 628. White accounted for **49.46% of classifiable automatically labelled ROIs**, not 49.46% of globally flowering native plant species. A 500-resample bootstrap by original observation gave a **descriptive** conditional 95% interval **48.04–50.74%**. It cannot correct photo visibility, taxonomic representativeness or colour-dependent missingness.

## Original geographic coverage and ROI classifiability

Latitude bands use geographic coordinates **only**: tropics from 23.44°S through 23.44°N; southern and northern extratropics otherwise. This is not a Köppen classification or a botanical realm analysis.

| Source-photo latitude band | Attempts | Auto-classifiable hues | Success among attempts | White among classifiable |
| --- | ---: | ---: | ---: | ---: |
| Northern extratropics | 14,756 | 5,911 | **40.1%** | 49.6% |
| Southern extratropics | 2,276 | 802 | **35.2%** | 48.5% |
| Tropics | 1,087 | 372 | **34.2%** | 50.0% |
| **All p0004 incremental photos** | **18,119** | **7,085** | **39.1%** | **49.5%** |

The **tropical photographs are only 6.0% of this p0004 incremental photo cohort**, and their ROIs are somewhat less frequently classified than in the northern extratropics. The same source has a strong monthly photo imbalance (April–July dominate). Do not infer that latitude causes ROI detectability: taxon mix, observation/photo selection, image composition, phenology, season and classifier suitability all differ.

### Consequence for the next acquisition allocation

**Preserve** the observation-hash-based, uncapped probability arm: this is the only arm currently offering consistent observable inclusion thresholds across the 93,209,351-photo frame. Its sequential compute batches are never valid global frequencies until a complete sample tier is finished.

A **separate enrichment arm** should preferentially supply independent plant taxa with few or no usable hue measurements, and neglected tropical/southern/seasonal/geographic cells. Freeze its selection algorithm and exact observation/photo IDs *before* seeing its image pixels. Record inclusion probabilities or explicit selection quotas, and label the enrichment arm so it is not naively pooled with the probability sample for white-flower prevalence. Retain all native/introduced/unresolved records in the metadata registry; nativeness should be an optional later analysis stratum, not a destructive acquisition filter.

The p0004 incremental cohort excludes some photos already captured in earlier FCP/MOIS sources and is not a census or unbiased complete p=.0004 global probability cohort on its own. The 42,111 source taxa include iNaturalist names and synonyms; taxon-ID coverage must not be interpreted as resolved accepted botanical species diversity.

### Reproducibility and evidence level

An independent Python audit read the saved signed source archive and 18 worker ZIPs, checked 71 source and terminal hash chains, source/terminal row identity and photo-ID uniqueness, and computed the region/taxon/month counts from all 18,119 real image terminal rows. The resulting JSON audit has the schema \`mois_complete_p0004_original_photo_coverage_diagnostic_v1\`.

All figures above concern true original-photo outcomes and their provenance; none are assumed photo-pixel colours from the 93m upstream metadata index, independently confirmed anthesis, or a global native-flora colour estimate. The denominator for a classifiable image share is the attempted source photos, whereas the white share is conditional on the automatically classifiable subset.
