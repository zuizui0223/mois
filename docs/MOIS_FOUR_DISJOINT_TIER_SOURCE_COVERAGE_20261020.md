# MOIS original-photo acquisition: four source-roster geographic and seasonal coverage audit

Date: 20 October 2026 (JST). This audit is a **photo-ID/observation-metadata check, not a flower-colour analysis**. Its purpose is to decide whether simply acquiring more randomly nested candidate photographs fixes global-seasonal sampling imbalance.

## Real source inputs and authentication

The original GitHub Actions pre-colour source ZIPs were downloaded and their original \`manifest.json\` bytes and **all source photo-shard CSV SHA-256 checksums** verified, then the exact original photo IDs and original observation UUIDs were compared across all four disjoint waves. No photo pixels were downloaded or classified in this audit.

| Original wave | Source run | Source manifest SHA-256 | Frozen new photo IDs | Original observation IDs | Source taxon IDs | Worker CSV shards |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| p0004 | [39353343403](https://github.com/zuizui0223/mois/actions/runs/39353343403) | \`c9e2cb8b80e696bf34c0e5f1cdba2bd3579e64c5126781c9aa75989222a437eb\` | 18,119 | 10,416 | 4,390 | 71 |
| p0008 | [39371199127](https://github.com/zuizui0223/mois/actions/runs/39371199127) | \`0038024ee40adfbc71aa427a112948ec5ce53cb53a7aae9879c1b32bc06eabde\` | 36,375 | 20,773 | 6,549 | 143 |
| p0016 | [39374811084](https://github.com/zuizui0223/mois/actions/runs/39374811084) | \`16baf617396dd7b388debc65fd25eb8f9041f49ca5e6f78d0b3e2898f6dfeed5\` | 73,872 | 41,745 | 9,421 | 289 |
| p0032 | [39401484641](https://github.com/zuizui0223/mois/actions/runs/39401484641) | \`3628794fbb70f9e1515dabea6db731862265c0ad6c251732cfd9ab4df4badce8\` | 148,071 | 83,835 | 12,656 | 579 |
| **Four distinct photo-roster waves** | | | **276,437** | **156,769** | *not additive across waves* | **1,082** |

All 276,437 original photo IDs are unique across the four sources. All 156,769 original observation UUIDs also belong to one source wave only. **No strict native-only filter or RGB screening filter was applied.** Source candidate IDs are not measured flower pixels or independent plant individuals.

## Geographic composition before flower detection

An observation's original latitude, not inferred botanical bioregion, determines the three latitude bands. The tropical band is **23.44°S to 23.44°N inclusive**. Every one of the 276,437 source photos had a valid point and observation date.

| New source wave | North extratropical photos | South extratropical photos | Tropical photos | Tropical fraction |
| --- | ---: | ---: | ---: | ---: |
| p0004 | 14,756 | 2,276 | 1,087 | 6.00% |
| p0008 | 29,708 | 4,564 | 2,103 | 5.78% |
| p0016 | 59,912 | 9,385 | 4,575 | 6.19% |
| p0032 | 120,315 | 18,596 | 9,160 | 6.19% |
| **Combined photo targets** | **224,691** | **34,821** | **16,925** | **6.12%** |

Of the 276,437 selected new source-photo targets, **81.28%** were in northern extratropical latitudes, **12.60%** in southern extratropics, and **6.12%** in the tropics. This near-constant geographical composition across increasing observation-inclusion probabilities reflects the **photographed source opportunity frame**, not any ecological scarcity of tropical flowering plants.

The p0032 source includes nominal taxon IDs of 7,820 in northern extratropics, 3,908 in southern extratropics and 2,343 in the tropics; these cannot be added to an accepted plant species count without Kew name and synonym resolution.

## Observation month and photographic pseudoreplication

Across the four disjoint source rosters, original observation months **April–July account for 152,211/276,437 = 55.06% of photo targets**. The 12 exact monthly source-photo counts are: Jan 8,100; Feb 9,148; Mar 17,380; Apr 38,533; May 38,705; Jun 37,758; Jul 37,215; Aug 30,139; Sep 24,290; Oct 17,071; Nov 10,369; Dec 7,729.

The 276,437 photographs belong to **156,769 observations**; more than one photo may come from the same original observation, so naive image-level confidence intervals exaggerate biological replication. Observer, taxon, observation and original calendar-day sampling must remain distinct in later white-colour / temperature / moisture / solar models.

## Next acquisition and inference design

1. **Do not interrupt or duplicate the already running p0008/p0016/p0032 full-resolution original RGB queues.** They are the comparably defined disjoint observation-hash probability arm. All completed pixel outcomes must be SHA-verified by the data-only checkpoint; queued candidates must not be counted as coloured flowers.
2. **Keep an independent, explicitly labelled geographic/taxon coverage-enrichment arm** when sourcing future images from the immutable 93m register. It should target neglected tropical and southern latitudes, winter/time-of-year coverage and underrepresented **accepted taxa**. Freeze inclusion rules and original photo/observation IDs *before* inspecting new image colours. Never let the enrichment arm contaminate the random arm's unweighted source estimation.
3. **Correct research estimates**, not raw data, for observer/photo/taxon/space/time selection, original observation clustering and hue-dependent ROI failure. If estimating a world-flora species-level prevalence, use one weighting unit per accepted taxon or an explicit hierarchical model; the observed photograph-level white fraction answers a different question.
4. WCVP/TDWG Level-3 native status is a later stratification field, not an upfront hard filter or proof that an individual photograph depicts a wild native plant. Keep introduction/unknown strata and original image-acquisition error classes.

## Status and scope

Audit outcome: **four source roster ZIP archives verified, every 1,082 source-photo shard checksum verified, zero repeated original photo IDs or observation UUIDs across these four new waves, exact region/month/taxon counts computed from source CSV rows**.

This is an **acquisition-coverage diagnostic**, not evidence of a white-flower environmental adaptation or the degree of native-flora whiteness. It cannot alone determine global biological abundance, phenology, pigment costs, flower detection recall or the true flowering status of a photographed organism.

The latest actual RGB terminal release remains [39413597007](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39413597007-v1), 667,011 genuine original-photo terminal IDs and 316,897 four-hue automatically classifiable photographs, with its own pixel/source model receipts.
