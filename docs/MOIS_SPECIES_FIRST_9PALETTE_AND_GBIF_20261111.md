# MOIS 2026-11-11: species-first nine-palette acquisition

## Immediate objective
Maximize the number of named plant taxa with at least one trustworthy 9-component, original-photo colour ROI instead of attempting all 93,209,351 source photos. Use the released 1,181,146 real original-photo outcomes (451,643 strict nine-palette ROI successes) against the entire 42,111-taxon source frame. Reuse measured original photo IDs without re-downloading images or inferring pixels.

## Operational status / stopping
- `no_strict9`: zero valid nine-palette photos; first acquisition priority.
- `provisional_single`: a clear nine-palette top component in the sampled images, even if just one photo; catalogue now, defer further work until zero-supported taxa receive attention.
- `sampling_consensus_single`: at least **3 distinct observation UUIDs** agree in a single clear palette maximum (component >=0.60 and top-vs-next margin >=0.20), without ambiguous observation-level vectors; stop ordinary new sampling for this taxon.
- `sampled_multicolour`: at least two distinct independent observation UUIDs support **each** of at least two clear hue states; stop ordinary new sampling, retain all observed states instead of forcing one colour.
- `multicolour_candidate`: conflicting colours but insufficient independent observations; acquire more if feasible.
- `mixed_image_unresolved`: ROI numeric vector valid but no dominant palette; keep technical/biological ambiguity rather than forcing a single hue.

Original observation UUIDs, rather than photo IDs, are independent evidence units. More than one original photo of the same observation is collapsed to one nine-component vector. If the observation UUID is missing, a measured photo can yield a *provisional* label, but never count as an independent replicate.

**The above labels are sampling decisions, not botanical proof of a species-wide fixed flower colour.** Nine-palette ROI could include non-petal regions, mixed colour parts, white balance, bronze artefacts; native/wild status is not independently verified. The 42,111 frame is nominal source taxa, not 42,111 taxonomically accepted GBIF species.

## Files and next data steps
1. The signed original-source three files are fetched in `.github/workflows/mois-species-first-palette9.yml`; their byte SHA values are verified before any analysis. Code and synthetic tests are in `scripts/analysis/mois_species_first_palette9.py` and `tests/test_mois_species_first_palette9.py`.
2. The script exports a full taxon catalogue, a GBIF name-matching roster and a signed calculation receipt. With `--candidates` it takes an external **existing source-only candidate photo roster** and selects at most one new original observation/photo per taxon per wave, removing all previously processed photo IDs/observations.
3. Only after classification: match each scientific species name to GBIF accepted speciesKey with doubtful/synonym/genus-only matches in separate categories. Download georeferenced occurrences as a proper GBIF export with DOI, quality filtering, and bias handling; join **species × grid presence** to the palette traits.
4. Never interpret a GBIF distribution coordinate as a photographed site-specific colour measurement. Mixed or spatially polymorphic taxa require explicit treatment, not a fixed single-colour map. The active MOIS manuscript focus remains global moisture-related colour structure.

No GBIF occurrence records have yet been downloaded by this branch; no additional original photos have been processed.
## 2026-11-11 signed full93m source backfill: ZERO-ROI9 taxa first

The previously existing FCP 100-observation/taxon candidate group00 was an unsuitable way to maximize previously unmeasured species: 124,500 candidate photos from 1,245 already-source-rich taxa yielded only **nine** taxa with zero previous strict9. A source-only queue of 818 distinct taxa was created but should not be used for the immediate *breadth* objective.

Switch to signed 93,209,351-photo-ID/42,111-taxonomy metadata tiers; process only taxa with `after_strict9==0`, preserving existing 1,181,146 attempted-original-photo IDs and observation UUIDs as exclusions. The source-only code is `scripts/acquisition/mois_species_first_zero9_source.py`. It checks original full-source archive and part SHA, chooses at most four independent observation UUIDs per zero-ROI9 taxon and prefers the smallest frozen observation-hash rank; all candidate colours remain unknown until image ROI is actually measured. Only CC-licensed and supported photo image extensions enter the queue.

**Signed group00 original metadata full 14-tier run #42846334381 SUCCESS**: 5,196,555 original photo metadata rows scanned, 691 zero-ROI9 nominal taxa with at least one source photo, 664 taxa with eligible new one-photo candidates and 2,285 candidate photos including fallbacks; 27 taxa had no eligible candidate after prior-photo/UUID/license/format gating. Images downloaded and new measured ROIs: **zero**.

Global remaining15 parallel run: `https://github.com/zuizui0223/mois/actions/runs/42846819580` with a separate signed all16 taxon/photo/observation union after successful group-level jobs. The final population counts are to be read from that job's verified receipt, not extrapolated from group00. Group00 artifact: `https://github.com/zuizui0223/mois/actions/runs/42846334381/artifacts/13779455206`.

**Data meaning:** These are source image opportunities, **not confirmed new nine-palette labels**. The next acquisition wave must pass these previously unseen originals through frozen actual flower-ROI9 inference and append terminal success/unknown by immutable image and observation IDs. An ROI failure leaves the taxon unresolved and triggers its next independent fallback. Once the global species colour catalogue stabilizes, name-resolve species to GBIF accepted keys before constructing the species-by-grid occurrence dataset. Do not spread a species-level hue across all GBIF locations for a suspected geographic polymorphism.

## 2026-11-11 07:20 JST — first real256 RGB9 result and exact 91-species append

The original **first256** actual-pixel worker run [42848145801](https://github.com/zuizui0223/mois/actions/runs/42848145801) succeeded. A completely new 256-taxon / 256-photograph / 256-observation source yielded **91** genuine measured strict9 numerical flower-region ROI colour vectors and **165** technical unknowns. It did not validate all specimens as flower petals.

The subsequent SHA-verified original-terminal-versus-derived-file replay in [run 42913236517](https://github.com/zuizui0223/mois/actions/runs/42913236517), source job **SUCCESS**, found:
- **72** newly measured taxa with an unambiguous dominant nine-palette hue (photo-level >=0.60 component and >=0.20 top-to-second margin); these now have `provisional_single` labels. Other **19** had valid nine-fraction ROI but no single confidently dominant palette colour and stay `mixed_image_unresolved`.
- The frozen complete four-arm ledger remains immutable: 1,181,146 attempted photos / 451,643 strict9. The signed 256-photo delta is separate, so the new exactly photo-deduplicated cumulative count is **1,181,402 attempts, 451,734 strict9 ROIs, 30,139 nominal source taxa with any strict9**. The exactly corresponding zero-ROI9 source taxa drop from 12,063 to **11,972**. Do not call 30,139 botanically confirmed species-wide fixed colours.
- The next measurement source has **1,024 distinct new photo and observation IDs**, one photo per source taxon: **167** previously attempted unresolved taxa using different observation/photo fallbacks plus **857** newly attempted zero-ROI9 taxa. The wave2 image inference job is part of run 42913236517 and is separate from the already-completed source audit.
- 42,111-taxon updated full catalogue, GBIF *name-match preparation only*, signed 256 append-only ROI9 delta, four 256-photo next-wave source shards, and calculation receipt are in the completed source job artifact [13807253594](https://github.com/zuizui0223/mois/actions/runs/42913236517/artifacts/13807253594).

**Critical validation fix:** The original frozen worker receipt SHA refers to `measured/shard_000/photos.csv`, not the derived `measured/zero9_first256_true_original_photo9_and_unknowns.csv`. The revised replay separately verifies the original terminal byte SHA, then checks all derived photo, observation, taxon, raw status, image identity, and successful nine numerical fraction strings against that immutable raw terminal. This guard passed on genuine pixels. Three synthetic incremental tests and three original species-first catalogue tests also passed. No historical original photographs were downloaded again.

No GBIF occurrence points have yet been retrieved. A species-wide photograph label attached to occurrence points is a *species-level trait overlay*, not observational evidence that a geographically polymorphic species has that colour at each GBIF location.

## 2026-11-11 08:38 JST — wave2 genuine 1,024 measured and appended; next 1,024 targeted

The completed wave2 1,024-original-photo inference run **42913236517** succeeded (2026-11-11 08:13 JST): **363 genuine strict numeric nine-palette ROI images, 278 clear dominant-hue labels, 661 nonresponse**. The four signed original worker receipts were recovered and checked against the derived all-1024 row set, without replaying any prior image.

The subsequent signature-audited true-photo append in [run 42925844766](https://github.com/zuizui0223/mois/actions/runs/42925844766) completed the source-processing job **SUCCESS**. It separated the **357 newly colour-supported nominal taxa** from six strict9 images measured in already colour-supported taxa. Full 42,111 nominal-source-taxon conservation, no repeated photo ID and no repeated independent observation ID were verified.

| Authenticated frozen and incremental result | Count |
|---|---:|
| Cumulative distinct original-photo processing attempts | **1,182,426** |
| Cumulative authentic original numerical ROI9 successes | **452,097** |
| Nominal taxa with at least one strict ROI9 photo | **30,496** |
| Nominal taxa with no strict ROI9 photo | **11,615** |
| Nominal taxa with any provisional dominant one-of-nine palette label | **28,312** |
| Wave2 newly first-covered strict9 taxa | **357** |
| Third-wave chosen distinct new nominal taxa with not-yet-attempted original photo | **1,024** |
| Third-wave retry-unclear taxa | **0** |

The third-wave selection deliberately prioritizes taxon breadth: 1,024 source taxa without recent one-photo measurement, each with one never-before-attempted photo and independent observation UUID, and no observed dominant nine-palette label. A provisional dominant hue is sufficient to stop **routine** extra photographs operationally, not to declare botanical flower-colour fixation. Original image colour measurement for the third wave is a separate job in run 42925844766. Do not count its pixels until the terminal and signed original numerical receipt are completed.

Reproducible complete 42,111-taxon catalogue after 1,280 new original photos, 1,024-new-taxon queue, 9-colour numeric wave2 delta and 42,111 GBIF-name matching roster are included in [run 42925844766 artifact](https://github.com/zuizui0223/mois/actions/runs/42925844766). GBIF occurrence records **have not** been downloaded or joined; accepted-name reconciliation and species-by-grid distribution remain downstream tasks, explicitly separated from photos. These data are photographed flower-region colour proxies, not botanical proof of true petals, native or globally fixed colour.

## 2026-11-11 12:04 JST — actual third wave 361 new strict9 taxa appended, fourth wave source frozen

Source-proven first256 + next1024 + third1024 each used distinct original photo IDs and independent observations, and the [signed third-wave worker run 42925844766](https://github.com/zuizui0223/mois/actions/runs/42925844766) ended SUCCESS: 1,024 authentic original photo ROI-v4 attempts, **361 strict9 numeric successes, 279 confident single-palette photo labels, 663 technical unknown**.

The third-wave exact-photo append and fourth wave queue [run 42950141749](https://github.com/zuizui0223/mois/actions/runs/42950141749) completed its first job SUCCESS after preserving all 42,111 original nominal taxa and rejecting any repeated image or observation UUID:

| Metric after all three new real image waves | Audited result |
|---|---:|
| Original photo attempts including original four-arm baseline | **1,183,450** |
| Real strict nine-palette numeric image ROIs | **452,458** |
| Original source taxa with at least one strict9 photo | **30,857** |
| Original source taxa with no strict9 photo | **11,254** |
| Original source taxa with any sampled dominant palette9 | **28,591** |
| Third-wave strictly new taxa with numerical nine-colour evidence | **361** |
| New fourth-wave first-photo opportunities | **1,024** |
| Fourth-wave retries to previously attempted taxa | **0** |

Updated immutable taxon atlas / genuine 1,024 image numerical delta / 1,024 next-first-photo queue and signed receipt: [GitHub Actions artifact 13824971351](https://github.com/zuizui0223/mois/actions/runs/42950141749/artifacts/13824971351). The **fourth** 1,024 original images were submitted to the separate ROI-v4 compute job in the same run; they cannot be added before successful original terminals are inspected.

Initial run failed due to incomplete small-fixture sample size; tests were updated to preserve the mandatory 1,024-photo live count. The next failure reflected a GH artifact ZIP root-directory mismatch; source paths were corrected. A subsequent failure was traced to pandas turning missing `dominant_hue9` strings into NaN, which excluded all unresolved species from the new queue. Only string-missing status fields are normalized; numerical strict ROI9 remains untouched. A synthetic regression now protects against this selection error. The first job subsequently completed SUCCESS with all real image/source SHA checks intact.

## 2026-11-11 — GBIF source identities and actual distribution capability

The all-nine-hue GBIF v2 COL XR calibrated name matching [run 42932546397](https://github.com/zuizui0223/mois/actions/runs/42932546397) queried **28,309** named source taxa with photographed dominant hues. Of these **26,852** were *exact, accepted, species rank, botanical kingdom* with unique COL XR keys; the remaining **1,457** names were withheld from automatic joining (synonym, higher rank, fuzzy or unmatched). Earlier 138 known accepted GET matches reproduced identically under the batch query. These counts refer to a named subset of the pre-third-wave 28,312 photographic hue labels, and exclude newly named wave3 taxa pending incremental matching.

The actual [26,852 accepted-species × nine continuous photo-colour fractions trait bridge](https://github.com/zuizui0223/mois/actions/runs/42934483448/artifacts/13816937878) completed SUCCESS: **14,435** taxa have a photographed single-hue evidence status and **12,417** have multiple hue records/candidates. The latter does **not** establish spatially variable genetics; both categories are photo-derived estimates, not botanical confirmation of petals.

For occurrences, an independently executed small accepted-species pilot read complete filtered GBIF pages for **nine** species into **213** EPSG:6933 100-km species-grid records; two further high-occurrence species were excluded by an explicit completeness cap [run 42931715566](https://github.com/zuizui0223/mois/actions/runs/42931715566). A second, colour-stratified 45 accepted species pilot was started in [run 42950078270](https://github.com/zuizui0223/mois/actions/runs/42950078270); retrieve its actual success receipt before making any claim about its outputs.

For full worldwide distribution coverage, the validated GBIF SQL spatial cube is **108 requests** × max 250 COL XR accepted species keys, ISEA3H resolution 8 (roughly 7,774 km² per cell). The public GBIF SQL validator returned HTTP 201 for the original COL XR syntax [run 42933021921](https://github.com/zuizui0223/mois/actions/runs/42933021921), but **the 108 bulk downloads have not been authenticated, submitted or assigned DOIs**. Hence full global species-by-grid occupancy, rainfall associations and causal ecological conclusions are not yet observed from GBIF. A species-level photographed palette overlaid on GBIF grid occupancy does *not* imply the flower colour is measured at each occurrence location.
