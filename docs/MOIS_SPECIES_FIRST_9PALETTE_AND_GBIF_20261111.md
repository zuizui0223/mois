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

### 2026-11-11 GBIF 45 COL XR accepted taxa: actual EPSG:6933 100 km occupancy

[Real 45 species pilot (run 42950078270)](https://github.com/zuizui0223/mois/actions/runs/42950078270) completed SUCCESS. An outcome-blind fixed hash of exact COL XR species keys selected **five species per nine photographed dominant hue** from the signed 26,852 accepted taxon roster, 45 species total. The GBIF `occurrence/search` API was actually queried, paging filtered georeferenced occurrence records to completion for species with <=1,200 qualifying records.

- **28** species had all filtered occurrence records enumerated; the resulting species × 100 km EPSG:6933 occupied grid table contains **1,296 rows**.
- **17** species exceeded the deliberate 1,200-record completeness cap, were not assigned truncated geographic ranges, and remain `above_enumeration_cap_no_grid`.
- A true accepted name/single nine-colour *photographic trait* appears as a species-level overlay; none of the GBIF occurrence points themselves was used to measure its local flower colour.
- Actual authenticated source roster, 45-species colour-balanced selection receipt, 28-species grid table and coverage/missingness statuses: [artifact 13825160751](https://github.com/zuizui0223/mois/actions/runs/42950078270/artifacts/13825160751).

The worldwide 26,852 COL XR accepted species distribution cube still requires GBIF **account-authenticated SQL downloads**. The 108 source-backed `SQL_TSV_ZIP` query JSON payloads are already prepared, and a representative query passed the GBIF public SQL validator (HTTP 201). Validation is *not* submission or execution. Official GBIF documentation: https://techdocs.gbif.org/en/data-use/api-sql-downloads .

## 2026-11-11 13:xx JST — Incremental GBIF COL XR naming for newly photo-coloured taxa

The previous 28,309 names were not re-queried. The verified third-wave original-photo atlas [artifact 13824971351](https://github.com/zuizui0223/mois/actions/runs/42950141749/artifacts/13824971351) contains 28,591 source nominal taxa with provisional dominant nine-colour labels (28,312 in the prior source atlas + 279 new hue-labeled taxa). These 279 alone were submitted to the same official GBIF v2 COL XR species match batch, keeping all unresolved names and source statuses separately.

New **real** exact accepted botanical COL XR species names: **267**, remaining **12** not auto-accepted; cumulative source-backed photo9 × accepted species name keys: **27,119**. The original nine continuous photographic colour fractions and previous 26,852 exact accepted species identities are preserved. All accepted keys remain unique, and nine-hue photographed polymorphism evidence states are kept. This is based on the [completed successful incremental run 42956460881](https://github.com/zuizui0223/mois/actions/runs/42956460881), [artifact 13827590735](https://github.com/zuizui0223/mois/actions/runs/42956460881/artifacts/13827590735). The new complete current accepted-species list plus all nine measured fractions is included there, alongside the 109 global ISEA3H-resolution-8 SQL download request JSONs.

**GBIF occurrence status:** The 109 SQL requests were prepared but **not submitted**; official GBIF authentication is needed to run the full 27,119-species occurrence cube. The smaller independent GBIF 45-species test had previously verified 28 completely enumerated species/1,296 100-km species-grid cells, and skipped 17 high-count species with explicit sampling completeness limits. Do not interpret the 27,119 accepted name keys as 27,119 fetched worldwide species distributions.

## Fourth-wave ROI9 append guarded against incomplete original images

The fourth 1,024-photo original ROI-v4 image inference is still in progress in [run 42950141749](https://github.com/zuizui0223/mois/actions/runs/42950141749) as of this update. A separate [source-authenticated incremental four-wave append](https://github.com/zuizui0223/mois/blob/analysis/mois-species-first-palette9-20261111/scripts/analysis/mois_species_first_wave4_incremental.py) and [GitHub Actions workflow](https://github.com/zuizui0223/mois/blob/analysis/mois-species-first-palette9-20261111/.github/workflows/mois-species-first-wave4-append-wave5.yml) were added. Synthetic regression tests passed, including preventing repeated original photo IDs/observation UUIDs, rejecting existing-colour taxon overwrites, and ensuring technical unknowns do not become measured colours.

The workflow intentionally checks the *entire real wave4 photo run's successful conclusion* before admitting a new numeric nine-palette success or generating a fifth-wave source queue. Its first source-free validation run completed successfully but **skipped image-result append** because wave4 had not yet completed. When genuine wave4 image terminal results are available, the code can verify four signed 256-photo receipts and produce a new 42,111-taxon catalogue and a disjoint next 1,024 *first-photo-of-previously-uncoloured-taxon* queue. The source data and software were not used to invent wave4 ROI success counts; no newer cumulative species-coverage count is claimed yet.

**Unchanged inference boundary:** All ROI values come from original observed image pixels with computer-vision flower-region masks. A given photo's dominant palette is an *operational photographed colour label*, not certainty of species fixation, pigment concentration, UV reflectance, local colour at GBIF occurrence points or native/wild plant abundance.


## 2026-11-11 fourth original-image wave signed and fifth source/GBIF COL XR updated

The completed original photo fourth-wave ROI-v4 inference [run 42950141749](https://github.com/zuizui0223/mois/actions/runs/42950141749) produced **1,024 never-before-attempted distinct species-photo/observation pairs**: **378** authentic strict original nine-colour fraction vectors, of which **290** had a clear dominant one-of-nine palette, and **646** were technical unknown. All four signed original worker terminal SHAs and original photo/observation/taxon identities were confirmed before adding any numerical colour to the catalogue.

A distinct, SHA-verified append job [run 42959597251](https://github.com/zuizui0223/mois/actions/runs/42959597251) completed SUCCESS, with an immutable delta plus new independent fifth-wave source. After the true fourth-wave delta:
- Original-photo processing attempts: **1,184,474** unique photos.
- Genuine nine-palette numerical ROI outcomes: **452,836** measured photos.
- Source nominal taxa having any strict9 original-image evidence: **31,235** / 42,111.
- Source nominal taxa with a recorded dominant photo nine-palette label: **28,881**.
- Source nominal taxa still without any genuine nine-palette ROI: **10,876**.
- The next batch contains **1,024** as-yet-uncoloured taxa, each with one new photographic/observation ID, zero ID overlaps or already-coloured species. It was signed into four 256-original-photo source shards.

The fifth wave's true original-pixel inference was **started** in [run 42959597251](https://github.com/zuizui0223/mois/actions/runs/42959597251). Its returned success count must not be assumed until all image workers finish and their signed terminal receipts are authenticated. The 42,111-taxon updated atlas, original strict9 delta, fifth-wave original source and audited counters: [artifact 13828163364](https://github.com/zuizui0223/mois/actions/runs/42959597251/artifacts/13828163364).

### Actual GBIF v2 COL XR accepted-species incremental name and trait linkage

Only the **290** new, photographed clear-hue source taxa added in fourth-wave photos were sent through current GBIF v2 COL XR batch matching, while preserving the exact prior **27,119** COL XR accepted species keys without re-querying them. [Run 42959945602](https://github.com/zuizui0223/mois/actions/runs/42959945602) completed SUCCESS with **275 additional strict exact accepted species** and **15** unresolved/review taxa, for **27,394** accepted unique COL XR species keys paired to all nine observed colour fractions and photographed multi/single-hue evidence states. These do not imply colour fixation or local flower colour at each GBIF occurrence location.

The complete current COL XR accepted-species × original nine-fraction trait table, unresolved names, and **110 prepared** ISEA3H-resolution-8 worldwide grid-download SQL request JSONs are [artifact 13828078781](https://github.com/zuizui0223/mois/actions/runs/42959945602/artifacts/13828078781). All SQL files remain **unsubmitted**: they require GBIF-authorized download credentials; no fictitious worldwide grid occupancy has been inferred from accepted names. Prior real pilot for nine colour classes queried 45 accepted species; 28 with fully enumerated presence data yielded 1,296 species × 100-km equal-area grid rows, whereas 17 high-occurrence species remained explicitly outside completeness cap.

These numbers represent **photo-based operational trait coverage**, not a complete census of botanical species, native flora, true petal hue, species-wide fixed colour, or globally unbiased prevalence estimates.


## 2026-11-11 current: signed fifth originals in progress; sixth 4,096 species-first source ready

**Current last completely integrated original photo baseline remains wave4:** 1,184,474 exact distinct original photo attempts, 452,836 authentic numerical nine-ROI photos, 31,235 nominal source taxa with strict9, 28,881 nominal source taxa with dominant nine-palette and 10,876 with no valid strict9. [Successful wave4 append and fifth source receipt](https://github.com/zuizui0223/mois/actions/runs/42959597251/artifacts/13828163364). **Fifth** genuine photo worker for exactly 1,024 different previously zero-ROI9 taxa remains in progress in [run 42959597251](https://github.com/zuizui0223/mois/actions/runs/42959597251); do not add unknown fifth success counts.

As of the signed fourth-wave COL XR incremental reconciliation, 290 newly photo-hue-labelled nominal taxa generated **275** strict exact accepted COL XR species keys, 15 remaining unresolved, and **27,394** total accepted species with nine authentic photo palette fractions, [completed run 42959945602](https://github.com/zuizui0223/mois/actions/runs/42959945602) and [artifact 13828078781](https://github.com/zuizui0223/mois/actions/runs/42959945602/artifacts/13828078781). There are now **110 prepared but unsubmitted** ISEA3H-resolution-8 GBIF occurrence download SQL files; accepted names are not worldwide species-grid occupancy.

For continued sampling efficiency, a [new source-only fifth→sixth incremental module](https://github.com/zuizui0223/mois/blob/analysis/mois-species-first-palette9-20261111/scripts/analysis/mois_species_first_wave5_sixth4096.py) and [four-parallel-worker sixth-wave workflow](https://github.com/zuizui0223/mois/blob/analysis/mois-species-first-palette9-20261111/.github/workflows/mois-species-first-wave5-append-sixth4096.yml) were added. The module requires all four SHA-certified fifth original pixel worker receipts, validates every old/new original photo and source observation UUID, appends only genuine nine-fraction numeric ROI successes to the 42,111-taxon atlas, then ranks the next **4,096 distinct still-uncoloured source taxa**, each with a new independent original photo/observation ID. Source partition design is 16 × 256-photo shards, four runner groups of 1,024 images, at most four GitHub Actions workers in parallel. This prevents repeated serial 1,024-image pipeline starts and repeated evaluation of already classified taxa.

The [initial sixth workflow check run 42960943494](https://github.com/zuizui0223/mois/actions/runs/42960943494) completed source-only unit tests **SUCCESS**. As the fifth real RGB inference had not completed, the source gate correctly skipped the append and sixth real image workers. **The 4,096-photo sixth work was prepared but has NOT yet downloaded or inferred RGB.** Do not count these candidate photos or taxa as measured until the gate passes and all 16 real FCP-v4 source-unit terminals are checked.

The image-derived nine palette labels are photographic flower-region traits, not proof of petal-only colour, physiological pigments, species-wide fixed colour, wild/native floral prevalence, or local colour at individual GBIF occurrence locations.

## 2026-11-11 15:12 JST — fifth actual 1024 authenticated; sixth 4096 real photo acquisition active

The [fifth actual-image worker run 42959597251](https://github.com/zuizui0223/mois/actions/runs/42959597251) finished **SUCCESS** at 14:31 JST, after all four signed 256-photo original terminal receipts were reconciled. Its **1,024 disjoint photos, 1,024 observations, and 1,024 previously zero-ROI taxa** gave **383 genuine nine-fraction numerical original flower-region ROI successes**, **287 clear sampled dominant nine-palette hues** and **641 technical unknowns**. Original worker receipt byte SHA-256: shard 0 `159330b6e56f33d9783cbd674ae1198301264f7107afebe5815de20474e8969f`, shard 1 `441d0be2054926e1088cf91e119cba043baffca8a714a4bef01616a262869b7d`, shard 2 `e78549c060bf3adcfe9ccd066761344d5d7d1ccbab8bb2a6431bd4151f0dfb27`, shard 3 `c32bda1ecdd8e42278c1250cbbf9ccb7d5bdebbaa5eb7e63091cb05e91d63e00`. Failed/mixed ROIs remain unknown and must not be imputed as measured colours.

[Run 42970701175](https://github.com/zuizui0223/mois/actions/runs/42970701175), fifth authenticated append source job **SUCCESS**, replayed these terminal receipt SHAs against each frozen source shard, checked original photo/observation/taxon IDs, finite nonnegative nine-palette fractions summing to 1 within tolerance, previous atlas and all five old genuine photo sources against pinned SHAs, and wrote [artifact 13832629260](https://github.com/zuizui0223/mois/actions/runs/42970701175/artifacts/13832629260). Signed fifth-append cumulative outputs: **1,185,498 distinct original photograph attempts**, **453,219 genuine nine-fraction original photo ROIs**, **31,618 of 42,111 nominal source taxa with at least one numerical ROI9**, **29,168 taxa with a clear sampled dominant palette**, **10,493 still with no ROI9**. Atlas SHA-256 `650fa24c27d972b7c67fac99734a873b489a2a51be1a0d34c59005d68dbe357c`; signed true fifth original delta `caa457e70b28a14477eb6727f53a273db74342bec625a4bf89ec74f9e28ecce7`.

**Sixth data acquisition is actually started, not just planned**: the same run has launched four parallel FCP ROI-v4 actual-pixel inference groups, covering **16 × 256 = 4,096 previously unattempted distinct source taxa**, one distinct original photo and observation per taxon. The sixth source has **zero retries**, retains all old photo/observation ID exclusions, and is frozen in the above signed artifact (sixth queue SHA-256 `22ad37aae2e3d48bce0578754bf6f3a921e0ef9610e8a01e24be09727ec7bd94`). Do not count any sixth ROI success before all sixteen terminal receipts have completed and are SHA-audited; the new `audit-sixth-4096` job in the same workflow enforces this requirement.

GBIF accepted-species bridge still has **27,394** strictly matched COL XR keys as of the fourth-wave reconciliation; the **287 new fifth-wave photo hue labels** are pending their *incremental* COL XR matches and must not be called accepted species yet. Previously prepared 110 worldwide SQL requests remain **unsubmitted**, without global GBIF occurrence download or DOI. The 9-colour photographic colour phenotype is not proof of species-wide fixation, local flower colour at any GBIF point, or petals-only native wild phenotype.


## 2026-11-11 15:08 JST onward — automatic next-wave watchdog + fifth numerical photo append

The user explicitly requested that each completed run start the next available original-photo acquisition automatically and that the sequence be monitored on a schedule. An enabled **hourly conditional automation** was created in Japan time, named **MOIS 9色取得・自動次波監視** (automation ID `Automation_11db7d3e5a7c8191a90ef2e889de709a`). It is instructed to check actual signed GitHub Actions image terminal receipts, guard against repeated original photo IDs and observation UUIDs or simultaneous identical source batches, append only verified 9-component ROI values, then progress to the next distinct source taxa, investigate errors, update PR #200, and notify only upon real changes. This is a scheduled orchestration and monitor, not a promise that inference succeeded without its terminal receipt.

At this turn the **fifth** original-photo ROI-v4 worker completed [run #42959597251](https://github.com/zuizui0223/mois/actions/runs/42959597251) **SUCCESS**. Its 1,024 independently selected new nominal source taxa produced **383 genuine strict numeric nine-palette ROI photos, 287 clear photo dominant palette labels and 641 technical unknowns**.

Following success, workflow [run #42970701175](https://github.com/zuizui0223/mois/actions/runs/42970701175) passed original-image and independent-observation ID checks. The finished fifth-wave append/frozen-sixth-source job reported:
- **1,185,498** cumulative exact original-photo attempts.
- **453,219** genuine strict nine-palette original-photo numerical vectors.
- **31,618** of 42,111 nominal photo-source taxa have at least one strict9 original photo.
- **29,168** nominal taxa have a provisional photographed dominant palette hue, and **10,493** do not yet have strict9 evidence.
- **4,096 different previously uncoloured source taxa**, each with one never-before-attempted original photograph/independent observation. The full queue was SHA-frozen in 16 × 256 photo-source shards, with four concurrent GitHub Actions runner groups of 1,024 photos each.
- Source appendix, fifth genuine original-photo numerical photo delta, receipts and source for sixth wave: [artifact 13832629260](https://github.com/zuizui0223/mois/actions/runs/42970701175/artifacts/13832629260). The four sixth image workers are now **running real original ROI-v4 inference**. The sixth-wave 4,096-photo successes are **not yet counted**; workflow includes a post-worker audit to verify each original SHA and preserve technical unknowns.

A distinct incremental GBIF v2 COL XR official name-matching [run #42971183953](https://github.com/zuizui0223/mois/actions/runs/42971183953) queried only the fifth-wave's **287 new clear photographed hue taxa**; **272** were accepted as exact botanical species and 15 are withheld for synonym/fuzzy/unmatched/uncertain status. Prior **27,394** accepted names were reused without requerying. Cumulative: **27,666** accepted COL XR species with nine real photographic ROI fractions, with **111** prepared but still *unsubmitted* global ISEA3H resolution-8 SQL occurrence download JSON files. [GBIF fifth incremental matching and nine fractions](https://github.com/zuizui0223/mois/actions/runs/42971183953/artifacts/13833417517).

**Scope** remains global species breadth of photo traits, never the full 93.2-million original image downloads. GBIF full 27,666-species distribution cube is not yet downloaded due to an authenticated GBIF account requirement. A species-level hue should not be shown as observed at each GBIF location or treated as proven taxon-wide fixed/genetic flower colour.

### 2026-11-11 15:17 JST — true fifth-only GBIF COL XR name match completed

The newly launched [GBIF fifth-only actual name match run 42971343937](https://github.com/zuizui0223/mois/actions/runs/42971343937) finished **SUCCESS** with original fifth 42,111-taxon atlas SHA and the prior wave4 COL XR crosswalk both independently pinned by their calculation receipts. It queried **only 287 newly hue-labelled photographed source taxa**, of which **272** were exact accepted rank-SPECIES Plantae in the current COL XR checklist and **15** were left unresolved/review. The previously accepted **27,394** original source species keys were reused rather than re-queried, yielding **27,666** uniquely accepted COL XR species keys paired with nine-component authentic photographed fractional traits. The [signed output artifact 13832704850](https://github.com/zuizui0223/mois/actions/runs/42971343937/artifacts/13832704850) also includes the complete accepted crosswalk and **111 prepared, unsubmitted** ISEA3H r8 worldwide GBIF SQL request files. No account-authorized global GBIF bulk download, DOI or complete species-by-grid occupancy has been produced. All photographic colour interpretations retain the existing limitations and unresolved names are not force-matched.


## 2026-11-11 16:24–16:35 JST — sixth real 4,096 originals audited; seventh real 4,096 started; sixth GBIF COL XR complete

[Sixth image run #42970701175](https://github.com/zuizui0223/mois/actions/runs/42970701175) is **SUCCESS**, including four actual original-photo ROI-v4 worker groups and the audit of all **16 × 256** photo source shards. [Audit artifact 13836871586](https://github.com/zuizui0223/mois/actions/runs/42970701175/artifacts/13836871586): **4,096** independent original photos / observations / nominal source taxa, **1,406 genuine strict nine-component original ROIs**, **1,087 clear sampled dominant nine-palette labels**, **2,690 technical unknowns**. Each source/receipt SHA, photo-ID, observation-ID and taxon-ID crosswalk was verified; strict numeric9 has nine finite nonnegative fractions summing to 1 ± 0.001, while unknowns retain no fabricated colours.

[Seventh acquisition run #42980109440](https://github.com/zuizui0223/mois/actions/runs/42980109440) passed all regression tests and signed sixth-append checks and **actually launched four simultaneous real-image seventh-wave workers**. No old original image bytes were downloaded or re-inferred. Cumulative authenticated append receipt:
- **1,189,594** actual original-photo attempts (1,185,498 + 4,096).
- **454,625** original nine-component numeric ROI observations (453,219 + 1,406).
- **33,024 of 42,111** source taxa have at least one valid numeric9 photo (31,618 + 1,406).
- **30,255** source taxa have at least one clear sampled photo-hue label (29,168 + 1,087).
- **9,087** taxa still lack a valid numeric9 observation.
- Signed complete sixth atlas SHA256: `df4d4ce1944182be7dbf30ea072168557aa9d3222540690f6cc56bda9deab640`; true sixth photo delta SHA256: `4cf743662443d54bedd86672182a561e391c1df38e8a0f793ceca4ae4eae73e9`. Audited sixth receipt SHA256: `04c27213ff381d1a99e1c58744377ede272acd9871edeb8a4dde957e46ed7e3a`.

**Seventh frozen photo source:** 4,096 distinct new original image/observation/taxon opportunities, including **3,353 never targeted source taxa first**, and **743 still unresolved taxa each receiving a *different previously unattempted original photograph***. No previously clear-labelled source taxon, photo ID, or observation ID is reused. It is split into sixteen signed 256-photo shards with four concurrent FCP YOLO11/EfficientSAM ROI-v4 jobs. Queue SHA256: `9bd88f8da5f782b07a580aff5b931a75083395368005b83b225b9e7293825389`. [Signed sixth atlas and seventh new photo work artifact 13836593480](https://github.com/zuizui0223/mois/actions/runs/42980109440/artifacts/13836593480). **Do not count the seventh numerical outputs or launch wave 8 before all workers and final seventh audit have completed successfully.**

Independent [GBIF v2 COL XR sixth-only name matching run #42980416767](https://github.com/zuizui0223/mois/actions/runs/42980416767) **SUCCESS**: the **1,087 newly clear genuine photographed hues** were matched *incrementally* without re-querying 27,666 earlier accepted names; **1,021 exact accepted species**, **66 newly queried names withheld**, cumulative **28,687** unique GBIF COL XR accepted botanical species linked to nine-component photographed hue traits. [Signed GBIF crosswalk and nine-fraction data artifact 13837940427](https://github.com/zuizui0223/mois/actions/runs/42980416767/artifacts/13837940427). Global ISEA3H resolution-8 SQL requests: **115 prepared, unsubmitted**. No account-authorized GBIF bulk occurrence download, DOI or locality-specific flower colour has been produced. A sampled photo hue is **not** evidence of taxon-wide genetic hue fixation or the photographed colour at each GBIF occurrence coordinate.
