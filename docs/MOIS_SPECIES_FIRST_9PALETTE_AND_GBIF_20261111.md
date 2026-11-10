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
