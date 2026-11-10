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