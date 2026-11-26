# MOIS 2026-11-19 — Priority 16 species: authentic original CC-photo contact visual gate

## Completed proof of actual source acquisition

[Repaired GitHub Actions #44547129368](https://github.com/zuizui0223/mois/actions/runs/44547129368) **SUCCESS**; signed artifact #14541906843, ZIP SHA256 `a6c96f5e2db118bef282688a237d720362d651a16a611ca01907ad12c5cb0998`.

- Source target: first **16 species / 16 genera** of the pinned 199-genus WCVP-native iNaturalist photo9 qualification queue. Queue exact Git blob `efc85eb9cc9fa842d630d69ab1962fa7a5029ab1`; source accepted 37,104 species iNaturalist taxon-ID crosswalk ZIP SHA256 `d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc`.
- Obtained **48 licensed original iNaturalist image bytes**, 3 per species from **48 independent observation UUIDs**. All photo IDs, observation UUIDs and original image byte SHA256 are unique in this wave; no image retrieval failures and no source API errors.
- Historical photo/observation exclusion: 1,231,337 prior photo IDs and 1,083,073 prior independent observation UUIDs, plus 93 source raw-image SHA. **The full set of later 2026 waves was not excluded**, so never call these guaranteed new to all MOIS sources.
- 16 original-image-derived contact JPGs produced. The original full photo bytes were used transiently for SHA and resized visual copies, not retained in this artifact. Photo licences: **34 CC BY-NC, 10 CC BY, 2 CC BY-SA, 1 CC0, 1 CC BY-NC-SA**. Each contact panel includes source photo ID and photographer login; the full original URLs, per-photo SHA and observation UUIDs are in `actual_CC_original_photo_UUID_SHA_metadata.csv`.
- **No botanical petal mask, independent flowering-stage proof or quantitative new 9-colour ROI computed**. None of the 16 species gains `sampling_consensus_single` status from this wave. PR #200's existing hue catalogue is unchanged.

## Image inspection: triage of photos, not confirmed petal colours

The three contact images per species were directly reviewed as a **rapid visual screen** (no blinded expert labels, no per-organ mask, no colourimetric calibrations). The question is *whether conspicuous flower corollas/perianth can be read at this image resolution*, not whether a plant bears flowers in general. Inconspicuous grass or sedge florets are not equivalent to a photographically measurable coloured corolla.

| Source species | Original source photo9 label (provisional) | Visual triage of 3 new licensed photographs | Status |
|---|---|---|---|
| *Aristida adscensionis* | pink | Grass habit / inflorescence, no visually measurable pigmented petals | No clear petal |
| *Lonicera caerulea* | white | Two images with pale cream flower corollas; another leafless branch image | Some clear corolla images |
| *Eleocharis uniglumis* | yellow | Sedge culms and spikelets; no visible coloured perianth | No clear petal |
| *Spondias mombin* | yellow | One image shows small pale flowers; two images show yellow/orange fruits | Limited / mixed organ |
| *Eugenia punicifolia* | red | Three images dominated by yellow/red berries and leaves | No clear petal |
| *Eremophila latrobei* | red | Red tubular flowers visible in two images; one shrub/habit image | Some clear corolla images |
| *Pedicularis lapponica* | white | Dark buds/spent structures versus one apparent pale open corolla | Limited / stage ambiguity |
| *Chamaecrista absus* | white | Mostly leafy shoots or developing structures; possible tiny flowers not readily resolved | Limited / ambiguous |
| *Hedysarum hedysaroides* | magenta | Dry pods/fruiting stems in all three | No clear petal |
| *Sporobolus pyramidalis* | yellow | Grass habit, panicles and dried seed structures | No clear petal |
| *Amphilophium paniculatum* | yellow | Violet/purple-toned floral structures and a visibly two-toned tubular corolla; source colour cannot be confirmed from these photos | Some clear corolla images |
| *Glycine canescens* | pink | Tiny lilac/pink flowers visible in one close-up, remaining images mostly vegetative | Limited / mixed stage |
| *Rumex fueginus* | yellow | Habitat, green inflorescences, mature structures, no discernible coloured perianth | No clear petal |
| *Carex glacialis* | yellow | Vegetative mat and sedge reproductive spikes, no visible coloured perianth | No clear petal |
| *Lythrum virgatum* | magenta | Dry structures in early images, distant purple bloom possible in one | Limited / ambiguous |
| *Banisteriopsis muricata* | red | Two photos visibly contain pale pink/white flowers, third is less resolved | Some clear corolla images |

**Visual triage totals**: four taxa with some clearly interpretable open corollas, five with limited or ambiguously resolved petal evidence, and seven with no photographically readable coloured petals in these 3-photo contact sheets. These counts are *not validated flower phenology labels or species-level colour frequencies*. All source labels are from a separate, prior collection of photo ROI proxies and cannot be directly called correct/incorrect solely by contrasting these new photographs.

**Why this matters:** nonflower organs (fruits, seeds, leaves, habitat images) survive even research-grade species-ID and CC licence screens. The species-first image process must explicitly gate **flowering observation** and then **actual visible petal/corolla ROI** before counting that image toward a plant's flower hue.

## Next operational gate

1. Retrieve extra species-matched independent observations **annotated as “Flowers” in iNaturalist** using `term_id=12&term_value_id=13` (documented official iNaturalist Help Center), and exclude all 48 photos/UUIDs from this wave. Lack of a flowering annotation means *unknown*, not confirmed absence of flowering.
2. Inspect photographs for visible flower organs, other tissues, colour calibration and blur/background. Preserve `no visible petal / uncertain / petal visible` plus per-photo reviewer decision.
3. Re-run the 9-hue ROI on **accepted floral-region source images only**. Do not infer flower-petal colour from red fruit or brown inflorescences. Cross-species LGM hue models remain exploratory until source colour classes are independently validated.
4. Later-stage strict 75-genera QA can be expanded only after a per-image botanical organ gate; adding nonfloral images does not increase effective valid species colour sample.

## Reproducibility

Code: `scripts/acquisition/mois_20261119_targeted_16_native_photo_contact_proof.py`. Workflow: `.github/workflows/mois-targeted-16-native-photo-contact-20261119.yml`. First source run [#44545919574](https://github.com/zuizui0223/mois/actions/runs/44545919574) was correctly diagnosed as a **false-success zero-image source** (7 source oversized JSON errors, 27 URL validation errors); fixed commit `8544fda27e44766b36be2e51fc3cb7a765affdad`. The second completed source run [#44547129368](https://github.com/zuizui0223/mois/actions/runs/44547129368) is **48 genuine original CC bytes**.

No changes to original WCVP-native classification, 199-genus prioritization, photos already measured, 9-hue labels, GBIF distribution or draft PR #206 climate-history model. Independent original petals, historical survival and floral pigment adaptation are still not established.
