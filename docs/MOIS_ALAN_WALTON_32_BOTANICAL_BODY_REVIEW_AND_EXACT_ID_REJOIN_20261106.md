# MOIS ALAN: botanical descriptions for the 32 unresolved Walton-associated species (2026-11-06)

## Completed work and what changed

The previous 32-row source-priority list was a preliminary literature triage, not a verified set of individual flower colours. This stage inspected species-level botanical descriptions and independently rejoined the named species to the SHA-frozen original MOIS photo/observation frame, WCVP accepted names and same-year VIIRS cells.

**32 primary species descriptions were inspected. 30 species have an applicable flower/corolla colour description; two need separate treatment because flower/petal colour is not comparable with the described anthers or catkins.** These 30 species correspond to **2,571 already measured strict9-success photo IDs**, including **71 same-year 2020 VIIRS photographs and 12 same-year 2024 VIIRS photographs**. This is newly connected descriptive evidence, NOT newly acquired photographs or botanically verified individual-image labels.

The full 32-species input independently reproduced **2,707 photo IDs / 2,707 observation IDs, 76 same-year 2020 VIIRS matches and 12 same-year 2024 matches**. Every species-specific photo/observation/VIIRS count agrees with the earlier priority ledger. No photos were added by joins or counted as independent twice.

## Reviewed descriptive states

| Description class | Species | Existing strict9 photos | 2020 same-year VIIRS photos | 2024 same-year VIIRS photos |
|---|---:|---:|---:|---:|
| Source describes nonwhite flowers; does not prove absence of rare white forms | 16 | 1,233 | 34 | 7 |
| White included in a reported colour range or variant | 11 | 1,142 | 32 | 5 |
| White corolla background with coloured pattern | 1 | 50 | 3 | 0 |
| Cream or off-white kept separate | 2 | 146 | 2 | 0 |
| Organ not comparable to the petal-colour target | 2 | 136 | 5 | 0 |
| Total | 32 | 2,707 | 76 | 12 |

**Retrieval-quality caveat:** the complete adopted description class is supported by a retrieved primary page body for 31 species. For *Arctium minus*, retrieved flora bodies describe purple/deep-pink flowers, whereas the additional white-colour record is supported only by the indexed Pladias species entry; direct retrieval of that page failed. This single modifier remains explicitly lower-tier, not fully body-verified. The other 29 flower-colour-applicable species, with their entire adopted class supported in a retrieved body, account for 2,538 existing photographs and 69/11 same-year VIIRS photographs. These are evidence-tier counts, not alternative ecological effects or sensitivity scenarios.

## Consequential source corrections

1. **Cardamine amara: genus-only evidence replaced by an explicit species description.** [InfoFlora / Flora Helvetica 2018](https://www.infoflora.ch/fr/flore/cardamine-amara.html), section Description de l’espèce, describes white petals, occasionally reddish, and usually purple-violet anthers. This supports a multistate species description for 52 existing photos, including three 2020 VIIRS photo matches. The prior RHS page's generic Cardamine description was not a species-level colour record.
2. **Lycopus europaeus: spatial floral pattern is not white/pink polymorphism.** [NZPCN detailed morphology](https://www.nzpcn.org.nz/flora/species/lycopus-europaeus/), attributed to Webb et al. (1988), describes a white corolla with purple markings on the lower lip and throat. The simplified description also allows pale-mauve appearance. Record the pattern and colour range rather than infer that every individual is uniformly white. [RHS](https://www.rhs.org.uk/plants/10613/lycopus-europaeus/details) independently describes small white flowers.
3. **Cirsium arvense: nonwhite is not invariant.** [E-Flora BC / Illustrated Flora of British Columbia](https://linnet.geog.ubc.ca/Atlas/Atlas.aspx?sciname=Cirsium+arvense), Flowers section, reports pink-purple disk flowers with rare white ones. This supports retaining the white exception. Earlier failure to retrieve the Flora of China page was resolved by an independently accessible primary flora, not by assuming its contents.
4. **Plantago lanceolata and Salix cinerea: keep organs explicit.** [RHS Plantago description](https://www.rhs.org.uk/plants/13217/plantago-lanceolata/details) distinguishes brown flowers and creamy-white anthers; [NZPCN Salix description](https://www.nzpcn.org.nz/flora/species/salix-cinerea/) distinguishes pollen-bearing male catkins and female catkins. Neither is converted into a white-petal phenotype. No willow pollination mode was inferred.

The remaining species descriptions were checked in RHS species summaries, named regional flora treatments hosted by Kew POWO, PlantNET NSW, InfoFlora, NZPCN, BSBI Fermanagh and the Missouri Department of Conservation. Species-specific sections were separated from generic descriptions, images, leaf colours and unrelated taxa. The local 32-species ledger retains URLs, section locators, short paraphrases, retrieval tier and uncertainty. It does not reproduce full source pages.

## Exact original-data reconstruction

Inputs were read from already downloaded numerical artifacts. Original photo/observation identities were required to be unique in the 1,032,882-row canonical frame; its strict9-success count remains 401,703. The target remains historical_main, successful original ROI9, photo dates 2019–2024, latitude −60° to +60° and valid longitude. WCVP joins use only unique exact nominal-taxonomic IDs resolving to an accepted species. Accepted names repeated across geographical WCVP distribution rows are collapsed before joining.

The 2020 VIIRS source has 1,059 candidate cells, of which 1,044 are complete; the 2024 source has 302 cells. A photo contributes to a year only if BOTH its observation year and original 0.05-degree cell match that year's satellite source. The 30 applicable species have 2020 support in **27 species / 71 photos / 60 cells**, and 2024 support in **10 species / 12 photos / 12 cells**. The two years are not pooled as a panel of identical species or sites.

Input member SHA256:

- Canonical: f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c
- WCVP source taxon resolution: dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2
- WCVP accepted names/distribution: 9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3
- VIIRS 2020 numeric cells: 3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3
- VIIRS 2024 numeric cells: 1c5f03e5a7c355a0f7045133140c5cf2e39a5802c548f89d1c2f66c8a235db5f

New output SHA256:

- Curated 32-species botanical-description CSV: 6c32bec228f791b3cb0e61d0faab32167ea227840f703e87d6664c5dc76b557d
- Public descriptive-state aggregate: ef275dfda15e2d0aab0bd02ea7f9bf64674c87ae8dbd5c106628de23ee2951e5
- Public same-year support / photo-ID-set hash table: 8ba6aed38f1ea05deba98f4f2f8287f297182532bb2fd0accb25ec3d40521e04
- Local verification receipt: 2d73dba584ba30fe8e81ee92bbbc3ed3862f81f2f64ec047faae3bd346ef6a9e

Source-record hashes cover the curated paraphrase and metadata, not raw website bytes. All five numerical input member hashes and all output hashes were checked. Separate standard-library CSV summation reproduced the pandas group totals. Eleven local unit tests passed for source identity, ambiguous mapping, genus-versus-species evidence, organ exclusion, unknown preservation, exact-year matching and prior-count agreement. These tests do not independently validate botanical interpretation. The locally executed script, tests, source descriptions and private count linkage are in the conversation reproducibility bundle, not in the public PR. No repository-wide CI success is claimed for this stage.

## Meaning for the ALAN analysis

This stage removes a concrete evidence-management bottleneck: the 32 previously unresolved species are no longer treated as one undifferentiated missing-colour block. Colour range, aging, pattern, cream and organ mismatch can now be retained separately. It does **not** identify the colour of every photographed individual, supply a prevalence of white morphs, establish effective nocturnal pollination, or create an independent diurnal-only control group.

Existing Island fixed single-state labels and the earlier five documented-white Walton species are unchanged. No artificial white/nonwhite assignment, original-image download, original ROI inference, MNAR tipping scenario, new ALAN coefficient or causal species-loss claim was made. Existing ecological results are preserved. PR #199 remains draft; main is untouched. Public files contain aggregates, methodology and hashes; the source-specific private support table and prior trait ledger are not publicly redistributed.
