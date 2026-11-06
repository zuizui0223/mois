# MOIS — nighttime moth pollen studies and 32 unclassified flower colours (2026-11-06)

## Status: data-source recovery, no new ALAN effect asserted

Previous source-verified Walton et al. (2020) nocturnal macro-moth pollen carriage source: 43 original GloBI botanical binomials, 40 WCVP-accepted MOIS photographed species, 3,481 existing successful original ROI9 photos in 2019–2024. Independent Island frozen species-direct single-colour evidence: **5 white, 3 nonwhite, 32 unresolved**. The 32 unresolved species have **2,707** previously successful original-photo ROI9 records, including 76 original 2020 same-year VIIRS-photo matches and 12 original 2024 same-year VIIRS-photo matches. No individual image recoloured or retrieved. Pollen carriage is not verified stigma deposition, seed output, or exclusive night pollination at unrelated photo sites.

## Two independent research datasets probed

Dryad public official API metadata and version numbers are confirmed in [run 42167793678](https://github.com/zuizui0223/mois/actions/runs/42167793678) SUCCESS: Diniz (2025) DOI 10.5061/dryad.sxksn03gp, version ID 405111; Atwater DOI 10.5061/dryad.gxd2547xh, version ID 412555. The actual Dryad source-file download returned **401 Unauthorized**, even after fixing a Dryad HAL file-ID parsing mistake. No rows from these supplementary files have been used or invented, and no attempts to bypass authorization were made. Diniz plant ASV taxonomic identifications may be higher-rank; do not turn genus- or family-rank ASVs into binomial species.

Independent publicly accessible Atwater et al. (2025) paper Table 2, DOI 10.1111/een.13446 (Florida sandhill moth body pollen and floral-visitation study): **7 named plant species**, 151 reported flower-visitor and 28 pollen-carriage records, different kinds of evidence. A new original WCVP/photo/VIIRS numeric-only [run 42168801373](https://github.com/zuizui0223/mois/actions/runs/42168801373) **SUCCESS**, 5 unit tests and SHA-locked original photo/observation IDs and 2020/2024 VIIRS sources. Source-seven overlap with MOIS 2019–2024: **2 accepted plant species / 66 original strict9-success photos / 1 independently documented white species**. Same-year satellite overlap: 2020 = 2 photos / 1 species / **0 documented white species**, 2024 = 3 photos / 1 species / **0 documented white species**. No new positive ALAN–white contrast is obtained. Source SHA artifact 13480242722 ZIP SHA256 8e3b50af502b9503cc9f0913743938dbdec905f004f6ac8c063a8d08c458689e.

## Private botanical literature triage: 32 unresolved Walton-linked species

Separate species-name-focused, publicly accessible botanical descriptions (Kew POWO, RHS, other flora) were associated with the 32 previously unresolved species, with photo counts and same-year VIIRS. **All are exploratory morphology/source notes, NOT accepted strict single-state labels** and not all 32 linked pages were individually HTTP-verified.

| External-source preliminary flower-state description | Species | Existing original strict9 photos | 2020 VIIRS photo matches | 2024 VIIRS photo matches |
|---|---:|---:|---:|---:|
| Typically nonwhite (does not exclude rare white individuals) | 17 | 1375 | 35 | 7 |
| Variable white / pink / purple morphs | 10 | 998 | 31 | 5 |
| Cream / off-white, requires measurement threshold | 2 | 146 | 2 | 0 |
| Anthers / flower organ ambiguous | 1 | 130 | 5 | 0 |
| Source only documents genus, not this species | 1 | 52 | 3 | 0 |
| Petal-specific phenotype inappropriate (willow catkins) | 1 | 6 | 0 | 0 |
| **Total** | **32** | **2707** | **76** | **12** |

Representative source-specific safeguards:
- [Kew: Convolvulus arvensis](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A265849-1/general-information) has white/pink corolla forms. [Kew: Alisma plantago-aquatica](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A30019518-2/general-information) describes white/pink petals; neither implies a species universally white.
- [RHS: Plantago lanceolata](https://www.rhs.org.uk/plants/13217/plantago-lanceolata/details) describes brown flowers with white anthers. White anthers cannot be treated as white petals or guaranteed animal-attractive floral colour; RHS also notes wind pollination.
- [RHS: Glechoma hederacea](https://www.rhs.org.uk/plants/8034/glechoma-hederacea/details) includes typically violet-blue flowers but pink/white variants. [Kew: Lycopus europaeus](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A449641-1/general-information) notes white corolla with dark purple spots.
- [RHS: Heracleum sphondylium](https://www.rhs.org.uk/plants/8656/heracleum-sphondylium/details) is off-white, [RHS: Sambucus nigra](https://www.rhs.org.uk/plants/16391/sambucus-nigra/details) is cream; do not silently call either a pure white single-state botanical species.
- [New Zealand Plant Conservation Network: Salix cinerea](https://www.nzpcn.org.nz/flora/species/salix-cinerea/) records yellow/green male/female catkins, not a directly comparable petal-ROI phenotype.

**New formally accepted independently documented single-state WHITE species added to original Island ledger = 0.** Prior 5 white species unchanged. The local private source-review CSV includes 32 public source URLs and measurement priorities; it was not published in the GitHub repository and contains no modified/redistributed Island species-by-trait master ledger. Private source priority CSV SHA256 fdac5c3c73d719ab2e162cf7c08ca940aa736ec18165e7e5418517ebbc4b14a2; aggregate summary CSV SHA256 8d254508fb46da0819d8180333b699c7b6f12b1350e6bc497aa6cc41b4620bcf.

## Next scientific gate

Species-level *measured* nocturnal visitation/pollen deposition (distinct from generic moth flower-visit associations), reliable plant corolla/petals including polymorphism, and independent diurnal comparison at exactly matching VIIRS years. Source-Walton UK moth pollen carriage ≠ local effective nighttime pollination at all MOIS photo sites. Source-only verification did not download original MOIS photographs, rerun ROI9, replace unknowns with white or nonwhite, run the previously rejected MNAR tipping sensitivity, or merge PR #199 into main.
