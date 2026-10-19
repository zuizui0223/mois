# MOIS × Kew WCVP: real 350,399 original photograph taxa matched to independent native/introduced botanical-country records

19 October 2026. **Completed authentic acquisition**, not an island-colour hypothesis test.
Source [successful original Kew/WCVP workflow](https://github.com/zuizui0223/mois/actions/runs/39225659840) and [source-verified taxon / native-region Release](https://github.com/zuizui0223/mois/releases/tag/mois-real350399-kew-wcvp-native-tdwg3-source-v1).

## Source integrity

The input is the existing, photo-ID-deduplicated **350,399 actual original-image ROI attempts** (164,297 machine four-hue labels). No new iNaturalist pixels or source-only metadata were labelled. The underlying five separate original RGB/ROI datasets and original checksums were validated before obtaining independent botanical source data.

Kew original source archive: [WCVP original download](https://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip), listed as updated **4 June 2026**.
- Original source ZIP SHA256: `d32ea2b3a85e489b14e83bcc9eae7274532e1d113753f7be290d4b2dfde573fa`; verified unchanged in repeated live acquisitions, now pinned in Actions.
- Original archive members: `README_WCVP.xlsx`, `wcvp_names.csv`, `wcvp_distribution.csv`. CSVs are **pipe-separated**, as verified against the original source headers. Unlike an R package's converted schema, this source is read exactly as published and actual column/schema choices are printed and retained as source provenance.
- WCVP accepted names were resolved from *exact* binomials or unambiguous synonyms to accepted WCVP species identifiers. Ambiguous synonym matches and unverified accepted species ranks are **not** forced to accepted species.
- Distribution status is retained for each accepted species × **TDWG botanical country (level 3)**; rows marked introduced, doubtful, extinct or contradictory cannot silently be called native.
- Only taxon-name and region-status subsets are archived; raw Kew data and all original copyrighted photographic pixels are not redistributed.

## Actually measured source matching output

| Original source concept | Actual measured-source result |
|---|---:|
| Unique iNaturalist taxon IDs with any original photograph attempt | **74,821** |
| Named iNaturalist source taxa with ≥1 machine four-hue label | **37,107** |
| MOIS taxon IDs with a verified unique WCVP accepted species | **73,335** |
| Distinct verified WCVP accepted species names/IDs after synonym collapsing | **73,235** |
| WCVP accepted species with ≥1 matching TDWG3 distribution record | **73,232** |
| Unique exact name→accepted candidate before accepted-species identity check | 73,746 |
| Names referring to unverified accepted species ranks or IDs | 411 |
| Multiple accepted species names/IDs from one source exact name | 318 |
| No accepted species match | 339 |
| Not a strict binomial | 7 |

**TDWG level 3 source-record counts** (accepted species × botanical-country region *pairs*, not independent taxa, population or photographs):
- 531,048 explicitly native.
- 124,724 explicitly introduced.
- 697 distribution doubtful.
- 622 extinct.

The data includes many taxa in the historic original FCP/deep photo acquisition as well as species-first extension arms. Its **74,821** source taxon IDs are not the same eligibility denominator as the separate **42,111-taxon September 2026 / 93,209,351 photograph metadata frame**. Those are *different sampling frames* and must not be equated in the paper.

## Important inferential gate: native in a botanical country is NOT native on that island

WCVP Level 3 is a *botanical country*, often larger than one island. A species marked native in a TDWG3 unit can nevertheless have an introduced subpopulation on a sampled island. The current data is **not mapped to each photo coordinate**, because that needs a separate original-source geo-boundary join with uncertainty buffers. Even after such a geographic join, the result is only `native_at_TDWG3`, not `native_on_specific_island` or `endemic_to_island`.

**H3 native-island original hypothesis remains NOT TESTED.** The next necessary steps are:

1. Join original source photo coordinates with an independently checked official/attributed WGSRPD3 boundary layer, flag coastal/uncertain/unmatched points, preserve source uncertainty (unknown is not precise).
2. Append *exact* taxon name accepted-ID to TDWG3 native/introduced status (only when unambiguous), and measure coverage separately for classifiable white and nonwhite photographs.
3. Overlay the previously audited GSHHG island/mainland landmass class **without treating WCVP native-to-botanical-country as proof of island nativeness**. Require island-specific authoritative botanical native range for an island-native/endemism confirmation, or report H3 underdefined.
4. Restrict any comparison to shared taxa/regions and sampling-support cells and maintain original photography-source cohorts, classification success and source observation effort. Only then estimate exploratory descriptive white composition on the support retained.

## Research interpretation

This is a large genuinely independent botanical range **metadata linkage** to a measured photo-colour source. It improves the macroecological data basis but is not a test of pigment protection, selection, taxonomic/canonical species turnover or local bee diversity. Botanical native status does not validate that the image depicts a physically open petal rather than a bract or fruit. New independent photo outcomes and botanical-stage calibration are still required for future confirmatory colour hypotheses.

Reproducibility:
- `scripts/acquisition/mois_wcvp_TDWG3_native_source_20261019.py`
- `tests/test_mois_wcvp_TDWG3_native_source_20261019.py`
- `.github/workflows/mois-wcvp-source-native-tdwg3.yml`
- [Signed actual local data CSVs and manifest](https://github.com/zuizui0223/mois/releases/tag/mois-real350399-kew-wcvp-native-tdwg3-source-v1)
