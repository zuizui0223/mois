# MOIS real 350,399 photographic pixels × WCVP regional native/introduced ranges

**Run:** 2026-10-19. **Status:** full 350,399-source-photo join **completed**, not a native-island adaptation hypothesis test.

- [Successful original source-SHA / original TDWG-3 boundary / WCVP join workflow](https://github.com/zuizui0223/mois/actions/runs/39229605151).
- [Immutable 350,399 original photo-ID-by-region table, summary CSV and receipt](https://github.com/zuizui0223/mois/releases/tag/mois-real350399-wcvp-tdwg3-photo-regional-status-v1).
- Script: `scripts/analysis/mois_real350399_wcvp_geographic_native_status_20261019.py`.
- [Previously completed official Kew WCVP accepted-name × native-range acquisition](https://github.com/zuizui0223/mois/releases/tag/mois-real350399-kew-wcvp-native-tdwg3-source-v1).
- Official TDWG [WGSRPD level-3 polygons](https://github.com/tdwg/wgsrpd/blob/52da7828aba9d461dd133c27b3bd7a4407161f54/geojson/level3.geojson), Git commit `52da7828aba9d461dd133c27b3bd7a4407161f54`, original Git blob SHA1 `91104e5159e31f88154833a51d3b0d1c9271083f`, original source SHA256 `c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840`. The original file bytes were not modified.

## Source identity and geographic coverage

Original image attempts: **350,399** with **164,297** genuine automatic four-hue ROIs. All five prior original image-photo releases and the accepted WCVP source tables were individually SHA-matched before any ecological classification. Technical unknown ROI colours retain their unknown state.

| Original geographic assignment state | Photographs |
|---|---:|
| Uniquely within one TDWG-3 botanical country polygon | 321,756 |
| Multiple TDWG-3 polygons at same original point | 1,131 |
| No TDWG-3 polygon interior or exactly on polygon boundary | 27,317 |
| Source latitude/longitude missing or invalid | 195 |
| **All true original image attempts** | **350,399** |

The original TDWG file has 369 botanical-country geometries; **19 invalid source rings** were repaired using GEOS `make_valid` with **polygon-only components**, while zero repairs remained invalid. No original source geometry file was replaced and no repaired dangling line was allowed to become a country interior. Every original coordinate with multiple, missing, offshore or coastline assignments remains unresolved rather than being snapped to a convenient mainland.

## Native and introduced `species-in-region` annotation, NOT a photographic individual status

| Original photograph's species × TDWG-3 region record | Original photos |
|---|---:|
| WCVP accepted species recorded **native in its assigned broad region** | **268,670** |
| WCVP accepted species recorded **introduced in its assigned broad region** | **30,960** |
| Source accepted species has no distribution row for this broad region | 17,460 |
| Original taxon exact accepted name unresolved | 4,618 |
| No uniquely assigned botanical country for original photo | 28,643 |
| WCVP broad region distribution doubtful | 25 |
| WCVP broad region distribution extinct | 23 |
| **All original photographs including unknown hues** | **350,399** |

The **268,670** in the first row is not a count of botanically validated **wild native flower individuals**. A photograph of a garden cultivar or an introduced population may occur inside the region where the species is generally native. Region-level distribution information also cannot distinguish an island population from its mainland relatives when both occur in a broad TDWG botanical country. Therefore:

- `native_at_TDWG3` is a **species-wide broad regional range annotation only**, never a confirmation that the pictured plant is native or wild at the source position.
- Neither introduced, uncertain, range-not-listed nor source-unclassifiable ROIs were silently recoded native.
- No source-only metadata entry was promoted to a photographic hue observation.
- **H3 (island-native white excess) remains untested**. To resolve H3, acquire island-level authoritative native/introduced/endemism distribution plus genuinely identified island landmass at the photograph coordinate and botanical organ validation. The old GSHHG-based exploratory island contrast cannot acquire genuine native-island status from broad WCVP alone.
- The above table gives *photo-attempt counts*, not verified-white photo counts. A separate colour analysis may use only the original four-hue classified subset, with the unsuccessful classifier/photo QC denominator retained.

## Intended next analysis and interpretation

For a descriptive sensitivity on this genuinely source-authenticated photographic sample, compare white-labelled original flowers among `species recorded native in region` versus `species recorded introduced in region`, with an unclassified-hue column and at least nominal genus × geographically shared cell/source-arm conditioning. Do not interpret this as a causal introduced-plant mechanism or evidence of biological selection without independent native locality, bloom validation and photographic effort data. Every native-island headline is withheld until island-specific native distributions can be obtained.

**Paper wording:** “We linked original photograph locations to broad WCVP species native and introduced TDWG level-3 ranges. This provides a region-level descriptor of the source photographed assemblage, not a determination of whether individual photographed plants or populations were wild and native, nor whether particular islands host native taxa.”
