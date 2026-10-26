# MOIS 2026-10-26 — bias reduction by **sample-composition-aware photo acquisition**, not species novelty

## Corrected core goal

Earlier attempts to add photographs only from species with zero accepted ROI9 colour were inefficient for the actual MOIS objective. In particular, a 64-photo pilot of previously uncaptured colour species does not directly correct north/south, tropical, polar or seasonal source imbalance.

**Current decision: keep all earlier photo colours and unknowns, but choose NEW photos to reduce source coverage gaps.**

We first audited **1,043,930** original photo ID attempts (historical 4-hue: **459,235**), then authenticated **400,710** strict actual nine-component floral ROI colour vectors in the same photo universe by recovering 10,169 real old-wave numeric source rows. We do not equate 93,209,351 metadata source photos with 93m real colour measurements.

## 1. Sampling units, colour success and what 'balance' means

The baseline is nonprobabilistic and comes from FCP fixed100/taxon, K1, legacy history and source-photo stratified FCP5. ROI9 success depends on region, photo quality, flower size and source cohort. The 400,710 strict photographs are NOT 400,710 independent plants, native populations, or a source-weighted representative global flora.

Observed strict original nine-colour successful photo counts **after historical recovery**:

| Previously observed broad latitude band | Strict9 successful photos | New deficit-targeted original photos |
|---|---:|---:|
| N temperate | 283,193 | 254 |
| N tropical | 23,705 | 347 |
| S temperate | 64,847 | 327 |
| S tropical | 9,529 | **410** |
| N subpolar | 18,219 | 310 |
| N Arctic | 1,209 | 142 |
| S subpolar | 4 | 2 |
| S Antarctic | 0 | 0 |

The four previous missing-coordinate strict9 rows remain in the 400,710 source denominator and are not assigned a geographic region. Original-source supports, local four meteorological seasons, source photo-ID sampling units and old per-region-season ROI9 success rates were audited before sampling; quotas are not based on the colour outcome of newly chosen photos.

**Be explicit**: 1,792 additional original images cannot numerically equalize 283,193 N temperate vs 9,529 S tropical. The benefit of the first batch is to obtain more independent species/local seasons in poorly covered supports, measure their failure rates on a common authentic ROI method and make later allocation source-backed, not to claim bias eliminated.

## 2. Actually completed pre-colour source experiment

[Signed pre-colour source release](https://github.com/zuizui0223/mois/releases/tag/mois-original93m-tier07-geoseason-balanced-2048-botharms-source-20261026-v1); [source proof and execution run](https://github.com/zuizui0223/mois/actions/runs/40288371055).

- Source universe: exactly **93,209,351** unique original photo IDs, 42,111 source taxa, original observation UUID.
- No repeated 93m scan: read already exported **tier07 1,193,829 photo metadata IDs**, verify all 16 original release SHA-signed part archives and individual photo/observation randomization.
- Source candidate pool after old photo AND observation UUID disjointness and format support: **658,509 independent original observations**, 1,154,190 source-supported unattempted photo candidates. The two numbers are different sample grains.
- **256 reference**: independently and uniformly hashed previously unattempted *observation* IDs; one independently ranked source photo from each.
- **1,792 targeted**: allocate quotas across hemisphere×latitude band×local recorded season. Target score depends on true previous strict9 support, previous photo-ROI nonresponse rate, available original source distinct taxa. Inside each region/season choose one independent observation/photograph per taxon via salted hash and then draw source taxa by a separate independent rank.
- New observed photos: **2,048 different photo IDs AND 2,048 different observation UUIDs**, all disjoint from the previous 1,043,930 attempts. Photo statuses and probabilities were frozen before any new photo pixels.
- Low-support southern Antarctic stratum has no candidate; no pseudo-photos fabricated. Missing-geopoint/month original source observations are retained in the source monitoring/independent reference and not silently classified as northern plants.
- The two arms MUST NEVER be pooled with each other or legacy 104m ROI9 and called an unbiased globally representative flora. The reference is conditional on the tier07 source and previous-photo exclusion. Targeted admission is informed by previous ROI success so unconditional 93m photo inclusion probabilities are NOT identified; per-source candidate and conditional stage factors, strata quota and detailed photo identity are saved.

## 3. Actual original image colour acquisition

[True new-photo FCP ROI-v4 measurement workflow](https://github.com/zuizui0223/mois/actions/runs/40289009421) uses **two pooled warm model processes**, each four 256-photo work units. Frozen image-side YOLO11 and EfficientSAM produce the same photo flower-ROI 9-palette proportions as the earlier source, not a new unvalidated bounding-box-derived Lab classification.

The collector records (for every actually attempted image) original photo ID, source observation UUID and taxon, latitude and observation date, the 256-reference vs 1792-targeted arm, extracted full original nine palette fractions if valid, original source 4 hue, detailed technical failures/unknowns, RGB file SHA, and *stratum-specific numerator AND denominator* of ROI9 success. No raw licensed photo bytes published.

The named source image is NOT an independent anatomical petal mask, and missing/unknown ROI is NOT absence of pigments. The 251-photograph device-Lab BOX pilot is not substituted for ROI9 true original RGB measurements. Same original method is essential for comparable regional missingness.

## 4. How to decide if bias actually decreased

Judge the completed data at **three grains**, not raw new photo count alone:

1. **Availability**: previously attempted/original source independent observation ID and strict9 success per latitude zone and season. Measure the change in region/season gaps in *new unique observed taxa and geography-season support*, not simply total 104m photos.
2. **Measurement response**: ROI9 success rates and technical reason by reference vs target region/season. The reference/control informs *observed photographic nonresponse* under a predeclared source subset, not arbitrary true native-flora frequency. Avoid quota refills based on whether newly decoded photos turned white or had no ROI.
3. **Analysis support**: number of unique taxa/accepted genera and taxon×5° locality×season comparisons, source availability and season overlap, effective genus clusters and independently heldout region; always retain the 400,710 historical strata and 617,165 original 9ROI-missing records. Do not claim associations have become causal through oversampling.

Only after these checks choose another **single source-anchored scaled batch** for the remaining weak-supported cells, not repeated blind 64/512-photo jobs or repeated 93m metadata scans. Maintain regional observational inclusion probability information in a separate design sidecar. To balance true dry/wet, altitude and island gradients, **source candidate** photos still need predictor/landmass annotation before the next focused wave; latitude alone is not a surrogate for source soil/moisture or island status.

## 5. Repository state and limitations

- [Current bias-targeted source selector](https://github.com/zuizui0223/mois/blob/feat/full93m-medium-pixel-firstpass/scripts/acquisition/mois_balance1043930_tier07_2048_source.py).
- [True original frozen full-ROI9 collector](https://github.com/zuizui0223/mois/blob/feat/full93m-medium-pixel-firstpass/scripts/analysis/mois_collect_balanced_real2048_original_roi9.py).
- The new image-measurement workflow may be **running** until source-signed outputs are complete; measure this status from the actual Actions receipt rather than counting selected photo IDs as measured data.
- Sampling quotas are tuned to observed geographic and seasonal support, not native island-biota, actual individual flowering status or human-verified petal organ colour. No protected image content is redistributed.
