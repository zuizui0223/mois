# Actual completed MOIS historical 100,543 original-photo date recovery — official 2026-09-27 iNaturalist Open Dataset

**Results verified 17 October 2026.** [Original source-photo/ROI colours](https://github.com/zuizui0223/mois/actions/runs/39002957961) remain exactly **267,425 real image photo IDs and 129,795 automated four-hue labels**. The previously exported **100,543 historic broad-taxonomy source photo IDs** were all original-date blank in the combined colour CSV. The [entire official metadata join succeeded](https://github.com/zuizui0223/mois/actions/runs/39009539731), not merely a simulation or 9,600-photo extrapolation.

## Complete provider source evidence

The official **2026-09-27** monthly metadata objects were fetched *once each*, directly streamed and checked through complete compressed byte count, end-of-gzip CRC and independently accumulated SHA256. This source scan processed **501,256,126 photographic metadata rows** and **282,933,221 observation metadata rows**. Their full original compressed SHA256 are `033b41cb9e0f28fa2015a66f9bc7e4e95d2ebb1bd96b174e310195eb241888be` (20,668,099,997-byte `photos.csv.gz`) and `e6c1706f2d585280b417d527465c2cbe41d85a7135b8bea3dfac7c4093e3dbf0` (13,420,521,411-byte `observations.csv.gz`). ETags and historical original photo table digest are source-frozen in [the numeric decision receipt](../results/mois_complete_historic100543_official_s3_original_date_recovery_20261017.json). This job did **not** redownload photographic pixels or recompute colour classification.

The exact provider join was **original `photo_id` → public `observation_uuid` → public snapshot `observed_on`**. Historical source numeric `observation_id` was preserved, but that ID does not appear in the iNaturalist Open Dataset observation table. The separate [original200 live public API source validation](https://github.com/zuizui0223/mois/actions/runs/39007591041) found matching public dates in **197 of 200** historic photos appearing in BOTH sources; there were **zero date disagreements among those 197** in that specific subset. This is cross-validation of the source chain, not a retrospective numeric observation-ID confirmation for all 100,543.

## Direct measured source outcome

| Real-source metric | Count |
|---|---:|
| Original distinct historical photo IDs, retained without loss | 100,543 |
| Historic original source photos present in public photo UUID map | 99,936 |
| Ambiguous/current moved UUID associations excluded from date evidence | 204 |
| Valid photo → current UUID → snapshot original `observed_on` | **99,732 (99.19%)** |
| No valid official-snapshot date, not imputed | **811** |
| Historic photographs carrying *both* public date and original automatic four-hue label | **45,591 of original 45,953** |
| Historic original white-classified photographs carrying snapshot date | 21,408 |
| Current snapshot taxon ID disagrees with frozen original source taxon ID | 77 |
| Entire genuine 267,425-photo dataset: date plus automatic hue **before → after** | **81,867 → 127,458** |

The total +45,591 newly date-resolved automatic hues are included in the 127,458; the original 81,867 source-date hues are not counted again. The historical exact 100,543 measured photo IDs and original 45,953 four-hue state labels are never changed. Original 811 missing dates are unknown, not absence of flowering. Full 100,543-photo detail and SHA/ETag receipt are held as the successful run artifact `mois-full-historic100543-real-photo-public20260927-S3-uuid-date-source-rejoin` (artifact ID 12092882545), not a public Git-tracked raw individual-photo table.

## First seasonal photo-composition evidence on the fully enriched historical source

The original historical broad/taxon-shallow photographic source *by itself* (45,591 classified and date-verified photo ROIs) gives:

| Original historical photographed hue | Local summer white photo/4-hue labels | Local winter white photo/4-hue labels | Summer–winter photo-white difference |
|---|---:|---:|---:|
| Northern extra-tropics | 5,303 / 10,902 = **48.64%** | 369 / 812 = **45.44%** | **+3.20 pp** |
| Southern extra-tropics | 1,369 / 2,796 = **48.96%** | 824 / 1,984 = **41.53%** | **+7.43 pp** |

When adding the separately measured high-depth FCP-derived source dates, the whole original 267,425-photo union is **48.84% summer vs 45.59% winter North** and **51.42% summer vs 38.93% winter South**. These are authentic **source-composition descriptive photo labels**, not verified blooming flowers, taxon-independent causal seasonal pigmentation changes, or independent externally sampled communities. Source taxonomy, observer effort, year/geography, current mutable provider metadata and whether a pictured organ is a true open petal may all affect ratios.

## Local source support after full dates: taxon and cell still confounded

Requiring **the same named original taxon × same 0.25° geographic cell × same meteorological season-year**, with at least **two different original observation IDs per season**, recovers only **two Northern source cell×taxon×years (eight genuine four-hue photographs) and zero Southern**. This is one additional Northern group relative to the earlier 9,600-photo date-enriched source; it remains wholly inadequate for any common *within-taxon* pigment–solar or pigment–rainfall season coefficient. Enlarging photo source taxa/date breadth is not equivalent to longitudinal observation of botanical morphs in one source locality.

**Necessary next distinction:** A *local photographic community assemblage* allows different source taxa to be photographed in winter and summer, but even that estimand has source effort and spatial-grid issues; it must not be called within-individual plasticity. Independent taxon-specific open-flower status, actual observation-year and month weather, and repeated source photographer/field effort remain unmeasured.

## Independence and provenance caveats

- The original RGB-v4 ROI automatic white/chromatic labels are **not** separately verified flower petals or reflectance/anthocyanin traits; the source's current public `Flowers` annotation may be absent and is not photo-specific.
- Monthly public `observed_on` is snapshot metadata, not necessarily camera-capture time and can change after 2026-09-27. Date/taxon disagreements are not rewritten silently; independently observed live API metadata may be newer.
- Original photo licensing/privacy and deleted/reassigned public photo IDs affect source coverage; original 811 source gaps remain explicitly missing in this versioned snapshot.
- FCP's within-species photographic `W/D/50-km` diversity/selection estimands are not imported; the historical MOIS 100,543-photo global water maxT slopes and already recomputed 267,425-photo annual precipitation/vapour estimates are **not** altered by metadata-only date restoration.
- An additional 9,600-ID current public API source was independently verified and may recover a small, provenance-distinct subset absent from this 2026-09-27 snapshot; it must be compared by exact literal photo ID, not added as though a separate set of images or independent ecological replication.

**Verified scientific inference:** It is now possible to study globally sampled photographic **seasonal flower-colour composition** across much more of the original source; current data still do not identify evolutionary colour balance, within-plant plasticity or a common solar/water causal mechanism.
