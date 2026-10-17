# MOIS: Full 149,653-date/hue original photo source seasonal taxon frequency decomposition

**Source completed 17 October 2026.** Research owner: `zuizui0223/mois`. Uses the same real source-backed global among-named-taxon floral-photo colour estimand as the MOIS integrated Research Article; retains independent FCP within-species pigment W/D/50-km phenotypic/spatial polymorphism study outside this analysis.

## Full genuinely measured source, photo denominator

After all 256 independent frozen FCP worker photo-pixel partitions completed, MOIS [full 309,439 original-photo source](https://github.com/zuizui0223/mois/actions/runs/39028400690) holds **309,439** unique original photograph IDs, **149,662** originally observed ROI-v4 automatic four-hue labels and **149,653** photos with both one verified original/public observation day and an automatically classifiable four-hue label. Exactly 19 original photos have no verified date, only nine of which have a four-hue label. Photo colour, source day, taxon ID and original photo identity are SHA-frozen: full 309,439-source photo GZIP SHA256 `546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e`. No new hue labels, downloaded image pixels, or inferred biological flowering states are introduced in this arithmetic analysis.

The independent ROI classifiability screen is not equivalent to botanical verification of a photographed *open flower*: bracts, nonopen petals, fruits and source visibility may bias a four-hue label. Seasons are localized meteorological seasons in extra-tropical original observation coordinates/dates (Southern summer in December–February; Northern summer June–August), excluding tropical source photos from the summer-winter comparison.

## Exact photographic composition identity, not field flowering biology

Let \(S\) and \(W\) be the source's original *classifiable, dated* photographed taxa in summer and winter. Within the taxa with at least one classified original photograph in **both** seasons, let \(p_{it}\) be the automatic white photo-label share for source taxon \(i\) in source season \(t\), and \(w_{it}\) the taxon's relative **photographic frequency** across photographed common taxa. The symmetric within-photo-rate term \(\sum_i \tfrac{w_{iS}+w_{iW}}{2}(p_{iS}-p_{iW})\) and photo-mixture term \(\sum_i \tfrac{p_{iS}+p_{iW}}{2}(w_{iS}-w_{iW})\) exactly sum to the common-taxon photographic contrast. To reconstruct the whole-source summer minus winter white-label difference, also add an explicit summer-only and winter-only photographed-taxon **support correction**, without assigning unphotographed taxa fictitious opposite-season hue values.

| Hemisphere | All-date/hue original photos, summer white/classified | Winter white/classified | Unpaired summer−winter photo-white percentage points |
|---|---|---|---:|
| Northern extra-tropics | **22,116 / 44,979** | **2,214 / 4,779** | **+2.842** |
| Southern extra-tropics | **3,461 / 6,832** | **1,880 / 4,760** | **+11.163** |

The global photographed seasonal taxon overlap is small: **651 of 7,209 Northern** taxa and **394 of 4,368 Southern** taxa with summer/winter original source photos occur in *both* sampled photographic seasons (Jaccard **0.0903** and **0.0902**). This is insufficient to claim that taxa without winter photos are biologically absent or do not bloom.

### Exact decomposition of the whole photographed source

| Symmetric original photo census component (percentage points) | North | South |
|---|---:|---:|
| White-photo label-rate change **within shared photographed taxa** | +3.1095 | +3.5748 |
| Relative original photo-frequency shift **among shared taxa** | −0.7728 | **+8.2125** |
| Taxa with classified photographs **in summer only** | +2.8345 | +3.9276 |
| Taxa with classified photographs **in winter only** | −2.3292 | −4.5522 |
| **Exactly reconstructed observed summer−winter white photograph difference** | **+2.8419** | **+11.1629** |

The single biggest positive term in Southern unpaired photographic summer-white enrichment is the **relative frequency of taxa photographed in both seasons (+8.21 of +11.16 points)**, not the smaller shared-taxon within-photo-rate term **+3.57 points**. The result survives adding 3,601 genuine new four-hue ROIs, but **is a photographic frequency accounting identity, not causal biological mediation, seasonal selection, pigment plasticity, flower availability or reproductive trait re-evolution**. Source photo frequencies are not normalized flowering abundances.

The separate spatially stricter comparison in the same 0.25-degree approximate site × meteorological season-year, at least five distinct source observations and two named taxa photographed with classifiable colour in both seasons, has only **27 Northern/28 Southern** matched communities. The [independent 999-geocell-cluster bootstrap source analysis](MOIS_FULL256_FINISHED_FCP_309439_ACTUAL_SOURCE_20261017.md) finds Northern **+6.70 points [−1.53,+14.92]**, Southern **+4.49 [−6.42,+13.18]**: both intervals include zero, and neither establishes a common local seasonal whiteness change among true open flowers. Coarsening to 1°/2° bins creates more support but pools different ecological localities.

## Reproducibility

`scripts/analysis/mois_full149653_original_photo_season_taxon_composition.py` source-locks the signed original 309,439-photo snapshot (photo IDs, raw four hue, original dates and taxa) before computing the exact four-term photographic identity for each hemisphere, plus 999 **exploratory taxon-resampled percentile intervals** for shared taxa, which are **not** independently geographic, field-population or 24-family maxT intervals. Synthetic regression tests check exact reconstruction and prevent assigning unphotographed taxa an inferred missing colour. The new [source-authentic GitHub workflow](../.github/workflows/mois-149653-full256-original-season-photo-taxon-mixture.yml) downloads ONLY the previous successful original photo IDs, regenerates the complete 149,653-date/hue photographic decomposition and uploads an independently verified source receipt and journal-figure CSV.

The original long-term genus-equal 309,439-photo global WorldClim water effects are a **different** estimator: rainfall **+1.988 pp per source SD**, maxT95 **[+0.552,+3.424]**, and vapour **+3.161 [+1.042,+5.279]**, both with simultaneously interval-nonconfirmed sign under genus×coarse geographic cell conditioning. No year-of-observation monthly solar, temperature, VPD or drought has been fitted from the original photo dates. Biometeorology must not be inferred from historical monthly climatological normals.
