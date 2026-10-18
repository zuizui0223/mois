# MOIS: island–mainland, elevation and Arctic-circle photographic colour composition

18 October 2026. Exploratory macroecological result, not FCP within-species selection.

## Source and successful execution

Frozen primary source: 313,535 genuinely image-processed photo IDs, 151,239 automatic four-hue labels, 71,427 white. No RGB reprocessing and no new photography in this analysis. Later differently recruited image cohorts were not silently pooled. The original colour/environment Parquet file hashes were rechecked.

Geography uses exactly island's GSHHG 2.3.7 high-resolution source archive, SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`. Mainland seed points identify continental sibling IDs 0,1,2,3,6; dateline-split pieces with the same sibling ID remain mainland. Greenland and large island landmasses remain islands, not assumed pollinator-poor. Exact coastal/boundary and offshore points are left unassigned; 6,511 attempts (2,935 classifiable) are not assigned to island/mainland. No country-code substitute or coastal snapping is used. Antarctica itself is outside the L1 layer and is not assigned by this first-pass analysis.

Successful workflow: https://github.com/zuizui0223/mois/actions/runs/39112342916
Artifact: `mois-real313535-island-mainland-macro-geography`, ID 12130981326.
Artifact ZIP SHA256: `0649622d9a8a577bca4ec7897c4d1065e2ef3234e46b3f3a9f2bc134b585ff33`.
Local checks independently verify the ZIP digest, all seven CSV/GZIP digests, and refit the adjusted coefficient from the exported 1,394-cell data to agreement within 1e-8 percentage points.

## 1. Island–mainland photographic composition

| Location class | Image attempts | Classifiable hue ROIs | White | Yellow/orange | Red/pink | Blue/purple |
|---|---:|---:|---:|---:|---:|---:|
| Islands | 41,828 | 18,864 | 50.933% | 29.591% | 10.623% | 8.853% |
| Mainlands | 265,196 | 129,440 | 46.710% | 31.218% | 10.701% | 11.371% |

The raw island–mainland white contrast is +4.223 percentage points. This is a composition of source-classifiable photographs, not the proportion of all real flowering individuals.

After equal weighting of photographed taxa within each 1-degree landmass cell, then equal weighting of cells:

| Minimum named taxa per cell | Island cells | Mainland cells | Island white | Mainland white |
|---|---:|---:|---:|---:|
| 1 | 1,906 | 6,125 | 46.729% | 46.589% |
| 5 | 610 | 3,102 | 50.324% | 46.892% |
| 10 | 352 | 2,133 | 49.186% | 47.289% |

This sensitivity is material: the island excess is not invariant to sparse-cell weighting and source inclusion. Do not foreground only the largest positive percentage.

## 2. Exploratory geographically conditioned comparison

Form cells from complete-case photo records first (heat, rain, solar and altitude all available); retain >=5 source taxa and >=10 classified photos per cell. Use only coarse geographic strata (10-degree latitude x30-degree longitude) with both island and mainland support. The matched source has 1,394 cells (177 island/1,217 mainland), 65,771 classified photo records, 35 coarse strata and 89 spatial clusters (10-degree cells). The response is the within-cell mean of per-taxon photographic white proportions, not a within-taxon colour change.

| Model | Island minus mainland white pp | Nominal spatial-cluster 95% interval | Nominal p |
|---|---:|---:|---:|
| Coarse geographic stratum effects | +3.740 | [+0.844,+6.635] | 0.01196 |
| Geography + warmest-month temperature + log(1+annual rainfall) + annual solar + altitude | +3.160 | [+0.276,+6.044] | 0.03212 |

Same exported complete-case cells, independently rerun locally:
- Minimum 10 taxa/cell: +3.081 pp [+0.118,+6.045].
- Exclude islands smaller than 5 km2: +3.095 [+0.258,+5.933].
- Exclude islands >=100,000 km2: +2.345 [-0.789,+5.478].
- Excluding the most represented island sibling 13 (17 cells) gives +2.964 [-0.040,+5.969]; excluding sibling 20 (14 cells) gives +2.593 [-0.167,+5.352]. These are post-hoc robustness checks, not independent confirmation.

The positive overall conditional association is promising but small-island generalization is not established. No multiplicity correction over the wider project's exploratory analyses has been applied. Shared coarse spatial strata are not matched biogeographic realms or identical species pools. Observer effects, native/introduced status and phylogeny remain unadjusted; no causal pollinator effect is measured.

## 3. Elevation and Arctic-circle source composition

| Elevation | Classifiable photo ROIs | White fraction | Blue/purple fraction |
|---|---:|---:|---:|
| 0–<500 m | 78,645 | 47.927% | 9.879% |
| 500–<1,500 m | 26,120 | 45.881% | 11.474% |
| 1,500–<2,500 m | 11,157 | 43.739% | 14.529% |
| >=2,500 m | 4,263 | 41.872% | 15.646% |

Only actually available original altitude fields are used. These are unadjusted global source mixtures, not same-mountain comparisons, a calibrated saturation gradient or formal alpine-biome membership. The high-elevation sample has a larger visible chromatic proportion, especially blue/purple, but region and lineage composition could contribute.

North of 66.562 N: 1,205 attempts, 632 classifiable ROIs, 257 photographed taxa; white 43.987%, yellow/orange 31.804%, red/pink 8.228%, blue/purple 15.981%. The raw Arctic-circle photos do not show a universal white excess. This is a small opportunistic source stratum and not a pollinator census or definition of the full Arctic biome.

## Scientific interpretation and next stage

The initial comparison supports studying colour composition as a geographical assemblage property: island photographs show a white excess in some weighting/conditioning schemes, whereas high-elevation and Arctic-circle photos have more blue/purple representation. Do not reduce this to a single 'pollinator-poor equals white' mechanism or a maintained 50:50 equilibrium. Species turnover is a valid MOIS response, but lineage independence needs independent analysis.

Next priority is matching actual biogeographic realms and available independent pollinator/plant reproductive data, plus robust high/low comparisons within mountain regions. Those enrichments can reuse the photo_ID-to-landmass cache and the stored species-cell table; no repeated image download, inference, or full coastline reconstruction is necessary. Seasonality remains SI. Visible white is not pigment absence; four RGB hue classes do not measure UV pigments, anthocyanin investment, calibrated colour saturation, reproductive fitness or biochemical costs.
