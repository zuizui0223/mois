# MOIS: environmental overlap, island composition and regional elevation

18 October 2026. New exploratory analyses of the existing, source-signed 313,535 original-photo dataset. No image redownload, RGB relabelling, or new photo acquisition was performed.

## Source and verification

Reuse the original geography output from workflow 39112342916, artifact 12130981326 (`mois-real313535-island-mainland-macro-geography`). Original ZIP SHA256: `0649622d9a8a577bca4ec7897c4d1065e2ef3234e46b3f3a9f2bc134b585ff33`.

All seven CSV/GZIP source-file SHA256 values were rechecked against `result.json`. The earlier +3.160 pp geography/climate-adjusted island coefficient was independently reproduced to within 1e-8 pp from the 1,394 exported complete-case cells. Six local unit tests passed. A separate statsmodels intercept-only, region-clustered fit reproduced the new paired estimate and interval to within 1e-10 pp.

The portable analysis package supplied with this result contains `macro_followup.py`, `test_macro.py`, the original source tables, every matched pair, all sensitivity outputs and a verification JSON. Executed analysis-script SHA256: `b66cbfdaeaf28fceef0cff037297501725cd14e53c7253aa9c0ebe05f6f5d029`. The numerical island summary is committed as `results/mois_macro_environment_matching_20261018.csv`.

## 1. Environment-overlap island–mainland comparison

The input consists of the earlier 1,394 photo-level-complete-case, species-equal 1-degree landmass cells: 177 island and 1,217 mainland cells. Each has >=5 photographed taxa and >=10 classifiable source photographs. Within each cell, a taxon's photographed white fraction receives equal weight. **The island and mainland cells do not have to contain the same taxa.** Assemblage turnover is part of MOIS's response, not FCP within-species polymorphism.

Within the SAME coarse 10-degree-latitude x30-degree-longitude stratum, perform optimal one-to-one matching WITHOUT replacement on standardized warmest-month temperature, log(1+annual rainfall), annual solar and elevation. Standardization is frozen on the original 1,394 cells. The primary caliper is Euclidean distance <=1 over these four standardized variables; 0.5 and 2 are exploratory tight/wide sensitivities. Matching does not read the hue outcomes. This is not true biogeographic-realm, observer or phylogenetic matching.

Primary result: **102 pairs**, **65 island landmasses**, **28 coarse regions**. Thus 75 of the initial 177 island cells do not obtain a unique mainland control within this caliper.

| Weighting on the matched sample | Island minus mainland white, pp | Nominal 95% interval |
|---|---:|---:|
| Equal matched pairs | +1.396 | [-1.553, +4.345], region-cluster t |
| Equal coarse regions | +0.092 | [-4.150, +4.335], t across region means |
| Equal island landmasses | +0.998 | Point estimate only |

Pair-equal white composition is 48.093% on islands versus 46.697% on mainlands. Tighter matching gives 55 pairs and -0.227 pp [-2.303,+1.848]; looser matching gives 141 pairs and +2.014 pp [-2.048,+6.077]. Removing islands >=100,000 km2 leaves 86 pairs and +0.704 pp [-2.411,+3.820]. These are alternative source supports and weighting rules, not independent confirmations or estimates of the same global population.

**Decision:** A universal island-white excess is not established by this common-support comparison. Keep the earlier all-cell +3.160 pp regression estimate as its own analysis, rather than silently replacing it. The change in result does not prove that the earlier difference was wholly caused by confounding: both the sampled population and its weighting changed.

## 2. Four hue classes on the same pairs

| Visible photo hue | Matched island minus mainland, pp | Nominal region-cluster 95% interval |
|---|---:|---:|
| White | +1.396 | [-1.553,+4.345] |
| Yellow/orange | +1.822 | [-0.934,+4.579] |
| Red/pink | -2.452 | [-3.872,-1.032] |
| Blue/purple | -0.766 | [-2.925,+1.394] |

Red/pink is a candidate compositional contrast, with the negative pair-equal direction retained in the caliper/area sensitivities. However, equal-region weighting gives -2.905 pp [-5.938,+0.129], which includes zero. All four changes sum to zero by composition. No multiplicity correction across the project's exploratory tests has been applied, and no universal named pollinator mechanism is inferred.

## 3. Nearby-region elevation composition diagnostic

Reuse the earlier species-equal 1-degree cell table, require >=5 taxa and >=10 photos, and compare **cell mean elevation 0–<500 m versus >=1,500 m**, within the same landmass and 5-degree x5-degree neighbourhood. This is an aggregate-cell diagnostic, not a reassignment of every original photograph to a high/low elevation band.

There are **44 neighbourhoods** with both categories: 191 low-elevation and 109 high-elevation cells, representing 16,144 and 8,425 classified photo records. Equal-neighbourhood white composition is **46.302% low versus 41.217% high**, a high-minus-low difference of **-5.085 pp**, nominal neighbourhood t95 **[-8.785,-1.386]**. The high-elevation side has lower white composition in **31/44 neighbourhoods**.

The four hue differences are white -5.085, yellow/orange +1.104, red/pink +2.318 and blue/purple +1.663 pp. The point difference is a larger chromatic share, **not greater calibrated saturation or a larger pigment concentration**.

Changing to 10-degree neighbourhoods gives 36 groups and -4.575 pp [-8.981,-0.169]. Leaving out each landmass in turn gives white point differences from -6.553 to -3.587 pp. The much higher >=2,500 m contrast has only eight 5-degree neighbourhoods (-9.528 pp [-18.750,-0.307]) or ten 10-degree neighbourhoods (-8.518 pp [-17.218,+0.182]); retain both rather than foregrounding one interval.

**Limitations:** These neighbourhoods are not identified mountain ranges and are not necessarily independent. The existing cell covariates have partial within-cell missingness, and mean altitude does not identify an alpine biome or the elevation of every flower. Temperature, rainfall, solar, observers and phylogeny have not all been controlled in this elevation diagnostic. It supports a geographically stratified compositional hypothesis, not universal alpine selection.

## Scientific scope

MOIS asks whether colour composition of different photographed floras is associated with geographic and environmental settings. Species turnover remains in the response. Island-white, Arctic-white and high-altitude-colour hypotheses must not be collapsed into a single claim that low pollinator availability inevitably produces white flowers. The present data do not measure pollinator scarcity, biochemical pigment cost, UV colour, botanical anthesis, flowering-individual abundance or reproductive fitness. Season remains SI.

The next useful data enrichment is an actual mountain/biogeographic-region crosswalk and independently measured pollinator or plant reproductive data. Those can reuse the source photo-to-landmass and species-cell caches rather than launching another image-measurement pipeline.
