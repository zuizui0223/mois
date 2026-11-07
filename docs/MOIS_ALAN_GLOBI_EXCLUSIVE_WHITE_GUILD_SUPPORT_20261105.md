# MOIS × VIIRS × GloBI: moth-visitor documented plants, excluding overlap (2026-11-05)

## Status: actual same-photo source audit DONE; triple-difference NOT identified

Two scientific tasks must remain distinct:

1. Species in an independently documented floral interaction with a **moth family** versus species with an independently documented interaction with a **butterfly family**, using GloBI partner taxonomy.
2. The hypothesis of **nighttime moth dependence** and differential loss of independently verified WHITE flowers across high ALAN exposure. GloBI visit records alone do not establish nocturnality, pollination effectiveness, exclusive pollination, abundance, plant community absence or local extirpation.

Frozen GloBI live API retrieval [run 41960275560](https://github.com/zuizui0223/mois/actions/runs/41960275560) succeeded. Numerical original flower-interaction claims 13,913 after strict insect-family/plant-lineage/binomial/study-title checks, 1,748 MOIS photographed species connected to some moth-or-butterfly evidence, 95,919 guild-tagged photo-record associations (NOT distinct photographs because 585 MOIS species had both guild records). Source GloBI live API is NOT stable version-pinned and multiple queries hit page caps. Missing records never mean absent pollination.

[New executed source-locked 2020 analysis run 41973792195](https://github.com/zuizui0223/mois/actions/runs/41973792195): **SUCCESS**, 8 unit tests, immutable 6 input SHAs and 1 numerical output SHA passed; [artifact 13389919241](https://github.com/zuizui0223/mois/actions/runs/41973792195/artifacts/13389919241), ZIP SHA256 349ed241ca4bc0aa9e5bb71c6bdd0750ab45ff07194065ff3b28b6c4c0124ac0.

### 1. Source group exclusivity

Recorded plant guild evidence, across ALL unique GloBI flower-visit plant names before MOIS restriction:

- Moth record but no butterfly record: **563** botanical names.
- Butterfly record but no moth record: **1,611** botanical names.
- Both types of record: **831** botanical names. **Excluded** from the mutually exclusive comparison, not randomly reassigned.
- No GloBI record is **not** a verified daytime plant.

Within the 2020 matched MOIS original source (historical_main only, real original ROI9 success), original 2020 VIIRS four-month radiance and WorldPop2020 source includes **1,044 valid 0.05° cells**. After accepted-species mapping, GloBI exact-binomial evidence and ROI9 successful photographs:

- Strict ROI9 successful **photo-to-GloBI evidence group records**: **1,786** (including species with both guild evidence).
- Aggregated exact 2020 0.05° cell × accepted species × guild rows: **1,765** (species equal within cell; duplicate photographs averaged).
- Only **108 cells** contain at least one moth-record-only species and at least one butterfly-record-only species with strict9 ROI.
- Of those, **only 4 cells** have at least 2 independently photographed species in each evidence group, **zero** have >=3 in each.
- Of the 1,044 available 2020 VIIRS cells, only **one** has independently documented species-direct WHITE flower evidence in both exclusive visitor groups.
- **ZERO cells** support all four distinct categories: moth record×white, moth record×nonwhite, butterfly record×white, butterfly record×nonwhite from independently described flora.

Thus the proposed biological WHITE×TRUE NOCTURNAL GUILD×ALAN triple interaction cannot be estimated validly from this matched source. Do not turn photographic original ROI9 white mass or recorded moth visits into proof of those categories.

### 2. What can be estimated without inventing biology?

Auxiliary outcome: at the SAME 0.05° cell, **mean original continuous ROI9 white component of moth-record-only species minus that of butterfly-record-only species**. All 108 supported cells were used; species are equally weighted within cell and each cell is one record. This is not a species-white fraction, not an interaction of actual night-vision or moth-dependent reproduction, and not the advertised causal triple difference.

Mean within-cell moth-minus-butterfly ROI9 white-mass difference = **+13.81 percentage points** among these 108 selected cells. Not a true difference in biologically verified white-flowered species prevalence or a between-treatment ecological effect.

Exploratory associations of this within-cell contrast against 2020 log1p VIIRS raw upward light radiance:

| Model | loglight coefficient (pp/log unit) | Nominal 95% spatial-block CI | nominal p |
|---|---:|---:|---:|
| Night light only | −3.287 | [−12.603, +6.029] | 0.481 |
| Night light + 2020 population + BIO5 + precipitation + elevation + photo effort, **without geographic FE** | −15.364 | [−30.556, −0.171] | 0.048 |
| Same predictors **plus 5° regional FE** | **−17.904** | **[−51.950, +16.141]** | **0.295** |

48 5° geographic blocks. Regression diagnostics cannot provide robust biological interpretation: **104/108 cells have fewer than two species in at least one guild**. The apparently narrower interval for a model without region fixed effects does NOT replicate with region fixed effects, and is not an independent significance result or confirmatory evidence of ecological sorting. These analyses were exploratory, no global familywise multiplicity claim. Within-species data, study effort, original photos with flower organs, sampling abundance, guild effectiveness, spatial exposure precision and recorded interaction source completeness remain unvalidated.

### 3. What the original proposed design can and cannot cancel

Restricting attention to WHITE-photo records can help reduce a *shared* white-colour observer preference. It does not remove differences in species/lineage photos, flower architecture, photography, landscaping, study intensity or genus-specific colour detection. Conditioning on documented flower colour can introduce differential selection if evidence availability depends on guild or urban setting.

A **cross-sectional 3-way interaction** (ALAN×documented white×independently confirmed nocturnal guild) is not a temporal triple-difference causal design. Actual temporal difference-in-differences would additionally require a treatment change and appropriate untreated controls and a plausible parallel-trends design; cell fixed effects alone do not guarantee it.

### 4. Immediate research decision

Do NOT advertise an estimate that white moth-dependent plants disappear under ALAN from current MOIS evidence. The next information bottleneck is **independent species-level pollinator time-of-day and effective pollen transfer**, alongside reliably documented white flowers. To recover adequate sample size, screen GloBI references by effect and visit time and retain signed provenance (not family-wide photoperiod inference). Original confirmed records from GloBI version 0.9 may be preferable to the currently truncated live API. Keep white, nonwhite, both-guild evidence and unknown as explicit categories, and report region/category positivity before estimating any ecological interaction.

No new photographic acquisition, ROI re-inference, unknown-to-colour promotion, forbidden MNAR tipping scenarios, or edits to main. Source only and all previous tested MOIS ALAN white results remain intact.

### Numeric provenance

- Canonical MOIS original-photo observation frame: SHA256 f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c.
- WCVP exact nominal-to-accepted source ID: dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2.
- WCVP accepted names: 9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3.
- Island direct independently described species colours: 7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577 (restricted individual species data NOT redistributed).
- Bounded live GloBI source claim cache: 6116acab80f2277af7ec109715e5815bf6c4461767ae61580da8d842c9e230e4.
- 2020 VIIRS and WorldPop exact original numeric cells: 3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3.
- 2020 GloBI-group spatial exploratory model CSV: c613189ed7cdc5f82f6cacb375ca0c6b1d6409d46c38d9dc6cd6dc796207881b.
