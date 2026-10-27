# MOIS GEB — Sampling bias, ROI missingness and MNAR sensitivity (SI)

**Source freeze: 27 October 2026. Real original photos, no new image acquisition.**

## Primary analytical decision

The manuscript targets **photographed, classifiable floral-colour composition**, not the true distribution of wild native floral pigmentation. Retain previously fitted world-wide annual moisture associations from **309,439 genuine original photos**, but disclose that positive global signs do not necessarily persist with genus × geographical-cell restriction and that colour classification is nonrandom. The updated four-arm **1,181,146-photo** strict nine-component dataset is a measurement-quality audit, not a refit of global climate associations.

## Three biases require three different protections

| Selection stage | What is over-represented? | What has actually been done | What remains not identified |
|---|---|---|---|
| **Observation and source acquisition** | Observers, regions, garden/introduced plants, conspicuous taxa, repeatedly photographed genera | Source photo-ID/observation dedup; photo-equal, taxon-equal and genus-equal denominators; genus × geography matching; separate source-reference and adaptive acquisition arms | True worldwide native-flora/flowering-individual probabilities; historical inclusion propensities cannot be reconstructed reliably |
| **Automatic floral ROI and four-hue/strict9 eligibility** | Large visible petals, well lit photos, favourable backgrounds, detector/segmentation/flip-gate suitability | Full attempted-photo denominators and success rates; regional/seasonal audits; unknowns remain missing; 512 original technical failures independently screened | True colour and whether an assessable flower is present in any failed photograph |
| **Taxonomic, climate and geographical support** | Wet and dry sites populated or photographed by different taxa and genera; seasonal taxon/photo weights | 999 whole-genus bootstrap/maxT on 309,439 original photos; source-restricted genus/cell comparison; matched island and elevation source; four-term original summer–winter composition identity | Universal local/individual pigment response, fitness cost/benefit and climatic causality |

**Inverse-probability weighting cannot fully fix these data.** It requires a defensible target population, known or identifiable selection probability and positivity. The older heterogeneous FCP/iNaturalist data cannot retroactively supply a full global flora propensity, and selection dependent on the unseen flower colour (MNAR) is not resolved by covariate weighting. Accordingly, no IPW-corrected world wild-flora white fraction is claimed.

## Real 313,535-photo original hue nonresponse sensitivity (executed)

The previously published [nine real photo rainfall × solar tercile cells](https://github.com/zuizui0223/mois/releases/tag/mois-original313535-four-hues-water-solar-genus-20261018-v1) were generated using all original attempted photo locations **before** selecting classifiable photos. Exactly **313,204 original photo attempts** have both rainfall and solar terciles; 151,107 have classifiable four-hue ROIs and **162,097 have no four-hue outcome**. The other 331 source photos lack these paired climate-bin data. This is a different eligibility definition from the 1.18m strict9 count.

For each wet and dry rainfall group at fixed solar level, write \(q=N_{\rm measured}/N_{\rm attempted}\), \(p=N_{\rm white}/N_{\rm measured}\) and \(u\) for an **assumed, unknown** white-label fraction among unclassifiable photo attempts. The counterfactual photo-label fraction is \(p^*=qp+(1-q)u\). This is a sensitivity relationship, not a measured hue or an imputed flower phenotype. It assumes, only for its counterfactual interpretation, that a missing photo could meaningfully receive a binary label; some missing photographs may have no flowers.

Under the scenario where the **dry-side unknown white fraction equals its observed dry white fraction**, the following wet-side hypothetical deviations erase the observed wet-minus-dry photographic contrast:

| Source solar tertile | Original dry/wet photo attempts | Complete-case wet−dry white-photo contrast | Wet unclassified white probability reduction needed to erase contrast |
|---|---:|---:|---:|
| Low | 32,271 / 30,348 | **+3.93 percentage points** | **7.68 percentage points** |
| Middle | 22,387 / 45,556 | **+5.33** | **9.74** |
| High | 50,084 / 28,496 | **+4.37** | **7.44** |

Weighting the three comparisons with the **number of all attempted wet/dry source photos in each solar stratum, before outcome classification**, yields a complete-case wet−dry white-photo difference of **+4.5504 percentage points**. Its direction would reach zero if the assumed wet-side unclassified white probability were **8.2506 percentage points below** the classified wet-side white fraction, keeping dry missingness MAR within its rain/solar cell. In the symmetrical scenario of a **5-point decrease** in wet missing white probability and a **5-point increase** in dry missing white probability, the contrast becomes **−0.6662 percentage points**. With wet missing −10 points and dry missing MAR, it becomes **−0.9648 points**.

**Actual numerical publication:** [Source-SHA audited MNAR tipping sensitivity release](https://github.com/zuizui0223/mois/releases/tag/mois-real313535-white-solar-wetdry-MNAR-tipping-source-20261027-v1), [successful test and execution run](https://github.com/zuizui0223/mois/actions/runs/40470760620). The source nine-cell CSV passed SHA256 verification; scenarios at hypothetical missing-value departures of 0, 5, 10, 15 and 20 percentage points and conservative full-unknown bounds are archived. **Zero original photos were downloaded, zero ROI labels added, no unseen biological pigment imputed.**

### What this does and does not say about the paper

The observed wet-minus-dry photographic white-label composition on a rainfall×solar source **can be erased by moderate colour-dependent nonresponse**. Thus its positive sign is not universally robust to missing-not-at-random scenarios. These descriptive comparisons do **not** control genus, observer, exact site and time simultaneously. Crucially, this is **NOT** an alternative estimate, a new confidence interval, or a mathematical correction of the different **309,439-photo genus-weighted +1.9884 pp/source SD annual rainfall slope**. The latter remains a supported *conditional photographic association* with a within-model simultaneous confidence interval excluding zero globally and including zero after genus/coarse geography restriction. It cannot be promoted to a universally causal true-flora pigmentation effect.

## Reviewer-facing reporting commitments

1. Show source image attempts and technically unknowns beside colour-success denominators; distinguish original four-hue and stricter nine-component measurement.
2. Report photo/taxon/genus weights and the limited matched genus-geography support. Do not pool outcome-targeted source arms with a probability-reference arm into an unweighted estimate of all flora.
3. Report the MNAR tipping above as explicit counterfactual scenarios, not recovered petal colour. Unknown is never automatically nonwhite, white or nonflowering.
4. Present globally positive environmental gradients alongside geographically conditioned intervals crossing zero and sparse same-species/season replication.
5. Avoid claiming a resolved biological mechanism or 50:50 evolutionarily maintained flower-colour equilibrium. The cost/fitness/pollinator pathways and actual wild-native status were not tested.

**Main paper:** report one short MNAR result and retain the original global/local moisture contrast. **SI:** this complete table, photo-classifiability by region/season/source cohort and all signed sensitivity outputs. This resolves the question of how sensitive observed associations are to missing-label assumptions **without pretending all photographic bias has been removed**.
