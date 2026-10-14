# Additional global species-colour × moisture evidence (mois)

**Original-source sensitivity results, 14 October 2026.** These are **independent analytical estimands from partly overlapping FCP photographic source records**, not independent botanical/photographic sample replicates. Source PRs and original SHA records retain priority over subsequent summaries. The primary worldwide interspecific colour–rainfall macro result is in [PR #193](https://github.com/zuizui0223/fcp/pull/193); paired spatial support in [PR #194](https://github.com/zuizui0223/fcp/pull/194).

| Test / scale | Original numeric observation | Required inference limit |
|---|---|---|
| Global 26-abiotic congeneric 4-colour/BIO12 | 20,694 species /1,934 genera; photo-white +2.255 pp/source SD, 24 contrast maxT95 [+0.319,+4.190] | Source-specific observational *interspecific* correlation |
| Strict photograph and coarse geography | +1.229 [−2.046,+4.504] and +0.808 [−3.556,+5.172] respectively | No independently established universal signed direction |
| Identical-species broad vs genus×cell | East 2,886 species: difference +3.178 pp/SD; 24-maxT95 [−0.258,+6.613] | Same membership controls taxon composition, not causal geography; familywise inconclusive |
| All36 full environment, unseen genus absolute composition | 20,786 species/4,853 genera; genus-equal Brier gain +0.004518 [0.002437,0.006602] | Mostly average taxonomic colour composition; not within-genus and within-species prediction |
| Within unseen genus species differences (all36) | 16,801 species/1,665 genera, incremental gain −0.000142 [−0.000748,+0.000433] | No supported universal among-congener difference transfer |
| Within-species photo deviations (all36) | 5,939 species/2,196 genera, gain −0.000670 [−0.001038,−0.000309] | Environmental **prediction model fails** for this target, not absence of physiological response |
| Add soil10 on identical 20,786 species | Incremental Brier +0.0000757 [−0.000516,+0.000653] | Added soil contribution unconfirmed in this fixed design |
| Unseen genus + another 90° sector + 500 km | All36 genus-equal Brier increment +0.000174 [−0.002629,+0.002846] | Transfer not robust across geographically isolated floras |
| Fixed nonlinear learner under distant heldout test | Nonlinear GEO→GEO+36 increment +0.004882 [+0.002320,+0.007496]; *strict photo* +0.000517 [−0.004122,+0.005363] | Nonlinear geographic baseline is worse than ridge; direct full-model comparison inconclusive, strict photo not supported |
| Original 26-variable four-hue Brier decomposition | White component +0.001517, 51.5% of the *same total* +0.002943 Brier improvement | An algebraic predictive-loss partition, not 51.5% pigmentation, variation or adaptation |
| Exactly the same 21,448 source photo IDs, 6,143 species | Annual rainfall–photo-white slope among congeners +2.722 pp/source SD, within named species −0.511, paired diff +3.233 [1.175,5.400] | **Same source photographs** directly show levels differ; posthoc association, not a causal scale-specific selection coefficient |
| Photographed ROI white fraction source sensitivity | 14,667 photos/3,351 species, between-congener +2.229, within-species −0.180 pp per *its own* SD; strict ROI quality between +1.116 [−1.169,+3.401] | Continuous photographed white pixels not spectrophotometry; strict quality unconfirmed |
| Original source unclassified photo colours | 54,590 of 100,543 photo colours unresolved | MNAR unknown-hue can erase wet–dry photo-white contrast; cannot call classifiable photo sample globally representative |
| Independent botanical NHM text hue | 670 eastern species/58 genera, genus-equal white–blue rainfall contrast +3.56 [−15.34,+22.60] | No independent phenotype direction; geographic sites still originally from FCP source |
| Global original classified white percentage | 21,574 /45,953 =46.95% | Classified-photo rate, **not** 50:50 genetic/floral morph fitness equilibrium |

### Upstream source and location of copied exact receipts

- [PR #168](https://github.com/zuizui0223/fcp/pull/168): `results/source_receipts/unseen_genus_colour_environment.json` — original root receipt and full original Actions run provenance.
- [PR #169](https://github.com/zuizui0223/fcp/pull/169): `results/source_receipts/unseen_genus_other_geography_500km.json` — stable folds, exact genus/sector exclusions and original independent replay.
- [PR #170](https://github.com/zuizui0223/fcp/pull/170): `results/source_receipts/nonlinear_geographic_transport.json`, `docs/supporting/NONLINEAR_GEOGRAPHIC_TRANSPORT_LIMIT.md` — nonlinear model, all-fold source hashes and no-rescue flags.
- [PR #167](https://github.com/zuizui0223/fcp/pull/167): `docs/supporting/FOUR_HUE_ENVIRONMENTAL_PREDICTION_DECOMPOSITION.md` — colour contribution to the same original four-colour predictive loss.
- FCP [PR #158](https://github.com/zuizui0223/fcp/pull/158) and [PR #159](https://github.com/zuizui0223/fcp/pull/159): `results/source_receipts/` matched categorical and continuous photo-colour records, missingness and photographed-white fraction.

These upstream FCP original branches remain in GitHub. This is a **migration**, not new estimation or a new substantive confirmatory result. Full original source executable code for PR #168–#170 remains referenced from those original commits; their signed source-verification receipts are copied intact.
