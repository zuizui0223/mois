# MOIS Hui–Walter type no-gold-standard image error audit (2026-11-05)

## Completed numeric source audit

[GitHub run 41942675493](https://github.com/zuizui0223/mois/actions/runs/41942675493): SUCCESS. [Artifact 13377340690](https://github.com/zuizui0223/mois/actions/runs/41942675493/artifacts/13377340690), ZIP SHA256 8b08267cabe0f7444fa920c1f34d6634ccf9505eebd7d99e339205b3662cb0aa.

The original frozen 131,072 one-row/photo numeric Parquet was SHA verified at 2d648c832929f8ba59301906d5ea8d100bdfec6449c6e2828f6b9a910ddd0a02. It contains 131,072 unique photo and observation IDs, original strict ROI9 success 48,784, technical ROI9 unknown 82,288. The 500px YOLO flower candidate detector and flower BOX Lab descriptor ARE NOT a second independent botanical white-versus-nonwhite PETAL test. An independent second white classifier was not present in these original numeric files.

Consequently: **0 valid paired independent white tests in this signed cohort; Hui–Walter fit NOT RUN; test Se/Sp and corrected environment effects NOT ESTIMATED.** The source-only result does not change original pixels, any strict9 counts, or original unknown status. No old workers restarted, no new images downloaded, no model inference or MNAR tipping scenarios.

## Structural identifiability

Classic 2 tests × 2 different-prevalence groups: 2*(2²−1)=6 observable degrees of freedom; 2 prevalence + 4 shared sensitivity/specificity parameters = 6. This is just identified in parameter count, subject to distinct true prevalences, equal Se/Sp between populations, conditional independence given the true state, a resolved positive-class orientation and regular conditions. A just-identified table has no overidentifying degrees of freedom to independently diagnose its assumptions.

Allowing BOTH tests' Se/Sp to vary separately across those 2 groups requires 2 group prevalences + 8 group-specific accuracy parameters = 10, versus 6 observable degrees: **UNDERIDENTIFIED**. More groups alone do not resolve freely varying accuracy for two tests: with G groups, observations 3G and parameters 5G.

Three distinct binary tests × 2 groups, if all accuracies vary by group, have 2*(2³−1)=14 observable degrees and 2+3*2*2=14 accuracy/prevalence parameters. Equal dimensions do not guarantee stable inference, conditional independence, biological validity or unique positive-class orientation. Shared image difficulty, background, morphology, algorithm training and common preprocessing can create conditional dependence even when model architectures differ.

## Define the biological target before assigning binary labels

A nonflowering image, fruit, bud, indeterminate flower organ, or failed mask cannot be converted to botanical nonwhite. The desired latent state must be explicitly defined, e.g. an **open flower's identifiable petal/tepal is predominantly white**. Detecting any flower BOX is a distinct latent task, not the same white test. Keep an explicit abstention/unknown state; conditioning on both tests' successful calls can create an outcome-selected population.

Independent VLM or another image model must receive the SAME source-photo ID and independently produce botanical white/nonwhite/abstain; it must not reuse original YOLO/SAM masks or labels. OWLv2 selecting a different photo or just producing detection boxes does not meet this requirement. Original-photo SHA, model SHA, test definitions, source arm and photo/observation IDs are separate provenance requirements.

## Ecological design and inference boundaries

Original acquisition mixes source-reference 32,768 and adaptive P1–P4 targeted 98,304 images. Do not pool to derive a global wild native floral white rate. Environment strata must be fixed by independent abiotic measurements before looking at test outputs; different observed positive rates alone do not prove different latent prevalences. Repeated observations, taxonomic composition, geography, season, observer and environment-linked image failures affect transportability. Corrected white-vs-moisture coefficients would additionally require exact sampling-frame weighting and explicit handling of abstention/nonflowering, not just test confusion-matrix estimates.

## Research novelty is conditional, not established

Hui & Walter (1980), Biometrics 36:167–171, PMID7370371 / DOI 10.2307/2530508 is the original method. Toft et al. (2005), Preventive Veterinary Medicine DOI 10.1016/j.prevetmed.2005.01.006 discuss violations. No-ground-data image classification accuracy via latent classes has prior art (Foody 2012, IEEE TGRS DOI 10.1109/TGRS.2011.2174156). A new MEE-grade contribution would require an identifiable and independently tested formulation for ecological image classifiers with correlated errors, environmental variation, source sampling probabilities and missing or nonflower targets. Transfer to new taxa/regions and credible adversarial checks are essential.

## Execution and next work

- Added scripts/analysis/mois_latent_error_identifiability.py; 13 focused unit tests passed; original 131,072-photo source passed SHA/ID/schema verification.
- Future same-photo source audit requires distinct independent model hashes, biological white-petal labels (white/nonwhite/unknown), preassigned environment strata, and exact photo/observation IDs.
- No estimator or ecological corrected prevalence has been run, because the second genuine white classifier data do not yet exist.
- Next concrete input: independently generated same-photo white/nonwhite/abstain outputs for a source-reference set, ideally a third architecturally independent classifier too. Human botanical truth can be held out as a later falsification check rather than treated as automatically existing.

Scientific outcome: **the proposed classic Hui–Walter idea is technically meaningful but NOT directly identifiable for environment-specific error with the present two-method MOIS source**. This stage certifies the gate, not the error estimate.
