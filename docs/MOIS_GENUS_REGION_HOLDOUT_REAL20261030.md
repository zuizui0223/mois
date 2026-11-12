# MOIS genus-turnover and out-of-region prediction — 2026-10-30

## Provenance / execution

- Frozen original ROI9 nominal-taxon regional analysis: [run 41018155739](https://github.com/zuizui0223/mois/actions/runs/41018155739), 30,634 regional source taxon rows, original canonical 1,032,882 observations / 401,703 strict9 success, adaptive arm not pooled.
- Previously independently retrieved Kew WCVP (June 2026) accepted-name source: [run 39225659840](https://github.com/zuizui0223/mois/actions/runs/39225659840). Accepted genus is extracted from WCVP accepted binomials; names with ambiguous synonyms/unverified accepted IDs are excluded; multiple MOIS source names sharing one accepted species ID are collapsed by their previously measured photo support.
- **New data-only genus run [41042569706](https://github.com/zuizui0223/mois/actions/runs/41042569706) SUCCESS.** Artifact [12970538161](https://github.com/zuizui0223/mois/actions/runs/41042569706/artifacts/12970538161), ZIP SHA256 `4b0371d57dec08f8127e1b5ebe39837acbb618b118926e56c808f7b46c6bd1d8`.
- SHA-pinned input GZIP: regional taxon `5f9dee3b0ab41e42d0b1b1dd06650eda2b954b89cb118ddc8bcdd1a813da1479`, WCVP taxon resolution `dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2`, WCVP accepted ID/name distribution `9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3`.

## Identity and matching

The six regional comparisons contain **30,634** sector×nominal taxon measured rows. Of these **30,158 (98.45%)** receive an exact unambiguous WCVP accepted binomial and genus. **476 rows/5,701 measured photo records** are left unresolved, never guessed. Synonym collapsing removes 104 double-counted sector×accepted species records, yielding **30,054 sector×accepted species** rows, with 417,111 measured photo records retained. These sector counts are NOT 30,054 distinct global species and are NOT incremental new photo acquisitions.

## Genus-level regional colour composition

Each accepted species has equal weight within genus; each accepted genus has equal weight within region. This differs from both photo-equal and nominal-iNaturalist-taxon-equal summaries in the preceding analysis. Differences are in original ROI9 *white colour mass* percentage points.

| Focal minus reference | Focal/reference genera | Shared genera | All genera equally | Only common genera equally |
|---|---:|---:|---:|---:|
| Arctic minus northern subpolar | 172 / 645 | 164 | **−2.664 pp** | **−3.766 pp** |
| South temperate minus north temperate | 2,705 / 3,029 | 1,038 | **−1.644 pp** | **−0.398 pp** |
| Island minus mainland (only previously geolocated source-photo cache) | 804 / 1,255 | 754 | **＋0.934 pp** | **＋0.707 pp** |

Four compositional terms (common-genus trait, common-genus relative weight, focal-only genera, reference-only genera) reconcile exactly to each all-genus contrast. However large and opposing exclusive genus terms frequently cancel, so the decomposition is **not** a numerical causal percentage attributed to lineage turnover. Shared-genus comparisons change taxon and region support. Island geography is the frozen 313,535-photo cache, not new full-cohort coastal point classification.

## Directional geographic holdout (no test species in training)

Train a region's equal-accepted-species mean *per genus* using only **reference region** samples. For the other (focal) region, predict only accepted species **absent** from the training region and requiring **at least 3 distinct training species in their genus**. Uncovered genera are not predicted. Compare MSE with the same supported focal species receiving a single reference-region equal-species mean; 1 − genusMSE/reference-meanMSE, not a predictive R² for all flora.

| Held-out focal ← training reference | Novel focal accepted species | With genus prediction | Prediction as % all focal species | MSE improvement over reference mean on supported species |
|---|---:|---:|---:|---:|
| Arctic ← subpolar | 132 | 99 | **27.7%** | ＋32.8% |
| Subpolar ← Arctic | 1,544 | 306 | 17.3% | **−5.6%** |
| South temperate ← north temperate | 8,492 | 2,030 | 20.6% | ＋2.0% |
| North temperate ← south temperate | 12,452 | 3,555 | 25.8% | ＋9.2% |
| Island ← mainland | 162 | 72 | **4.9%** | ＋25.4% |
| Mainland ← island | 1,538 | 446 | 15.7% | ＋17.1% |

**Interpretation:** Some genus-level transfer exists **only where an out-of-region training genus has sufficient independently photographed species**. A 25.4% relative MSE improvement for 72 island-novel species does **not** support broad island-flora predictability; 95.1% of all focal accepted island species were outside the novel-species prediction target. South→north and north→south support is low, and north→south improvement is close to zero. The northern subpolar←Arctic direction is worse than the plain reference mean.

All six estimates are descriptive without independent confidence intervals. Test species have not entered training; regions are broad groups, not truly independent source-observer or biogeographic realm holdouts. Original source species are photographically classified, not verified wild native flowering individuals. Shared accepted species were **excluded** from the novel-species predictive accuracy calculations.

## Beyond white

The same species and genus source tables preserve blue/purple/bronze component fractions plus normalized original nine-palette entropy. In genus-equal contrast, bronze Arctic−subpolar = ＋1.983 pp (common genera ＋2.765 pp); south−north temperate = ＋1.894 pp (common genera ＋0.541 pp); island−mainland = −0.313 pp. **Bronze is not botanically proven brown or low conspicuousness**; the palette contains no independent green, corrected chroma, background contrast, floral display size, or perianth reduction. No scientific 'dull flower rate' has been inferred.

## Verification / next gated task

10 focused Python tests passed locally and in Actions; GitHub results and local outputs have byte-identical SHA256 for all five numeric CSV files. An independent pandas recalculation of six heldout MSE skills agreed within 1.84×10⁻¹⁶ and all 15 genus decompositions reconciled numerically. No photo downloads, re-inference, unknown-colour promotion, MNAR sensitivity or original manuscript edits.

The WCVP archived **derived** species tables include accepted name and ID but do not include family. Before claiming family-level composition, a separate SHA-authenticated original WCVP accepted-species → family lookup is required. No pollination guild, floral display or anatomical labels exist in this run.
