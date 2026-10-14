# mois — global interspecific floral colour × moisture

**Scientific scope:** the association between photographed floral colour composition **among plant species** and moisture-related environment across the world. This repository is the destination for the cross-species global climate/soil analyses originally developed in [zuizui0223/fcp](https://github.com/zuizui0223/fcp). The separate **fcp** repository focuses on *within-species* photographed flower-colour variation, recurrent white–chromatic colour axes and geographic partitioning.

## The main ecological question

**Do photographed white versus coloured floral states covary with precipitation and atmospheric moisture among congeneric species, and how well does this relationship persist under photographic, spatial and alternative biological-measurement controls?**

The worldwide source contains **100,543 distinct original photo IDs** across **42,111 named taxa**, of which **45,953** received a classifiable visible flower hue and **54,590** remain unclassified. Main climate26 hue–environment contrasts involve **20,694 named species /1,934 qualifying genera**; soil36 uses a *different* complete-case population of **16,801 /1,665** for the comparable congeneric hue slopes. The global genus-equal white-colour coefficient for annual precipitation is **+2.255 percentage points per original source environmental SD** (within-scope 24-contrast maxT 95% [+0.319,+4.190]); atmospheric vapour pressure is **+3.952** ([+1.199,+6.704]). These are retrospective photographic *between-species* associations, not causal water-response, genotype-frequency or flower-pigmentation benefit/cost effects.

**Nonconfirmation matters:** strict-quality precipitation +1.229 pp (95% [−2.046,+4.504]), within-genus×coarse-cell +0.808 (95% [−3.556,+5.172]); the same-taxonomic-support paired spatial contrast has a simultaneous interval spanning zero. Independent NHM-described floral hue also fails to independently establish a positive geographic rainfall direction. Positive global and negative geographic/source checks belong in the same article.

## Source files transferred byte-for-byte from fcp

- [Global colour × 36-abiotic direction methods/report](docs/FCP_GLOBAL_COLOUR_ENVIRONMENT_DIRECTION_AND_CAUSAL_GATES_20261013.md); [source script](scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py); [1,280-colour×environment result ledger](results/fcp_macro_direction_20261013/direction_ledger.csv); [decision receipt](results/fcp_macro_direction_20261013/decision_receipt.json)
- [Same-species-matched macrogeographic support report](docs/FCP_MACRO_HUE_MATCHED_SPATIAL_SUPPORT_20261014.md); [source script](scripts/analysis/fcp_macro_spatial_support_decomposition_20261014.py); [576-row spatial support/difference ledger](results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv); [original source result](results/fcp_m1_matched_spatial_support_20261014/result.json)
- [Standalone main paper draft](docs/MOIS_GLOBAL_MOISTURE_MANUSCRIPT_20261014.md) and [full source/transfer manifest](docs/MOIS_ORIGINAL_FCP_PROVENANCE_AND_TRANSFER_20261014.md)
- [Supporting original-photo and alternative-source receipts](results/source_receipts/) retain missingness, matched within/between photo slope, ROI validity, botanical-colour and genus/spatial transport limits.

Original source identity: GitHub FCP [PR #193](https://github.com/zuizui0223/fcp/pull/193), [PR #194](https://github.com/zuizui0223/fcp/pull/194), source code commits `0adef37c75c0e5bc135c447a861cd77eed34dbb2` and `111eea149f5ac970036f82d7b646d79fcd378d68`. **Historical FCP files were copied, not deleted**; their blob hashes are checked in the migration manifest. Large original zipped image/abiotic archives live in the preserved source-run workflow artifacts/recovery routes, **not** in this repository: photo ZIP SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`, environment ZIP SHA256 `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`.

## Reproducibility and inference

The transferred analysis scripts and source-verifier scripts preserve the original model definitions and original source-archive checks. The original [source-verification run #38502660130](https://github.com/zuizui0223/fcp/actions/runs/38502660130) documents archived full source replay. The lightweight `mois` CI separately guards committed source tables and statements without pretending to redownload original photographs. Original source scripts may still require their historical raw ZIP/CSV environment and specific dependencies before a fresh end-to-end run.

**Not in scope:** prospective within-species white–chromatic test `W`, 1,499-species high-depth D and within-species 50-km local decomposition, which remain in [fcp](https://github.com/zuizui0223/fcp). Neither a photographic 46.95% white fraction nor the worldwide moisture coefficient proves a 50:50 evolutionary balance.
