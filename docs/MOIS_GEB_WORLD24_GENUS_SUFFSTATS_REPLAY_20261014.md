# Complete source-matched 24-contrast genus-bootstrap replay — 2026-10-14

**Status:** PASS for the original 24-comparison *world / climate26 / classifiable all-photos / genus_geo / genus_equal* model family. Independent re-implementation from the two exact original measurement-table ZIP archives, using the source frozen resampling method and random seed. **No new biological photo/plant data, hypothesis family, search or coefficient selection.** This replication is the original model and original source—not independent scientific replication.

## What was independently computed

Starting with the original photo-colour measurement table ZIP SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0` and the original environment table ZIP SHA256 `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`, an independent source reader recovered 100,543 unique photo IDs, 45,836 classified climate26-complete photo records, and 25,028 source taxa. The ≥3-species nominal-genus gate retained **20,694 species and 1,934 genera**. Photo response means and site predictors were averaged **per named source species**, centred within its original nominal genus, standardized using all 45,836 climate26 source photo environmental means/SD, and adjusted for the three original broad geographic terms (absolute latitude, sin and cos longitude).

Rather than copying individual source photo or species rows, this independent reconstruction formed for each nominal genus its **genus-equal weighted crossproducts** `X'X/n_genus` and `X'Y/n_genus`, with `X` containing three geographic nuisance columns plus the six original biological climate axes (BIO5, BIO12, BIO14, BIO15, annual radiation and mean vapour pressure), and `Y` containing the four mutually exclusive photographic flower-colour fractions. The genus names were replaced by a stable **anonymous positional index**, not exported in the shared internal sufficient-statistic file.

For each of six focal environmental axes, solve the original four-covariate normal equations `[abs_lat, lon_sin, lon_cos, focal_axis]`, using the original **999 multinomial genus resamples** with frozen seed **2026101430**. Compute source unadjusted 2.5–97.5% bootstrap intervals; standardize every replicated error by its original cluster bootstrap SD; take the 95th percentile of the **maximum absolute standardized error across all 24 axis×hue contrasts**; reconstruct all 24 original scope-simultaneous maxT intervals.

| Independent original-source recheck | Result |
|---|---:|
| Named original source species admitted | 20,694 |
| Original genus clusters | 1,934 |
| Original genus bootstrap resamples | 999 |
| Original scope family size | 6 climate anchors × 4 photographic hues = 24 |
| Reconstructed maxT critical value | **3.0109108050453743** |
| Original photo-source CSV point coefficients matched | **24 /24** |
| Original percentile 95% intervals matched | **24 /24** |
| Original maxT 95% intervals and sign flags matched | **24 /24** |
| Maximum absolute point/interval difference versus original frozen 1,280-row CSV | **3.55 ×10⁻¹⁵ percentage points** |

Only the original **white–annual-precipitation** and **white–annual-vapour-pressure** individual contrasts have scope-specific maxT intervals excluding zero. The nominal confidence intervals of any individual nonwhite category do not replace their scope-corrected intervals. All original source results are read, not modified, after independent reconstruction.

**The original source reviewer-visible result record** (the 24 selected original ledger lines plus exact replay seed, maxT value and original-source comparison error) is committed as `results/mois_world24_genus_moment_source_replay_audit_20261014.json`. This receipt contains no photo ID, original individual plant position or source genus name and is checked in CI. The reproducible **replayer program** is `scripts/analysis/replay_mois_world24_from_genus_moments_20261014.py`.

## Restricted anonymous-group intermediate — NOT a public data publication

A separate **internal-only** package was generated with genus-indexed 1934×9×9 normal-equation moments, 1934×9×4 hue responses and minimum-three-species cluster counts (original-source `.npz` SHA256 `77dff4acc5756610e463bec3e5877c711f1695d8984be0b675e86f7d8afb4aae`), the original source 1,280-row coefficient table, a replay script, a provenance manifest and the complete original-source equality report. The intermediate **has not been committed to public GitHub** or deposited with a reviewer/stable DOI service. It has no original source taxon/photo identifiers, named genera or per-site coordinates. **It does still contain aggregated geographic-predictor covariance information**, so privacy and photo/data-licensing approval is required before external redistribution. Removal of labels is *not* a legal or statistical guarantee of anonymization.

The replayer expects the authorised internal moments file:

```bash
python scripts/analysis/replay_mois_world24_from_genus_moments_20261014.py \
  --suffstats path/to/mois_anonymised_genus_sufficient_stats_INTERNAL_20261014.npz \
  --manifest path/to/mois_anonymised_genus_suffstats_manifest_20261014.json \
  --ledger results/fcp_macro_direction_20261013/direction_ledger.csv \
  --output world24_recheck.json
```

The ZIP and the derived moments are **only handed to the study owner for internal verification**, not published in the repository and not an anonymous GEB reviewer link. Current CI verifies the committed 24-row receipt against the original 1,280-row scientific ledger, and compiles the replayer; it **cannot rerun these 999 worlds without the privately held moments**. The original source photos are still not included in either recovered ZIP, so photometric classifier remeasurement remains unavailable.

## Inference and submission status

This check upgrades numerical reproducibility from **8 main point coefficients** to **all 24 same-source world-climate26 family point coefficients and their original 999-bootstrap percentile/maxT uncertainty**. It does **not** independently recheck every original source family across four sectors/strict image quality and all other models, or replace quality, geographic, independent botanical-colour or source selection nonconfirmations. The 24-comparison maxT adjustment is **within one original model scope**, not adjustment for the entire post-outcome macroecological programme. Neither original photographed flower-colour source, source trait-classification missingness, environmental covariate confounding, morph-relative fitness nor selection causal inference is resolved by numerical equality.

**GEB submission still blocked:** author-approved individual-record/re-identification risk review; photography and licensed derived data permission; stable reviewer-anonymous data/code deposit with persistent identifier; final figures, title page, conflicts, permissions and coauthor approval. A private internal reproducibility file is **not a journal-compliant data deposition**.
