# Original-source independent matched-spatial 576-row reconstruction — 14 October 2026

**Status: PASS, source-data statistical verification in this working environment.** This recheck builds on the full original-source climate26 24-contrast bootstrap verification. It targets a **different inferential question**: whether changes in photo-white–moisture slopes after altering geographic conditioning on the *same named species* exceed the original paired uncertainty.

## Design preserved verbatim

Two checksum-frozen original source-table archives (photo-colour/site SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`; all36 original environment SHA256 `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`) were opened locally and joined by photo ID, observation ID and taxon ID. The 26 nonsoil climate fields were complete for **45,836 classified photographs**, and the source species means were created separately for all-world and East 90–180° photo opportunities. A source-derived 18×9 equal-area grid assigned each species-mean photo site to a coarse cell.

Within each scope, Stage A fitted species from eligible genera with broad geographic coordinates; B retained only eligible same-genus×coarse-cell species while using the same broad coordinate specification; C reused **exactly the same species/photo outcomes/environment values as B**, but changed to within-genus×coarse-cell conditional response and predictors. Neither the original observation source nor the sampled plant populations were experimentally altered. Four original weighting/region combinations were checked: world/East × genus-equal/species-equal. Every condition comprises three stage estimates and three paired stage differences, six source ecological axes and four mutually exclusive visible photographic flower hues.

A **second independently coded source reader and weighted least-squares/multinomial bootstrap implementation**, separate from the original FCP #194 script, reprocessed the exact original table ZIPs using the original 999 paired genus resamples and original region-dependent seeds. Every source comparison used identical genus bootstrap worlds for A/B/C. Original paired difference maxT intervals are adjusted for **24 source hue×environment contrasts within each region/weight/paired-difference family** only.

## Reproduced comparisons, without a new alternative choice

| Source region/weight | B & C identical species | Broad B white/BIO12 pp/SD | Local C white/BIO12 pp/SD | B−C pp/SD | Original 24-maxT 95% interval for B−C |
|---|---:|---:|---:|---:|---:|
| World; genus-equal | 9,370 species /845 genera /14,268 photo opportunities | +1.32586 | +0.80798 | **+0.51788** | **[−1.25351,+2.28927]** |
| East 90–180°; genus-equal | 2,886 /280 /4,149 | +1.10200 | −2.07568 | **+3.17768** | **[−0.25804,+6.61340]** |
| World; species-equal | 9,370 /845 /14,268 | +1.95017 | +0.70883 | **+1.24135** | **[−1.24305,+3.72574]** |
| East 90–180°; species-equal | 2,886 /280 /4,149 | +4.59150 | +1.24768 | **+3.34382** | **[−0.51029,+7.19794]** |

All four B−C 24-maxT 95% intervals include zero. The East *nominal* (uncorrected) interval can exclude zero even though the source predeclared 24-maxT interval includes zero. A large visual point difference is not a familywise-confirmed local adjustment contrast.

| Technical source reproduction check | Independently verified |
|---|---:|
| Original macro stage/difference coefficient rows | **576/576** |
| Original point estimates | **576 values** |
| Original ordinary bootstrap confidence interval endpoints | **1,152 values** |
| Original paired-difference maxT confidence interval endpoints | **576 values** |
| Original identical-genus source resamples per scope | **999** |
| Maximum absolute difference from the original 576-row frozen CSV | **1.60×10⁻¹⁴ percentage points** |

The counts total **2,304 independently compared numeric fields**. Reconciliation source: `results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv` and the historic PR #194 decision receipt. The public **aggregate-only** recheck is `results/mois_original_matched_spatial_full576_replay_20261014.json`. The newly independent optional source-table replay is `scripts/analysis/replay_mois_original_matched_spatial_576_20261014.py`; it accepts original source ZIP paths and the frozen 576-row source CSV, preserving model/data identity. Source ZIP files are **not** committed. Public CI compares the aggregate receipt with all relevant saved original source anchors but cannot rerun the full program without private raw source table access.

## Ecological interpretation and submission safety

- The broad world association remains an observed positive **source-conditional between-species photographic colour–water relationship**, not a universal physiological response or genotype adaptation.
- Conditioning on same-genus×same-large-cell changes the statistical estimand. Holding named species constant eliminates one simple source-composition explanation for changes in fitted slope, but a **24-maxT interval spanning zero** means the direction of the *difference* is not established at the source-family threshold; this does not prove the two target effects are equal.
- The East sign change after tighter spatial conditioning is an **exploratory point estimate**; it is **not** a regionally confirmed opposite water response.
- No matched biological individual, seed/pollen relative fitness, pigment chemistry, genetic flower-colour polymorphism or plant-level water treatment is present.
- The original source ZIPs contain third-party image/observation IDs, photo licence fields and site coordinates (not the original camera image pixels). Their use for independent internal checking does not authorise publication of raw records. **An external licensed/anonymized reviewer-accessible archive remains unresolved.**

**No retroactive preregistration claim:** the paired spatial comparison was constructed after the original broad source was inspected. Its historical 24-maxT pertains to specified within-scope contrasts, **not** the entire historical family of selected regions/models.
