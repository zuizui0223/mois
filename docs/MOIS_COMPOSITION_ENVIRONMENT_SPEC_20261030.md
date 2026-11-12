# MOIS: taxon composition, nine ROI components and environmental associations

User-approved direction: separate regional taxon composition from shared-taxon differences, extend beyond binary white using the measured nine ROI components, and compare elevation associations before/after accounting for temperature, rainfall and solar radiation. This is an exploratory extension of already viewed data, not a newly sealed confirmatory test.

## Frozen populations and source protection

Reuse the canonical outcome-blind one-photo-per-observation frame (1,032,882 observation UUIDs; 401,703 strict9 successes) from run 41007444495, the original georestored-v2 four-arm palette from run 40440441203 (1,181,146 photos; 451,643 strict9 successes), and only exact photo-ID matches to the original 313,535-photo landmass cache from run 39112342916. Three SHA256 values are pinned in the script. Photo/observation/taxon identity and strict9 outcomes must agree before joining. No spatial extrapolation of the old island cache.

Primary composition and environment analyses use only `historical_K1_FCP_source_heterogeneous`; subsequent adaptive refill arms stay separate in source/region coverage tables. This remains a heterogeneous photographic source, not a probability sample of native flora. Source-subarm differences and unmeasured photographs remain explicit limitations.

## Regional composition, not removal of a nuisance

Contrasts: Arctic >=66.5622N vs 55--66.5622N; southern vs northern temperate bands with absolute latitude 23.5--55 degrees; islands vs mainlands on the exact cached-ID subset. Equal nominal `inat_taxon_id` weights replace photograph-frequency weights. Each equal-taxon regional difference is reconstructed exactly from shared-taxon outcome differences, shared-taxon relative-weight differences, focal-only taxa and reference-only taxa. Show the number/proportion of shared taxa and separately report the common-taxon-only estimand. The common subset is not the same target population as the total assemblage. No algebraic component is a causal, evolutionary or physiological attribution. Genus/family/phylogeny and functional-guild attribution are not yet fitted.

## Beyond white, without inventing dull flowers

Analyse nine recorded fractions, normalized nine-palette entropy and dominance. These are automatically segmented ROI composition descriptors. White is not low conspicuousness; bronze is not botanically verified brown; palette entropy is not calibrated chroma, flower size or floral display. Green as an independent colour, calibrated Lab, flower-background contrast, display area and perianth reduction are not present. Record their absence explicitly; do not manufacture an inconspicuousness index. Unknowns remain unknown.

## Altitude and measured environmental pathways

Use the SAME climate-complete historical strict9 rows, SAME equal-taxon total weights and SAME source-subarm x 5-degree cell x local-season fixed effects in both models:

- M0: elevation only.
- M1: elevation + warmest-month maximum temperature BIO5 + log1p annual rainfall + annual solar radiation.

Predictor SDs are fixed using all attempted eligible observations, not successful ROIs alone. Outcomes are all nine components plus entropy and dominance. Nominal 95% intervals use 10-degree spatial cluster-robust CR1 covariance; repeated taxa across blocks, adjacent spatial dependence, missingness and botanical error remain outside these intervals. Multiple endpoints are descriptive without project-wide multiplicity correction. Coefficient attenuation or sign change is not a causal mediation estimate. Report collinearity, same-sample support and the zero-sum compositional coefficient check.

## Executed implementation plan

1. Tests first: ID drift/duplication, unknown preservation, palette closure, exact decomposition, empty common support, geographic boundaries and regression parity with explicit weighted dummy variables.
2. Minimal standalone numerical module; tests passed locally (10 tests). No inference-library imports.
3. New isolated branch and a new data-only workflow. Read only existing numeric artifacts/releases. Archive outputs with input and output hashes, then independently verify the white outcome and aggregation locally.
4. Open a review PR without auto-merging, altering the old manuscript or restarting any old image worker. No MNAR scenarios, image downloading, recolouring or human-review requirement is added.
