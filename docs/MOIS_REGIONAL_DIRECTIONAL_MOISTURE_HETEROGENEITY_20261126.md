# MOIS — Regional heterogeneity of directional moisture–colour composition slopes (2026-11-26 JST)

## Executed analysis: only the existing MOIS database

Follow-up to [directional moisture analysis](https://github.com/zuizui0223/mois/blob/d3eb1e9e873ec2809c07c34c9417f7caa4ed45d1/docs/MOIS_DIRECTIONAL_NATIVE_HUE_MOISTURE_NULL_20261126.md). No Island data, new photos, new occurrences or external traits. Reused the prior reconstructed native NPZ with exact SHA256 b290f13acdca2264b17c378cdafc7af5d56464952788dbad8f0e710da9d18693, 37,104 catalogue species, 34,781 native donors, 2,271,852 strict-native pairs. Rechecked all 18,577 x 9 signed cell counts and no duplicate species-cell pairs. This is not a new original GEOS overlay or an independent species-membership reconstruction. Existing native reference, geography and modern climate ZIPs and members were reverified against the prior pinned SHA256 and receipts.

Replayed the same 999 full9 and 999 white-fixed species-label draws (seed 20261125). Every one of the previous global slope vectors exactly reproduced (maximum absolute difference 0), including observed coefficients. All native donor strata defined by genus x native footprint bins x six-sector native affinity as before. Labels globally synchronized across cells and sectors, not separate independent permutations by region. Full9 mobility 19,111 species in 3,638 strata; white-fixed8 7,124 in 1,632 strata.

## Design

Primary previous 5,491 within-500km edges, each endpoint >=100 native species; 5,473 white-fixed sensitivity edges. For six-region comparison retain both endpoints in the same fixed coordinate sector and at least two edges per 500km block within sector: 5,488 full9 edges, 5,470 white-fixed edges (three excluded in each). Same wet-minus-dry colour response and ten covariates as prior directional analysis. Fit each sector's own within-block centered regression, including all the previous nuisance terms. Rescale all rainfall slopes to the same global full9 SD of log1p annual precipitation contrast (0.16146474296483554), avoiding sector-specific units. Minimum 100 edges and 10 blocks required for each sector. All six pass. These are coarse coordinate sectors, not ecological biomes or formal biogeographic realms.

Primary test: for each observed or globally synchronized null draw, subtract the six-sector mean slope vector from each sector, flatten the contrast matrix, center with the common mean of the 1,000 observed-plus-null vectors, and calculate squared Mahalanobis distance using a shared covariance pseudoinverse. Symmetric treatment of all draws and observed vector, rank p including observed. Rank 40 (5 region contrasts x 8 compositional dimensions) in full9, rank 35 in white-fixed8. Secondary: all-sector joint vector test without subtracting regional mean; Euclidean squared regional contrast departure sensitivity; within-sector joint vector tests with Holm6. Individual hue x sector slopes are descriptive diagnostics; post-result multiplicity audit across 54 full9 or 48 whitefixed coefficients. The analysis plan was specified before this numerical run but prior global results were already known; not independently preregistered.

## Results

| Fixed coordinate sector | Full9 edges | Within-sector 500km blocks | Joint vector rank p | Holm6 p |
|---|---:|---:|---:|---:|
| Americas_N | 1,445 | 103 | .338 | 1.000 |
| Americas_S | 837 | 57 | .268 | 1.000 |
| Africa_Eurasia_W_N | 971 | 88 | .026 | .156 |
| Africa_Eurasia_W_S | 298 | 37 | .234 | 1.000 |
| Asia_Oceania_E_N | 1,151 | 132 | .900 | 1.000 |
| Asia_Oceania_E_S | 786 | 46 | .797 | 1.000 |

Blocks can straddle sector boundaries; counts per sector are not necessarily independent or additive. All permutations are synchronized. Primary full9 regional heterogeneity rank p=.291 (score 44.17228145; rank40; null reference 95% [24.7874,57.8120]); all-sector joint p=.248, Euclidean sensitivity p=.340. Whitefixed8 heterogeneity p=.551 (score 33.62019447, rank35); all-sector p=.455, Euclidean p=.681. Minimum sector joint p=.123 in whitefixed8, Holm6=.738. No corrected support for exceptional regional directional response.

For transparency, full9 excess coefficients in pp / common rainfall SD (observed minus mean null):

| Sector | White | Yellow | Orange | Red | Pink | Magenta | Purple | Blue | Bronze |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Americas_N | +.0418 | +.0381 | +.0124 | -.0020 | -.0188 | +.0090 | +.0107 | -.0034 | -.0877 |
| Americas_S | -.0352 | -.0770 | -.0051 | +.0421 | +.0287 | +.0106 | +.0502 | -.0367 | +.0224 |
| Africa_Eurasia_W_N | -.1799 | +.0679 | +.0599 | +.0089 | +.0338 | +.0212 | -.0486 | +.0106 | +.0261 |
| Africa_Eurasia_W_S | -.2101 | +.2531 | -.0658 | +.0227 | -.0695 | +.0303 | +.0785 | -.0375 | -.0018 |
| Asia_Oceania_E_N | +.0902 | -.0274 | +.0203 | -.0032 | -.0140 | -.0028 | -.0328 | -.0238 | -.0065 |
| Asia_Oceania_E_S | -.0014 | -.0270 | -.0063 | +.0398 | -.0168 | -.0262 | -.0084 | +.0562 | -.0099 |

The smallest individual unadjusted p=.002 (orange, Africa_Eurasia_W_N, observed +.062897 vs null +.003022) becomes .108 after Holm54; smallest white-fixed p=.018 becomes .864 after Holm48. These are post-result diagnostics, not a new orange-selection hypothesis confirmed by this sample. Regional directions can differ descriptively; no evidence that differences exceed the strict conditional null. Non-significance is not equivalence, and conditioning may absorb genuine historic ecological structure.

## Conditional prediction tests

Outcome: vector of observed wet-minus-dry hue proportions (pp) minus **analytical** conditional null expectation, calculated from mean label proportions inside each fixed donor stratum; fixed nonwhite denominator for white-fixed sensitivity. This uses the known full frozen catalogue; held-out predictions are conditional spatial diagnostics, not unseen-species or independent botanical validation.

Known-sector 5-fold CV: 177 projected 1000km blocks for full9, 176 for whitefixed8, greedily balanced across sectors, each block kept intact; zero shared cell endpoints between train and test. Same common nuisance slopes and sector intercepts in each model; add common rain slope or six separate sector rain slopes. Training-only covariate scaling. Leave-one-sector-out: no sector intercept or region-specific slope extrapolation, pooled nuisance model vs pooled model plus common rain. All models least squares; loss = mean squared error across edge and hue proportions. Bootstrap 20,000 1000km blocks conditional on fixed OOF predictions; does not propagate model refit or new-region uncertainty.

| Model | Full9 MSE (known-sector blocked CV) | Whitefixed8 MSE |
|---|---:|---:|
| Baseline without rain contrast | .7972141167 | 1.7938679210 |
| Common rain contrast slope | .7972494711 | 1.7939046735 |
| Sector-specific rain contrast slopes | .7977932888 | 1.7959840445 |

Sector-specific vs common rain full9 MSE **worsens by 0.06821%**; absolute reduction -.0005438177, conditional bootstrap interval [-.0015585968,+.0004846989], positive in 2/6 sectors. Whitefixed8 worsens by 0.11591%, interval [-.00497426,+.00050039]. Common rain vs baseline also no improvement. Leave-one-sector-out full9 baseline .7972594631 -> common rain .7974981195, deterioration 0.02993%; interval for MSE reduction [-.0006798395,+.0001585273], positive in 4/6 sectors. Whitefixed8 baseline 1.7947708314 -> 1.7953653092, deterioration 0.03312%. Thus regional flexibility does not demonstrate additional predictive value even within known regions, nor does common slope generalize with a measurable benefit.

## Verification and deliverables

Actual local computation completed both 999-draw models, 102 sector-hue coefficients, four multivariate prediction comparisons per model. Independent SVD replay of omnibus rank p/score, direct per-coefficient mean and p recomputation, independent prediction MSE and fold endpoint leakage verification passed. Maximum absolute recalculation error 2.842171e-13, compositional closure error <=3.22e-15. Initial container execution request with unsupported streaming mode did not launch; reran normally and completed. No GitHub Actions execution is claimed.

Downloadable local reproducibility archive MOIS_regional_moisture_20261126.zip (5,740,034 bytes; SHA256 3de5218d542284d25d001ba4f2c0d95319b25e3a85be2a3e19f4e344951117aa) contains actual region_moisture.py, verification script, reused dependencies, 102-row CSV, observed and all null vectors, joint scores, full OOF predictions and responses, summary JSON and manifest. Actual code SHA256 ffeb6af08f76a04f6d8797a3ba942a09c11f993c965aba5d41228547b95d58e6. Both NPZ draw SHA256 cc4f3241764f5c7081da06751b28d295b7be9d2ee6ccecf62dbb448ff99109d7 and cac98baea3dc8204936c9b300b7c11025bb507ab59ed4c863c9068115a8cef8c. All ZIP members validated against manifest. Raw species-cell rows and photos not redistributed. This GitHub commit saves the result report; archive contains code and numeric arrays. No main merge, PR #206 remains draft.

## Scientific decision

No corrected support for directional moisture-response heterogeneity at the tested coarse sector scale, conditional on the specified genus/native footprint/geographic affinity null. Do not turn descriptive sign reversals into proven historical biogeography or conclude that environment never matters. Stop expanding individual moisture x colour x region tests unless new independent evidence motivates them. A different and potentially more informative MOIS-only question is whether geographic flower-colour composition can be reconstructed from taxonomic turnover using training-region genus colour distributions (hold out entire regions and species, compare with global frequency baseline, measure incremental predictive value and coverage); use actual documented occurrence assemblages and avoid target-species colour leakage. This proposed next analysis is distinct from the present directional null and is not yet executed.
