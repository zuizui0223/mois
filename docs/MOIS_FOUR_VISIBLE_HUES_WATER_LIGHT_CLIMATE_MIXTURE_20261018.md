# MOIS four visible photographic hues: which colours occupy the changing white fraction?

**18 October 2026, main article complementary ecology; observation-month seasons belong in SI.**

The main source holds **313,535 actual original processed RGB photo attempts, 151,239 classifiable visible four-hue flower-like ROIs**, **71,427 white**. Its apparently near-half white proportion is NOT a maintained exact 50:50 evolutionary frequency or fitness optimum. A contrast of white against one undifferentiated "pigmented" class obscures the biological differences between yellow/orange, red/pink and blue/purple signals. Visible nonwhite is NOT equivalent to anthocyanin and visible white need not lack ultraviolet or other pigments.

**Question:** When photographed white hue becomes less common in dry or high solar sites, does yellow/orange, red/pink, blue/purple or some mixture become more frequent? Are the differences a direct source composition pattern only, or does climate improve prediction of four mutually exclusive source hues for nominal taxa withheld by genus/geography?

Original-source-locked design:

1. Read the exact previous verified GitHub Releases, checking 313,535 actual photographed photo IDs, 151,239 colour-known states, source SHA and full original BIO1–19/VAPR/solar/terrain site Parquet. Do NOT count metadata-only provider IDs or the separately season-gap-targeted 256 image source. No new RGB processing or hue relabelling is necessary.
2. Create annual precipitation tertiles and solar tertiles using **all photo attempts before hue filtering**, so classifiability selection is explicit; report the actual 9-cell attempted-photo counts, known four-hue counts and their multinomial proportions. Report all four "wet minus dry" photographic hue differences within fixed solar tertiles; the four changes mathematically sum to zero within each solar level.
3. Create **four real photo-label fractions per named source taxon** only after including photos with all climatic predictors, and restrict to nominal genera with at least two taxa. 5-fold heldout-genera and separate 5-fold approximate 10° source-taxon-centroid geography test generalization. In each training fold standardise only the training sites and PCA-compress **temperature 11**, **moisture BIO12–19 + annual VAPR 9**, and **solar 1** using 90% retained environmental variance (max six components per block) in addition to nuisance geography and source recruitment; fit genus-weighted **multi-output Ridge** to the four class proportions. Normalize heldout predictions to a four-class simplex (nonnegative; sums to 1). Compare weighted multiclass Brier score for the full climate model, geography only and each full-minus-one-climate-block fit on the **same held-out taxa**. Loss improvements are predictive, not causal biological fitness evidence.
4. Preserve source confidence boundaries: automated ROI has not independently confirmed open petals, genetic pigment pathways or flower fitness; taxonomic/observer/geographic/source-classifiability biases remain. Annual WorldClim 1970–2000 is a historical climate normal, not actual observation-month weather.

The main 33 abiotic variable study already showed **water and temperature are different leading conditional prediction candidates** depending on sample soil completeness and whether holdout groups are genera or geographic blocks. The new four-hue perspective asks whether whole-source white-versus-other averaging conceals ecologically distinct colour-category associations.

Source-locked workflow emits a nine-cell visible-hue table, 12 wet–dry colour contrasts (four hue classes × three solar bins), genus/geography-heldout climate block multiclass Brier table, and a checksum receipt. On successful trusted-main CI, a versioned GitHub Release retains the actual data/model derivatives. Raw copyrighted original photo bytes are never re-hosted.

Scripts:
- scripts/analysis/mois_original313535_four_visible_hues_climate_mixture.py
- tests/test_mois_original313535_four_visible_hues_climate_mixture.py
- .github/workflows/mois-313535-four-hue-climate-genus-geography.yml
