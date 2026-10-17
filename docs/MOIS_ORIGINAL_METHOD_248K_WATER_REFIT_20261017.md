# MOIS original six-anchor moisture model, now on 248,312 original photo IDs

**Priority:** increase measured photo count to the source/provider/compute limit, then update **the same global photographed four-colour–moisture question** as the historical 100,543-photo study. Source images may be shared with FCP, but its intraspecific D/W/50km polymorphism analyses are not copied.

The fixed current data contains 248,312 original photograph IDs and 120,807 four-hue-classified photographed floral regions. A second, independent new-photo measurement workflow is actually processing the next 4,096 original photo IDs, preselected without inspecting colour from a 106,376-ID sample of 1,000 plants. The source archive has 3,430,502 IDs and a larger 42,111-taxon upstream opportunity frame; none is incorrectly called 'already classified'.

## Original method, not a new favourable-test invention

This PR imports **the existing**, byte-identical MOIS historical estimator functions from `scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py`. Do not substitute the exploratory genus-cluster/BH 248k model for this main moisture analysis.

- The **same six original ecological anchors** are selected independent of outcomes: warmest-month temperature BIO5; annual rainfall BIO12; driest-month rainfall BIO14; rainfall seasonality BIO15; solar radiation; atmospheric vapour pressure.
- Four fixed coarse hues: white, yellow/orange, red/pink and blue/purple. Species photo fractions, **genus-centred and genus-equal**, controlling `|latitude|`, sine/cosine of longitude.
- Same **999 shared-genus bootstrap draws**, four hues × six preselected anchors, 24-contrast familywise **maxT simultaneous 95% CIs**, along with species-equal and same-genus/coarse-cell sensitivity.
- Re-estimate the genuine historical source-only sample as a comparator **on the same source/normalisation definitions**, then refit the newly enlarged photo data without double-counting photo IDs or equating the two partially shared sources to independent biological replication.
- WorldClim v2.1 original photo-site 1970–2000 climatology is used; annual vapour pressure is sampled from the authentic official 10-minute 12 monthly TIFFs, not derived from a made-up weather value. Existing WorldClim site BIO5/BIO12/BIO14/BIO15 and 12-month solar are source-verified. WorldClim solar raster values are in kJ/m²/day; a legacy raw column label `Wm2` is inaccurate and will not be presented as a valid source unit. Coefficient standardisation is unaffected by unit multiplication.

The historical original [MOIS main source](docs/MOIS_GLOBAL_MOISTURE_MANUSCRIPT_20261014.md) reported a global congeneric white–rainfall slope +2.255 percentage points per original source SD, and a vapour slope +3.952. This comparison **does not claim identical 26-variable complete-case source support**: the photographic sample has expanded, the original strict-image filtering columns have not yet been reconstituted for every new image, and the new vapour source is sampled separately. Unlike the historical manuscript, this stage is explicitly a refit on a broadened shared-source cohort, not a new independent confirmation.

The upcoming report must foreground **bootstrap familywise intervals and the same coarse-cell geographic restriction**, not nominal p values or unadjusted positive contrasts. A positive photographic association does not demonstrate water-dependent pigment fitness, heritability, experimental causality or FCP's genetic colour polymorphism.

The complete photo-environment frame and SHA for the 248k records are in [MOIS original WorldClim workflow](https://github.com/zuizui0223/mois/actions/runs/38992341847). This PR's source-scoped actions produce machine-readable expanded water slopes, actual classifiable support and complete-source SHA receipts. **No unclassified source photo is included as a white/coloured denominator.** 
