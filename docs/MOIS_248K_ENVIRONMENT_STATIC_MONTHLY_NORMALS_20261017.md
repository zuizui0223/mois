# MOIS: independent 248,312-photo static environment and observation-month climate-normal join

**Repository:** `zuizui0223/mois`. This is an input-enrichment layer for global **between-species floral photo-colour ecology**, not an extension of FCP's within-species D/W and 50-km geographical polymorphism models.

## Actual measured image frame

Original source data are unchanged: **248,312 distinct photo IDs, 120,807 classifier-assigned four-hue photo labels**. Historical MOIS 100,543 photos and their interspecific original precipitation/vapour result remain frozen; FCP high-depth/photos and six genuine already-measured ROI pilot groups were brought in at the **photo/colour-data level only**, with 3,411 original literal overlaps removed as duplicate photo records, not as evidence to exclude source reuse.

### Real data enrichment this execution

The script `scripts/analysis/mois_sample_248312_worldclim_static_and_monthnormals.py` extracts raster values at each original photographic latitude and longitude. Inputs are:

- **WorldClim v2.1**, global 10 arcminute 19 BIO indices [official documentation](https://www.worldclim.org/data/worldclim21.html); 1970–2000 long-term baseline (not temperature in the 2023 or 2025 photograph month).
- Original **frozen FCP published 10m bio and srad ZIPs**, SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf` (see original exact source checks and GH Actions), and `c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026`. The workflow checks the exact official release SHA, and refuses source drift.
- WorldClim 10m monthly precipitation, monthly mean temperature and monthly solar radiation normals at the **authentic photographed calendar month**, only where a real observed_on date exists. These are average *seasonal climate conditions*, not observed-year monthly weather values.
- WorldClim 2.5 arcminute elevation (~4.6km latitudinal cell width), sampled at original photo sites. The program estimates **coarse 2.5 arcminute-grid** slope, local relief/roughness, TRI, TPI, eastness and northness from surrounding DEM cells using geodesic kilometer adjustments and refuses cells without valid neighbors. This is a landscape-scale terrain index, **not** site-specific fine topography or field slope.
- Original historical MOIS SoilGrids variables (pH, organic carbon, nitrogen, clay, available-water proxy) are preserved unchanged. **No new soil values are claimed** for the FCP high-depth records until a verified new SoilGrids spatial data import is completed.

### Source heterogeneity and time hierarchy

The photo-source descriptor distinguishes the original globally broad/shallow MOIS sample from high-depth FCP source species and other technical ROI groups; source and species/genus weighting must be explicit in subsequent environmental models. Original photo IDs may be reused between repositories without treating MOIS's environment slopes as independent confirmation of within-species FCP polymorphism.

The temporal variables are **three separate blocks**: (i) immutable geographic/soil/terrain backgrounds, (ii) long-term climate gradients (BIO and monthly mean climatologies), (iii) actual changing observation-year/month temperature/precipitation/radiation/VPD/water deficit, which require a separate real-year TerraClimate join. A photographed month's WorldClim normals are **not** the actual monthly weather experienced in the observed year. Photo time-of-day **must stay unknown** until an authoritative timestamp with timezone has been recovered; hour cannot be inferred from date alone.

### Reproducibility

`.github/workflows/mois-original-248k-worldclim-site-time.yml` retrieves the real compressed MOIS measured-photo table from original successful [MOIS run 38990701473](https://github.com/zuizui0223/mois/actions/runs/38990701473), downloads the actual raster ZIPs (pinned 10m bio/srad originals and documented official other variables), runs point sampling for all 248,312 original photos, retains nodata explicitly, and publishes a photo-level compressed combined environmental dataset, source-normal comparison diagnostics, photo- versus species-equal quartile descriptive summaries and a SHA-signed receipt.

Next model: at the **between-species** level, estimate whether new year-round photographed hue composition tracks BIO temperature, precipitation, radiation, elevation, terrain and soil after accounting for genus, source, observer effort and geography. Any date-resolved *within-site* climate-response extension remains a separate MOIS seasonal model, not an FCP D/W polymorphism model.
