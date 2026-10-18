# MOIS: why source-aligned reuse is much faster than another 4096-image run

**18 October 2026.** The current original photographic study has 313,535
genuine image-pixel attempts, 151,239 classifiable four-hue ROI outcomes and
71,427 white. A new fifth 4096 original photo ROI wave is being run as
16 expensive 256-image workers with only 2 workers in parallel; ~63% of the
first 512 genuine photo attempts had no interpretable four-hue ROI.

A separate [FCP 49,800-photo genuine full-wave run](https://github.com/zuizui0223/fcp/actions/runs/38992543115)
is already **successfully completed** with every original image actually
processed using frozen ROI-v4/EfficientSAM. Its final measured-photo artifact
includes exactly 49,800 original source photos from 498 named taxa.
Rather than paying the download, image decode, YOLO inference and SAM
segmentation cost again, this workflow REUSES those prior real pixels'
complete original photo/colour/observation/source signatures.

Older MOIS FCP reuse workflow was gated on a historical 256,504-source
MOIS archive and its source import could not attach directly to the
new 313,535 photo-feature release. This workflow reads the actual
313,535 immutable original photo-ID/Zstd Parquet source instead and
validates the entire FCP original 49,800 measured result receipt
before even counting a single new original hue.

## Do not confuse 49,800 *measured in FCP* with 49,800 *new to MOIS*

A photo can occur in both repos; repeated original photo IDs are counted
once, historical MOIS hue/day/taxon values are never overwritten,
and discordant original classifier labels are explicitly audited.
The independent ongoing fifth4096 source photo-ID reservation and the
separate 4096 SI-targeted season reservations are also excluded from
this new source arm so their future actual pixels won't be counted twice.

Only the verified genuinely novel photo IDs from complete FCP
measurement are appended. The exact count and new classified-hue yield
must be computed by CI, **not guessed from FCP's total 49,800 attempts**.
Unclassifiable FCP photo ROIs remain UNKNOWN, not "nonwhite."
The resulting full photo table and new-source-delta parquet plus audit
receipt are published in a versioned GitHub Release; no raw licensed
original RGB image bytes are copied.

## Critically preserve the ecological meaning

FCP wave4 samples **498 deliberately high-depth taxa**, not the full
world flora. The additional visible-hue results must be identified as a
separate source-arm, and must NOT be naively combined as an unbiased
new global white/coloured prevalence or treated as independent
confirmation of FCP polymorphism and biological selection. Instead
compare original 313,535 core/photo/genus baseline and this new
high-depth source arm separately, with species/genus/source weights,
photographic missingness and geographic uncertainty.

**Next pipeline step:** Extract annual WorldClim 19 BIO/VAPR/solar and
coarse DEM terrain ONLY for genuinely novel 49,800-source photo IDs
not already in MOIS; reuse original annual climate and terrain Parquet
values for the original 313,535 photos. Soil gaps are not imputed.
This is much more efficient than rerunning original YOLO/segmentation
on identical photographs, and permits later snapshot-level climate
models after source-selection sensitivity.

Code: scripts/analysis/mois_instant_reuse_fcp49800_after313535.py
Workflow: .github/workflows/mois-reuse-complete-fcp49800-after313535-without-redownload.yml
