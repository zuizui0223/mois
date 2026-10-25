# MOIS: flower-region dominant CIELAB Stage1 v2 (2026-10-25)

## Research decision

Do not require anatomical classification into individual petals, calyces,
stamens, bracts, florets or semantic flower organ masks for the 93,209,351-photo
Stage1 pass. Treat the **dominant colour of a YOLO-detected flower region**
as the operational *estimated major petal colour*. It is a photo-ROI proxy,
not a histologically identified petal or calibrated spectral reflectance.

Keep **white** as the focal ecological gradient, but jointly retain **lightness
L*, chroma C*, hue h°, runner-up colour and nine historic palette proportions**
to detect greenish/brownish or otherwise inconspicuous colour responses on
islands. Brown, green and white are never merged into one category called
"drab". A genuine green flower and green foliage inside the detection box
cannot be separated by this method alone; preserve that uncertainty.

## Actual implementation

New `scripts/acquisition/mois_flower_box_lab_colour.py`:

- 500px original image already decoded for full93m Stage1, a single YOLO
  flower-detection box per photograph; sample 80×80 sRGB box pixels once.
- Convert sRGB (D65) to CIELAB with standard piecewise gamma linearisation
  and XYZ-D65 whitepoint. Fix bins 12.5 L* × 20 a* × 20 b*; select most
  populous Lab bin, report **median L*/a*/b*** of its original pixels as
  principal modal colour. Keep the modal fraction, second-bin fraction,
  and RGB display swatch separately. This is more precise than blindly
  converting a nine-colour label to a fake CIELAB point.
- Report (C^*=\sqrt{a^{*2}+b^{*2}}\), hue (h=\mathrm{atan2}(b^*,a^*)\)
  in [0,360). The mode hue is **missing** below C*=12, rather than claiming
  0° for a white or grey flower. Circular chroma-weighted mean hue uses only
  pixels C*≥12; diffuse/inconsistent directional signal remains undefined.
- Report region L* and C* 10th, median, 90th percentiles; area fractions
  below C*=12, nominal white-like (L*≥80 and C*≤15), nominal green-like
  (hue 105–190° with C*≥20), nominal brown-like (hue 20–100° with
  L*<55, C*≥18). **These fixed thresholds are descriptive QA fields,
  not calibrated botanical class labels.** Brown has a warm hue and low
  L*; it is not equivalent to every low-saturation or dark patch.
- Center-vs-periphery median Lab ΔE76 is retained as a possible foliage/
  background/multicolour diagnostic, **never as a hard automatic exclusion**.

Signed Stage1 `floral_box_lab_v2` mode appends these features to the
original five HSV box proxies (backward-compatible historical data), pins
new photo packet schema
`mois_all_original_photo_medium_pixel_stage1_CIELAB_dominant_v2` and
model identity + `mois_flower_box_srgb_cielab_mode_v1`; the existing v1
packet SHA/receipt cannot silently be resumed under v2.

`mois_93m_persistent_worker.py --colour-feature-mode floral_box_lab_v2`
implements it in the persistent same-fetch source pipeline; the huge archive
job count or original photo-ID sampling probability is unchanged.

## Executed validation

1. [Offline CI](https://github.com/zuizui0223/mois/actions/runs/40197634252)
   passed synthetic canonical D65 RGB, flat white (hue undefined), pure
   green and brown, 50/50 mixed-colour, signed full-source packet and
   v1-to-v2 resume-refusal tests.
2. [256 ORIGINAL K1 photos downloaded and measured with warm pinned
   YOLO11](https://github.com/zuizui0223/mois/actions/runs/40197802202)
   completed: original source 128 previously four-colour-successful and
   128 original ROI-failed photographic images, exactly 256/256 decoded;
   251 YOLO flower-box candidate images received actual device-Lab major
   ROI measurements. Median of the 251 dominant Lab modes: L*=**68.43**,
   C*=**20.42**. Photos with ≥25% nominal green-like BOX pixels: **92**;
   ≥25% brown-like BOX pixels: **33**; ≥25% white-like BOX pixels: **37**.
   **These are deliberately outcome-balanced technical calibration images.
   They are NOT counts/rates of botanically green, brown or white flowering
   plants, and cannot establish island prevalence.** The boxes may include
   leaves, sepals and background; there is no independent organ truth.

## Scientific analysis

- For existing 1,043,930 fully attempted original photos, 390,541 strict
  source-recovered nine-part ROI vectors allow original flower-ROI dominant
  palette and second-largest component; 36,224 ambiguous are separate
  exploratory candidates. [Actual signed 1m-photo main-colour proxy
  computation](https://github.com/zuizui0223/mois/actions/runs/40197966700)
  does not require re-downloading images.
- **Lab cannot be recovered from nine palette fractions alone.** Older
  source ROI Lab hard-mask summaries, if retained, should be reused by
  exact photo ID. Otherwise use new actual sRGB image access, not a
  palette swatch-to-Lab substitution.
- Multivariate island analysis should first keep white percentage as a
  core hypothesis, test separately adjusted L*/C*, and use circular
  hue models or sin/cos, not straight degree differences. Include
  moisture, high elevation, source taxon/genus, geographic support,
  season, observation process and ROI/nonresponse probabilities.
- More geographically uniform source photo sampling is required before
  inferring true green/brown floral prevalence. A modest *blind human
  image audit* is desirable to distinguish green/brown flowers from leaves
  or damaged flowers, **not a mandatory organ segmenter for 93m images**.

## Status

Feature **implemented** with offline tests and a real 256-image Lab pilot.
Global 93,209,351 original photo pixels have not been scanned with Lab v2,
and no large billed worker fleet was launched. No botanical petal-specific
ground truth or independent RGB colour-chart calibration has been obtained.
