# MOIS: actual 569-photo low-resolution screen comparison, 20 October 2026

## Decision

**Do not enable the 320px screen as an image-exclusion filter. Continue the ordinary full-resolution ROI acquisition.** Even the least restrictive tested threshold (0.05) would lose 7 of 199 automatically classifiable photos, including 5 of the 93 white-labelled ROIs, for only 5.20% projected ROI-call time savings in this cohort. This is a decision from this pilot, not a universal rejection of all possible screening models.

## Completed real-image data

Source: run [39369517433](https://github.com/zuizui0223/mois/actions/runs/39369517433), disjoint tier03, compute bucket 0/128. There were **569 new original-photo attempts**, after excluding 3 previously processed photo IDs. **567 images were acquired**, **2 acquisitions failed**, and **199 photos were automatically classifiable into four hues**. The other **368 acquired images remain unclassifiable**, not confirmed non-flowering and not chromatic by default.

The 199 classified ROIs contain **93 white**, **81 yellow-orange**, **8 red-pink**, and **17 blue-purple** labels. All 567 acquired images received the ordinary full ROI-v4 baseline AND the separate 320px shadow screen. Zero screen runtime errors were recorded. No photo was discarded by the screen.

## Hypothetical screening effects recomputed from image-level CSV

| Detector confidence threshold | Images that would be skipped / 567 | Baseline-classifiable photos lost / 199 | White-labelled ROIs lost / 93 | Projected ROI-call time saved |
|---|---:|---:|---:|---:|
| 0.05 | 95 | 7 (3.52%) | 5 (5.38%) | 5.20% |
| 0.10 | 140 | 12 (6.03%) | 6 (6.45%) | 10.23% |
| 0.25 | 225 | 31 (15.58%) | 13 (13.98%) | 22.68% |

Actual baseline ROI calls totalled 1,978.196271434 seconds. The extra shadow screen took 19.672176363 seconds. Net projected savings subtract screen time on all acquired photos from baseline ROI-call time spent on the hypothetical rejected subset. These are **not realized production speedups**, do not include download/model-initialization/palette-classification time, and must not be extrapolated as an end-to-end 93m-image acceleration. The pilot itself did extra work.

The reference labels are the existing automated ROI classifier, **not independent botanical ground truth**. The small red-pink and blue-purple groups do not establish robust hue-specific recall. This single selected compute bucket is not a world-flora colour-frequency estimate. Preserve the data and use no automatic image exclusion.

## Source evidence and independent arithmetic check

Completed artifact: `mois-tier03-actual-pixels-and-shadow-comparison`, ID **12242257994**; ZIP SHA256 `248ee0a1dc6065ee2e1206362a1afc107c66bdf47b4940aa87f361237bc51e21`.

- `measured_photos.parquet`: SHA256 `ea918906b65fccb9c6fac8bda209a32aee3fd61b1b275d6f81fde0920440adfd`.
- `screen_shadow.csv`: SHA256 `adf6187e209098820e8ce2d883aa603a5f5610b6dd7342fec49748a43a778846`.
- `source.json`: SHA256 `285e9f86e13f4d9330c71ee3719898268338069e3db59f639c99bd88b7e84fd1`.

The ZIP and these three source hashes were checked after retrieval. All 567 unique measurement IDs, hue counts, threshold-specific discarded/classifiable counts and projected timing sums were independently recomputed from `screen_shadow.csv` and reconciled to `screen_shadow_summary.json`. Raw licensed photo images are not included in the result archive. The genuine new measured rows can be included in the completed-pixel checkpoint, with original unknowns and acquisition failures retained.
