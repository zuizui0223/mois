# MOIS: complete two real original-photo wave additions and refit original water model

This source workflow is a continuation of the user-requested primary data expansion, **not a new taxon-polymorphism paper** and not a human botanical review. Photo IDs may be shared with `fcp`; only FCP's photo-geometry D/W/50km polymorphism estimands are disallowed in MOIS.

## Current real input

The original MOIS+FCP photo-level union was **248,312 unique original photo IDs / 120,807 true ROI-v4 four-hue photographic labels**. The *first* newly measured original photo wave [run 38994341776](https://github.com/zuizui0223/mois/actions/runs/38994341776) completed **4,096 actual downloaded original photos**, with **1,576 new four-hue classifications**; the complete SHA-enforced reassembly was independently replayed in [run 39000213324](https://github.com/zuizui0223/mois/actions/runs/39000213324). Result: **252,408 original photos / 122,383 four-hue classified photos**.

A *second* original source wave was independently frozen from all **3,430,502** authentic iNaturalist original source photo IDs (and explicitly outside the original 106,376 photo screen and current 248,312 measured photos) in [run 38995068942](https://github.com/zuizui0223/mois/actions/runs/38995068942). The real second wave is [run 38995240689](https://github.com/zuizui0223/mois/actions/runs/38995240689), 16 × 256 original photo-ROI source work units. **No source-only image is counted as a new colour measurement.**

## Source- and outcome-safe completion logic

The workflow `.github/workflows/mois-final-two-real4096-waves-water.yml` is triggered only on successful completion of the original second-wave source-specific GitHub Actions workflow (or on explicit manual invocation), and it rejects unfinished, failed or partial source-photo shards. It downloads the original signed photo-ID source receipts and **all 32 actual ROI source-photo terminal artifacts**, checks unchanged algorithm/model file hashes and every original photo-ID and source cohort membership, refuses to count the same original photo ID more than once, and publishes exactly **256,504 unique original photos**. The exact number of additional four-hue classifications comes **only from completed original image pixels**, not an estimate from the 39% first-wave yield.

The following job then joins genuine WorldClim site BIO5/BIO12/BIO14/BIO15, long-term solar normals and annual vapour pressure to **only the new original photo sites**, keeping original 248k site climate values intact. It fits the exact same 2026-10-13 transferred six-anchor × four-hue genus-centred/genus-equal 999-bootstrap and 24-contrast maxT procedure originally applied to the 100,543-photo core. A subset fit for the historical original photo source is retained as a nonindependent reference; same-genus/coarse geographic-cell conditioning is a critical secondary limit.

**This is not proof of heritable flower-colour polymorphism, climate-dependent fitness or local adaptation.** The original photo classifier has technical missingness; original WorldClim baselines are not the actual weather in the photographed year/month. The 3,420,814 unmeasured source photo IDs in the 1,000-taxon original upstream frame are **capacity for more real image processing**, not already measured data, and the additional 42,111-taxon candidate frame is a distinct larger data source problem.

This linked workflow makes the user-desired continuation photo measurement → real photo-ID union → original moisture refit executable without requiring another sequence of manual botanical review handoffs.
