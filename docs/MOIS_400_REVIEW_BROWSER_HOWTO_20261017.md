# MOIS: Browser-assisted, independent botanical open-flower review

**Owner:** `mois` seasonal/global floral colour analysis. This is an *interface only* for two real independent botanical reviewers; it **does not create a biological observation or colour classification**. The FCP within-species D/W paper and original MOIS global between-species moisture coefficient are unchanged.

## How the reviewer works

1. Open [the single-file review browser](../tools/review/mois_400_botanical_review_browser_20261017.html) locally in a web browser (desktop is recommended). No installation, account, server or upload is required.
2. Reviewer A loads only `mois_400_botanical_open_flower_reviewer_A_UNFILLED.csv`; reviewer B separately loads only the **B** original CSV. Both originals are SHA-frozen in the handoff documented in [MOIS_400_BOTANICAL_REVIEW_HANDOFF_20261017.md](MOIS_400_BOTANICAL_REVIEW_HANDOFF_20261017.md).
3. Enter a distinct `reviewer_id` and record the actual visible photographed flower stage and the tissue inside/at the apparent ROI. Check the **original photo**, not an observation-level iNaturalist flowering/fruiting annotation. The dropdowns preserve the exact authorized stage/organ values.
4. Move through 400 photos (different random orders for A and B), saving a `*_REVIEWED_DRAFT.csv` regularly. The same exported CSV can be reimported to continue. The original 400 rows and their license, photo URL, attribution and blind tokens are preserved by the UI.
5. Return A and B CSVs separately to the study coordinator. Do **not** open the other reviewer's answers or the coordinator-only source/date/season mapping.

**What the browser shows:** original iNaturalist photographic URLs, taxon name, photographer attribution, and the two classification questions. **What it must not show:** photo observation date/year, location, climatic data, colour classifier output, white-rate denominator, or previous reviewer answers.

The browser does **not** automatically save or upload review data. Its source CSV and in-progress human decisions remain in browser memory until downloaded to a CSV; open pictures are fetched from their normal iNaturalist S3 image server, which requires internet. A reviewer may select `uncertain` when the photo is inaccessible, unrecognisable or lacks an identifiable open flower. None of these photos is silently replaced by another after outcome review.

## Independent decision QC and third adjudication

The validator `scripts/analysis/mois_validate_400_blinded_botanical_reviews.py` compares both filled forms with the exact original SHA-pinned blank forms and a private coordinator mapping. It checks unchanged 400-token/photo/URL/licence membership, distinct reviewer identities, complete stage+organ votes, and source-level agreement. The generated `THIRD_REVIEW_PENDING_BLINDED_ORIGINAL_PHOTOS.csv` contains only disagreements and remains blind to the A/B decisions and climatic/temporal information. A third independent botanist can fill that file; the validator accepts it with `--adjudication`.

Only independently agreed or adjudicated **open flower plus actual petal/tepal** observations can become candidates for *future* original image flower-colour ROI measurement. This does not convert the original human stage answers to white/nonwhite colour labels or prove an adaptive seasonal response.

## Verification

`.github/workflows/mois-400-browser-review-safeguard.yml` downloads both pristine, original A and B 400-photo sheets, checks the exact SHA256 values, tests JavaScript syntax and CSV quotation/newline round-tripping with every real review row, and ensures the page contains no photo-month, climate, or automatic hue/source coordinate fields. The historical sample and all biological outcomes remain unchanged.

The standalone HTML is an additional convenience artefact; reviewers may use a spreadsheet instead, as long as they return the *same* source fields and controlled stage/organ values. User-supplied adjudicated human answers should not be uploaded into public GitHub Actions logs or public source repos. Only redacted aggregate counts can be committed after the independent human judgments.
