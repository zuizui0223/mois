# MOIS: signed complete second p0004 original RGB union — 2026-10-20

## Completed deliverable

The 18,119 newly targeted p0004 original photographs have all reached terminal states and were archived in the immutable [second-wave photo-pixel release](https://github.com/zuizui0223/mois/releases/tag/mois-full93m-observation-second-original-RGB-p0004-v1). The existing 612,483 genuine FCP/MOIS/first-wave original image outcomes were then SHA-verified and photo-ID-deduplicated together with the second wave, yielding the **fully published, signed 630,602-photo p0004 union**:

**[mois-genuine630602-original-pixels-first-two-93m-waves-v1](https://github.com/zuizui0223/mois/releases/tag/mois-genuine630602-original-pixels-first-two-93m-waves-v1)**.

| Result scope | Original photo attempts / unique photo IDs | Genuine frozen four-hue ROI labels |
| --- | ---: | ---: |
| Original verified FCP/MOIS + first global source draw | 612,483 | 295,801 |
| Source-signed new p0004 photos | 18,119 | 7,085 |
| **Full p0004 union, exactly one record per photo ID** | **630,602** | **302,886** |

Within the new 18,119 terminal photos: **18,065 actual image acquisitions**, **54 acquisition failures**, **7,085 valid four-hue labels**, **10,980 acquired images without a classifiable four-hue flower**. The 7,085 newly classified photos comprise white **3,504**, yellow-orange **2,316**, red-pink **637**, blue-purple **628**. The second source wave added **253 nominal iNaturalist taxon IDs** with classified photos that had not been classified in the historical 612,483-photo source. A nominal taxon ID is not an independent accepted botanical species and this count must not be interpreted as such.

This is a complete **nested observation-probability p0004 source wave** with genuinely measured image pixels, *not* a 93,209,351-image flower-colour census. The global upstream ID index remains metadata, and photo detection/classification failures are not inferred to be non-flowering or nonwhite. No WCVP native-status filter, geo/climate model, or raw licensed image redistribution was applied.

## Incident, root cause, code correction and production validation

Full original image acquisition [run 39353343403](https://github.com/zuizui0223/mois/actions/runs/39353343403) succeeded, but the two original signed append attempts 39381982184 and 39382015130 failed at Parquet serialization with an original geographic field schema mismatch:

`pyarrow.lib.ArrowInvalid: Could not convert '29.7127033333' with type str: tried to convert to double (column latitude)`.

The immutable historical Parquet source stored geographic/climate metadata as numeric columns, while the new CSV-based photo source retained those values as text. The original concatenation also added missing columns to the historical DataFrame in memory. Neither error justifies redownloading or recolouring images.

[PR #131](https://github.com/zuizui0223/mois/pull/131) adds defensive numeric-metadata alignment to the historical source column schema, rejecting unexpected nonnumeric values instead of silently erasing them, and copies the original historical photo table before column alignment so that input records remain immutable. Unit regressions first failed (2 failed, 5 passed), then the corrected append CI and repository checks succeeded. A `main` push trigger runs the signed append once after merge, instead of requiring a new original photo-acquisition run.

**Production append run [39383128360](https://github.com/zuizui0223/mois/actions/runs/39383128360): SUCCESS**, original pixels reused (zero photo downloads), output 630,602 distinct original photo records with no duplicate IDs.

## Independently checked release files

Original data-only combined union:

- `mois_FCP_and_first_second_93m_actual_original_photo_ROIs.parquet`, SHA-256 `5d8580b291079b6c689433bc84e3b7446c177fd177bc5ad84daafad6129f47d2`.
- `mois_p0004_actual_photo_source_arm_hue_quality.csv`, SHA-256 `398dc89181a112400fa50ee642e42c1716b7684e3b5a8b418560cf4b9c9766f9`.
- `mois_p0004_actual_photo_append_receipt.json`, SHA-256 `4737a52e0c826a31cc34393379e51bbdb3e51d9316396ab3476fa527b00e6588`.

After production publication, [workflow artifact 12246062868](https://github.com/zuizui0223/mois/actions/runs/39383128360/artifacts/12246062868) was downloaded independently as a 23,076,151-byte ZIP (SHA-256 `6eb4f8dd914af56bf0fbf5bf160381d2a7bfa5d5b5a6c87d0226c4a346f155b3`). All Parquet/CSV hashes matched the published receipt; independently re-summing the eight source-arm rows yielded exactly 630,602 attempts and 302,886 classified results. The source/Parquet row-level and original-image lineage gates also passed in production; local container verification checked SHA and CSV/JSON arithmetic, not every Parquet row.

## Relation to the newer live incremental checkpoint

The later source-authenticated [partial checkpoint 39381982167](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39381982167-v1) includes *the same full* 18,119 p0004 photo IDs **plus** sealed p0008, tier03 shadow and p0016 subsets. Its **639,363 attempted / 306,186 classified** values are not additive to the complete 630,602-photo p0004 union. Do not sum the two releases.

- Complete p0004 source: DONE and signed.
- Live p0008 source: 36,375 new planned photo attempts, some sealed, remainder processing.
- Live p0016 source: 73,872 new planned photo attempts, some sealed, remainder processing.
- Low-resolution screening: disabled as an image exclusion due to observed white-hue false negatives.

The next genuine cumulative checkpoint must de-duplicate all actual source photo IDs, and never treat in-flight photo targets as measured.
