# MOIS — full p0008 original-photo RGB completion and six-arm pixel checkpoint

**Verified 2026-10-21 (JST)**

## 1. Genuinely completed original-photo outcomes

The *complete* third nested original-photo increment, p0008 run [39371199127](https://github.com/zuizui0223/mois/actions/runs/39371199127), finished all **36 of 36 source-signed RGB groups** (143 complete original photo-shard terminals). Its immutable official [full-wave release](https://github.com/zuizui0223/mois/releases/tag/mois-full93m-observation-third-original-RGB-p0008-v1) contains **36,375 unique original photo attempts** and **13,874 classifiable frozen four-hue images**, with **36,331 actually acquired originals**, **44 image-acquisition failures** and **22,457 acquired-but-unclassifiable original photos**. The original source `source.json` SHA-256 is `0038024ee40adfbc71aa427a112948ec5ce53a7aae9879c1b32bc06eabde`. None of these 36,375 photo IDs are metadata-only outcomes.

The next source-authenticated [six-arm checkpoint 39446764842](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39446764842-v1) was **successfully produced in GitHub Actions** and independently downloaded and checked by SHA-256 and CSV arithmetic.

| Source wave | Actual unique attempted originals | Original 4-hue classifications | Still pending new photo attempts |
| --- | ---: | ---: | ---: |
| Immutable FCP/MOIS / first p0002 source baseline | 612,483 | 295,801 | n/a |
| Full p0004 observation-hash source | 18,119 | 7,085 | **0** |
| **Full p0008 observation-hash source** | **36,375** | **13,874** | **0** |
| Partial p0016 observation-hash source | 17,408 | 6,698 | 56,464 |
| Partial p0032 observation-hash source | 11,264 | 4,358 | 136,807 |
| Partial independently sampled tropical/southern geo05 arm | 3,072 | 1,057 | 19,206 |
| Completed 569-photo 320px shadow study (no photo exclusion) | 569 | 199 | **0** |
| **All distinct genuine original-photo outcomes** | **699,290** | **329,072** | **212,477** |

Compared with the last checkpoint in the preceding conversation [39423567174](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39423567174-v1), this is **+21,015 actually attempted original photo IDs / +7,999 four-hue classifications**. This increase comes from p0008 +7,703/+2,923; p0016 +5,120/+1,931; p0032 +6,144/+2,443; and geo05 +2,048/+702. The complete p0004 and independent shadow original photo data were not double-counted.

## 2. Measured four-colour data versus source metadata

| Complete original photo wave | White | Yellow/orange | Red/pink | Blue/purple | Auto-classifiable fraction of photo attempts | White among auto-classifiable |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| p0004 | 3,504 | 2,316 | 637 | 628 | 39.10% | 49.46% |
| p0008 | 6,762 | 4,489 | 1,345 | 1,278 | 38.14% | 48.74% |

The two disjoint waves show a similar **photograph-conditional automatically classified hue composition**, but this is not evidence of 50/50 species-level evolutionary equilibrium or a global native-flora white fraction. The source opportunity frame overrepresents northern extratropical latitudes, some flora is species-filtered before upstream photo indexing, and original photos from the same observation are not independent. Flower-stage/ROI failures can also differ by habitat, taxonomy and colour. All estimates for species or environment need explicit observation/taxon/region/time conditioning and missingness audits.

The original **93,209,351 image IDs** still denote source photo opportunities, *not* 93 million photographs already measured for flower colours. The geo05 arm's original observation-level region-specific `selection_probability` is kept separate from the nested uniform probability arm and is now retained in the compact cumulative Parquet, too. Never mix geo05 rows unweighted into a global white-flower share.

## 3. Immutable evidence and automated continuation

Original six-wave data-only checkpoint: [run 39446764842](https://github.com/zuizui0223/mois/actions/runs/39446764842). Both tests and recovery succeeded. The ZIP artifact ID is **12276332510**, byte size **39,068,561**, SHA-256 `2fedd721d82377cabf97e134800d077a4132666b9309cec891cdacd596b6d86a`. Independently inspecting its four members verified the 3 output-file SHA256 digests in `checkpoint_receipt.json`, summed the per-wave × hue CSV rows and reconciled **699,290** original photo attempts and **329,072** classifiable hues. The production collector additionally checks each original source-shard SHA, photos.csv/terminal checksums, original photo IDs and unchanged FCP full-resolution ROI method before publishing.

This is a **partial compute-completion checkpoint**, *not* an unbiased probability-complete sample from the underlying entire 93m source index, even though p0004 and p0008 source increments individually finished. The 212,477 remaining planned selected originals are not counted as white, nonwhite, nonflower, or otherwise measured.

The ongoing p0016, p0032 and geo05 original photo workers have not been cancelled, duplicated, modified or forced through a lossy 320px screen. The authenticated source-aware 3-hour GitHub Actions watcher is configured to initiate data-only checkpoints **only when new successfully sealed original photo worker artifacts are present**. Do not start a competing wave until the current outstanding processing capacity is considered.
