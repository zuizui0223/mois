# Exact source-input recovery and licence gate — global mois analysis

**Last audited:** 14 October 2026. **Action status:** the original GitHub Actions archive *identifiers* and source SHA256 values are known and historically verified; **the two ZIP files have not been uploaded to mois** by this documentation work. No current source permissions or new model estimation is implied.

## Original full-data artifact identifiers

| Artifact | Verified originating FCP reference | Source SHA256 / artifact identity | Availability |
|---|---|---|---|
| Deduplicated original photo/site/colour ZIP | [FCP Actions run 37861582446](https://github.com/zuizui0223/fcp/actions/runs/37861582446), artifact **11587115220**, named `fcp-global42111-original-photo-real-site-climate-soil-20261009`; original archive ~34.3 MB | `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0` | GitHub metadata checked 14 October 2026: **not expired**, scheduled expiry **2027-01-06T23:49:56Z**. Download before expiry; historical access can change. |
| Site-joined all36 abiotic ZIP | Artifact **11653183273**, required as `FCP_SOURCE_ALL36_ZIP_ID` in FCP PR #193 original source-verification workflow | `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca` | Exact byte SHA was verified by successful full-source FCP run **38502660130**, but current artifact expiry has **not** been independently retrieved in this audit. Do not presume it is permanent. |

Authoritative original source verification is [FCP Actions 38502660130](https://github.com/zuizui0223/fcp/actions/runs/38502660130) and FCP code `.github/workflows/fcp-macro-abiotic-hue-direction-original-20261013.yml` at source commit `0adef37c75c0e5bc135c447a861cd77eed34dbb2`. That workflow downloaded **both exact GitHub artifact ZIPs**, verified their SHA256, regenerated the 20 macro cohorts and 1,280 direction rows, independently reconstructed 352 original source WLS estimates, and checked original positive/negative results. Later FCP #194 matched 576 row results and 288 direct WLS results remain source-verified. These are computational replay, not independent photography or experimental replication.

## How to recover, without accidental Git redistribution

From an authenticated checkout of the **mois** repository with GitHub CLI (`gh`):

```bash
bash scripts/recover_mois_original_inputs.sh
# optionally use an external directory:
MOIS_SOURCE_DIR=/secure/local/flower_colour_inputs bash scripts/recover_mois_original_inputs.sh
```

The script reads the original two GitHub Actions artifact URLs using only the exact IDs above; validates each ZIP against its original SHA256 and tests ZIP integrity if `unzip` is installed. The files are kept in `raw_inputs/` (ignored by Git) by default or in the explicit path chosen by the user. Failure of an artifact endpoint or checksum stops the operation. It does **not** reclassify source photos or refit the 999 original bootstraps. **The script has been committed but has not been executed as part of the light manuscript CI**, so its successful download is not yet claimed.

## Required steps before journal reviewer access or source archival deposit

1. **Recover** both ZIP archives into a secure local directory while original GitHub Actions artifacts remain available, then confirm the two source SHA256 values and exact source photo/taxon/observation identity join.
2. **Inspect photo rights and privacy**. iNaturalist photographs and locality fields can have varying licences, contributor attribution and coordinate visibility. Do not publicly redistribute unreviewed images/metadata merely because GitHub Actions could download them; check per-photo licence and geoprivacy and prefer minimally necessary derivatives for stable publication.
3. **Package analysis and data with an anonymous reviewer interface** (Dryad private peer-review link or equivalent genuinely stable authorized repository). Include data dictionary, photo-hue code, WorldClim/SoilGrids version/depth/surface provenance, all source-environment join rules, source scripts and independent verifiers. Record the private reviewer link in the blinded manuscript **only after it exists**; the current placeholder is intentional.
4. **Archive publicly with a persistent identifier** after permission review and journal policy confirmation. GitHub and expiring Actions artifact URLs do **not** satisfy Global Ecology and Biogeography's long-term archive rule by themselves. Document the eventual DOI, licences, checksums and embargo/review access appropriately.
5. **Run end-to-end re-estimation on recovered inputs** before claiming independent current-repository raw-source replay. The previous FCP 2026-10-14 Actions run already verified source output once; the transferred mois tests still only verify copied source bytes, numerical receipts and synthetic invariants.

### Study boundary

Only global **between-species** moisture/photo-colour results are required for this manuscript. Do not copy the FCP 149,900-photo high-depth H1/H2 within-species set or its prospectively unopened species into the macro data archive. The same-original-photo between/within BIO12 descriptive sensitivity may be archived with the global macro source, provided source image/observer rights are respected.
