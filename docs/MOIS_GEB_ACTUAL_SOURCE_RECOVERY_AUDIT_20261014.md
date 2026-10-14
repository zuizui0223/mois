# Real source ZIP revalidation and photographic data-release boundary (14 October 2026)

**Status: source bytes recovered and audited in this working execution environment; not deposited or public-released.** The full original input archives were obtained from historical FCP GitHub Actions artifacts **11587115220** and **11653183273**. This step is no longer just an artifact-ID lookup: *both actual ZIP files were downloaded*, their complete SHA256 verified and their compressed member integrity tested. This is a computational recheck on the *same historical source records*, not a new biological source cohort.

## What was recovered and what was not

| Historical source | Original artifact | Exact SHA256 check | Package contents |
|---|---:|---|---|
| Original flower-colour/site measurement-table archive | 11587115220 (34,282,732 bytes) | `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0` — **PASS** | 8 compressed-table/summary members; 42,111 breadth original photo records and 85,337 original species×cell photo records |
| Expanded 36-environment/site measurement-table archive | 11653183273 (38,759,029 bytes) | `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca` — **PASS** | 4 compressed-table/summary members; 100,543 source photo-site matched environment records |

**Material distinction:** these two archives contain measured photo-colour and environmental tabular data, including source identifiers, licensing/attribution and photographed coordinates. **They do not contain the original photographs' image pixels**. A true independent image-classifier remeasurement would additionally require original photographs, ROI-v4 model and segmentation pipeline weights/versions. Full statistical source-table replay is a different requirement from pixel-level remeasurement.

### Executed sample-identity audit

An independent Python implementation re-opened both original archives directly, deduplicated source photo IDs, joined `inat_taxon_id`, `observation_id` and `photo_id` one-to-one, restricted exact source classifications and photo-site environmental completeness, and reconstructed the original genus-equal world climate26 BIO12 and mean vapour-pressure slopes. It **did not** import the original estimator or rerun the 999 whole-genus bootstrap worlds.

| Recomputed denominator / numeric | Independent real-ZIP result |
|---|---:|
| Unique historical photos | **100,543** |
| Source four-hue classifiable | **45,953** |
| Climate26-complete classifiable photos / species | **45,836 / 25,028** |
| Soil36-complete classifiable photos / species | **31,925 / 20,786** |
| Worldwide nominal genus-eligible species / genera | **20,694 / 1,934** |
| World genus-equal original white–BIO12 photo-colour slope | **+2.254545235123577 pp/source SD** |
| World genus-equal original white–vapor-pressure slope | **+3.951850325152856 pp/source SD** |
| Independent four-hue coefficients, 2 axes | **8/8** agree with the original frozen 1,280-row coefficient CSV |
| Largest absolute coefficient difference | **3.51×10⁻¹⁴ percentage points** |

The independent raw-source checker is `scripts/analysis/verify_mois_original_source_two_axis_wls_20261014.py` and the **aggregate-only**, machine-readable executed receipt is `results/mois_original_source_independent_recovery_audit_20261014.json`. It stores **no** photo IDs, full latitudes/longitudes, individual photographer attributions or plant locations. The CI checks exact source-receipt and stored coefficient consistency and compiles the new script, but the raw ZIPs are deliberately absent from GitHub CI. This is an **independent source-table point-estimate reconstruction**, not independent source replication or bootstrap CI recalculation.

## Original photograph licences and locality status (aggregated only)

The recovered original source photo rows include the following recorded licence categories, **not verification of legal rights to republish**:

| Original photo licence | Unique photo rows |
|---|---:|
| CC BY-NC | 78,367 |
| CC BY | 13,795 |
| CC0 | 3,229 |
| CC BY-SA | 2,715 |
| CC BY-NC-SA | 2,433 |
| CC BY-NC-ND | 2 |
| Missing/unspecified | 2 |
| **Total** | **100,543** |

The original recovered source status flags were **100,348** `VALID_PUBLIC_ORIGINAL_PHOTO_POINT` and 195 other statuses: 91 taxon-ID inconsistencies, 53 cultivated/not-wild flags, 31 without unobscured public positions, 14 source observations not returned, five missing exact photo identities and one historical-cell mismatch. They are **source status fields**, not independent reassessments of modern geoprivacy. No individual photo/location information is published in this aggregate audit. The presence of noncommercial, no-derivative or unspecified photo licences, attribution fields and precise locations means that archiving raw records or original pictures needs explicit human permissions and privacy/licence review first.

## Executable recovery/replay for an authorised reader

```bash
# From an authenticated local clone of mois
bash scripts/recover_mois_original_inputs.sh

# Requires Python with numpy and pandas, plus the already committed 1,280-row CSV.
python scripts/analysis/verify_mois_original_source_two_axis_wls_20261014.py \
  --photo raw_inputs/original_100543_photo_sites.zip \
  --environment raw_inputs/original_100543_photo_sites_all36_environment.zip
```

The recovery helper itself remains **unexecuted** in this code/repository CI; the 2026-10-14 retrieval was performed through the connected GitHub artifact download and the verification in a separate working environment. A local clone of mois, GitHub Actions and the journal's independent reviewer data archive **do not automatically inherit those recovered files**. The archival path and reviewer link are therefore **not yet established**.

## Remaining nonoptional publication gates

1. **Rights and original-site privacy:** inspect licensing and attribution per retained original photo, geoprivacy and source-site coordinate policy; separate permitted derivative statistics from images or person/source-level records. Do not upload the raw ZIPs directly into public GitHub.
2. **A reviewer-accessible, anonymous durable deposit:** create a stable, licence-compliant archive containing sufficient input derivatives and complete model code with a private reviewer-access path and permanent DOI/publication link. The present SHA values and temporary Actions artifacts cannot replace it.
3. **Extended verification:** reproducibly execute the full 20-source-scope bootstrap pipeline with 999 whole-genus resamplings and the exact spatial-support source pipeline against the newly recovered local ZIPs; the source FCP original historical run already achieved this but the new `mois` manuscript workflow currently has only eight independently reverified point estimates and receipt CI.
4. **Human GEB submission sign-off:** complete authorship, ORCID, title-page, cover-letter PDF, data/photographer licences, acknowledgement, funding and anonymization checks, final figure accessibility/resolution and permissions.

**Scientific nonclaim:** Verified positive source-wide white–BIO12 and vapour-pressure point coefficients cannot replace strict-photo/source-locality null sensitivity, external botanical nonconfirmation or demonstrate pigment/fitness adaptation. The independent eight-coefficient check does not turn post hoc ecology into a new confirmatory biological experiment.
