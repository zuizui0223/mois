# Upstream FCP source capacity audit: 93.2 million photo opportunities, not 100,543

**16 October 2026.** The earlier Mois answer incorrectly treated the **100,543 original measured-photo union** as the limit of the underlying photographic source. The user correctly requested a check of the upstream **FCP** source rather than a reanalysis of only the derived climate/soil ZIPs.

## Verified upstream evidence (real original provider metadata, not new colour measurements)

The upstream [FCP PR #200](https://github.com/zuizui0223/fcp/pull/200) and its [full-source count run #38749681793](https://github.com/zuizui0223/fcp/actions/runs/38749681793) used a pinned **2026-09-27 iNaturalist Open Data metadata snapshot** (35,966,549,880-byte archive, SHA256 `bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e`). The target is the **original frozen 42,111-species frame**, independently verified against frame SHA256 `5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc`. An independent audit [Actions #38818718877](https://github.com/zuizui0223/fcp/actions/runs/38818718877) passed and emitted the following exact **source metadata opportunity** counts:

| Original upstream source count | Confirmed |
| --- | ---: |
| Target plant species in the frozen frame | **42,111** |
| Species with at least one source photo candidate | **42,095** |
| Species with no candidate photos | **16** |
| **Photo opportunity records** over these target species | **93,209,351** |
| Species with ≥40 source photo opportunities | **32,794** |
| Species with ≥100 source photo opportunities | **26,510** |
| Species with ≥80 source photo opportunities | **28,066** |
| Median photo opportunities per species | **212** |
| 99th percentile | **37,472** |
| Maximum recorded for a species | **301,529** |

This is approximately **927× the earlier 100,543-photo source union**, but the numerator is **provider metadata opportunity records** and the denominator is **earlier selected original photo IDs**. They are **not equivalent measures of independently measured floral colour**. Do not call the 93.2M unique licensed downloaded flower images or actual floral phenotypes: source photo-ID manifests and collision/rights checks remain incomplete. Some photos may lack visible/evaluable flowers; photo quality, colour classification and environmental joins must still be performed.

A separate [128-partition capacity audit #38819205800](https://github.com/zuizui0223/fcp/actions/runs/38819205800) validated planning the entire source into **128 source partitions**, with **1,456,459 nominal 64-photo compute units**. This is **capacity planning only**, not completed image retrieval, phenotype measurement or a valid worldwide environment model.

## Upstream processing status, not to be overstated

FCP PR #200 already contains original ROI-v4 detector + EfficientSAM photo-processing infrastructure, hash-bound original data acquisition, checkpointed distributed work queues and an independent limited original-pixel technical pilot. Its PR notes a 96-image technical run where **96/96** pixels were downloaded/processed, **39/96** four-state classifiable, **53/96** ROI/flip-gate unclassifiable and **4/96** ambiguous. This is a *small quarantined technical pilot*; it is **not** a representative estimate of worldwide photo classification success and not a full 42,111-species colour result. The first full metadata source manifest attempt failed on disk capacity, and the later counts-only source scan was successful. Full source photo-ID freeze and 93M image-pixel measurement are **not completed**.

The earlier Mois 100,543-photo original measured source and 45,953 classifiable hue results remain valid **for their own selected observational cohort**, as do the 16,801-soil36/13-variable original analysis. They cannot be re-labelled as the complete available iNaturalist photographic universe.

## Practical path for a single floral colour–climate paper

1. **Freeze the full 42,111-species source** and validate all original photo IDs, observation IDs, taxon mappings, photo rights and nonoverlap with previous samples. Avoid claiming 93.2M *unique* photo IDs until that manifest is complete.
2. **Design a feasible predeclared phenotype sampling stage** across the whole source frame (balanced genus/climate/geography and photo opportunity, with a frozen maximum per taxon for *measurement workload*, not silently selecting species on observed hue). A 93M-photo complete pixel run is vastly larger than the original 100,543-photo measured cohort. If all 93M are truly measured, use distributed checkpointed workers and verify their full results; otherwise disclose the prospective sampling design rather than imply an all-photo census.
3. **Run original ROI-v4 floral visibility/segmentation and four-colour measurement** on genuinely acquired, legally usable pixels, retaining explicit unclassifiable outcomes, selection probabilities and observer/source metadata. No unclassified photo is forced to be white or nonwhite.
4. **Join the selected/validated new colour outcomes to photo-site climate, altitude and soil** with missingness and source-support accounting. Only after an actual new measured-colour result exists, estimate pooled white/13-variable and genus-within/between associations with the same predeclared denominators and familywise uncertainty.
5. **For the existing single GEB draft**, keep the verified original 100,543/45,953 and 16,801-source results as the present paper. Describe the 93.2M metadata as **a much larger upstream available sampling frame, not 93.2M measured phenotypes or an ecological result**. Do not delay honest source-supported interpretation solely to claim complete measurement of all 93M photos; if a new broad cohort is measured, document it as a planned study extension or prespecified validation.

**Scientific distinction:** upstream metadata count, retrievable image count, classifiable flower-image count, classifiable plant species count, and climate/soil-complete model sample are five different denominators. All publication tables must explicitly label which one they use.
