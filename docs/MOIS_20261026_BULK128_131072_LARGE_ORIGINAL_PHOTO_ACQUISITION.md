# MOIS: one-shot 128-shard large original-photo acquisition (2026-10-26)

## Executed source-selection result

**Real signed source-only work completed:** [GitHub Actions source run 40338267506](https://github.com/zuizui0223/mois/actions/runs/40338267506) / [frozen original 131,072-source-photo ID roster](https://github.com/zuizui0223/mois/releases/tag/mois-full93m-bulk128-tier10-131072-original-photo-source-v1). All **128 immutable original taxon-modulo source partitions** were examined using exactly **one of the 14 observation-hash tiers** (tier10), not an additional scan of all 93,209,351 source-photo IDs.

Actual pre-colour counts:

- **9,543,472** original source-photo metadata records read once; all 16 tier10 archives and all 128 part gzip checksums matched the previous full-93m SHA registry.
- **131,072 actual previously unattempted original photo-ID source opportunities selected**, precisely **1,024 per each of 128 large compute shards**, each internally feeding four 256-photo warm FCP ROI-v4 work units (512 small photo units total).
- Source frame: 42,111 nominal iNaturalist plant taxa; **32,768** source-reference photo targets and **98,304** adaptive gap-targeted photo targets, provenance and probabilities separated (NOT a combined unbiased global-flora estimator).
- The source priority receipt proves **zero P5-targeted photos**. Within the combined 131,072 source roster, priority labels are **P1=520, P2=39,624, P3=61,902, P4=2,301**; remaining source-reference shots are largely P5. The combined P1–P4 counts include reference entries and are not all targeted.
- Latitude metadata sources include **19,453 northern tropical** and **13,004 southern tropical** photos (also including reference entries). Both hemispheres' recorded local winter, spring, summer and autumn are preserved.
- Source acquisition deduped against **1,043,930** existing original image attempts, **2,048** subsequent true ROI-v4 measurements, **64** method-distinct original 500px-Lab pilot photos, and **4,096** P2-first separately reserved source photographs, using original photo IDs and observation UUIDs. No already measured or reserved source photo IDs were newly selected.
- Old 93m opportunity photos remain immutable and available for other studies; source temporal/geographic taxon sampling weights cannot be interpreted as plant abundance or verified native/wild petals.

## Efficiency and actual worker execution

- Source metadata cost: one scan of the existing tier10 group only (16 archives, approximately 411 MB before subsetting), not repeated 93m full-index scans.
- [Signed source bundle](https://github.com/zuizui0223/mois/releases/tag/mois-full93m-bulk128-tier10-131072-original-photo-source-v1) is **19,163,860 bytes compressed**; [verified 16-group repack release](https://github.com/zuizui0223/mois/releases/tag/mois-bulk128-tier10-source16-efficient-shards-v1) contains **16 small approximately 1.2 MB source packets**. Any one photo worker retrieves only its own 8-shard source bundle, NOT 19MB of unrelated source photo IDs.
- Every **1,024-photo large worker** uses the previous source-pinned YOLO11 + EfficientSAM ROI-v4 and reuses its model process for four consecutive <=256-photo units. Bounded original S3 downloads are prefetched. Matrix has **128 requested photo shards, at most 16 concurrent workers**, no synthetic colour values, no change to the frozen original ROI thresholds.
- [Actual 131,072-photo acquisition Action](https://github.com/zuizui0223/mois/actions/runs/40340338433) launched with all 128 source shards and initial 16 workers. A queued worker is NOT a completed photo download; genuine original-image success and 9-component ROI status are recorded only in signed completed numeric terminals.
- Each completed numeric shard is saved as a workflow artifact and a SHA-verified immutable Release (raw licensed photo files are never redistributed). A completed shard is checked and reused on retry, avoiding repeat image inference where a sealed receipt exists.
- [Source-verified collector source](../scripts/analysis/mois_bulk128_collect_original_roi9.py) retains strict ROI9, original 4-hue labels and unclassified/missing cases. The separate [data-only collector](../.github/workflows/mois-bulk128-data-only-collector.yml) can aggregate fully completed shards or produce a count-only partial receipt without opening raw original images.

## Scientific and operational limits

1. **131,072 is currently a verified newly selected source-photo work queue, not yet the number of fully measured flower colours.** A real photo only counts after its 4-unit worker terminal has completed and is independently SHA-audited.
2. **Gap-targeted selection depends on previously measured classification outcomes**, so its observed colour fraction cannot be treated as a population probability sample. Reference rows are separated and report conditional source probabilities; the technical worker probability on adaptive rows omits the unknown second-stage inclusion probability.
3. Frozen original flower ROI-v4 frequently returns unclassifiable photos. All technical unknowns stay unknown. This bulk-acquisition action does not certify floral organ identity, wild/native provenance, true petal colour spectra, pollen interactions or individual fitness.
4. 128-way source sharding is for data-parallel throughput. The rate of successfully classified genuinely unrepresented taxa, not the number of queued photographs, determines ecological coverage improvement. Source-priority and region/season audits must be rerun on the completed actual pixels before claiming a better ecological sampling frame.
5. Large-scale GitHub-hosted processing uses substantial concurrent runner minutes, network bandwidth and provider requests; **max-parallel=16** is an intentional cap, not 128 simultaneous S3 download processes. This respects memory/CPU limits and avoids repeated model bootstrap per small unit.

**Current task scope:** source queue and optimized worker orchestration have been executed; final actual-photo ROI9 measured count and exact new classifiable taxa remain dependent on the running workers, not presumed successful.
