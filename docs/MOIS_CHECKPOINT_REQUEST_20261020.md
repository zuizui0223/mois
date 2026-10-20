# MOIS: verified completed-pixel checkpoint request (2026-10-20)

This source-backed milestone intentionally triggers the existing **data-only** checkpoint workflow on main. It reads **previously completed and SHA-verified original-photo work units**; it neither downloads raw iNaturalist images nor reruns flower-colour inference. The ongoing p0004, p0008 and p0016 acquisition jobs are left untouched.

## Source snapshot that justifies a new immutable checkpoint

At the source check around 06:20 UTC / 15:20 JST on 2026-10-20:

| Fixed original-colour acquisition run | Completed sealed RGB worker groups | Total planned photos |
| --- | ---: | ---: |
| p0004 `39353343403` | 16 of 18 | 18,119 |
| p0008 `39371199127` | 5 completed | 36,375 |
| p0016 `39374811084` | 1 completed | 73,872 |
| shadow `39369517433` | 1 completed, full cohort | 569 |

The preceding immutable checkpoint is [39376435276](https://github.com/zuizui0223/mois/releases/tag/mois-real-pixel-checkpoint-39376435276-v1), with **628,412 unique original-photo attempts** and **302,026 machine-classified four-hue outcomes**. More worker artifacts have arrived since that checkpoint, but their new photo counts/classifiable labels MUST be recomputed by the validated production collector, not projected from the scheduled source rows.

The checkpoint source registry pins these exact four run IDs, their original photo roster SHA256 and the immutable FCP/MOIS prior baseline. New and previous photo IDs must remain disjoint, missing units must stay *pending*, acquired-but-unclassifiable images cannot become coloured flowers, and legacy unknown labels are preserved exactly.

## Operational result criterion

A successful [checkpoint workflow](https://github.com/zuizui0223/mois/actions/workflows/mois-completed-pixel-checkpoint.yml) will publish `mois-real-pixel-checkpoint-<run-id>-v1` with `checkpoint_receipt.json`, both Parquet outputs and photo/hue counts. Only after the run is successful and the new receipt matches all file digests should its totals be described as confirmed.

This milestone path is an explicit source-backed way to request an additional checkpoint during active long-running image acquisition without modifying or restarting the measurement workers. It does not schedule repeated artificial snapshots.
