# MOIS third independently measured original photo wave and same original water inference

**Scientific objective:** Increase the number of original iNaturalist photos actually scored by the original ROI-v4 four-colour method, and then refit the **same** genus-equal, geographically adjusted white–rainfall and white–vapour observational model as the initial 100,543-photo MOIS study. Sharing original photograph sources with FCP is allowed; importing FCP's within-species D/W/50-km polymorphism *statistics* is not.

## Authentic expanding photographic sample

The old 100,543-photo source with 45,953 four-hue classifier assignments was enlarged to 248,312 original source photo IDs with 120,807 classified images by reusing independently completed FCP source-photo colour tables. The first newly computed 4,096 iNaturalist image sources were then physically downloaded, segmented and classified in 16 original photo-ROI jobs. All original SHA results were materialized into **252,408 original photo IDs / 122,383 four-hue classifications**.

The **second independent 4,096-photo source**, deliberately chosen from outside the original 106,376-photo precolour screen, is processing at [run 38995240689](https://github.com/zuizui0223/mois/actions/runs/38995240689). It is admitted to the complete 256,504-photo source and identical water analysis **only after its 16 authentic image-ROI shards have all succeeded**, as implemented in merged PR #28.

A **third further disjoint 4,096-photo source** has now been screened against all old measured photos, the entire first screened pool, and all second-wave photo IDs and selected from the full 3,430,502 original photo-ID data frame, with actual original imaging running in [the third-wave workflow](../.github/workflows/mois-original-third4096-real-ROI.yml). Third-wave IDs are not treated as measured or classifiable until their original ROI terminal records exist.

## No invented photograph outcomes and original 999-bootstrap method

`scripts/analysis/mois_append_real_third4096_photo_wave.py` accepts the SHA-verified **complete 256,504-photo first-two-wave dataset** plus the signed original third-wave 4,096-photo source and all 16 actual ROI-v4 terminal files. It requires exact photo-ID uniqueness across all three waves, source photo/species/observation/observer parity, unchanged fixed YOLO/EfficientSAM method hashes, genuine four-state hue or explicit unclassifiable outcome and complete 4,096-photo terminals. It refuses partial/failed new-photo shards.

The resulting source must have **260,600 unique original photo IDs**, and its actual final number of classifiable photographic hues is determined by image pixels rather than extrapolated from the earlier 38.5% success rate. It preserves all past 248,312 original photo IDs and already published observational MOIS moisture results.

Then the original method `mois_refit_original_water_with_actual_new_original_photos.py` is invoked on the actual enlarged original photos with `--number-of-new-4096-waves 3`. It sources genuine WorldClim BIO5/BIO12/BIO14/BIO15, solar and atmospheric vapour at the new original photo coordinates while retaining original prior photo/site values. The original **four hues × six preselected climate/water axes**, genus conditioning and genus-equal weighting, 999 shared-genus bootstrap replicates and 24-contrast maxT simultaneous 95% familywise intervals are unchanged. The same-genus coarse geographic-cell sensitivity remains an essential limitation.

`.github/workflows/mois-third-original-real4096-auto-water.yml` can run at most once per source completion (source-run and source-artifact gates). Its first CI/pull-request stage never counts photos before all 16 ROI terminals are available. Once the real photo wave and earlier two-wave union exist, it assembles and refits the original moisture model and publishes a source SHA-bearing finished analysis. A bounded schedule is used to avoid losing the dependency if GitHub completes one source before the other; this does not assume that either source has already finished.

**Causal constraints:** WorldClim annual averages are not observation-year local weather; four-color photograph classifiers are not proof of field-observed petal colour, genotypes, heritability or relative fitness. Photo reuse across FCP and MOIS is not independent-provider replication.
