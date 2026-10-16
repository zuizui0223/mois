# 全元データから「観察日あり＋実際に開花している写真」を抽出する

**研究方針（2026-10-16）：42,111種・93,209,351枚の元写真候補をすべて対象にする。** 旧100,543枚の測定済み標本、1種40観察の上限、旧来の写真数上限はいずれも本解析の母集団に使わない。

季節解析に入れる花色写真の採用条件は、(1) 元の2026-09-27 iNaturalist観察表に有効な観察日があり、(2) その写真自体に開花中の花が写り、元ROI-v4画像測定により4状態の花色を分類できること。この二つを混同せず、元データから順番に確認する。

## 1. 元データの列と開花状態の制限

元のiNaturalist Open Datasetには taxa・observations・observers・photos の4表がある。observations表の **observed_on と observation_uuid** を、photos表の **observation_uuid と photo_id** に結合する。既存のFCP全件原本インポータがこの照合を実装しており、種別の元写真マニフェストには観察日が含まれる。2026年10月時点のAPIで観察日を推測して埋めない。投稿日時・写真アップロード日時も代用しない。

重要：この4表には **各写真が本当に開花中の花を写していることを保証する注釈がない**。現在のiNaturalistのFlowers注釈は観察記録単位であり、複数写真のすべてが花を写していることを意味しない。また注釈がないことは非開花を意味しない。したがってFlowers注釈だけで元の写真候補を先に削除しない。

## 2. 観察日が判明している元写真を上限なしで列挙する

[全元写真の原本観察日監査コード](../scripts/acquisition/mois_full_source_dated_flower_candidates_20261016.py) は、FCPの凍結済み原本source shardを1区画ずつ走査する。入力はSHA256で確定した元種別photo_manifest.ndjsonと、独立監査済み42,111種の全写真数CSV。元のYYYY-MM-DDの日付が実在し、原本スナップショット日（2026-09-27）以前なら採用する。

- 種ごとの抽出上限なし：日付のある全元写真IDを候補とする。
- 元写真全件数・観察日有効件数・日付不明／形式不正／スナップショット後の日付件数を区別して記録。
- 元の観察UUID、写真ID、撮影地点、観察者ID、URLは保護されたローカルの候補キューのみへ保存。
- 北半球と南半球の季節を逆転して整理。熱帯は暦月を維持し、雨季・乾季を決めつけない。
- この段階の写真は **開花未確認**。公表結果に白花率や開花写真数を捏造しない。
- 原本の128 source shardsすべての確定前に、1 shardの集計を全球の完成値と呼ばない。

実行例：

    python -m scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 \
      --sealed-original-fcp-shard "$FCP_SOURCE_SHARD" \
      --independent-42111-photo-counts "$FCP_COUNT_CSV" \
      --out-private-dated-photo-queue "$MOIS_PRIVATE_ROOT/dated_candidates.ndjson" \
      --out-public-source-receipt "$MOIS_PRIVATE_ROOT/dated_candidates_receipt.json"

## 3. 実際に咲いている写真だけを花色集計する

[FCP全件画像測定済み原本から開花・日付両方を確認するコード](../scripts/analysis/mois_original_full_fcp_dated_flower_ROI_20261016.py) は、元のFCP ROI-v4/YOLO/EfficientSAMの正規測定結果を使用する。花の見える写真で measurement_status が classified_four_state_morph、かつ元観察日が有効という **二条件を満たす写真だけ** が季節別白花率の分母になる。白花の分子には morph=white を追加要求する。

それ以外の写真は偽の有色花・白花に割り当てず、日付不明数・花色評価不能数として残す。元FCPの完全な実画像計測証明・オリジナル測定方法SHA・全写真結果SHAを要求し、旧10万枚の既測定結果を新しい9,320万件の花色結果に混ぜない。

**月別に2種類の結果を保存する。** まず、元の分類可能な写真のうち白花の割合。次に、同じ月の各種の白花率を計算してから種を等重みで平均した白花率。後者が開花種の季節的な入れ替わりを表す主解析である。いずれも「写真で確認できた開花花色」であり、野外に存在する全開花個体・全開花種の無作為母集団ではない。

## 4. 月別気候と比較する

[真正な48枚のWorldClim月別ラスタ原本](../results/mois_worldclim21_dynamic_48_real_source_receipt_20261016.json) は気温・日射・降水量・水蒸気圧の12か月平年値である。後続の[地理セル固定効果モデル](../scripts/analysis/mois_monthly_white_species_equal_cell_FE_20261016.py)では、同一地域の月別花色と各月の気候を照合する。緯度・経度、標高、土壌、地形の比較的安定したセル間の違いはセル固定効果で調整する。これらの独立した地域間係数は別モデルで評価する。

1970–2000年の月別平年値は、観察年の実測天候ではない。花色ROI分類成功率、日付欠測、開花写真の観察努力、種・観察者・地点の反復を別々に示す。

## 5. 今の実績と未完了工程

独立監査済みの元写真候補は **93,209,351枚**。これは写真IDの機会数であり、開花確認写真数ではない。大容量の元36GBアーカイブ、永続SQLite、128区画の全写真ID確定と元画像大量ROI測定は未完了。以前のGitHub-hosted実行はディスク容量／実行上限により停止した。全件を実際に処理するには永続ディスクのある適切な計算基盤と、写真提供者の利用条件に沿った取得・レート管理が必要。

現時点で実際の「観察日判明写真の全件数」「開花画像として分類可能な全件数」「月別の全球白花率」「白花率×日射・気温・水分の全球係数」は**未算出**。追加した2つのソース監査・花色採用コードは合成データでテスト済みで、未測定の画像を開花扱いしない。

[FCP元全写真IDの永続化手順](https://github.com/zuizui0223/fcp/blob/analysis/fcp-fresh-uncapped-zero-20261015/docs/FCP_REAL_42111_PHOTO_ID_FREEZE_HANDOFF_20261016.md) を参照する。


## 6. 128 original-source shards to actual monthly T/S/P/VPD seasonal white inference (implemented)

The source workflow is now runnable in distinct, SHA-verified durable phases. Its workers MUST have access to the frozen original provider archive and the completed FCP durable source; no GitHub-hosted short CI run has been claimed to complete these 93-million-source stages.

### Recover original observation dates in bounded, idempotent source units

[Safe original full-source shard range runner](../scripts/acquisition/mois_run_resumable_full_original_date_shard_range_20261016.py) takes at most four consecutive FCP source shards **per COMPUTE job**. Four is NOT a photo/taxon sampling limit; it simply bounds one resumable work unit. Previously completed source date queues are independently SHA-verified and skipped, and partial/unsealed outputs are rejected for inspection instead of overwritten.

    python -m scripts.acquisition.mois_run_resumable_full_original_date_shard_range_20261016 \
      --fcp-source-shards-root "$FCP_DURABLE_ROOT/complete_42111_metadata_source" \
      --original-42111-count-ledger "$FCP_ORIGINAL_COUNT_CSV" \
      --private-dated-photos-root "$MOIS_PRIVATE_ROOT/dated_original_42111" \
      --start-shard 0 --stop-shard-exclusive 4 \
      --out-progress-receipt "$MOIS_PRIVATE_ROOT/progress_000_004.json"

Run remaining disjoint ranges 4:8, 8:12, etc through 124:128, with unique progress receipt filenames, only AFTER their original FCP FCP-SHA-sealed source shards exist. Each result is specifically labelled **PARTIAL DATE-KNOWN PLANT PHOTO CANDIDATES, ZERO FLOWERING CONFIRMED**.

Once and only once all 128 original dated source queues are sealed, execute the independent [128-shard full metadata observed-on reconciliation](../scripts/acquisition/mois_aggregate_all_original_dated_shards_20261016.py):

    python -m scripts.acquisition.mois_aggregate_all_original_dated_shards_20261016 \
      --fcp-sealed-shards-root "$FCP_DURABLE_ROOT/complete_42111_metadata_source" \
      --mois-private-dated-shards-root "$MOIS_PRIVATE_ROOT/dated_original_42111" \
      --original-42111-count-csv "$FCP_ORIGINAL_COUNT_CSV" \
      --out-public-all-original-dated-photo-receipt "$MOIS_PRIVATE_ROOT/all_original_dates_audit.json"

This independently checks each sealed shard source and date-queue SHA, disjoint source taxon membership, all 42,111 original taxon identities, exactly 93,209,351 total photo opportunities, snapshot-date validity/exclusion conservation and geographically localized original-photo month counts. Missing one shard/one taxon/one photo opportunity is not silently called complete. The output is **original observed_on candidate incidence only**, not actual flowers, blooming photographs or measured white flowers. Source protected IDs, observation UUIDs, users and coordinates remain outside public reports.

### Real flower-positive photos: observation-first, species-equal, then climate comparison

The completed *original FCP whole-photo* fresh ROI-v4 reassembly supplies an exact date-verified, **visibly classifiable flower photo** manifest through [positive bloom and original date gate](../scripts/analysis/mois_original_full_fcp_dated_flower_ROI_20261016.py). This gate requires the original ROI-v4 source method fingerprints; it does NOT reuse old 100,543 measured photos or infer absent Flowers annotations as nonblooming. A classifiable flower photograph is evidence of visible floral hue; precise anther/stigma-stage anthesis and open-bud phenology are **not independently classified**, so avoid claims about exact flowering developmental state.

The newly implemented [full original bloom-photo → 48 monthly WorldClim → seasonal species-equal FE analysis](../scripts/analysis/mois_full_bloom_monthly_TSPVPD_observation_species_FE_20261016.py) now directly accepts the SHA-sealed complete photographed bloom positives (private NDJSON), the receipt from the complete 93M-photo ROI run, and the [real source-verified 48-band WorldClim rasters](../results/mois_worldclim21_dynamic_48_real_source_receipt_20261016.json). It validates EVERY positive bloom's original date, taxon, observer, month, hemisphere and provider-attested ROI hue; the private output keeps only taxon × month × climate grid-cell aggregates. Temperature (°C), solar radiation (kJ m-2 day-1), precipitation (mm/month) and derived average-month VPD (kPa) are extracted from the photo's *original* location and observed month. Source 1970–2000 normals are not observations from the photographed year.

    python -m scripts.analysis.mois_full_bloom_monthly_TSPVPD_observation_species_FE_20261016 \
      --original-full-date-and-bloom-positive-photo-ndjson "$MOIS_PRIVATE_ROOT/actual_dated_blooms.ndjson" \
      --original-full-positive-bloom-receipt "$MOIS_PRIVATE_ROOT/global_month_white_summary.json" \
      --verified-original-worldclim48-root "$MOIS_CLIMATE_REAL_ROOT" \
      --verified-original-worldclim48-manifest-sha256 9d0cf37a64534542d7788debef28c13f106a15c3bf7bf91be0653f3ce106baae \
      --out-private-observation-species-site-month-climate "$MOIS_PRIVATE_ROOT/original_species_cell_month.csv" \
      --out-public-species-equal-month-white-climate "$MOIS_PRIVATE_ROOT/season_white_TSPVPD_site_FE.json"

**Distinct photographs are NOT necessarily distinct observations:** the full original FCP source can have several original pictures belonging to the same iNaturalist observation_uuid. All photo information is preserved, but the source UUID is first collapsed to an **equal-weight observation-level white fraction**, then averaged among original named taxa within source 10-minute climate cell × calendar month; otherwise a 20-photo single observer event falsely counts as 20 independent flowering observations. For the worldwide zone × observed-month report, a named taxon has one vote after aggregating its independent observed photographs. Standard photo-count-weighted white is retained as a clearly separate descriptive denominator. These observed flowering-community means deliberately preserve species turnover; the same-genus/species-only analysis remains secondary.

Three exploratory conditional models are calculated only when supported: exact geographic-cell fixed effects; cell + common month fixed effects; and cell + Northern/Southern/tropical climate-zone × month effects. The original stable latitude/longitude/soil/altitude/terrain terms at the climate cell scale are absorbed, not separately identified. Full-rank 4-axis T/S/P/VPD support and >=30 independent repeated-month climate cells gate preliminary one-way geographic-cell cluster intervals. These do NOT fix observation detectability, shared genus/observer correlation, source photo ROI nonresponse, or identify evolutionary floral fitness.

### Source status and scientific claims

The genuine provider-level 42,111-taxon 93,209,351-photo **candidate counts** are independently audited. The full original date-known candidate extraction, open-flower/colour ROI measurement across all date-qualified photographs, complete original positive-photo-to-monthly-climate join and worldwide adjusted seasonal T/S/P/VPD flower-white model remain **UNEXECUTED on the entire real source**: FCP's 36GB original archive, permanent SQLite and all 128 source shard identities still need source-host durable execution. These scripts passed synthetic full-chain/duplicate-observation/partial-shard tests; no real worldwide seasonal white rate or climate coefficient has been manufactured.


## 7. Date-first photo image processing: all valid-day original pictures, no undated-image waste

**Optimization confirmed on 16 October 2026:** The main Mois seasonal study does not require reprocessing the pictures which lack the original snapshot observed_on day at all. It *does* require processing **EVERY date-valid original photo**, not just the old 100,543 pictures, any 40-photos/taxon sample, or only positively annotated iNaturalist Flowers observations. The independently attested 93,209,351 source PHOTO opportunity frame and the ORIGINAL DATE EXCLUSION denominators remain intact.

**The DATE-FIRST pipeline is implemented as its own separate whole-source Mois scientific scope.** It is not a false claim that FCP's separate all-93M-photo polymorphism study has completed: all original FCP source-ID shards are still required, but the subsequent pixel ROI work is restricted by the user's explicit original-date eligibility condition.

After one of the complete, source-SHA-pinned original FCP date-shard queues has been generated:

1. [Full date-known photo → location-blind ROI packet dispatch](../scripts/acquisition/mois_prepare_all_dated_original_roi_packets_20261016.py) consumes the COMPLETE per-shard private date-known photo queue, validates all source UUID/taxon/observer/calendar/site/photo URL and SHA ledger, then emits protected worker/acquisition/private join files in **at most 4,096-photo compute packets**. The source photo roster is not sampled, and all original photo IDs with usable observed_on are dispatched. A **4,097-photo single-species/observation synthetic regression** verifies the source overflows into the next compute packet instead of dropping picture 4,097. The worker sees only anonymized method-unique IDs, image filename, image licence and image bytes, **never genus, source date, latitude/longitude or historical flower-colour results**.

2. [Original frozen FCP ROI-v4 work unit](../scripts/acquisition/mois_execute_all_dated_original_roi_workunit_20261016.py) runs **64 real original pictures maximum per machine batch** from those packets, enforcing the original source-measurement Git blob hash, YOLO and EfficientSAM file hashes, photo-source licensing, and host-global image-provider pacing. This is a bounded *compute unit* and never a seasonal photo cap. Temporary original image pixels are discarded. Transient provider errors abort the whole unit; explicit permanent source-photo acquisition failures are recorded **unclassifiable**, never as absent/coloured/white flowers.

3. [Source-hash sealed per-shard photo-result reassembly](../scripts/analysis/mois_reassemble_dated_original_roi_shard_20261016.py) requires **EVERY original date-valid photo** to have an attested ROI terminal (including documented failures) before the shard may be marked completed. It joins private original photo ID, observed_on, original observer/taxon/site and actual source four-hue morphology only AFTER the location-blind image evaluator. Only positive photo ROI-v4 four-hue classifications are exported to protected actual-bloom NDJSON. The public result separately counts date-valid image opportunities, image acquisition failure, ROI missing/ambiguous, classified flower colour and white. An absent work unit or an altered photo/method/receipt hard-fails without global success.

4. [All-128-shard **global DATE-FIRST blooming-photo census**](../scripts/analysis/mois_complete_original_date_first_blooming_photo_census_20261016.py) requires the complete independently verified original metadata date census, PLUS every separately sealed source-shard photo ROI result. It reconciles the exact full 42,111-taxon and 93,209,351-photo source denominator, the number with an original snapshot date, the exact number of date-known pictures on which image ROI was attempted, classified positive floral-hue photos, their white subset, distinct original observation UUIDs, taxa and genera, and photo-classification missingness. Original duplicate photo IDs across shards are prohibited. Species-equal white fractions average *observation-level* white composition first, not raw photo repetitions. **Even ONE incomplete date-shard or image-work-unit stops the worldwide claim.**

5. [Full photo-positive original month × 48-source-layer climate and spatial controls](../scripts/analysis/mois_full_bloom_monthly_TSPVPD_observation_species_FE_20261016.py) accepts this distinct **DATE-FIRST full 128-shard positive-flower census** OR FCP's separately completed full-photo remeasurement census, with different explicit metadata/image-work denominators. The prior version demanded FCP image measurement of **every single** source photo, even without original observation day; now the Mois seasonal method can skip those uninformative photo pixels while retaining their counts. It extracts actual WorldClim original-site climatological month temperature, solar radiation, precipitation and mean-month-derived VPD; multi-photo observations collapse before species-equal cell/month community white statistics. The three exploratory models still control exact climate cell, cell+month, and cell+latitude-zone/month, subject to real repeated-classified-site and collinearity/cluster uncertainty gates.

**Protected host execution skeleton**, for each 128-source-shard whose frozen original dated queue is already complete (using one packet number and image work-batch example here):

    python -m scripts.acquisition.mois_prepare_all_dated_original_roi_packets_20261016 \
      --original-private-dated-photo-queue "$MOIS_PRIVATE_ROOT/dated_original_42111/dated_candidates_shard_000.ndjson" \
      --original-source-date-receipt "$MOIS_PRIVATE_ROOT/dated_original_42111/dated_candidates_shard_000_receipt.json" \
      --out-private-roi-packets "$MOIS_PRIVATE_ROOT/dated_roi_shards/shard_000_of_128"

    # Original provider rights and global rate limits MUST be reviewed before downloading.
    # FCP_PROVIDER_GLOBAL_RATE_FILE must be configured for coordinated permitted access.
    python -m scripts.acquisition.mois_execute_all_dated_original_roi_workunit_20261016 \
      --private-original-photo-roi-packet-dir "$MOIS_PRIVATE_ROOT/dated_roi_shards/shard_000_of_128/packet_000000" \
      --private-original-date-dispatch-root "$MOIS_PRIVATE_ROOT/dated_roi_shards/shard_000_of_128" \
      --private-original-roi-terminals-root "$MOIS_PRIVATE_ROOT/sealed_roi_shard_000" \
      --packet-number 0 --batch-number 0 \
      --fcp-source-repository "$FCP_SOURCE_REPO" \
      --original-fcp-roi-measure-script "$ORIGINAL_FCP_ROI_SCRIPT" \
      --original-fcp-yolo-detector "$ORIGINAL_FCP_DETECTOR" \
      --original-fcp-efficient-sam "$ORIGINAL_FCP_EFFICIENT_SAM" \
      --provider-access-reviewed

    # Run EVERY bounded work unit, then and only then:
    python -m scripts.analysis.mois_reassemble_dated_original_roi_shard_20261016 \
      --protected-date-first-packet-root "$MOIS_PRIVATE_ROOT/dated_roi_shards/shard_000_of_128" \
      --protected-sealed-original-roi-terminal-root "$MOIS_PRIVATE_ROOT/sealed_roi_shard_000" \
      --out-private-original-actual-bloom-photo-ndjson "$MOIS_PRIVATE_ROOT/roi_positive_shards/shard_000_of_128/positive_actual_bloom_photos.ndjson" \
      --out-public-original-date-roi-shard-receipt "$MOIS_PRIVATE_ROOT/roi_positive_shards/shard_000_of_128/roi_shard_receipt.json"

After all 128 such ROI shards **and** the independently sealed all-128 original-date metadata census:

    python -m scripts.analysis.mois_complete_original_date_first_blooming_photo_census_20261016 \
      --all-original-date-census-receipt "$MOIS_PRIVATE_ROOT/all_original_dates_audit.json" \
      --all-128-complete-original-roi-shard-root "$MOIS_PRIVATE_ROOT/roi_positive_shards" \
      --out-private-global-photo-bloom-records "$MOIS_PRIVATE_ROOT/all_original_date_verified_bloom_photos.ndjson" \
      --out-public-date-verified-global-photo-white "$MOIS_PRIVATE_ROOT/all_source_date_first_month_white.json"

That final private positive photo file then feeds the climate model documented in section 6, using the **date-first completed global ROI receipt** instead of FCP's separate full-image ROI remeasurement receipt. The pipeline **cannot yet be executed worldwide** because the real 128 original source-ID manifests are not durably completed. Synthetic tests of complete source-date → blinded photo ROI → shard/global photo outcomes → month WorldClim, including multiple original photos per observation and a missing packet/shard, have passed; they are *engineering tests only*, not new world phenological/climatic results.


## 8. Early original observed-on audit from genuine permanent FCP SQLite

**New, shorter route to the first true 93-million PHOTO observation-date result.** Even before 128 source-shard photo-ID files have been emitted, FCP's byte-verified, globally deduplicated permanent original source SQLite already contains the full licensed FCP photographic table and original observation UUID/date/latitude rows. The new [full-original-source SQLite observed-on opportunity auditor](../scripts/acquisition/mois_real_fcp_staged_sqlite_original_observed_on_counts_20261016.py) directly joins those tables on the original observation UUID, counts every nonambiguous photo and classifies the immutable original 2026-09-27 observed-on using a strict Gregorian day parser. It independently demands the original source SQLite SHA256, full 42,111 original taxa, exact 93,209,351 unique photo opportunities and photo-count ledger SHA256. It gives date-known / invalid / absent original photo opportunity counts, original observed year and month×hemisphere/source-photographic taxon support, without retaining personal photo IDs in the aggregate and **without opening any original photo pixels**. It is NOT yet a blooming-flower count or white-flower estimate.

Run only on an operator-provisioned persistent FCP stage that truly passed the original full 36-GB provider metadata ingest. For example, from the Mois repository root:

    python -m scripts.acquisition.mois_real_fcp_staged_sqlite_original_observed_on_counts_20261016 \
      --fcp-complete-original-staged-sqlite-root "$FCP_DURABLE_ROOT/source_20260927/complete_42111_metadata_source" \
      --independent-42111-original-photo-count-csv "$FCP_DURABLE_ROOT/references/all_42111_species_photo_opportunity_counts.csv" \
      --out-public-original-observed-on-photo-opportunities "$FCP_DURABLE_ROOT/analysis/mois_full_source_observed_on_months.json"

FCP's manual source stage requires a persistent self-hosted runner, legitimate original S3 archive/permission, byte-pinned metadata and large durable storage. The archived FCP **photo-count-only** artifact does NOT contain the entire photo/date table or stage SQLite. As of this writing the 93M original SQLite has not been made available through GitHub, so the actual date-known count is **NOT YET COMPUTED**; the new code passed source-schema synthetic CI, not production global data. It is nonetheless an earlier result gate than all 128 photo-ID manifests or pixel processing.

## 9. Actual original FCP ROI-v4 code compatibility (not just mock tests)

[Actual independent original-FCP-backend integration Actions #38881206709](https://github.com/zuizui0223/mois/actions/runs/38881206709) succeeded. It checked out the immutable FCP PR200 source commit c00e3db156ec2f214b13e83f3e91dd8d863d83dc and byte-verified its true ROI worker Git blob 45de671444c0fac4fb4611375d19ed7e5a1608f3. The Mois adapter imported the *real* source and exercised its original nine-colour/photo-hue terminal validator, function ABI, permanent versus transient failure classes and original S3 URI rejection, not a fake injected backend. The original frozen scientific measurement code commit is 9fae6ccdf684a46026f72ba12e98de2c5c54bf2a. This proves production source **interface compatibility** only; original model weights, permissions and actual provider photograph pixels were not opened in this integration CI. Thus the real-world full image measurement remains uncompleted.
