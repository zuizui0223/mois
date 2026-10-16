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
