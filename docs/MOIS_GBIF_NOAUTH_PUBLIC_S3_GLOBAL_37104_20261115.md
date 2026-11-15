# MOIS: GBIF認証不要の全球37,104種×等面積格子データ取得（2026-11-15）

## 主結論

**GBIF SQL Download APIを使うには認証が必要だが、GBIF全出現記録の月次公開スナップショットはAWS S3 / GCSから認証なしで取得できる。** 前の回答で「認証が唯一の障害」と述べたのは誤り。これは公式GBIFによるPublic Cloud Datasetの利用で、GBIFユーザー名・AWSアカウント・アクセスキーが要らない。

出典: [GBIF cloud services](https://techdocs.gbif.org/en/cloud-services/)／[AWS public Parquet data](https://github.com/gbif/occurrence/blob/dev/aws-public-data.md)。

## 実確認したスナップショット

- AWS公開bucket: gbif-open-data-us-east-1
- Snapshot: occurrence/2026-11-01/occurrence.parquet/
- 完全Parquet object数: **10,346**
- source総圧縮バイト: **271,777,586,222 bytes = 271.778 GB**
- source citation.txt: **GBIF.org (1 November 2026) GBIF Occurrence Download https://doi.org/10.15468/dl.jeesns**
- 参照source author citationのSHA256: d9770dd7ef38ee8fb0e461f460e5c26abed60ed4a13319747b7b7276f8c20fc0
- [最初の真正Parquet読込 Actions #43721976307](https://github.com/zuizui0223/mois/actions/runs/43721976307) SUCCESS。409,863記録（うち植物94,741）、COL XR accepted species key一致4,512、specieskey型string。匿名でParquet実行できることを確認。
- [全10,346 objectの匿名インベントリ Actions #43722246618](https://github.com/zuizui0223/mois/actions/runs/43722246618) SUCCESS。11ページのS3一覧を完全読込し、16個の非重複partitionへファイルサイズを均等割当。artifact #14164808564 SHA256 473c9eff45abba0cafa589494973c48af918a0251fcc3d7fb857d8f0b090cc4c。

## 本物のGBIF分布座標の匿名取得と格子化

[少数ファイルを真正に処理した Actions #43722722347](https://github.com/zuizui0223/mois/actions/runs/43722722347) SUCCESS。artifact #14165001819 SHA256 e52ca6c9673b48a5fa80b6c535b391bcc416373122a94541db24c71e6c2bce27。

| 6ファイルの実行 | count |
|---|---:|
| Parquetの全GBIF行 | 2,641,098 |
| 37,104 COL XR受容種キー一致かつPlantae/PRESENT | 90,809 |
| source specieskeyに対して学名binomialが異なる行 | 2 |
| 生物地理品質フラグ・座標の追加QA除外 | 349 |
| **利用できた位置情報付きGBIF行** | **90,458** |
| **species×100 km等面積セル（EPSG:6933）** | **69,881** |
| species×country code pairs | 26,457 |

この6ファイルは分析パイプラインの検証専用。全球coverage・37,104種すべての有効GBIF分布は**まだ**完成していない。AWS/GBIFの認証値を空にして実行している。

## 全球全件：16分散処理

- 起動した分散ワークフロー：[run #43723276949](https://github.com/zuizui0223/mois/actions/runs/43723276949)
- Workflow: .github/workflows/mois-gbif-anon-full-37104-species-16shard-grid-20261115.yml
- 計16独立jobのうち、最初に8jobが約16.99GBのS3 objectsを各担当（残り8jobは後続）。対象のParquetは必要な列のみリモート読込、元271.8GBを一時的に全ダウンロードしない。
- 受容種はCOL XR 37,104の完全ID台帳に由来し、specieskeyと学名のbinomialを二重照合。kingdom Plantae／PRESENT、緯度経度の数値範囲・0/0除外、4つの重大geospatial issue排除、不確実性50km以下（欠損は別扱い）でフィルター。
- 座標は **EPSG:6933の100 km×100 km等面積の矩形グリッド**へ変換。以前のGBIF SQLで提案した **ISEA3H resolution8の六角形グリッドとは別の格子定義**であり互換と扱わない。
- 16jobがすべて成功した後のみ、manifestの10,346個のpath, bytes, ETagと実行ログをexactly-once照合。合致した場合にだけ全球37,104種×格子の表と受容種×9色品質表を作る。失敗・未取得を分布不在に補完しない。
- 新たな個別GBIF download DOIは発行しないが、上記公開GBIFスナップショットDOI・citation.txtと独自抽出原本のSHAを保持する。

### 生態学上の解釈の境界

- 写真9色は花被の真値保証ではなく、以前の5写真mask QAでは葉・背景の誤検出を実際に確認。
- 種レベルの撮影色はGBIF出現地点における局所花色ではない。GBIF出現地点は在来／外来区分も含まない。
- まずGBIF species×grid presence、photo hue evidence status、GBIFの観測努力と分類群の系統構成を区別してから、地理・気候・土壌・島嶼環境と関連付ける。

### 提供可能な出典と実行コード

1. scripts/analysis/mois_gbif_anonymous_snapshot_probe_20261115.py
2. scripts/analysis/mois_gbif_anonymous_snapshot_inventory.py
3. scripts/analysis/mois_gbif_anon_species_grid_scan_20261115.py
4. scripts/analysis/mois_gbif_anon_merge_16shards_20261115.py
5. .github/workflows/mois-gbif-anon-full-37104-species-16shard-grid-20261115.yml

従来の149分割SQLや認証用PR #201は**必須の入口ではない**。公開スナップショットの実行が失敗した場合の代替・別検証として保存する。
