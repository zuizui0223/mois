# MOIS — GBIF 37,104受容種を1本のSQLにまとめる実験（2026-11-15）

## 結論：公式バリデーターで37,104種を一括受理

[GitHub Actions #43712422647](https://github.com/zuizui0223/mois/actions/runs/43712422647) がSUCCESS、GBIF公式公開SQLバリデーターが **HTTP 201** を返し、1本のSQLに指定した **37,104のCOL XR受容種キーを全部そのまま返した**。欠落0、重複0。

| 対象・技術検査 | 実数 |
|---|---:|
| 厳密に名寄せしたCOL XR Plantae species | 37,104 |
| 1本のSQLで指定したspeciesKey | 37,104 |
| バリデーターの応答で保持されたspeciesKey | 37,104 |
| SQLテキストの長さ | 328,247 bytes |
| GBIF POST本文JSONの長さ | 328,309 bytes |
| APIバリデーター | **HTTP 201、受理** |
| 実際の認証付きSQL送信 | **0** |
| GBIF downloadKey/DOI/種×格子データ | **未取得** |

GBIF SQLは、georeferenced/PRESENT/非geospatial-issueという基本品質フィルタで、COL XR specieskey、ISEA3HCode解像度8（等面積セル）、記録件数をGROUP BYする。37,104個すべてのspecieskeyをIN句で限定し、オフターゲット種が混入しない設計。

**従来の250種×149 SQL分割は、構文／POSTサイズのために必要ではなかった。** いまは一括SQL 1本が本命。サーバー実行負荷・出力サイズが原因で実ジョブが失敗する場合のみ、149分割をフォールバックとして使う。

### 重要な境界

SQL **validation** が通ることは、GBIF非同期実行がSUCCEEDEDになる保証ではない。実際の認証付きGBIF request、ZIP、DOIはまだ存在しない。37,104種の全球分布を「取得済み」と主張しない。

## 固定sourceと再現情報

- 元COL XR受容名artifact #14157705123、ZIP SHA256 abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a、内部CSV SHA256 5a77f3c968cf221287cd98a855559b3c3e8a8cf1322c08add28878adeb293df8。
- 公式validator [output artifact #14161421174](https://github.com/zuizui0223/mois/actions/runs/43712422647/artifacts/14161421174)、ZIP SHA256 ecbf0d142fbe1be9670e70907006eb5321db1d1adf12a571f973cb0a6155d690。
- 1本のPOST JSON payload SHA256 7b760453322c79d901dee392558d71cedfed0b1d5d6b93c3b71781e700fc4892。公開validator HTTP 201の戻りデータSHA256 4dbc7dc47b18a4c62a6421183a7a377345a1d7ca100c1cfd288de6757125a191。
- 実装コード：scripts/analysis/mois_gbif_37104_one_sql_validate.py、テスト：tests/test_mois_gbif_37104_one_sql_validate.py、CI：.github/workflows/mois-gbif-37104-one-SQL-size-validation-20261115.yml。

## 認証付き一括取得の安全な入口も追加済み

Workflow: .github/workflows/mois-gbif-37104-one-SQL-auth-submit-20261115.yml。

- pushの間はバリデーター・署名済み元SQLとCredentials存在の検査だけであり、**自動POST 0**。
- [実行 #43712761737](https://github.com/zuizui0223/mois/actions/runs/43712761737) SUCCESS時点でGBIF_USERNAME（またはGBIF_USER）とGBIF_PASSWORDのGitHub Repository Secretsは **未設定**。実ダウンロード送信0、GBIF downloadKey0、DOI0。
- 安全な実送信には、GitHub Secretsの設定と明示的なmanual dispatch（YES指定）が必要。Workflow_dispatchはこのworkflowがdefault branch上にないと使えないため、現時点の未マージPR #200だけでは手動起動できない。
- パスワードをチャット、Git、JSONに書かない。認証付きPOSTで201が出た場合も、その後の非同期ジョブ完了とDOI・格子ファイルを別途検証する。

GBIF公式技術文書も数万種の大量リクエストは個別searchの連打でなくDownload APIを使うよう推奨している。AWS/GCS/BigQueryの月次公開スナップショットも代替手段だが、別途クラウド処理資源が必要。

次の解析は、GBIFの真正な37,104種×ISEA3H8格子の取得・checksum・DOI確認の後。写真の種レベル色はGBIF地点の局所花色ではなく、花弁の植物学的真値にも未確認の部分があるため、十分なROI QCなしに適応的花色地理差とは解釈しない。
