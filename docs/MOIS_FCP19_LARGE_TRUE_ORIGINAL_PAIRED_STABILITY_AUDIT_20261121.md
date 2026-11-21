# MOIS — 2026-11-21：19枚の同一写真を`large`／真正`original`で比較したFCP ROI9完全監査

## 研究上の判定

**花被が見えると予備判定された19枚で、入力画像を`large`から真正`original`へ変えても、9色ROIの技術的成功枚数は6/19のまま。** 前回分類不能13枚から1枚が回復する一方、前回成功6枚から1枚が失敗へ変わる。厳格な9色単色ラベルは**5→4枚**。フル原画像を使えば花色測定問題が解決するという結論は支持しない。

これは**写真・モデルを固定した画像版差**の実際の再実行結果。19枚は16種・全世界の非無作為な開花注記付きCC画像から440px予備目視で選んだ便宜標本。19件は独立した19種でも、植物の真正花弁色の検証でもない。花被器官の手動マスクは**全件未完了**、種別の確定花色0種、全球37,104種の旧ROI写真色プロキシ・GBIF・LGM係数は変更していない。

## 本物の原画像バイトを使った対比較

| 写真9色の数値ROI | `/original`数値有効 | `/original`技術不明 | 合計 |
|---|---:|---:|---:|
| `/large`数値有効 | 5 | 1 | **6** |
| `/large`技術不明 | 1 | 12 | **13** |
| **合計** | **6** | **13** | **19** |

FCP既定の固定`image_loader`はソース契約上`/large`だけを許可している。そこで凍結FCP`run_batch(loader=...)`を使い、写真ID・観察UUID・ライセンス・`/original` S3ホスト・原画像SHA256・画像デコード・12 MiB上限を確認する別ローダーを付与。画像の取得方法だけを変え、**FCP YOLO11検出器・EfficientSAM encoder/decoder・計測ソースGit blob・閾値は両群同一**。個々の原画像の生バイトSHAは先行46枚のSHAと**全19枚で一致**、取得失敗0件。

先行画像：開花注記`12|13`付きの真正原画像46枚（別の前48枚を除外）と花被候補19枚の目視。前回[19枚を`large`で計測したActions #44753624667](https://github.com/zuizui0223/mois/actions/runs/44753624667)、source artifact #14636423108、ZIP SHA256 `92315ef6679ced282fd822dcb08336a541add546afa949fb657cda76f6fd90ea`。

### 1. `large`で分類不能だった13枚：フル原画像で1枚のみ数値救済

- 前の花マスク検出なし4枚→フル原画像でも0枚が厳格数値ROI成功。
- 前の向き・反転不安定6枚→1枚のみ成功、5枚は依然失敗。
- 前の色パレット曖昧3枚→0枚成功、全3枚が依然失敗（うち1枚は反転色不安定に失敗理由が変更）。
- **真正`/original`画像SHA一致13/13**。9成分数値ROIが出た1枚は*Hedysarum hedysaroides*／photo ID `699083444`。FCPの4状態は`yellow_orange`に通過しても、9成分はbronze=.5474、white=.3583、pink=.0674、yellow=.0269で単色閾値には達せず、**花被を黄や橙色と確定できない**。
- [元原画像13の真正CI #44838300475](https://github.com/zuizui0223/mois/actions/runs/44838300475) SUCCESS、artifact #14679264010、ZIP SHA256 `042cfe9e3be17d5b214c65659d65c98d53a519fcb152260e685c0cd989e71f0e`。個票は[13ペアCSV](../analysis_specs/MOIS_FCP13_TRUE_ORIGINAL_vs_LARGE_PAIRED_RETEST_20261121.csv)。

### 2. `large`で数値有効だった6枚：真正原画像では5枚のみ維持

- 6枚のうち5枚はフル原画像でも9成分数値が有効。
- うち4枚で厳格単色9色ラベルが**完全一致**：*Lonicera caerulea*写真675214739（白）、*Eremophila latrobei* 731316888（赤）、*Amphilophium paniculatum* 590141358（白）、*Banisteriopsis muricata* 100176025（白）。**これらはあくまでFCPの写真ROI色、花弁真値とは異なる。**
- *Lythrum virgatum*写真722033950は両画像サイズで9成分が数値有効だが、ピンク・マゼンタの混色で単色閾値を通らない。
- *Banisteriopsis muricata*写真199487935は`large`で暫定白の数値ROIを取得できたが、真正`original`では花ROI／反転ゲートに失敗。これは技術成功でも**画像サイズによって消失しうる**具体例。
- [6枚の真正CI #44839432470](https://github.com/zuizui0223/mois/actions/runs/44839432470) SUCCESS、artifact #14680281175、ZIP SHA256 `c9cb7d46d98269503ffcaefe9ed7af76703e69e9fb89e06445b48ad4a273e10b`。6枚の原画像生バイトSHAは先行46枚と6/6一致、モデル証明も旧FCP19と完全同一。

## 3. 19枚の統合監査記録

[真正source統合CSV](../analysis_specs/MOIS_FCP19_LARGE_VS_FULL_ORIGINAL_IDENTICAL_PHOTO_COMPARE_20261121.csv)（SHA256 `7b2b34eb10bbcd77823d56134d92c8754cfec575250be87a353058d09b0778cc`）では、19行全件にphoto ID・observation UUID・フル画像SHA・FCP画像SHA・`large`と`original`の技術状態・マスク失敗理由・9成分数値・9色暫定ラベルを保持する。技術不明の9色分率は欠損値として残し、白=0に強制しない。

- `large`数値ROI有効6枚／厳格単色9色5枚。
- 真正`original`数値ROI有効6枚／厳格単色9色4枚。
- 数値結果が両方に存在する5枚、片方だけの数値結果2枚（改善1／悪化1）。
- photo ID／observation UUIDは同じ画像19件をペア比較したもので、追加の生物学的な独立観察数は**0**。
- 解析対象写真の**真正花弁器官の独立承認数0**、種レベルで確定した花色数0。

## 4. 生態学的結論の境界と次の行動

**現在の一次課題は画像数・解像度だけでなく、何の器官を「花色」として測定するかという同定問題。** 特に雄しべや萼、苞、果実、草本の目立たない小花を含む植物群では、種の同定済み観察・Flowers注記があるというだけでは真正花弁のピクセルを保証しない。未確認のROI結果を全球解析の種別9色に昇格すれば、系統・器官・開花ステージの観察偏りを植物地理的な適応と誤認しうる。

次の実証ステップは、19枚のフル原画像を確認できる**同定可能な花冠／花被の器官マスク（複数人または独立アノテーション）**を確保し、マスクとFCP自動抽出のIoU・色差を比較すること。4枚のマスクゼロ、反転失敗群、曖昧パレット群を分けて失敗モードを再評価し、安易な閾値緩和による色の水増しを禁止。手動正解無しでさらに9色ラベルを増やすのは主目標を前進させない。

実装：
- `scripts/acquisition/mois_FCP13_original_resolution_retest_20261121.py`＋`.github/workflows/mois-20261121-fcp13-true-original-resolution-gate.yml`
- `scripts/acquisition/mois_FCP6_full_original_positive_controls_20261121.py`＋`.github/workflows/mois-20261121-fcp6-full-original-palette9-controls.yml`
- 原本数値・固定SHA・修正理由はPR #200のブランチに隔離。PR #206本論の全球種色・LGM・WCVP在来結果に自動投入しない。
