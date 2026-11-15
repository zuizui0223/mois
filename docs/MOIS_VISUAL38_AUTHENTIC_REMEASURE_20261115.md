# MOIS：38既存写真の真正FCP-v4再計測と花被色QA（2026-11-15）

## 署名付き実行結果

[Actions #43673586307](https://github.com/zuizui0223/mois/actions/runs/43673586307) の3ジョブ：source SHA検証・固定モデル実測・独立監査がすべてSUCCESS。従来strict9=falseだった38写真を再測定し、**13枚で有限非負かつ和1の9色分率**を得た。うち**10分類群で明瞭な代表色候補**、3分類群は9数値が存在するが単色判定の閾値を満たさない。残りはROI/flipゲート23枚、曖昧な9色組成2枚。画像取得失敗0枚。

| 指標 | 数 |
|---|---:|
| 対象 | 38旧photo ID、38分類群、38独立観察UUID |
| 再測定後9成分数値成功 | **13** |
| そのうち明瞭な9色代表候補 | **10** |
| 9成分成功のみ・明瞭な単色なし | **3** |
| ROI/反転品質ゲート失敗 | 23 |
| 色組成の曖昧さで失敗 | 2 |
| 新規ユニーク写真ID | **0** |
| 本体atlasへの数値・種名採用 | **0** |

利用画像はiNaturalistオープンデータのprovider **large派生画像**であり、original解像度とは呼ばない。別の写真を追加取得したものではない。対象は可視花を選んだ38件なので13/38を全球失敗群への回収率と見なさない。

最新参照atlasの1,231,341ユニーク処理試行・467,565数値ROI写真・39,911数値分類群・39,172代表色分類群・2,200未数値分類群は変更なし。主論文の白花率や環境係数の再推定もしていない。

## 再現とデータゲート

- Source準備は従来の90枚目視QA、118万枚のsource写真ID台帳、2026-11-15署名済み42,111分類群atlasを照合し、38分類群すべてが最新でも9色数値未取得であることを確認した。
- 38行source SHA256：6c035941f6b7b77e0369ec8bbb62f6f1abf8046ad0187634a284d9aabfb9ce4b。凍結入力の検査6件SUCCESS。
- Source [artifact #14145521186](https://github.com/zuizui0223/mois/actions/runs/43673586307/artifacts/14145521186) ZIP SHA256 b35e963d208a82da580dc6d093b95b000530ce9f2204791102c3a70e49e8608f。
- Fixed method commit 9fae6ccdf684a46026f72ba12e98de2c5c54bf2a、measurement Git blob b1d7091db44682559bbe4b48a9f3b8db648501a6。YOLO11 detector SHA256 f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5。EfficientSAM encoder/decoderも署名固定。
- Worker [artifact #14145555692](https://github.com/zuizui0223/mois/actions/runs/43673586307/artifacts/14145555692) ZIP SHA256 8d49aeb4b9362f76a3b2b290850b6e7b0a42847052c2e4dc5e9bb7e25d79bce4。元terminal CSV SHA256 6b3f7ea4981af8ec9546eafad4e5421765d8e674e4b64c4218194bddbed85e8a。Worker receipt SHA256 3f59132aa9ef408b4b32bb6ecfa76ce94f21c120351990016545691cc6882684。
- Independent audit [artifact #14144984441](https://github.com/zuizui0223/mois/actions/runs/43673586307/artifacts/14144984441) ZIP SHA256 a4d76f529046749b5e48a96b375c118cd662e912cd42a6837d7584f1b2e8f4a5。全38写真のmeasurement status・9分率・candidate_taxon_patch_NOT_PROMOTEDを出力し、本体atlasを変更しない。

## 新規13件に対する写真内の色・器官QA（AI一次判定）

13枚の既存medium縮小画像を再表示し、実測9色と写真中の花被の見た目が方向的に一致するかだけを確認した。**FCPが出したマスク自体は保存されていない**ので、以下は優先監査フラグであり植物学的真値・測定器の誤検出確定ではない。

| 種名 | 自動代表色候補 | 画像の見た目からの要点 |
|---|---|---|
| Aristolochia macroura | white | 茶色〜ベージュ主体の花被に対し白95.7%。**要確認** |
| Aristolochia macrocarpa | white | 暗赤褐色の花に対し白99.8%。**要確認** |
| Aldama ghiesbreghtii | yellow | 黄色い花と方向整合。ただし花マスク未確認 |
| Pleurothallis lopezii | なし | 数値9色成功、明瞭単色なし |
| Phyllarthron arenicola | white | 淡紫〜白の花、単色白100%は不確実 |
| Aganisia pulchella | white | 白花と方向整合、マスク未確認 |
| Artabotrys zeylanicus | white | **要確認**：緑色花被に対して白100%。前回の緑色花被パッチでは9色分母0 |
| Bignonia sciuripabulum | white | **要確認**：紫花に対して白99.5% |
| Eremophila linearis | なし | 数値9色成功、明瞭単色なし |
| Libidibia parvifolia | bronze | **要確認**：黄花に対してbronze67.0% |
| Kalanchoe waterbergensis | なし | 数値9色成功、明瞭単色なし |
| Droogmansia pteropus | purple | 紫花と方向整合、マスク未確認 |
| Daviesia audax | red | 赤橙花と方向整合、マスク未確認 |

代表色候補10件のうち、**視覚的な食い違い要確認5件・方向整合4件・境界曖昧1件**。これは非無作為13件のAI視覚スクリーニングであり、測定精度推定ではない。

考えられる原因は複数：非花被ROIへの検出、花被と苞/雄しべの混在、緑・茶・黒の12アンカーからの除外後に残った9色画素だけで再正規化される仕様、provider large画像と元写真の処理差。現状は元マスク未保存であり、相違を単一原因に断定しない。

## 次の限定主タスク

**要確認5件を優先し、モデル検出マスクの位置・花被器官の正しさ・12アンカー生画素数・9色残存率を監査する。** 現行の凍結9色法を勝手に変更しない。数値が出たからといって13分類群を本体atlasへ自動昇格しない。署名付き元画像と植物学的QAがそろうまでは13件すべて「再測定候補」のままとする。

関連：scripts/analysis/mois_visual38_existing_v4_20261115.py、tests/test_mois_visual38_existing_v4_20261115.py、.github/workflows/mois-visual38-existing-v4-20261115.yml。いずれもPR #200の未マージ解析ブランチ上。
