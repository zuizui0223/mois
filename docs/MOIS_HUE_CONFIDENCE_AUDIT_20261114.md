# MOIS 9色・独立観察信頼度の実データ監査（2026-11-14）

**本監査は元画像の再推定ではない。** GitHub Actionsの署名付き原本を取り出し、ZIP SHA256・内部atlas SHA256・COL XR bridge source SHAを独立照合した後、42,111分類群の既存写真ROI9列を再集計した。写真ROI9は開花器官の独立確認済み花弁色、種全域で固定した花色、遺伝多型の証拠ではない。

## 1. 署名付き原本に対する総数

| 単位 | 実測値 |
|---|---:|
| 全原画像の処理試行 | 1,231,337枚 |
| genuine strict numerical 9-component ROI | 467,565枚 |
| 42,111 source分類群に帰属する処理試行 | 1,198,627枚 |
| 42,111 source分類群に帰属する9-component ROI | 450,006枚 |
| 9色数値ROIを1件以上持つ分類群 | 39,911 |
| 明瞭な代表色を持つ分類群 | 39,172 |
| 数値ROIがまだ0件の分類群 | 2,200 |
| exact COL XR accepted species × numeric9 trait | 37,103 |

**母数の注意**：総写真試行と42,111 source分類群内の試行の差は32,710枚、数値ROIの差は17,559枚。source分類群外の既計測結果を42,111分類群の色割合の分母に混ぜない。

## 2. 元atlasの独立観察に基づく分類状態

| atlas `hue_status` | 分類群 | COL XR exact accepted species |
|---|---:|---:|
| sampling_consensus_single | 2,876 | 2,686 |
| sampled_multicolour | 5,014 | 4,810 |
| multicolour_candidate | 8,026 | 7,607 |
| provisional_single | 23,256 | 22,000 |
| mixed_image_unresolved | 739 | 0 |
| no_strict9 | 2,200 | 0 |

- `sampling_consensus_single`: 明瞭な代表色を持つ独立観察3件以上、他の曖昧な数値ROI観察なし。**これは標本内での一致であり種全域固定の証明ではない**。
- `sampled_multicolour`: 少なくとも2色の各々を独立観察2件以上が支持。種内の遺伝多型とはまだ断言しない。
- `multicolour_candidate`: 独立写真の色が一致しないが、各色への反復支持は不足。
- `provisional_single`: 標本内の代表色を暫定的に割り当てたが、単色の独立観察コンセンサスは未成立。

## 3. 9色代表色の実測内訳

| 9 palette | 代表色分類群 | 代表色39,172中の割合 | 単色コンセンサス群 |
|---|---:|---:|---:|
| white | 22,121 | 56.47% | 2,204 |
| yellow | 6,837 | 17.45% | 596 |
| orange | 260 | 0.66% | 0 |
| red | 1,269 | 3.24% | 29 |
| pink | 1,676 | 4.28% | 11 |
| magenta | 561 | 1.43% | 1 |
| purple | 1,321 | 3.37% | 10 |
| blue | 1,332 | 3.40% | 14 |
| bronze | 3,795 | 9.69% | 11 |

**重要な選択バイアス**：白の比率は代表色39,172分類群では56.47%、単色コンセンサス2,876分類群では76.63%（2,204/2,876）。撮影機会・ROI検出成功・曖昧さ・色別分類群構成が異なるため、厳しい品質フィルター群を全球の白花率と同一視しない。

## 4. 追加QAで検出した台帳の計数問題

- **3,812分類群**で `independent_roi9_observations > after_strict9`。内訳は `provisional_single` 3,602分類群、`mixed_image_unresolved` 210分類群で、既存の多色・単色コンセンサスのstatus群には発生しなかった。現時点ではどちらが異なる母集団を数えたのか原画像ごとの台帳で確認できていない。原因未確定のため単なる「誤記」とは扱わず、source/計測/分類カウントを再照合する。
- ROI9ゼロの**2,200分類群すべてに1件以上の写真処理試行**があり、うち**1,139分類群で3回以上**、**96分類群で10回以上**が失敗している。次は無差別な元画像追加より、この繰り返し失敗群の花器官検出・開花段階・ライセンス・低解像度画像混入の監査を優先する。
- 既存ROI9の`bronze`を含む9色はあくまで**写真内検出領域の色**であり、花弁に属するかどうかの独立検証が必要。複数色状態は必ずしも種内多型ではない。
- 37,103 COL XR accepted keysとatlas 42,111 taxonomy IDの照合は1対1かつ色・status・計数列で不一致0。だが全球GBIF occurrence SQL 149件は準備のみ（認証送信0・全球occurrenceダウンロード0）。

## 5. 再現に使用した固定ソース

- 源atlas: [Actions artifact 14102611429](https://github.com/zuizui0223/mois/actions/runs/43566369683/artifacts/14102611429)
- 源COL XR: [Actions artifact 14102518161](https://github.com/zuizui0223/mois/actions/runs/43566923701/artifacts/14102518161)
- atlas ZIP SHA256: `17da6386cc2970f8f55f6ee5373e0565027cac7d6c9f2f4d2fc84c532f91df6d`
- atlas内部CSV SHA256: `588f74124ff2982af5ada2fe8c0258b0e28f9d91010e0451b0b3023e7c67fbab`
- COL XR ZIP SHA256: `7b8d366d803c9f21d328495bc2517526e1ade6b448bb0c66fa731062e6ae4425`
- COL XRの`receipt.json`に記録された源atlas内部CSV SHAと上記の内部CSV SHAが一致。受容名37,103行に対応する全てのatlas色/status/数値・処理試行カウントが一致。

**判定**：種数優先の写真ROI9カタログ作成は進展。今後の主要ボトルネックは(1)2200ゼロROI分類群の器官検出失敗と(2)39,172代表色の信頼性、(3)3,812台帳計数不整合。GBIF全域分布を写真上の固定花色として塗る解析は保留。


## 追補：3,812件の計数不整合は原因特定・修復済み（2026-11-14）

上記「追加QAで検出した台帳の計数問題」の未確定状態は後続の原本比較で**解決**した。第8〜57波の追加計測更新処理が、既に数値ROIを持つ未分類群の 'after_strict9' を1に上書きしていた。色の数値測定そのものではなく、分類群別累積写真成功枚数が不正確だった。

初期atlas（11月10日、原本署名固定）との独立照合から、全42,111分類群の数値ROI枚数の合計は旧 450,006 から**457,415**へ修復、3,812分類群の欠落**7,409枚**を復元した。全原画像成功枚数467,565、代表色39,172分類群、37,103 exact COL XR受容種の数・色値・分類statusに変更はない。受容種のうち3,427件で枚数欄のみ修正した。過去の原本は書換えていない。

修復根拠・式：[MOIS_STRICT9_COUNTER_RECONCILIATION_20261114.md](MOIS_STRICT9_COUNTER_RECONCILIATION_20261114.md)。[GitHub Actions #43586480972](https://github.com/zuizui0223/mois/actions/runs/43586480972) は署名付き3原本からの訂正台帳再生成・回帰テスト・3成果物保存まで**SUCCESS**。訂正済み解析には旧 'after_strict9' 列ではなく、新しい訂正済み台帳を用いる。
