# MOIS：固定FCP ROI-v4の誤測色5写真・実マスク位置検証（2026-11-15）

## 実行したこと

2026-11-15 [Actions #43694635524](https://github.com/zuizui0223/mois/actions/runs/43694635524) **SUCCESS**。以前38件を真正再計測した[元run #43673586307](https://github.com/zuizui0223/mois/actions/runs/43673586307)の署名付きsource ZIPと数値terminal ZIPを再利用し、視覚的に花色と測色値が不一致だった5写真について、同じ固定**FCP-v4 YOLO11 detector＋EfficientSAM**を同じprovider large画像で実行した。

今回の画像byte SHA、マスク面積、水平反転mask IoUと色差ΔE、**12アンカーの全画素カウント**は5件すべてで過去の署名付きterminalと一致。計算条件を変えて新しい花色を捏造せず、検出器が**実際に選んだ画素の位置**だけを別のQA手続きで取り出した。原RGB画像はGitやGitHub Actions成果物へ保持せず、診断用2値マスクPNG5枚と元画像上の検出bbox位置、mask重心・面積を別artifactに保存した。

固定run [artifact #14153673735](https://github.com/zuizui0223/mois/actions/runs/43694635524/artifacts/14153673735)、ZIP SHA256 **648e56389d6dc466483190a9449d794c05d54ffb69758d414ac9feb49df34fa8**。5件とも元photo-byte hashと12色画素数完全一致、再現テスト5件成功。これは**写真全域／花の部位の真値ではなく、モデルが測った場所の再現監査**。

## 写真と凍結ROIマスクの重ね合わせ結果

実画像medium派生版に、large測定画素から得たmaskを最近傍で縮小重ね合わせ。**AIによる写真内器官の一次確認**で、植物学専門家の独立アノテーションではない。したがってピクセル1個単位の花被IoUの真値評価とは区別する。

| 既存写真 | 元測色候補 | 画像面積に占めるmask | 位置・機構の一次判定 |
|---|---|---:|---|
| [Aristolochia macroura / 164665576](https://www.inaturalist.org/photos/164665576) | white 95.7% | 1.60% | 花被の隣にある**葉**を選択。白の計数は花被に属しない |
| [Aristolochia macrocarpa / 111641760](https://www.inaturalist.org/photos/111641760) | white 99.8% | 0.16% | **画像最上端**の狭い背景／周辺部分を2か所選択。赤褐色の花被を取れていない |
| [Artabotrys zeylanicus / 255188296](https://www.inaturalist.org/photos/255188296) | white 100% | 2.55% | 画像**上端の左右の明るい背景・ボケ**を取得し、中央の緑色花被を完全に外している |
| [Bignonia sciuripabulum / 446129531](https://www.inaturalist.org/photos/446129531) | white 99.5% | 0.16% | **右端の2か所だけ**に微小な検出。主要な紫色花冠が含まれず、白ラベルを代表花色にできない |
| [Libidibia parvifolia / 487988985](https://www.inaturalist.org/photos/487988985) | bronze 67.0% | 6.98% | 18検出枠が**実際に黄色の花房にも重なる**。こちらは完全な背景false positiveではなく、花房内の暗い組織・器官混合と固定アンカーパレットの問題を区別する |

**5件の機構分類：主な花被から外れた検出が4件、花房には重なるが色解釈が不確実なもの1件。** 5件とも花被の信頼できる9色真値としては採用しない。

前回の緑の花被の手選パッチでは9色の分母が0となったArtabotrysについて、今回の固定v4 maskの20,028画素は白11,834／brown7,879／green42／black273だった。この食い違いは**緑アンカーを除外することだけでは説明できず、花被以外の明るい背景のマスク化が主要な原因**と実マスクから判明した。

Bignoniaでもmask1,281画素にはwhite788／green489／yellow4、**紫・青0**。花は紫でも実測ピクセルはほぼ白／緑だった。Bignoniaのmask最小外接枠は写真縦方向に長いが、これは上下に離れた2つの小領域をまとめた外接矩形であり、縦長全領域を検出したわけではない。

Libidibiaではmask48,784画素から**bronze23,625／yellow11,550／brown10,681／green2,265／black582**が数えられた。**green、brown、blackを除外した9色の再正規化**によってbronze67.0%、yellow32.8%となる。ただし撮影条件・蕾・花糸・花梗などの混入、影、組織の区別はまだ定量的に検証されておらず、「パレット自体のバグ」と一意に断定しない。

## データと検証境界

- 元38枚の取得・測色source ZIP SHA256：b35e963d208a82da580dc6d093b95b000530ce9f2204791102c3a70e49e8608f。
- 元38枚worker ZIP SHA256：8d49aeb4b9362f76a3b2b290850b6e7b0a42847052c2e4dc5e9bb7e25d79bce4。terminal内部SHA256：6b3f7ea4981af8ec9546eafad4e5421765d8e674e4b64c4218194bddbed85e8a。
- 凍結FCP method commit：9fae6ccdf684a46026f72ba12e98de2c5c54bf2a、detector SHA256：f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5、EfficientSAM encoder/decoderは元receipt署名と一致。
- 今回のマスクPNGは**診断用の別成果物**。主43,000近傍分類群カタログの写真計数・形質データ・9色結果には戻していない。主42,111分類群の累計写真1,231,341／数値ROI467,565／明瞭な代表色39,172は今回変更0。
- 画像のsourceは全部 CC BY-NC。オーバーレイによる二次配布は著作者名・写真ID・ライセンスを併記し、元データの帰属条件を守る。GitHubには原RGBを追加しない。
- **これは明らかなdiscordanceを事後選別した5写真**であり、全467,565枚の誤検出率・白色比率の補正値は算出できない。観察数を追加したという主張もしない。

## 次に必要な検証

**優先すべきは9色閾値の緩和ではなく、花被と葉・背景を区別できる独立した器官領域の正確性評価。** 現行の固定法は再現性のため維持し、改善方法を試すなら別版の花被ROIと事前固定した独立holdoutで、(1)petal IoU (2)背景誤検出 (3)色別失敗率 (4)緑花被の欠測を評価する。

- コード：scripts/analysis/mois_visual5_frozen_mask_qa_20261115.py
- 回帰テスト：tests/test_mois_visual5_frozen_mask_qa_20261115.py
- ワークフロー：.github/workflows/mois-visual5-frozen-mask-qa-20261115.yml
- 写真別の公開可能な非RGB数値QA表：results/mois_visual5_true_mask_position_qa_20261115.csv
- 5件の直接オーバーレイPNG・著作者表記は研究者向け私的QAパッケージで管理（Git追跡外）。
