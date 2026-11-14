# MOIS 2026-11-14｜9色ROI未取得2,200分類群の写真原本IDに基づく失敗監査

## 完了したこと

- [修復済み42,111分類群台帳](https://github.com/zuizui0223/mois/actions/runs/43586480972)（SHA固定）で `after_strict9==0` を抽出：**2,200分類群**、全て処理試行1件以上、合計 **7,838写真試行**、数値ROI9成功0。
- 失敗試行の多さで3階層：**10回以上96分類群／1,278回**、**3〜9回1,043分類群／4,759回**、**1〜2回1,061分類群／1,801回**。
- 96分類群には2026-11-10初期atlasの署名付き元写真 **256件** の試行記録がある。2016年などの撮影日や原写真URLを画像内の花とみなさず、[1,181,146枚の元写真ID台帳](https://github.com/zuizui0223/mois/releases/tag/mois-original-1181146-exact-photo-crossarm-451643-strict9-taxongap-20261027-v1)（Parquet SHA256 `3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835`）と突合。
- 元写真IDと観察UUIDで既存の実測が厳密に9色ROI失敗だったことを確認し、独立観察を優先して **96分類群・220枚（観察UUIDも220件すべて存在し重複0）**をレビュー候補化した。内訳は分類群ごとの主レビュー写真**96枚**、予備**124枚**。初期256件のうち36件は最大3件/分類群の制限により今回のキューから除いた。
- [Actions #43589424434](https://github.com/zuizui0223/mois/actions/runs/43589424434) は元source SHA、42,111分類群台帳SHA、初期atlas SHAを検証し**SUCCESS**。[artifact #14111328732](https://github.com/zuizui0223/mois/actions/runs/43589424434/artifacts/14111328732) に全2,200分類群優先表、220枚の実原写真ID/観察UUID/URL、96枚のレビューHTML、技術的4色陽性だが9色ROI失敗の29分類群別表を出力。
- **この監査で元画像ピクセルのダウンロード・再推定0、新規ROI9成功0**。花や花弁の存在を自動判定したり、失敗写真へ花色を補完したりしていない。

## 見つかった技術的な再検査候補

220枚の履歴写真の中で、**31枚・29分類群**は旧 `original_roi_four_hue` が `mixed_uncertain` 以外の4色技術ラベルを持つのに、同写真の `strict9=false` だった。

| 旧技術4色ラベル | strict9失敗写真 |
|---|---:|
| yellow_orange | 16 |
| white | 10 |
| red_pink | 3 |
| blue_purple | 2 |
| **合計** | **31** |

31枚は全て `source_subarm=historical_main`。これは**測定パイプラインの応答条件が異なることを示す技術的手掛かり**で、独立した植物学的証拠ではない。29分類群の優先検証表では、旧4色ラベルを9色ROI成功に「救済」していない。

写真照合例：[Peyritschia virletii / photo 107788481](https://www.inaturalist.org/photos/107788481)（旧white）、[Piper guahamense / photo 673663811](https://www.inaturalist.org/photos/673663811)（旧yellow_orange）。これらのリンクは実際に元台帳へ一致するphoto IDであり、**ここでは写真内の花器官有無を未判定**とする。

## 目視レビューで分ける判定（未入力）

1. **open_flower**：開花した花器官が画像内に視認できる。
2. **bud_only / fruit_or_vegetative**：蕾、果実、栄養器官等はあるが開花花器官を判断できない。
3. **nonplant / image_unavailable**：写真対象外またはリンク欠落。
4. **uncertain**：低解像度・部分遮蔽・植物同定等で判定不能。

開花した器官が確認された画像のみに、`open_flower_but_ROI9_failed` を付与する。未確認のまま `flower_absent` などと決めない。

## 次に必要な検証と境界

最初の96枚のblind HTMLレビュー結果を、旧4色技術情報を見ずに記録し、その後29分類群の技術的4色陽性との対応を確認する。同じ分類群内で主写真と予備写真の相違や観察UUIDを保持する。花が実際に写っている場合のみ、適切な花器官ROI/9色分率の再計測を別の原画像監査手続きとして検討する。**旧4色ラベルは数値9色の代用ではない**。

原写真参照リンクは写真のホスト・ライセンスを変更せず、写真データやサムネイルを再配布しない。全96分類群について「花の有無」判定完了とはまだ主張しない。元の37,103 COL XR受容名と全球GBIF occurrence取得状況（全量SQLは未認証提出）は変更なし。

### 再現コードと結果

- [抽出・SHA検査・blindレビューHTML生成](../scripts/analysis/mois_zero9_failure_review.py)
- [回帰テスト](../tests/test_mois_zero9_failure_review.py)
- [GitHub Actions source-only QC](../.github/workflows/mois-zero9-96-failure-review.yml)
- [220枚/29分類群データと署名付きreceipt](https://github.com/zuizui0223/mois/actions/runs/43589424434/artifacts/14111328732)
