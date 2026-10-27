# MOIS — 128シャード131,072枚：元画像9色ROI × 500px花候補Lab 統合済み（2026-10-27 JST）

## 現在地

- **131,072 original photo IDs × 131,072 observation IDs を重複なく1行/写真で結合**。元の9,320万枚のソースを再走査せず、既存tier10の128シャードから選定済みの写真を使用。
- 131,072枚は、500px高速画像デコードとフル解像度FCP-v4厳格9色ROIの**2方式を同一写真で実測した数**であり、262,144枚の独立写真とは数えない。
- [元データ確定キュー](https://github.com/zuizui0223/mois/releases/tag/mois-full93m-bulk128-tier10-131072-original-photo-source-v1)：source-reference **32,768**、色データ不足を補うP1-P4 targeted **98,304**、P5 targeted **0**。予備選定時に他の計測済み・予約済み写真/観察UUIDを除去。収集は自然界の野生在来植物相に対する確率的母集団標本ではない。
- [フルサイズROI9元画像実測run](https://github.com/zuizui0223/mois/actions/runs/40340338433)：128/128画像推論成功、128 worker jobsはGitHub Release書込403で**技術的失敗扱い**。全128 SHA付き数値計測Actions artifactsは保存済み。
- [500px花検出元画像実測run](https://github.com/zuizui0223/mois/actions/runs/40341341397)：128/128画像計測ジョブ成功、最終collectは`numpy`未インストールで**技術的失敗**。全128 SHA付き画像計測Actions artifactsは保存済み。
- **[データ専用回収・同一写真統合run 40431588179 成功](https://github.com/zuizui0223/mois/actions/runs/40431588179)**：全256画像worker成果物、全512元画像FCP terminalとソースSHAを検証、同一photo/observation/taxonを照合。画像の再ダウンロード・画像推論は**0**。
- **[公開済み統合済みParquet・厳格ROI9＋Lab・地域季節/分類群優先度別集計](https://github.com/zuizui0223/mois/releases/tag/mois-128shards-131072-fullROI9-plus-500px-sourceID-joined-20261027-v1)**。データ処理は統合済み。実装を含む[PR #189は親PR #188へマージ](https://github.com/zuizui0223/mois/pull/189)済み（merge SHA `540fe4f15a81cb096d2027313ba3d64fe91511cb`）。**mainへの科学解析コード全面反映はまだ別工程**。

## 確定した量的成果

| 指標 | 比較用 reference | 不足補完 targeted | 合計 |
|---|---:|---:|---:|
| 元写真・元観察ID（重複無し） | 32,768 | 98,304 | **131,072** |
| 500px画像を取得・デコード | 32,699 | 98,028 | **130,727** |
| 500pxでYOLO花候補を検出 | 31,873 | 95,727 | **127,600** |
| 元画像で真の9成分ROI成功 | 13,243 | 35,541 | **48,784** |
| 元画像9成分ROI判定不可 | 19,525 | 62,763 | **82,288** |
| 500px陰性だが元画像9色成功 | 26 | 69 | **95** |

取得時点の既知9色ゼロ分類群に新たに厳格9色測定が得られたのは、**targeted P1：184分類群、targeted P2：7,010分類群**、合わせて**7,194名目分類群**（targeted内でtaxon ID重複除去済み）。これは凍結済みソース設計での**新規カバー分類群**であり、他の並行進行中FCPやmainのcross-arm台帳を照合した**世界全体の最終純増数ではない**。

## 方法論上の新しい結論

- 500pxの検出候補があるのに**78,911枚で真の9色ROIが成功しなかった**。単純な花検出では、ROI/flip-gateの失敗構造は解消しない。
- YOLO500が花候補なしとした写真のうち**95枚には真の9色ROIがあった**。500px陰性の全件除外は、実際に測定可能な花色を落とす。次の大規模パイプラインでは、陰性の独立校正標本（優先的に地域/季節/分類群で層化）を必ず保持する。
- 500pxの`floral_box_` CIELABは花候補BOX内の**sRGB連続色指標**であり、植物学的に切り出した花弁色や分光反射率ではない。フルサイズ9色ROIにも開花/自生種/花器官の真偽は未検証のものがある。技術的unknownを勝手に白/有色へ置き換えない。
- source-referenceと色データ不足補完targetedは、抽出時の条件・確率が異なる。**2群を無補正で足して世界白花率に使わない**。
- 今回2方式を同じ131,072枚に別々に実行したため、S3アクセスの重複コストが発生した。次の改善は**各写真を1度だけ取得・デコードし、高速BOX検出と必要な本格ROIを同じワーカー内で実行、SHA完了パケットを再利用すること**。ただし必要な原本の解像度差とfalse-negative検証に対応する。

## 次工程・統合の境界

1. この統合済み本物の原画像9色・500px連続色は **写真ID単位では統合完了**。既存の古い131,072枚画像ワーカーを再実行しない。
2. その後、[PR #188](https://github.com/zuizui0223/mois/pull/188)で引き継いだ処理コードをmainへ昇格させるには、mainの並行収集中photo_id/observation_id/taxon_idのcross-arm台帳とSHAで正確な重複除去を実行すること。未照合の最新main総数と`48,784`や`7,194`を機械的に加算しない。
3. 未分類82,288枚の視覚的監査、否定/陽性の真値・生態学的欠測選択性、および地域/季節別の厳格なバイアス補正が、今後の世界花色分布に対する生物学的検証の主目標。

本ドキュメントは確認済みリリースと成功したデータ専用復旧runの数値に限る。ワーカーの`failure`表記を画像計測失敗と取り違えず、実際の公開物とSHA/ID照合を基準とする。


## 2026-10-27 JST — Signed full four-arm 1,181,146-source-photo union

**This supersedes the earlier “main cross-arm dedup pending” status above for the specified four frozen original-photo arms.** It does **not** claim all future main-branch acquisitions are automatically included.

- **[Data-only global 4-arm SHA join: SUCCESS](https://github.com/zuizui0223/mois/actions/runs/40438094879)**, **[immutable one-photo-per-ID Parquet / 42,111-taxon before-after CSV / source-SHA receipt](https://github.com/zuizui0223/mois/releases/tag/mois-original-1181146-exact-photo-crossarm-451643-strict9-taxongap-20261027-v1)**, [PR #190 merged to `main`](https://github.com/zuizui0223/mois/pull/190), merge commit `4b0ee7fb32efa9f0cc808d2b74bde2bd8389054e`.
- Exact original photo-ID source arms (all with true frozen YOLO11/EfficientSAM original image measurement): `1,043,930/400,710` (previous main original source baseline) + `2,048/682` (geo-season reference+targeted) + `4,096/1,467` (first-colour P2) + `131,072/48,784` (bulk128). The numerator in each pair is attempted genuine source photos, the denominator after the slash genuine strict original **nine-component ROI colour**.
- Exact full four-arm **1,181,146 unique photo IDs, 451,643 strict9 images, ZERO reused photo IDs across arms**. All original photo IDs and nine components saved in the new Parquet; no previously unknown photos silently turned into pigmented or white.
- In the frozen 93,209,351-upstream-image / **42,111 nominal-source-taxon** frame, strict9-coloured taxa **21,235 → 30,048 (+8,813 true net distinct taxon IDs)**; without strict9 **20,876 → 12,063**. This is the taxon-ID cross-arm pure gain, not `85 + 1,467 + 7,194` from independently reported cohorts. Among 8,813 newfound nominal taxa: **186 P1** previously source-rare zero-colour, and **8,627 P2** previously source-rich zero-colour. Remaining zero-colour: **16 P0** without any source photos, **2,862 P1**, **9,185 P2**.
- Supplemental 2,048, 4,096 and 131,072 photo/observation selections each have all independent source-observation UUIDs, zero exact photo-ID overlaps, and **zero new observation UUIDs reused across supplementary source arms or against the older 1,043,930 baseline**. The older historical master itself includes **148,264 extra photographs belonging to already represented observation UUIDs**. These are genuine separate images but not 148,264 independent ecological observation events; hierarchical/multilevel models must account for observer/observation, species, geography and sampling design.
- **64 separate genuinely new 500px-box Lab pilot photos** remain outside the strict original 9ROI estimate; no double-counting or “Lab proxy is flower petal” promotion. The 131,072 bulk128 images' medium500 BOX descriptor is an auxiliary *same-photo* signal, **not another 131,072 original photos**.
- **Scientific boundary:** `1,181,146/451,643/30,048` is a **full SHA-certified join of these four frozen acquisition designs**, not a global independent census of all wild-native flowering plants or a continuously updated union of any later main/FCP newly completed image waves. Preserve separate probability-reference and adaptively selected taxon-gap cohorts. Do not estimate worldwide white-flower frequency unweighted.
- Primary future execution target: **9,185 source-rich nominal P2 taxa** still with zero strict9 and **2,862 rare P1**; improve organ-specific ROI classification and taxon/season/geography sampling support, and independently audit the 95 medium500 false negatives discovered in the same bulk128 dual-method images before excluding detector negatives.

The cross-arm **data reassembly and publication work is complete for the four source-frozen arms**; subsequent main additions require independent source-photo/observation/taxon SHA delta verification, not a blind addition of new ROI or source target counts.
