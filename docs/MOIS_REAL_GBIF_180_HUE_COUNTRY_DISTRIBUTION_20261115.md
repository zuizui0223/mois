# MOIS 2026-11-15：写真9色×実GBIF分布の180受容種解析、および全球37,104種の空間格子取得ゲート

## 1．再現性のある観測単位と区別

2026-11-15の署名付きiNaturalist花領域9色atlasは42,111名目分類群のうち39,912群に数値ROI、39,173群に代表撮影色を割り当てた。GBIF v2 COL XR受容学名と厳密に一致するのは37,104種。色は**写真から算出した種レベルの特徴**であり、GBIF座標で測定された個体色・野生在来色ではない。過去の5写真mask独立検証では葉・背景等の誤検出があるため、39,173群の花弁色を全面確認済みと主張しない。全分類群・学名・写真9分率・過去の累積写真カウンタの上書き問題については [全数監査](MOIS_COMPLETE_HUE_QA_AND_COLXR_GBIF_20261115.md) を参照する。

## 2．GBIF公式APIでの180種実分布の取得

[Actions run #43707058183](https://github.com/zuizui0223/mois/actions/runs/43707058183) SUCCESS、[artifact #14158822487](https://github.com/zuizui0223/mois/actions/runs/43707058183/artifacts/14158822487)、ZIP SHA256 **2181898c7428172cd039b208a5e3fc5169c12d3a8d324bef05a84f3a6a89a446**。

対象は37,104 exact COL XR受容種から、**9 photographed modal hues×各20種**をGBIF分布データ閲覧前に固定したspeciesKeyのSHA256順位から選定した180種（seed `MOIS-COLXR-20261115-COUNTRY20-EACH9-HUE-BLIND`）。GBIF公式 occurrence/search の georeferenced facet: `hasCoordinate=true`、`hasGeospatialIssue=false`、`occurrenceStatus=PRESENT`。GBIF APIレスポンス180件の生JSON、参照URL・SHA256、種×国・種×大陸の取得カウントを保存した。

| 集計 | 数 |
|---|---:|
| GBIF facet APIレスポンス / 実GBIF記録のある種 | 180 / 180 |
| APIが返した座標付きレコード数の合計 | 2,167,354 |
| 種×国（出現記録1件以上） | **1,913** |
| 種×大陸（出現記録1件以上） | **400** |
| 元の全数37,104種のGBIF occurrence取得完了 | **NO** |

2,167,354は**GBIFが返した件数の合計**であり、216万行の個別occurrence CSVを自分で取得した数ではない。国別の記録数は国土面積や観測努力に影響され、現地の定着や在来性を確認したものではない。

## 3．分布国数：色だけで明瞭な差は検出されず

9色均等抽出180種の国数のKruskal–Wallis: H=4.6102, **p=0.7983**。写真modal白 vs ほかの色、応答log(1+記録国数)、GBIF座標付き記録数log(1+n)調整、HC3: 白色係数 **+0.06308**、95%区間 **−0.3021～+0.4282**、p=**0.73494**、adjusted R²=0.3854。GBIF記録数と記録国数のSpearman rho=+0.49665。努力量のproxyを入れた後の色の分布幅係数は現状では支持されない。

## 4．大陸×写真色（探索、未確証）

各色20種のうち南米にGBIF記録がある種：
白9、黄8、橙14、赤8、桃10、赤紫6、紫2、青3、銅9。種を交換単位とする9色×南米presenceの10,000回 permutation **p=0.0039**、7大陸でのBenjamini–Hochberg q=**0.027**。南米のGBIF記録10件以上を要求した場合にも、橙10/20 vs 紫1/20で方向は維持。

ただし大陸の選択は探索的・事後であり、仮説検証の独立標本にはなっていない。写真色とGBIF経緯度が同じ個体・地点で測られていないこと、分類群の系統差、移入外来種、国別観測努力、source taxon bias、オフターゲット花領域測定の混入を理由として**地理環境における花色選択・適応の証拠とは解釈しない**。

写真色の証拠は同数・同強度ではない。今回180種の内訳：
`provisional_single`100、`multicolour_candidate`49、`sampled_multicolour`29、`sampling_consensus_single`2。このため全180種で遺伝的に固定した花色を示したものではない。

南米の探索的な差の再現には、同じハッシュ順位で**各色の21～40番目にあたる別180種**を取る。初回180種とのCOL XR key重複0を保つ。独立holdoutジョブ [workflow](../.github/workflows/mois-gbif-180-colour-country-holdout-20261115.yml) の結果を本記録へ追補する。

## 5．全37,104種のGBIF全球ISEA3H空間格子への制約

[GBIF SQL validation + auth-gate run #43708178817](https://github.com/zuizui0223/mois/actions/runs/43708178817) **SUCCESS**、[artifact #14159013597](https://github.com/zuizui0223/mois/actions/runs/43708178817/artifacts/14159013597)、ZIP SHA256 **8c642ab261e6b0498e33f7f305946a89af610d9c2c75c07cfda7c81207e3571b**。

受容37,104 specieskeysを**250種×148件＋104種×1件＝149件**の非重複GBIF ISEA3HCode(8,latitude,longitude,0)分布格子集計SQLに変換。149件内のspeciesKey総数・重複0・格子関数を独立検証、GBIF公式SQL構文validator HTTP **201** 受理を確認。

しかし実行環境のGBIF_USERNAME（またはGBIF_USER）とGBIF_PASSWORD Secretsは**未設定**だったため、認証SQL送信0、正式download key/DOI 0、global species×cell取得0。アップロード先のGitHub Actions workflowは [こちら](../.github/workflows/mois-gbif-colxr-37104-global-cube-auth-gate-20261115.yml)。秘密情報を本文・CSVへ記録していない。

**国×種facetの180種と、全球ISEA3H種×セル37,104種は別のデータ層**。前者で緯度・標高・降水・日射・土壌の全球セル回帰が完了したとは主張しない。global cubeを実際に取得できた後にspecies×equal-area gridのpresenceへ変換し、観察努力、季節、在来性、系統構成を別に考慮する。

## 6．再現コードと監査成果

- scripts/analysis/mois_gbif_colour_country_facets_20261115.py：受容名のみ、色均等ハッシュ事前固定、GBIF country/continent実facet。
- tests/test_mois_gbif_colour_country_facets_20261115.py：欠測やfacet集計のエラー停止、ハッシュ順位の重複防止。
- .github/workflows/mois-gbif-180-colour-country-20261115.yml：180実APIレスポンス取得・記録数努力調整OLS。
- .github/workflows/mois-gbif-180-colour-country-holdout-20261115.yml：発見標本と重複しない180種の検証。
- .github/workflows/mois-gbif-colxr-37104-global-cube-auth-gate-20261115.yml：149件全GBIF ISEA3H8 SQLリストと認証ゲート。
- 全42111分類群技術QA＋37,104種の撮影9色traits ZIP SHA256 eff371f8da3a70d23371a1b78236e089da1206b82ec9fb68ef860926e2965cf4。GBIF実取得180種ZIP SHA256 2181898c7428172cd039b208a5e3fc5169c12d3a8d324bef05a84f3a6a89a446。


## 追補（2026-11-15）：独立180種で南米の色別差は再現せず

上記の「独立holdoutジョブ待ち」は終了した。[run #43708230511](https://github.com/zuizui0223/mois/actions/runs/43708230511) が正式SUCCESS、source ZIP/artifact #14159523021 digest ec1df9eb35e856b1020e92a755aff8be8c121b640602e8e3191d528cc531df1aを確認。先の180種とCOL XR accepted keyの重複0、9色×20種ずつ。

南米のorange対purple presenceは探索14/20対2/20 → holdout **8/20対8/20**。独立9色×南米10,000 permutation p=0.8569、GBIF記録数とphoto statusを調整した9色GLM-LRT p=0.8255、orange>purple Fisher p=0.6262。**元の南米色別パターンは支持されない。**

全360受容種のphoto evidence statusは provisional_single 205、multicolour_candidate 93、sampled_multicolour 59、sampling_consensus_single僅か3種。写真modalとmean9最大色が不一致59種。GBIF country breadthの白写真係数はdiscovery p=0.735、holdout p=0.0095、統合 p=0.0439 と標本間で異なり、一般則や種全域の白花固定色の証拠に採用しない。

統合した完全なソース・判定・GBIF空間格子認証ゲートは [docs/MOIS_REAL_GBIF_360_HOLDOUT_PHOTOHUE_QA_20261115.md](MOIS_REAL_GBIF_360_HOLDOUT_PHOTOHUE_QA_20261115.md)。全球37,104種のGBIF SQL149件は公式validator HTTP201を通過したが、GitHub GBIF Secrets未設定のため正式なspecies×gridとDOIは未取得である。
