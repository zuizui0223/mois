# MOIS 2026-11-15: 全42,111分類群の写真9色品質監査とGBIF分布解析への接続

## 技術的な9色監査の結論

署名付き source atlas artifact 14157022969 (ZIP SHA256 df8151f0343ebfd32a0c45181e0e05c3d910b29407195f3a99ecb8f28f52f2f1) と COL XR accepted-names＋numeric trait artifact 14157705123 (ZIP SHA256 abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a) を全件照合した。COL XR checklistKey は 7ddf754f-d193-4cc9-b351-99906754a03b。

| layer | count |
|---|---:|
| nominal source taxa | 42,111 |
| genuine photo attempts | 1,231,342 |
| genuine 9-fraction photo ROI successes | 467,566 |
| source taxa with >=1 numeric ROI | 39,912 |
| source taxa with sampled photo modal hue | 39,173 |
| source taxa without numeric ROI | 2,199 |
| numeric ROI but no dominant modal hue | 739 |
| exact accepted COL XR Plantae species linked to photographic traits | **37,104** |
| taxa whose modal photographic hue differs from max of mean 9-fractions | **2,594** |
| accepted species with this modal vs mean9 distinction | **2,466** |
| technical cumulative-photo counter reset still present in source | **3,812 taxa / 7,409 photos** |
| corrected 9-component successes assigned inside nominal 42,111 taxa | **457,416** |
| original numeric successes outside nominal frame | 10,150 |

全39,912数値分類群の9成分は有限・非負・合計1±0.001。ROIゼロ群2,199の9色値は欠測。37,104の受容種キーは重複0、Plantae・rank species・exact accepted、色モード・状態・9分率・元の分類群キーに不一致0。写真ごとの多数決modalとmean9で最大の色は異なる要約統計量であり、2,594差を計算エラーと決めつけない。旧累積写真計数の上書き問題は元台帳に残っていたため、2026-11-14の署名付き派生修復表から差分+1を足して別表に再修復。原本の写真色や受容学名を改変しない。

## photo evidence status / 42,111 taxonomy

| status | n |
|---|---:|
| sampling_consensus_single | 2,876 |
| sampled_multicolour | 5,014 |
| multicolour_candidate | 8,026 |
| provisional_single | 23,257 |
| mixed_image_unresolved | 739 |
| no_strict9 | 2,199 |

単色一致は写真標本内の停止規則であり、種全域で色が固定した証拠ではない。反復支持の複数色は写真色の差であり、遺伝多型の直接証拠ではない。最近の5写真の実マスク位置QAでは4写真で葉・背景・花被以外を検出しており、**39,173分類群の真の花被色が全て確認済みとは主張できない**。bronzeは器官混入に敏感。native/wild statusも未確認。

## Photo modal 9-hue — exact accepted 37,104 species

| hue | accepted species | sampling_consensus_single |
|---|---:|---:|
| white | 20,911 | 2,051 |
| yellow | 6,491 | 564 |
| orange | 243 | **0** |
| red | 1,214 | 28 |
| pink | 1,591 | 10 |
| magenta | 528 | 1 |
| purple | 1,258 | 9 |
| blue | 1,259 | 13 |
| bronze | 3,609 | 10 |
| **total** | **37,104** | **2,686** |

受容種の写真modal白比率56.36%は、全球在来植物相の花被色比率ではない。厳しい単色一致群の白比率76.36%、橙0種という明瞭な色別技術・撮影選択の偏りがある。コンセンサス群のみの9色×環境の全域推定は不可能。

## GBIF distribution analysis next

全球 ISEA3H resolution8 species×grid occurrence SQL は37,103種に対し149件準備済みだが、認証提出0、DOI0。新たな1種 Linzia rosenii を追加して37,104受容種。過去の実GBIF45種100kmグリッドパイロットは28種のみ完全列挙で、17種は取得上限超過としてグリッド不明（欠測）を維持。今回の新しいcountry/continent facet pilotは9色から各20受容種を**GBIF結果を参照せずsource-keyハッシュで事前抽出**し、座標あり・座標問題なし・PRESENTに絞った実GBIF出現国／大陸を取得する。研究目標は範囲・地域差と写真色構成の関連を検討することだが、GBIF occurrence地点の局所花色・野生在来分布・生物学的因果機構と同一視しない。

最終的な全球気候・標高・土壌・島嶼性との species×grid 推定には、認証されたGBIF SQL/DOIと全37,104受容種の真の空間セルデータが別途必要。国単位facetの標本で全球格子解析が済んだとは主張しない。

## Reproducibility

2026-11-15の完全監査Pythonおよび全42,111分類群の修復済みCSV.gz、37,104受容種の品質フラグ付きCSV.gz、2,594のmodal/mean不一致表とreceiptは会話添付「MOIS_全42111分類群_9色QC_COLXR37104_20261115.zip」に保存。SHA256 eff371f8da3a70d23371a1b78236e089da1206b82ec9fb68ef860926e2965cf4。実マスク位置の独立検証は docs/MOIS_FROZEN_V4_FIVE_MASK_LOCATION_QA_20261115.md。