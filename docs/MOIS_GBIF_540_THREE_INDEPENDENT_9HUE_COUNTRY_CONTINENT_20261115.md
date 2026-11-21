# MOIS 2026-11-15：写真9色×実GBIF分布540種・3独立標本の定量再現判定

## 実際に完了したこと

- source＝GBIF COL XR exact accepted 37,104受容植物種。**9写真modal huesごとに20種を事前ハッシュから抽出した180種×3分離標本、計540受容種**。GBIF accepted speciesKey重複0。
- GBIF公式 occurrence/search country/continent facets を3標本で180件ずつ本当に取得。GBIF APIが返した**座標付きoccurrence数合計6,023,946**。これは約602万行をダウンロードしたという意味ではなく、540回のAPI検索に含まれるレコード集計値。
- GBIF original facet APIの全540件は各run ZIPにJSON/URL/SHA256で保存。署名されたsource ZIP3つ：discovery [artifact #14158822487](https://github.com/zuizui0223/mois/actions/runs/43707058183/artifacts/14158822487)、holdout1 [#14159523021](https://github.com/zuizui0223/mois/actions/runs/43708230511/artifacts/14159523021)、holdout2 [#14160537231](https://github.com/zuizui0223/mois/actions/runs/43710488238/artifacts/14160537231)。
- full 540×3-cohortの再現CI [run #43711365162](https://github.com/zuizui0223/mois/actions/runs/43711365162) **SUCCESS**、[artifact #14160269115](https://github.com/zuizui0223/mois/actions/runs/43711365162/artifacts/14160269115)。ZIP SHA256 **c60d2a25a6daf846fbc883b93e973a7cffa60f753667c8974480e4e2133f937a**、内部3CSVのSHA256全一致を別途確認。再現コードは [.github/workflows/mois-gbif-540-three-cohort-QA-20261115.yml](../.github/workflows/mois-gbif-540-three-cohort-QA-20261115.yml) 内に保存。

## 主要結果：グローバルに安定した写真9色×GBIF記録国数の差は得られない

| sample | n | 白写真色 vs その他8色: β (log(1+国数)) | GBIF記録数調整後HC3 p | 9色全体の非調整Kruskal p |
|---|---:|---:|---:|---:|
| discovery | 180 | +0.0631 | 0.7349 | 0.7983 |
| holdout1 | 180 | +0.4941 | 0.00949 | 0.1549 |
| holdout2 | 180 | +0.00822 | 0.9730 | 0.3510 |
| **3標本統合** | **540** | **+0.18246** | **0.12461** | — |

統合モデルはGBIF座標付きレコード数log(1+n)と標本cohortを調整、HC3 robust CI **[−0.05040, +0.41532]**。9色カテゴリー全体を同様にlog GBIF記録数、3 cohortと写真色evidence statusで調整したjoint Wald p=**0.64001**。独立holdout1でのみ白色の記録国数係数が有意だが、holdout2ではほぼ0。**「白花植物ほど広い分布」という確証的な一般則にはできない。**

## 6大陸×9色の独立再現

6大陸（アフリカ・アジア・欧州・北米・オセアニア・南米）のpresenceをそれぞれ10,000回の写真色ラベル置換で評価、各標本内6大陸BH補正。**三標本すべてで有意となる地域は0。** 南米の橙写真色の記録あり数は14/20→8/20→7/20、紫は2/20→8/20→8/20。探索で見えた橙偏在は新しい二標本では消えた。

アフリカは1・2標本で一時的な色差が見えたが、詳細な色方向も一致しない。たとえばpinkのアフリカ記録は1/20→12/20（最初の2標本）、purpleは11/20→2/20。三標本目のアフリカ9色全体の置換p=**0.05639**、6大陸BH q=**0.33837**。分布記録と写真ラベルの共起は**同じ局所写真の色**や在来分布ではない。

## 色分類の品質制約が最重要

540 COL XR受容種に対する写真由来9色状態の内訳：

| photographic evidence status | accepted species |
|---|---:|
| provisional_single | 320 |
| multicolour_candidate | 130 |
| sampled_multicolour | 83 |
| sampling_consensus_single (3+ observations) | **7** |
| total | **540** |

540種のうちmodal photograph hueと平均9分率の最大色が一致しない種は**79種**。**この540種について外部から専門家検証済み花被ROIは0**。この技術的な9色ラベルは花弁の真値・遺伝的に固定した色・局所植物群集の色ではない。GBIFはgeoreferenced・no geospatial issue・PRESENTを指定しているが、**在来か導入かを種・国ごとに識別できていない**。記録国数はGBIF記録量に依存し、真の分布面積/適応度の指標ではない。

従って本時点で科学的に言えるのは「**撮影標本の技術的9色と種レベルGBIF記録国数に、3独立パネルで安定した一般関連は検出できなかった**」までである。有意差なし＝真の効果ゼロの証明ではない。気候・土壌・標高・島嶼性の全球種×格子推定は、**正式GBIF ISEA3H r8 occurrence cubeの認証付き取得（37,104種149 SQL・未送信）**と花器官QAの別評価が必要。

## 正確な再現資産

- CI source：[GBIF 3-cohort signed workflow](../.github/workflows/mois-gbif-540-three-cohort-QA-20261115.yml)
- 6大陸2-cohort test implementation：[scripts/analysis/mois_gbif_continent_two_cohort_qa_20261115.py](../scripts/analysis/mois_gbif_continent_two_cohort_qa_20261115.py)、[回帰テスト](../tests/test_mois_gbif_continent_two_cohort_qa_20261115.py)
- CI output artifact #14160269115：
  - 540_species_exact_accepted_country_photo_hue_QA.csv
  - 540_country_range_white_vs_nine_hues_by_cohort.csv
  - 18_three_cohort_sixcontinent_permutation_QA.csv
  - receipt.json
- 参考：より厳密な540種の写真trait-QA／mean9最大色不一致は、全37,104種の署名付き写真9色表（ZIP SHA256 eff371f8da3a70d23371a1b78236e089da1206b82ec9fb68ef860926e2965cf4）へ厳密に結合した別パッケージで監査済み。ソースは改変しない。
