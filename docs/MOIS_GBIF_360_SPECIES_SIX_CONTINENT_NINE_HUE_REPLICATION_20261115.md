# MOIS 2026-11-15：全球6大陸×写真9色・GBIF実分布の360種独立再現

## 実行結果

GBIF公式 occurrence/search の実際に取得済み国・大陸facetだけを再集計した。探索180受容種と独立holdout180受容種、GBIF COL XR specieskey重複0。**各標本の9写真色状態×20種、同じ6大陸を別々に10,000回置換検定し、6大陸のBH多重検定を適用**した。

[実行 Actions #43710777192](https://github.com/zuizui0223/mois/actions/runs/43710777192) は元ZIP2本SHA検証、4件の回帰テスト、実source計数と置換集計、結果artifact保存がすべてSUCCESS。[signed artifact #14160322719](https://github.com/zuizui0223/mois/actions/runs/43710777192/artifacts/14160322719) ZIP SHA256 585013f8626bb429e068276c92f785483994c8f7649c7f1d3dc7c9afb420aa65。

| 大陸 | 探索180（p / 6大陸BH q） | 独立180（p / BH q） |
|---|---|---|
| アフリカ | 0.02210 / 0.06629 | 0.00140 / 0.00840 |
| アジア | 0.82842 / 0.92751 | 0.20818 / 0.57894 |
| 欧州 | 0.63694 / 0.92751 | 0.93331 / 0.93331 |
| 北米 | 0.92751 / 0.92751 | 0.54165 / 0.81247 |
| オセアニア | 0.55774 / 0.92751 | 0.28947 / 0.57894 |
| 南米 | 0.00390 / 0.02340 | 0.85551 / 0.93331 |

**両標本で6大陸BH後にq<0.05となった地域はゼロ。** 最初の南米の橙色偏在は再現せず、アフリカは二標本とも未補正p<0.05だが、初回 q=0.066で不支持。

さらに、アフリカでの具体的色別の**方向自体が反転**した：
- pink 写真色がアフリカにGBIF記録あり：探索1/20→holdout12/20。
- purple 写真色がアフリカにGBIF記録あり：探索11/20→holdout2/20。
- orange も探索10/20→holdout4/20。

したがって「同じ大陸で総体的に9色の偏りが見える」ことと、「同じ花色が同じ大陸に多い」は同義でない。探索有意な大陸を事後選別して進化・選択・生物地理シンドロームに格上げしてはならない。**pが有意でないことは効果ゼロの証明でもない**。

## 分析と出典の境界

- 2本の凍結 GBIF artifact ZIP SHA256：
  - discovery: 2181898c7428172cd039b208a5e3fc5169c12d3a8d324bef05a84f3a6a89a446。
  - disjoint holdout: ec1df9eb35e856b1020e92a755aff8be8c121b640602e8e3191d528cc531df1a。
- 両標本合算のGBIF返却座標付き occurrence 件数4,935,142。これはAPI集計値であり全個別occurrence行をダウンロードしたことを意味しない。
- 分析単位はaccepted種×大陸にPRESENTかどうか。南極大陸は2種のみの記録のため今回の仮説検定から除いた。
- 分類群の受容種名と写真9色代表は一致しているが**真の花被ROIが全写真について植物学的に確認できたわけではない**。写真色はGBIF記録地点の局所色ではなく、種内固定色でもない。GBIF記録は在来／移入statusが未確認で、努力偏り・属/系統構成の交絡がある。
- 6地域・2独立標本の12検定全体のBH補正も別カラムに保存。両標本で同じ色方向が支持された地域は存在しない。
- 全37,104受容種のGBIF ISEA3H8 species×gridは未取得。認証付きGBIF SQLが必要であり、この360種の大陸facetを全球高解像度格子と取り違えない。

## 追跡可能なコードと成果物

- scripts/analysis/mois_gbif_continent_two_cohort_qa_20261115.py
- tests/test_mois_gbif_continent_two_cohort_qa_20261115.py
- .github/workflows/mois-gbif-continent-two-cohort-replication-20261115.yml
- CI artifact #14160322719：source_360_taxa_x_six_GBIF_continent_presence.csv、12_region_cohort_9hue_tests_permutation_BH.csv、receipt.json。
