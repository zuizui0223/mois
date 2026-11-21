# MOIS 2026-11-16 — 全37,104種の写真9色 × 真正GBIF分布 × 気温・降水

## 1. 全球GBIF取得・監査は完了

[GBIF匿名公開S3全量の16区画＋merge（run #43723276949）](https://github.com/zuizui0223/mois/actions/runs/43723276949)：全17ジョブSUCCESS。元10,346 Parquet（圧縮271.778GB）をobject path/ETag/bytesでexactly-once確認した。GBIF・AWSアカウントは不要。出典：[GBIF Occurrence Snapshot 2026-11-01](https://doi.org/10.15468/dl.jeesns)。

| 実数 | 値 |
|---|---:|
| 全GBIF occurrence走査 | 3,887,364,403行 |
| COL XR受容37,104種に対応し座標QA後に残った出現記録 | **396,142,825行** |
| species×EPSG:6933 100km等面積セル | **5,493,931組** |
| 少なくとも1セルの有効記録がある受容種 | **36,758種** |
| 今回の品質条件で記録0（自然分布不在ではない） | 346種 |

元S3全量から得たZIP：[artifact #14166701658](https://github.com/zuizui0223/mois/actions/runs/43723276949/artifacts/14166701658)、SHA256 422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2。格子はEPSG6933の矩形で、以前提案したISEA3H r8の六角格子とは異なる。

## 2. 全球の「白色種ほど広域分布」仮説を検定

[globally signed source & genus FE run #43790371870](https://github.com/zuizui0223/mois/actions/runs/43790371870) SUCCESS、[artifact #14192037667](https://github.com/zuizui0223/mois/actions/runs/43790371870/artifacts/14192037667)、ZIP SHA256 63a260c0a9c68d5c7198c6f4a821bd391eae25954a31cf060e6b1fbac544e313。

- 白写真modal種の占有100kmセル中央値48、非白37。
- log(1+セル数) ~ 写真白 vs 非白＋log GBIF記録数＋写真evidence status、属単位cluster：β +0.03745、p=0.0147。
- **同属固定効果を加えるとβ −0.00558、95%CI [−0.02639,+0.01524]、p=0.5996**。白・非白の種を含む2,376属・29,544種を等属重み：β −0.00898、p=0.506。
- 写真modal hue とmean9最大成分の不一致種を除いても、同属等重み：β +0.00136、p=0.9236。
- 9色全体の属内Wald χ²=22.34、8df、p=0.00432。ただし白対ほかの8色を別々に比べ、BH補正q<0.05の色は**0**。細かな有意性を適応と誤認しない。
- 種の記録された100kmセルの平均絶対緯度：同属白効果 −0.029度、p=0.8317、同属等重み+0.062度、p=0.6189。138個の高緯度セルでは中心がprojection外で、緯度要約のみに含めず元分布セルは維持。

**観測上の結論：白写真色の単純な記録分布幅差は同属種比較で確認できない。** 属FEは完全な系統補正ではない。GBIF記録数は努力量と実分布の両方を反映する。

## 3. WorldClim 2.1（1970〜2000年）を全球GBIF分布へ接続

[WorldClim公式BIO1/12](https://worldclim.org/data/worldclim21.html) の10 arc-minute基準ラスタを取得し、各100km矩形セル内25点（5×5）の陸上気候値から温度・降水の平均値を推定した。世界の異なる年のGBIF出現記録・写真日付と同時の気象ではない。

[世界全種×WorldClim run #43791207244](https://github.com/zuizui0223/mois/actions/runs/43791207244) SUCCESS、[artifact #14192392368](https://github.com/zuizui0223/mois/actions/runs/43791207244/artifacts/14192392368)、ZIP SHA256 82f06f69043a2080a961d5586ab4ecc2342b7d7d739c4ac819147971f4d8e069。WorldClim原ZIP SHA256 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf。

- 全格子18,577、気候25点すべて欠測（海上/沿岸など）2,859。
- 少なくとも有効気候を得た種36,755、格子の75%以上で有効な種**36,704**。
- 同属固定効果・属クラスタ分散、GBIF記録数、平均絶対緯度、写真evidence status、および他方の気候変数で調整。

| 白写真色vsその他：GBIF記録セルの気候 | β | 95% CI | p |
|---|---:|---|---:|
| BIO1年平均気温（°C） | **+0.09445** | +0.00676〜+0.18214 | 0.0348 |
| log(1+BIO12年降水mm) | **+0.01519** | +0.00011〜+0.03027 | 0.0483 |
| BIO1、属を等重みにした感度 | +0.05268 | −0.03018〜+0.13554 | **0.2128** |
| log BIO12、属等重み | +0.01185 | −0.00183〜+0.02553 | **0.0896** |

最初の2項目はBH q約0.0483だが、**推定差は非常に小さく、属を等しく数えると有意性が消える。** 花色適応の証拠としては不十分。

[16条件のWorldClim感度検定 run #43791652492](https://github.com/zuizui0223/mois/actions/runs/43791652492) SUCCESS、artifact #14192482999、ZIP SHA256 d0ff92fbee5f8965d8ed38ca5a40c21e3de9a62980a4c3eafe6dbe6b4847eedb。被覆率閾値50/75/90/100%、写真modalとmean9矛盾種の除外の有無で再検定した。

- 温度の差は+0.09〜0.12°Cの小さな方向を維持。ただし属等重みでは非有意。
- 降水log係数+0.013〜+0.015は被覆率100%でp=0.110、写真色不一致種除外の75%被覆でp=0.067へ変化。**水分に関するp<0.05は感度条件に頑健でない。**

## 4. 非常に重要な制約

今回の37,104種はGBIF COL XR exact accepted nameであり、**全部の本当の花弁色が認証された種ではない**。実際のFCP ROI-v4マスク点検では葉・背景を花色と誤測定した事例がある。受容種の写真色単色コンセンサスは2,686種だけで、その他は暫定単色・多色候補・反復複数色。写真modal9と平均9色最大成分不一致は2,466種。

GBIF出現記録は自生・外来・栽培を区別していない。撮影種の色をGBIF記録各セルに投影することは、**セル内の実測花色を意味しない**。WorldClimの1970〜2000年の平年値とGBIF記録・写真撮影年は一致しない。海岸・小島の100kmセルへの気候割り付けは特に偏り得る。

ここまでの信頼できる成果は、**真正GBIF全球分布×写真色技術分類×気候を全種へ接続し、同属比較・空間品質・分類色の不確かさを保って解析した**こと。選択・適応・花色進化の因果的結論はまだ支持されない。

## 再現ソース

- scripts/analysis/mois_global_photo9_grid_latitude_20261116.py
- tests/test_mois_global_photo9_grid_latitude_20261116.py
- .github/workflows/mois-global-37104-photo9-gbif-range-latitude-20261116.yml
- scripts/analysis/mois_global_worldclim_hue_20261116.py
- tests/test_mois_global_worldclim_hue_20261116.py
- .github/workflows/mois-global-37104-GBIF-WorldClim-BIO1-BIO12-20261116.yml
- scripts/analysis/mois_worldclim_species_hue_sensitivity_20261116.py
- .github/workflows/mois-37104-worldclim-coverage-sensitivity-20261116.yml
- 全公開GBIF匿名処理：docs/MOIS_GBIF_NOAUTH_PUBLIC_S3_GLOBAL_37104_20261115.md
