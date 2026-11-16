# MOIS 2026-11-16: LGM気候履歴・レフュージア・生産力の事前解析設計

## 現在地と結論（結果ではなく検証設計）

MOISの現在の対象は全球の写真ベース9色組成と水分環境。既存の真正GBIF公開スナップショットに対応する受容種37,104種、種×EPSG:6933の100 km格子5,493,931組（有効出現約3.961億件）を継続利用する。11/16のAGB/NDVI/WorldClimの完全結合は36,753種。主たる解析単位は**種1行**であり、種×格子を独立種として複製しない。

既存の木質バイオマス解析で白写真色vs非白の環境AGB中央値は43.384 vs 37.298 Mg/haだが、属等重み・気候調整の白効果はβ=-0.00564（95%CI -0.02229〜0.01102, p=0.5071）。RFの全群AUC 0.55581、AGB除去0.55489。**木質AGBが強い独立予測因子である証拠はない。** AGBは木質地上部量であり、草本/花の実測バイオマスや繁殖投資ではない。

**LGM・レフュージアと花色との関係は未計算。ここに有意係数を報告してはならない。**

### 事前固定した問い

- H1（主検定、無方向）：現在の気候・木質バイオマス・NDVIと記録努力を調整しても、現在のGBIF分布域におけるLGMからの気温差・降水差は写真9色の種ラベルと条件付きで関連するか。
- H2（補助）：地域内／属内で、NPP（純一次生産量）が木質AGB/NDVIと異なる情報を提供するか。NPPはフラックス、AGBはストック。NPPも個体の花への配分の証拠ではない。
- H3（探索）：長期の気候安定性とNPPに交互作用があるか。方向・効果量は未指定。
- H4（第二段階、現データだけでは未検証）：気候安定地、氷期残存域、氷期後の再定着域で花色・系統組成は異なるか。「現在地点の低いΔ気候」≠「その種のレフュージア」。

## Stage 1: 同じ全球格子への古気候抽出

1. 現在とLGMの整合性を優先し、**WorldClim v1.4の1960–1990年基準**（10 arcmin BIO1/BIO12）と、WorldClim v1.4にバイアス補正されたLGM（約22 ka）CCSM4, MIROC-ESM, MPI-ESM-P の10 arcmin BIO1/BIO12を対応させる。現代WorldClim v2.1 1970–2000の値とLGM v1.4を直接差し引かない。
2. 既存の重複なし**18,577個のGBIF 100 km格子**に25点（5×5）の陸上サブサンプルを適用。現在・過去で対応する陸上点のみに基づく、気温差（°C）、|気温差|、log1p降水差、|log1p降水差|を出力。温度デシ℃のスケール（×0.1）、降水mm、nodata、CRS、GeoTIFF/BILを原ファイルで検査。
3. GCMごとに差分を独立保存し、そのモデル間中央値／範囲を計算。観測値ゼロと無効値（氷床・海・欠測）を混同しない。LGM時の氷床・旧海岸線・当時の地上部が未確定のセルは**レフュージア判定不可**のまま保存。
4. 現在の種×格子出現を元に、各受容種の「現在占有する格子における気候変化の平均」を記録。分布範囲の大きい種が反復測定として過重になるのを防ぎ、各種につき1行を保つ。これは祖先生育地や移住ルートの復元ではない。
5. レシートに原ファイルURL・SHA256・気候モデル・時代・単位・欠測理由・カバー率・プログラムのcommitを記録する。

参照元：
- https://worldclim.org/data/v1.4/worldclim14.html
- https://worldclim.org/data/v1.4/paleo1.4.html
- https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip
- https://geodata.ucdavis.edu/climate/cmip5/lgm/cclgmbi_10m.zip
- https://geodata.ucdavis.edu/climate/cmip5/lgm/mrlgmbi_10m.zip
- https://geodata.ucdavis.edu/climate/cmip5/lgm/melgmbi_10m.zip

## Stage 2: 生産力と時系列気候安定性

- MODIS MOD17A3HGF v6.1、年間NPP（kg C m^-2 yr^-1）を年別で同じ格子に収集し、2001–2024など品質が同等な年を固定して平均と経年変動を集約。異なる生物群系／土地利用で感度分析。年NPPは植生生産量フラックスで、写真の花のバイオマスや受粉成功ではない。
- ESA CCI v7 AGB（Mg dry woody mass ha^-1）、NDVI、NPPを同じ「生産力」係数に混ぜない。相関、群アブレーション、森林/非森林の層別で独立情報を比較。
- CHELSA-TraCE21kの約21 ka〜現在100年刻みのBIO1/BIO12は、LGM-現在の2時点差では測れない時間的安定性を要約する第二段階で使用（モデル・時代・季節データの整合確認が必須）。
- 氷期レフュージアの主張は、氷床履歴、当時の海陸配置、種/系統ごとの過去生育適地や化石・花粉・遺伝地理などを別途照合するまで保留する。

参照元：
- https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD17A3HGF
- https://climate.esa.int/en/catalogue/6429d1aafe1e43b9b414e4a5a7f8b903/
- https://www.chelsa-climate.org/datasets/chelsa-trace21k-centennial-bioclim
- https://pubmed.ncbi.nlm.nih.gov/32714360/ （*Campanula americana*、温度と氷期後の集団史に応じた花色変異の先例）
- https://pubmed.ncbi.nlm.nih.gov/21979937/ （LGM以来の気候速度・固有性の先行研究）

## 解析単位・モデル・採否条件

主モデル：写真white vs other の種レベル二項応答、同属内比較（属で中心化し属クラスタSE）、現在BIO1、BIO12、緯度、NDVI、AGB、log(1+GBIF記録数)、色evidence statusを調整。説明変数群に3GCMの歴史気温・降水差を加え、条件付き増分を検定。属等重み、属ホールドアウト、地域ホールドアウトを区別。9色の全構成を別途多項応答で評価し、「白だけ」から普遍的適応則を導かない。

**事前ゲート：**
- 36,753種を母集団とした完全ケース数と欠測数を明示。被覆率閾値50/75/90%を感度分析。
- 3つのGCMで符号と効果量を独立評価し、恣意的な有利モデル選択を避ける。
- 高品質単色写真根拠（標本内3独立観測一致等）による部分集合でも再現性を検査。真正な花器官ROI検証なしに植物学的花色と呼ばない。
- GBIFの自生・外来・栽培混在を可能な範囲で感度解析。現在分布への種色投影は地点固有の花色測定ではない。
- 結果を報告する際は、係数・区間・実際のN・QA除外数・複数検定処理・原データSHA・現在と過去の期間を添える。
- 気候史の有意性が得られても、選択・定着フィルター・単なる系統的保存の因果判別にはならない。

## 継続位置

既存のsigned全GBIF・9色・生態環境は
- `docs/MOIS_GLOBAL_37104_GBIF_PHOTO9_RANGE_CLIMATE_20261116.md`
- `docs/MOIS_GLOBAL_PHOTO9_MOISTURE_NDVI_BIOMASS_MULTIMETHOD_20261116.md`
- `scripts/analysis/mois_global_worldclim_hue_20261116.py`
にある。現時点では**設計の保存のみ**。LGM/NPP取得・モデル実行が完了したと書いてはならない。
