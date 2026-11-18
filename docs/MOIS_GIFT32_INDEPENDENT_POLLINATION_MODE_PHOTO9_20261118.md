# MOIS：独立送粉様式GIFT 3.2×全球写真9色データの生態学的照合（2026-11-18）

## 何を実行したか

既存island研究の**独立したGIFT v3.2 direct documented pollen-vector-mode表（5,771種、biotic 4,538／abiotic_wind 1,233）**を、MOISの真正**37,104 COLXR accepted種×写真modal9色×現在GBIF分布域環境**へ、受容学名の**完全一致のみ**で照合した。異名・Fuzzy name・分類群幅の違う結果を無理に転用していない。原MOIS署名付きCSV adapterとGIFTの元ZIP・内部CSVをSHA256で照合してから解析。入力種ラベルの変更・写真原画像の再計測はなし。

- 元island独立送粉モード： [Actions run #34803493307](https://github.com/zuizui0223/island/actions/runs/34803493307)、artifact #10332810251、ZIP SHA256 `074a52323ba835201043a3a5980cc2b2d0ff79878e83bc9be19511ae83a1c0e8`。内部mode CSV SHA256 `cddb8bb02feede8fbb45e80bc62cb25384a8899d05312843eef0d56dee57ab65`。
- 元MOIS source 37,104 accepted speciesのParquet ZIP SHA256 `d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc`、内部Parquet SHA256 `e772092266ebe603c26414e30a9f31f09990c2cfc608b590bb80a1e5c7ff306b`。
- 今回のエクスポート： [MOIS Actions #44267638254](https://github.com/zuizui0223/mois/actions/runs/44267638254) **SUCCESS**、artifact #14411131896、ZIP SHA256 `d024e4f497c67c0eb3073244dc84189e50fe1c22deec2231f27adef8ab5d7035`、内部CSV SHA256 `af74536f54317f3fc923b8f5f00eb61f5cba3803a19b9be021b073a7a1ea2d46`。

## 確認された実データの対応率

| 指標 | 真正値 |
|---|---:|
| 元GIFT pollen-vector mode記載種 | **5,771** |
| MOIS写真9色分類済み受容種 | **37,104** |
| Exact binomialで双方に存在する種 | **3,196** |
| そのうち動物媒 `biotic` | **2,756** |
| そのうち風媒 `abiotic_wind` | **440** |
| 実際に重複する属 | **1,328** |
| **同じ属に動物媒/風媒の両方が載る属** | **11属** |
| その11属に含まれる一致種 | **54種** |

一致した種を1種1票にすると：

| GIFTの独立送粉様式 | MOIS写真白ラベル | MOIS写真褐銅ラベル |
|---|---:|---:|
| biotic | 1,711/2,756 = **62.08%** | 162/2,756 = **5.88%** |
| abiotic_wind | 314/440 = **71.36%** | 80/440 = **18.18%** |

風媒−動物媒の写真白ラベル差は+9.28ポイント。**属を単位に再標本化**する3,000回bootstrapでも動物媒−風媒の差の95% percentile区間は **−15.55〜−3.20 pp**、写真褐銅ラベルは **−17.54〜−6.81 pp**。しかしこれは、**両様式が異なる属に集中している事実を統計的に取り除いた区間ではない**。同属内の両様式比較は11属・54種しかないので、属固定効果での高精度な推定・機構証明をしていない。単純な9×2の種独立χ²検定も計算されているが、**属依存を無視したp値を科学的主張に使用しない**。

## 生態学的に何が言えるか

**写真花色の全球的な地理構造の一部が、独立の動物媒・風媒植物種プールの違いと整合する**という技術的な出発点にはなる。しかし、
1. 元GIFTの`biotic`は**動物媒全般**。マルハナバチ・蝶・鳥などの**独立した訪花ギルド区分ではない**。
2. 風媒植物の写真白/褐銅ラベルは、実際の花冠色でなく葉・樹皮・花序などの器官差を拾った可能性が高く、写真9色は真の花被色未検証。
3. 送粉様式と属は非常に共変し、**11属のみが属内比較可能**。全属へ汎化する独立の送粉モード効果は識別できない。
4. 送粉様式は種の属性であり、**現在の同一100km地域における送粉者の種組成・来訪頻度・送粉成功率は計測していない**。
5. GIFTは全37,104種のうち約8.61%しか一致せず、GIFT調査の対象選択が非無作為である可能性が高い。多様な写真色・地域への欠損補完や、記載がない場合を風媒とみなすことを行わない。

この独立送粉様式とMOIS色の対応は、ユーザーの生物間相互作用仮説を検証する第一段だが、**現状では分類群構成の独立した因果識別には届いていない**。次の主線には花から**実際に訪花・花粉運搬が記録された動物ギルド**の一次資料（例：DoPI v2の678kレコード、欧州限定）を別結合する必要がある。

## 再現性

エクスポート元は今PRの`scripts/analysis/mois_export_photo9_species_environment_for_independent_pollination_20261118.py`とActions workflow。実クロスウォーク解析スクリプト`mois_gift_photo9_mode_analysis.py`は会話に同梱した検証可能な別ZIPに保存（第三者GIFTの個別種ラベル自体はこのPRで再配布しない）。3,196種の全データを隠したまま、aggregate countsと11属の検証表・計算記録を配布。元GIFT版の識別はIslandの元artifactで担保する。

出典：Weigelt, König & Kreft (2020), *Journal of Biogeography*, doi:10.1111/jbi.13623。Denelle, Weigelt & Kreft (2023), *Methods in Ecology and Evolution*, doi:10.1111/2041-210X.14213。GIFTデータソースは公開・非導出参考文献のみ利用する原`island` H5cが定めた条件に従い、個別元文献を論文化時に引用する。

**このbranchはMOISのPR #206（地理勾配・花色集中）から独立。mainへ自動マージしていない。**
