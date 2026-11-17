# MOIS｜花色を主線とする Fourth-corner JSDM：真正GBIFでの地域別事前監査（2026-11-17）

## 結論・いま進めた工程

JSDMは、MOISの**植物種の現在GBIF記録分布（地域×種）**を応答として、**花色（種の形質）×現在環境／NPP／LGM古気候**のFourth-corner係数を推定するために利用できる。花色を直接予測するRF・GAMとは推定対象が異なるが、**生物地理区ごとの花色組成**を問う本研究に役立つ。

まず巨大な37,104種×18,577格子の一括モデルを行わず、事前に花色を見ずに指定した3地域の真正GBIF記録行列と、WCVP地域自生可否・環境データ・種別花色ラベルの3表を作った。

**GitHub Actions [run #44002905583](https://github.com/zuizui0223/mois/actions/runs/44002905583) SUCCESS**、[artifact #14293695159](https://github.com/zuizui0223/mois/actions/runs/44002905583/artifacts/14293695159)、署名ZIP SHA256 `c07239b2890811b8400984a5657563f07ff6545d26ffdb719205779216bfa586`。追加の**[run #44003290846](https://github.com/zuizui0223/mois/actions/runs/44003290846) SUCCESS**（[artifact #14292573380](https://github.com/zuizui0223/mois/actions/runs/44003290846/artifacts/14292573380)、ZIP SHA256 `e9136c8c79a55986ebb6034273027554c02aaf391c6cf60e182942a792e8edb6`）では、GBIF記録努力の補助量（その格子内の全対象種GBIF occurrence計数・記録種数）をX表へ付与して3地域の結果を再実行した。どちらも画像色ラベルによる対象選択は行わない。解析の主線は**花色**であり、属・Jaccardを主目的変数に変更しない。属の固定効果・重み・選抜は今回のパイロットに使わない。

実装 `scripts/analysis/mois_jsdm_fourthcorner_gbif_native_pilot_20261117.py`、小規模テスト `tests/test_mois_jsdm_fourthcorner_gbif_native_pilot_20261117.py`、Actions `.github/workflows/mois-JSDM-GBIF-fourthcorner-native-region-preflight-20261117.yml`。

## 実データ前処理結果

| 花色を見る前に固定した地域 | 100 km格子 | 植物種 | 500 km空間ブロック | WCVP自生かつGBIF記録ありの種×格子 | 自生範囲内だがGBIF未記録の種×格子 |
|---|---:|---:|---:|---:|---:|
| 北米温帯（25–55°N, 125–65°W） | 300 | 153 | 55 | 5,683 | 9,125 |
| 欧州温帯（35–65°N, 15°W–40°E） | 300 | 170 | 50 | 11,168 | 11,610 |
| 東アジア温帯（25–55°N, 95–145°E） | 300 | 173 | 47 | 5,328 | 16,914 |

追加でWCVPその地域の自生範囲にないと記載された種×格子の組み合わせはY=NAにして、**勝手に不在0へ変換しない**。地域は事前の経緯度で固定し、サイトは固定SHA順、種は**白／有色ではなく自生GBIF記録支持数＋固定SHA順**で選択。学習対象の種内に白・有色の割合の偏りは残る。種別色が複数種に共通したり、反復観察が不足する点も保持する。

| 事前パネル | 写真whiteラベル種 | 3観察以上の単色一致 status |
|---|---:|---:|
| 北米 | 93／153 | 13 |
| 欧州 | 93／170 | 14 |
| 東アジア | 106／173 | 16 |

従って**独立した花被色の真正な外部検証が済んだ種は0**。写真白ラベルは種に与えたproxyであり、GBIF記録地点で撮影された花の実測値ではない。

## JSDM本体で何を推定するか

解析対象は、サイトi、植物種jのGBIF出現**記録**の有無。花色t_jは種別ラベル。現在気候・水分・日射・地形（ソースがある場合）・NPP・木質AGB・LGM現在差をx_i、空間・未測定共変動の潜在因子をu_iとする。

`logit[P(Y_ij recorded)] = alpha_j + beta^T x_i + gamma^T(x_i × colour_j) + u_i^T lambda_j + bias_i`

ここで主に見る `gamma_{white×NPP}`、`gamma_{white×LGM warming}` は、環境によって**白色ラベル種が相対的にどのくらいGBIFに記録されるか**を表す。真の花弁色の選択強度や真の占有確率ではない。モデル自体の応答は花色ではなく記録出現であるため、並行して花色を目的にしたRF/GAM/GBIF格子白率モデルを残す。

- 主候補：R `gllvm` Fourth-corner GLLVM。Xは現在気候・NPP・LGM履歴とサイトの記録努力代理変数。TRはwhite vs nonwhite、9色のone-hot（orangeの反復不足で過度な係数解釈は不可）。低次潜在因子1–2とrow effectを試す。公式ドキュメント：<https://cloud.r-project.org/web/packages/gllvm/vignettes/vignette1.html>。
- 補助：`Hmsc`（種形質の階層効果・必要なら系統情報）、`sjSDM`（大規模計算）。それぞれの主効果を同じものとして解釈しない。
- ベースライン：花色を目的とする既存RF/HGB/GAM／正則化ロジスティック回帰。JSDMのみが花色を直接測定できるわけではない。

## 実モデルの最初の診断：gllvmのTMB初期化エラーとNA保持

- 原本照合と同一の300格子・170候補種の欧州パネルからGBIF記録JSDMを試験する[GitHub Actions run #44003711557](https://github.com/zuizui0223/mois/actions/runs/44003711557)を実行した。公式gllvm 2.xのインストールと元Y/X/TR・GBIF記録努力のSHA検査は**PASS**。
- モデルの最初の`gllvm(..., TR = white_proxy, num.lv = 1, family = "binomial", method = "VA")`実行時、TMBの`MakeADFun`内で**Segmentation fault（exit code 139）**。**係数・予測・適合度は得ていない**。GBIFの真正な不在データがないまま非自生域NAを0で埋める処置はしていない。
- エラーが入力のNAパターンか、パッケージ／R/TMBの組合せか、第四コーナーTMBの実装かは未確定。候補ごとに切り分ける必要がある。まず**NA応答を明示的に扱うHmsc**へ移し、`scripts/analysis/mois_hmsc_fourthcorner_NA_mask_smoketest_20261117.R`で少数反復の動作のみ検証する。5サンプル＋5burn-inで済んでも**収束や独立した花色×環境係数の検証とは呼ばない**。

## 検証を進めるゲート

1. **観測過程**：GBIFで記録がないことは、現地で調査して不在と確認されたことではない。XへGBIF記録数・記録種数を観測努力の粗いproxyとして加えても、本来の検出確率や非記録の理由は識別できない。したがって本パイロットでは「GBIF**記録**JSDM」と表記し、占有・定着JSDMとして断定しない。
2. **同じ格子／種の事前固定パネル**で現在気候のみ → +NPP/AGB/NDVI → +LGMとの差を比較。標準化と補完は学習分のみで計算する。LGMは現在分布格子の約22ka推定気候であり実際の古代分布／レフュージアではない。
3. **500km空間ブロック単位holdout**で、未学習格子のGBIF「記録あり vs 可能域内未記録」の予測log loss、PR-AUC、Brier scoreを比較する。未記録0はpseudoabsenceであることを毎回明記。サイト潜在因子を検証時のYから推定してしまわず、新サイト予測は潜在因子0／母集団平均で評価してリークを防ぐ。
4. 花色の誤ラベル・写真evidenceのstatus群で感度分析。先行の3観察以上コンセンサスは3パネル各13–16種しかなく、**その小標本だけのJSDMは不安定**。新しい花被直接画像認証標本か外部の種形質資料が要る。
5. 検証結果がRFと一致しなくても誤りとは限らない。RFは花色の種別予測、JSDMは出現記録の種別多変量予測と形質×環境の相関で、推定対象が異なる。残差相関を種間相互作用の因果証拠にしない。

## リソースを節約するルール

37,104種を全18,577格子へゼロ充填した上で標準gllvm/Hmscへ投入しない。生データの5,493,931種×格子の**記録あり**は8.57億でなく約550万であり、密な種×格子配列の多数の0は真の非出現ではない。まず固定3パネルの小規模・明示的なpseudoabsence挙動を検証する。次に生物地理区ごとの学習・空間holdoutを拡張し、個別の草本・島嶼性・送粉シンドロームの問いは後段で処理する。

**現在実際に達成したのは真のsourceからのJSDM入出力行列の生成と交絡チェック用の前処理**。次段として`scripts/analysis/mois_gllvm_photo9_fourthcorner_recording_pilot_20261117.R` と `.github/workflows/mois-GLLVM-real-fourthcorner-GBIF-native-pilot-20261117.yml`を追加。欧州の300格子・500km空間ブロックで、公式gllvmの1潜在因子モデル（現在環境のみ vs NPP＋LGMと白色交互作用）を試験実行する。**別runのSUCCESSが確認されるまでは、JSDMを適合して係数が出たという主張はしない**。
