# MOIS：写真9色×GBIF実分布360種の独立再現検証（2026-11-15）

## 結論

**南米では橙色写真ラベルの種が紫色より多く記録されるという探索的なパターンは、別の180種では再現しなかった。** 検証標本は最初の180種とは受容speciesKey重複0で、同一の事前ハッシュ順位から9色各20種を独立抽出した。

| 南米にGBIF記録あり | 最初の180種（探索・事後） | 独立180種（事前検証） |
|---|---:|---:|
| orange photo hue | 14/20 | **8/20** |
| purple photo hue | 2/20 | **8/20** |
| 9色×南米presence置換検定 | p=0.0039 | **p=0.8569** |
| orange＞purpleの片側Fisher検定 | 事後 | **p=0.6262** |
| GBIF南米記録10件以上の片側Fisher | 感度探索 | **p=0.3705** |
| GBIF記録数と写真色状態調整9色GLM-LRT | 感度探索 | **p=0.8255** |

従って、最初の南米の方向性を普遍的な花色適応・生物地理法則とは主張しない。これは明確な否定的再現結果。

## 実データと完全照合の範囲

- GBIF discovery 180種、Actions #43707058183、artifact #14158822487、ZIP SHA256 2181898c7428172cd039b208a5e3fc5169c12d3a8d324bef05a84f3a6a89a446。
- GBIF holdout 180種、Actions #43708230511、artifact #14159523021、ZIP SHA256 ec1df9eb35e856b1020e92a755aff8be8c121b640602e8e3191d528cc531df1a。
- 全写真9色・受容37,104種の厳密対応ZIP SHA256 eff371f8da3a70d23371a1b78236e089da1206b82ec9fb68ef860926e2965cf4。再結合で全360種のspeciesKey・inat分類群ID・modal9色・photo evidence status一致、重複0。
- 対象は9色×各40種、計360種。GBIF公式 occurrence/search の公開country/continent facet、hasCoordinate true / hasGeospatialIssue false / occurrenceStatus PRESENT。**APIから得た位置情報付きoccurrence件数の合計4,935,142、種×国組合せ3,424**。これらはGBIF側の集計値であり、約493万件の個別occurrence行をダウンロードしていない。
- sample中のdistinct generaは297、複数種を含む属に所属する種は107。系統構造は未制御、在来／導入を区別していない。受容学名がGBIFで正確に一致しても、写真由来の種レベル色をGBIF点の実測色へ置換してはならない。

## 360種の写真色の証拠品質

| photo evidence status | Species |
|---|---:|
| provisional_single | 205 |
| multicolour_candidate | 93 |
| sampled_multicolour | 59 |
| sampling_consensus_single | **3** |
| **total** | 360 |

写真modal9色とmean9最大アンカーの不一致**59種**。全360種における独立した真の花弁色認証は未実施。多数決と平均アンカーは異なる要約統計量なので不一致を一律のバグと断定しない。**360種中で厳密な単色写真コンセンサスは3種だけ**で、遺伝的に固定した花色に基づく進化解析と断言できない。

## 記録国数と写真白色（探索的分析）

| 標本 | N | white vs nonwhite係数 | 95% HC3区間 | p |
|---|---:|---:|---|---:|
| discovery | 180 | +0.0631 | −0.3021〜+0.4282 | 0.7349 |
| holdout | 180 | +0.4941 | +0.1207〜+0.8675 | 0.00949 |
| 360種（cohort固定効果追加） | 360 | +0.2675 | +0.0073〜+0.5277 | 0.0439 |

応答はlog(1+GBIF記録国数)、共変量はlog(1+GBIF座標付記録数)。九色全体のKruskal–Wallisは discovery p=0.7983、holdout p=0.1549、統合 p=0.1665。白色係数は標本間で一致しておらず、統合p<0.05のみをもって「全球白花植物は広域分布」とは結論できない。属単位クラスタ頑健標準誤差でも統合p=0.037だが、植物の真の花色・在来分布・系統とGBIF努力偏りが未検証であることは変わらない。

## 全37,104受容種の全球種×格子解析の未解決条件

正式GBIF SQL集計（ISEA3H解像度8）は、37,104種を149分割（250種×148件＋104種×1件）で準備し、speciesKey重複0を再確認。GBIF公式SQL validator HTTP201を得た：[Actions #43708178817](https://github.com/zuizui0223/mois/actions/runs/43708178817) artifact #14159013597、ZIP SHA256 8c642ab261e6b0498e33f7f305946a89af610d9c2c75c07cfda7c81207e3571b。

**GBIF_USERNAME（またはGBIF_USER）とGBIF_PASSWORDのGitHub Secretsが未設定**。従って認証付きSQL送信0、download key0、DOI0、全球37,104種×セル記録取得0。国別facetの360種は、全球のセル×気候・標高・土壌解析の代用ではない。本人が認証情報を安全に設定した後、初回250種を一度だけ提出して真正GBIFダウンロードとDOIを監査し、残り148件の取得を段階的に拡大する。パスワードは報告書やチャットに記録しない。

## 再現性とプログラム

- scripts/analysis/mois_gbif_colour_country_facets_20261115.py
- tests/test_mois_gbif_colour_country_facets_20261115.py
- .github/workflows/mois-gbif-180-colour-country-20261115.yml
- .github/workflows/mois-gbif-180-colour-country-holdout-20261115.yml
- .github/workflows/mois-gbif-colxr-37104-global-cube-auth-gate-20261115.yml
- 対話に添付した360種統合QAは、二つの真正180-species source ZIPと全37,104 species photo trait ZIPを各SHA256で照合し再集計。結果には全360種の学名・色・国数・生態的品質フラグ、南米の両標本色別集計、白色係数の感度解析と再実行コードを保持。

**主要な再現判定：南米のorange対purpleの地理パターンは支持されず、他の写真9色×GBIF出現国数の関連も全球の生態学的機構や在来分布の証明には未到達。**


### 認証設定後の初回250種の安全な入口

149 SQLの完全な署名付き原本から先頭250受容種を、**ユーザーがGitHub ActionsのGBIF_USERNAME（もしくはGBIF_USER）とGBIF_PASSWORDをRepository Secretsへ設定した後にだけ**一回提出する専用のone-shot workflowを、未マージPR #200へ追加した：

- .github/workflows/mois-gbif-37104-first250-on-demand-after-secrets.yml
- 直前の2026-11-15認証ゲートのreceiptはGBIF_CREDENTIALS_NOT_CONFIGURED、正式送信0、downloadKey0、DOI0。パスワード値はGitHubのSecretsでのみ受け取る。直接チャット・Git・結果CSVへ保存しない。
- **GitHubのworkflow_dispatchは対象ワークフローがdefault branchに存在しないと起動できない**。本workflowは現時点で未マージ解析ブランチにあるので、実際の手動実行には適切なレビュー／default branch反映が別途必要。今回も実送信は一切していない。
- SQL POSTが受理されてもDOIとgrid ZIPは非同期で後から作られる。HTTP 201とキーの受信のみで全球格子取得「完了」と判定しない。
