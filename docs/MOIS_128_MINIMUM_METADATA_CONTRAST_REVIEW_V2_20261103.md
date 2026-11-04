# MOIS：128写真の花被・背景コントラスト評価v2（2026-11-03）

## 実行済み

[Actions run 41562692689](https://github.com/zuizui0223/mois/actions/runs/41562692689) SUCCESS。元v1と新v2の合計18テスト、両HTMLのJavaScript構文、128件の元写真IDと空欄人手評価のSHAを検証した。

[評価者A/BのHTMLと監査台帳（artifact 13197751320）](https://github.com/zuizui0223/mois/actions/runs/41562692689/artifacts/13197751320)。ZIP SHA256：5a4beaaaf137444b1b06c7c480356f384a301d5c3fcd75af39de3c8b2d058598。

元128写真roster SHA256：5fef860e4d2d0ff4b2dfec850253e4a08c8dc87878647d6712bc61e4b75dcb56。A HTML SHA256：60db1819c12156d6dffd84cfb5af3a6391affa258c066e112a27c82ea13d96c9。B HTML SHA256：1ec8c493c883670ce5ba27eee2617fca1d2009df97737b2bdcd4a0f2b0ad6291。

v2では評価者HTML内の埋込rosterを写真IDのみとした。元観察UUID、分類群ID、地域、現地季節、取得系列、検出BOXやROI失敗の事前情報はHTMLと書き出すCSVに含めない。評価者に元の写真ID自体を隠せないため完全盲検ではない。既存のv1は保存し、GitHub main・ROI9結果を変更しない。

## 実際の評価対象分母（128枚）

緯度帯：N_arctic 8、N_subpolar 14、N_temperate 21、N_tropical 26、S_temperate 21、S_tropical 30、S_subpolar 4、unknown 4。現地季節：summer 43、autumn 30、spring 26、winter 25、unknown 4。取得設計：historical 51、P2-first 32、bulk128 25、bias-refill 20。技術BOX事前診断：候補 111、陰性 14、技術エラー 3。

このrosterには島／大陸の区分がない。北極は8写真しかなく、128全体も意図的に元strict9失敗から選択された。地域差・世界植物相の低顕示率をここから推定しない。

## private人手評価・再検証

2名が独立にローカルHTMLを使い、手動で一致する元写真ファイルを読み込む。実際の開花・花被・局所背景を確認し、全て肯定なら各3点を選択。人手CSV・原画像は公開GitHubへアップロードせず、研究調整者側でのみ保管する。

コーディネーターのローカル実行コマンド：

    python -m scripts.analysis.mois_petal_background_blind_gate_v2 validate-private --roster ORIGINAL_128_FROZEN_ROSTER.csv --reviewer-a reviewer_A.csv --reviewer-b reviewer_B.csv --output private_contrast_v2

ソース台帳と完全に一致する128写真ID、異なる評価者ID、双方の開花・花被・背景判定一致、同じ画像SHA、花被・背景各3点のRGB記録から再計算可能な色差のみをprivate計測候補とする。画像SHA一致のみでは元写真IDの真正性は証明されないため、実写真との同一性は人手確認も必要。全てのtechnical unknownと原strict9成功数は変更しない。

v2の追加ブラウザ操作スモークは、この環境ではfileとlocalhostの両方で ERR_BLOCKED_BY_ADMINISTRATOR のため実施不能だった。JavaScript構文およびPython18テスト成功と区別する。v1の合成画像に対するブラウザ操作は過去に別途確認されている。

**真の植物学的確認：0/128、真の花被・背景コントラスト実測：0/128。** 新規写真取得・再推論・仮定的MNAR感度分析は行っていない。測定が実際に完了しても、ΔE76は非校正sRGB画像の局所可視色差にすぎず、UVや送粉者の視覚コントラストではない。
