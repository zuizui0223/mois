# MOIS QA180：元写真・観察・ライセンスの独立監査（2026-11-17）

## 結論
**従来「180種・538枚・538独立観察」とされていた盲検標本は、538の別写真ID・180受容種・現在の独立観察534件である。** iNaturalist元台帳に数値観察IDとUUID観察IDが混在し、同じ現観察IDに収束する4組があった。標本の再抽出は行っていない。

元の538写真とそれぞれの現在のiNaturalist観察の**写真ID対応は538/538確認**。写真ライセンスと著作者表示の公開APIメタデータも538/538で取得。旧写真台帳の分類群IDと現行分類群IDの相違は19写真・10盲検種。分類変更・誤同定をここで原因断定しない。

**画像配信URLの到達確認は538/538（100%）**。最初の一律 `medium.jpg` 試験では315成功/223失敗、うち222は拡張子を `.jpeg` に変えて確認。最後の1枚 `MOIS-QA-032-P01`、photo ID `2963387` は公式APIでの実URLが大文字 `medium.JPG` であり、同じURLのHEADが200・`image/jpeg`。**315×jpg＋222×jpeg＋1×JPG＝538**。画像取得・色計測・花被検証の成功率ではない。

## 原本・処理を変えていないこと
- **真正原本台帳**：1231244行の元写真ID・観察ID・iNaturalist分類群ID。GitHub Actions artifact #14093884178、原ZIP SHA256 `ac558f06e5a0f75cf0fab2d0a37e3f17b2fe3428f788cccb801b5a01dd2b81ac`、内部gzip `f945779b760b46ace93bcfa10763d2ad56f251dbc9c8146d3e44c30b55dca710`。
- **受容種・写真9色属性**：37,104 COL XR exact accepted species。artifact #14157705123、ZIP SHA256 `abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a`。
- **固定済み盲検評価票**：`qa/mois_QA180_20261117/blinded_original_photo_review_sheet.csv`、canonical CSV SHA256 `422d7d0a2a66cf83e3fffcb693007a31a60852a8ee4ae36453fd1a6b171d7442`。9写真色×modal/mean9一致・不一致18層に各10受容種、合計180種。モデル色・証拠statusの回答欄公開なし。写真ID重複ゼロ。2独立の植物学的花被色レビューは**未実施**（0/538）。
- iNaturalist v1数値ID観察126と、v2 `fields=id,uuid,taxon,photos` によるUUID観察412を照合。v2で得た現数値IDからv1詳細写真メタデータに復元し、ライセンス・写真帰属・写真リストを確認。v2 `fields=all` は既知API500になるため使わない。
- 写真ピクセルの新ダウンロードは0。源写真、花被マスク、種花色のモデルラベルを変更せず、QAの回答に自動色をリークしていない。

## 独立性と観測過程

| QA項目 | 実数 |
|---|---:|
| 固定blind species/photo | 180 / 538 |
| 台帳の写真・観察ID文字列の各ユニーク数 | 538 / 538 |
| 現在の**独立数値観察ID** | **534** |
| ID別名が現在同じ観察に収束した組数 | 4 |
| 影響種 | `MOIS-QA-009`、`MOIS-QA-074`、`MOIS-QA-076`、`MOIS-QA-129` |
| 正規化後に独立観察が3件ある種 | **174** |
| 正規化後に独立観察が2件ある種 | **6**（旧々2種＋ID別名4種） |
| 今の観察と元写真IDが対応する行 | **538** |
| 現在写真ライセンス/著作者メタデータを得た行 | **538 / 538** |
| 元と現在のiNat taxon ID不一致 | **19写真/10種** |
| 公開観察の「開花」注記（観察単位のみ） | **221写真の記録**（写真別開花・花被の証拠ではない） |
| provider medium HEAD到達確認（拡張子は元URLに従う） | **538/538** |
| 花被色・マスクの独立専門家認証 | **0/538** |

元観察IDの別名4組は、異なる写真IDを再利用したという意味ではない。**同じ撮影イベントの複数写真が別ID表現で独立と誤認された**。今後の色コンセンサスを計算するとき、写真でなく現観察IDでまとめ、2件観察しかない6種を3観察コンセンサスに自動格上げしない。盲検抽出標本は一切差し替えない。

## 現在写真ライセンス
- `cc-by-nc` 427、`cc-by` 72、`cc-by-sa` 17、`cc0` 14、`cc-by-nc-sa` 8。いずれも現APIによるメタデータで、過去にダウンロードした際の許諾条件の記録とは区別する。
- 現時点で原RGB画像はGitHubへ再配布しない。後で画像自体を表示・頒布するときは写真ごとの帰属、CC-SA/NC条件を個別に確認する。

## 誤認・因果解釈の防止
- **S3 HEAD=200**は元写真の内容・花被・正解マスク・色が正しいことを保証しない。
- **iNaturalist observation-level Flowers label** は538画像のそれぞれが花を示すかどうかの真値ではない。
- **COL XR accepted species** は旧iNaturalist taxon IDが今日も不変である保証ではない。taxon変更19写真はレビュー時に現在帰属と原台帳の分類群を併記し、無断改名しない。
- **9色ROI strict** は花被器官の領域検出真値ではない。旧5例のマスク再現テストでは4例が花被外を選択。過去の写真から見える白率やNPP/LGM相関はこの観測誤差を免れない。
- 過去の全球RF（36,581種）は現代気候→NPP/植生→LGMの属外予測増分がゼロを跨いだ。今回の538写真のsource-QAは**花色真値検証ではない**ので、RFを再適合したことや生物学的効果が強まったことを主張しない。

## 最優先の次工程
1. 元写真ID=2963387の現在の公式写真URL・拡張子・ライセンスを検査済み：`https://inaturalist-open-data.s3.amazonaws.com/photos/2963387/medium.JPG` はHEAD=200、現行ライセンス`cc-by-nc`。画像ピクセルは未取得。
2. URL・現行分類群・ライセンス・現在正規化観察IDを内部レビューメタデータとして残し、2名の独立盲検評価者で花被可視性と可視9色を別々に記録。画像の撮影者・著作権表示を保持する。
3. 「花被可視／正しいROI／色9分類」の条件付き誤分類を区別。4 ID別名の種はクラス間の疑似反復を除外する。層別選抜の比率で集計し、単純538画像率を全球率に外挿しない。
4. QAに通った形質で初めて、同一標本・同一属外ホールドアウトの水分・NPP・AGB・LGM-RFを再評価する。独立レフュージア・化石なしに過去分布適応とは言わない。

## 実行と監査リンク
- [Source registry 538照合・HEAD537の前段（Actions #44089586261）](https://github.com/zuizui0223/mois/actions/runs/44089586261)。
- [JPG→JPEG 222件回復（Actions #44091032539）](https://github.com/zuizui0223/mois/actions/runs/44091032539)、[signed artifact #14331891757](https://github.com/zuizui0223/mois/actions/runs/44091032539/artifacts/14331891757)、SHA256 `8544466b2b8c941706c283327814fece64d23b4f9c8062b13d63d2f4bb4e1bbe`。
- [UUID＋数値IDを現在IDで統合した監査（Actions #44092370121）](https://github.com/zuizui0223/mois/actions/runs/44092370121)、[signed artifact #14331878149](https://github.com/zuizui0223/mois/actions/runs/44092370121/artifacts/14331878149)、SHA256 `46f05942ebd9a07190a984a68d275edfb7820c9f990e3e25a1f8539372108170`。
- [残り1写真の公式API拡張子・HEAD検証（Actions #44093077272）](https://github.com/zuizui0223/mois/actions/runs/44093077272)、[artifact #14332916293](https://github.com/zuizui0223/mois/actions/runs/44093077272/artifacts/14332916293)、SHA256 `238c25e717f2212d9252bf3621363385d8be7f8e79e33edf378b6d2bb0b6b981`。
- GitHub上の数値原本・ライセンス付き画像の二次配布は行っていない。独立QAの生写真ライセンス・撮影者の利用条件は別管理。
