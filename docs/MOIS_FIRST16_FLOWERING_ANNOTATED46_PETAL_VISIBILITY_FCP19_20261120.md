# MOIS 2026-11-20：開花注記付き真正原画像46枚と花被の可視性ゲート

## 判定・現在地

**16種の優先WCVP在来候補に対し、先行取得の48原画像（5枚のみiNaturalist開花注記）を写真ID／観察UUIDで除外し、開花注記`term_id=12&term_value_id=13`を元観察JSONでも確認した真正CC原画像を46枚追加取得した。** 15種で独立観察3枚、`Carex glacialis`のみ1枚。原画像バイトSHA256は46枚すべて確認。API失敗0、画像取得失敗0。これは**開花が観察に注記されている**証拠であり、当該写真に種別の真正花被が明瞭に写ることは保証しない。

- 真正取得：GitHub Actions [#44751415736](https://github.com/zuizui0223/mois/actions/runs/44751415736) **SUCCESS**、artifact **#14635796309**、ZIP SHA256 `c2ae31d6d492ba08b07457ff68dfc39c76870723974383219700d7ec7ed5395e`。
- 既存48枚原本：Actions [#44547129368](https://github.com/zuizui0223/mois/actions/runs/44547129368)、artifact #14541906843、ZIP SHA256 `a6c96f5e2db118bef282688a237d720362d651a16a611ca01907ad12c5cb0998`。過去の1,231,337写真ID・1,083,073観察UUIDも照合済み（ただしそれ以降の全MOIS取得波との非重複は未証明）。
- 真正写真46枚の観察UUIDは相互にすべて異なる。各写真の実画像SHA・CCライセンス・元画像URL・注記JSONをZIP内CSVに保存。フル原画像バイト列は一時的に使用、長期成果物にはサムネイルとSHA・URLを保持。
- iNaturalist公式注記条件：[How to use iNaturalist Search URLs](https://help.inaturalist.org/en/support/solutions/articles/151000170738-how-to-use-inaturalist-search-urls)。「Flowers」は`12|13`。元観察JSONを再検査した上で受け入れ、APIパラメータだけを証拠にしない。

## 440px縮小画像の予備目視による花被可視性：46枚

**この判定は写真からの予備トリアージで、植物学専門家による確認・完全解像度の器官マスク・CIELAB色量化はまだ行っていない。** 

| 候補種 | 花被色を読めそうな候補 | 原寸確認必要 | 今回は色計測不適 |
|---|---:|---:|---:|
| *Aristida adscensionis* | 0 | 0 | 3 |
| *Lonicera caerulea* | 2 | 1 | 0 |
| *Eleocharis uniglumis* | 0 | 0 | 3 |
| *Spondias mombin* | 0 | 1 | 2 |
| *Eugenia punicifolia* | 0 | 2 | 1 |
| *Eremophila latrobei* | 3 | 0 | 0 |
| *Pedicularis lapponica* | 2 | 1 | 0 |
| *Chamaecrista absus* | 2 | 0 | 1 |
| *Hedysarum hedysaroides* | 3 | 0 | 0 |
| *Sporobolus pyramidalis* | 0 | 0 | 3 |
| *Amphilophium paniculatum* | 1 | 1 | 1 |
| *Glycine canescens* | 2 | 1 | 0 |
| *Rumex fueginus* | 0 | 0 | 3 |
| *Carex glacialis* | 0 | 0 | 1 |
| *Lythrum virgatum* | 2 | 1 | 0 |
| *Banisteriopsis muricata* | 2 | 0 | 1 |
| **計** | **19枚（9種）** | **8枚** | **19枚** |

### 注記が特に有効だった例と限界

- *Hedysarum hedysaroides* は先行3写真では果実・莢が中心だったが、今回は紫系の開花器官を視認できる候補が3写真あった。*Lythrum*・*Chamaecrista*・*Glycine*にも花被候補がある。
- *Eugenia* は1枚に白い雄しべが見えるが、花弁の色領域を雄しべから分離できず`uncertain`。雄しべや果実の色を花被色と混同しない。
- *Aristida*・*Eleocharis*・*Sporobolus*・*Carex* は開花注記が存在しても、草本の穂・小穂などが中心で**9色の花被器官として計測できる画素が確認できない**。これは種が花をつけないという判断ではない。
- 先行48枚で「花冠が見える写真あり」は4種だったが、今回の独立46枚で「花被色の読取候補」は9種。非ランダムな注記付き探索・別画像・非盲検の判定なので、改善率を普遍的な検出確率の因果効果としない。

## 証拠ファイルと計測入口

[全46写真の予備目視トリアージCSV](../analysis_specs/MOIS_FIRST16_46_FLOWERING_ANNOTATED_THUMBNAIL_PETAL_TRIAGE_20261120.csv) は、真正photo ID・観察UUID・原画像SHA・CCライセンス・原画像URLと、各写真の`clear_candidate` / `uncertain` / `not_measurable` を1行1写真で保持。CSV SHA256 `5f580ace30f8b359e3b9b7efdc22bdd90f57fb810a09d4d517aca7459399e572`。

**19 `clear_candidate` だけ**を元のFCP YOLO11/EfficientSAM ROI-v4に入力する診断コード `scripts/acquisition/mois_20261120_annotated_19_FCP_ROI9_provisional.py` を追加。真正ZIP・トリアージCSV・photo IDをSHAで固定し、元の`mois_full_index_real_original_photo_measure`のFCP一括ワーカーで9成分数値と技術上の不明を分ける。Workflow `.github/workflows/mois-20261120-FCP19-flowering-ROI9-provisional.yml` （[Actions実行一覧](https://github.com/zuizui0223/mois/actions/workflows/mois-20261120-FCP19-flowering-ROI9-provisional.yml)）。

**重要：技術ROI9が出ても元の植物種写真色ラベルを自動変更しない。** 19写真は人手で見える候補にすぎず、器官マスクは未認証。数値ROI結果は別の探索的レイヤーに隔離し、元37,104種・全球GBIF分布・LGM関連の効果推定量に自動追加せず、バイアスへの示唆として使う。

## 次の生態学的判断

新しい注記付き画像からの取得は効率を改善しうるが、**花が「開花中」と注記される傾向が種によって異なる**。草本小花、花冠が目立たない種、果実/葉が多い写真、開花注記の欠損は同じ標本枠でない。この検出過程を制御しない限り「白花種が多い／少ない」「LGMレフュジアで特定色が残存した」等の生態学的結論には結び付けない。

PR #200のデータ取得側への分離実装。主系統へのマージ、写真元ラベルの上書き、種内固定色判定・LGMレフュジア因果は行っていない。
