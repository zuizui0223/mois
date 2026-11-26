#!/usr/bin/env python3
"""Render exact saved same-cell statistics as shareable scientific figures."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd


def source_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paired', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    source, out = Path(args.paired), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    font = Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        family = font_manager.FontProperties(fname=str(font)).get_name()
    else:
        family = 'sans-serif'
    plt.rcParams.update({'font.family': family, 'font.size': 10.5,
                         'axes.unicode_minus': False, 'svg.fonttype': 'path',
                         'axes.spines.top': False, 'axes.spines.right': False})
    stat_path = source / 'native_intro_samecell_statistics.csv'
    comp_path = source / 'native_intro_samecell_colour_composition.csv'
    stats = pd.read_csv(stat_path, float_precision='round_trip')
    comp = pd.read_csv(comp_path, float_precision='round_trip')
    selected = []
    definitions = [
        ('both_richness_ge30_equal_cells', '30種以上・セル等重み'),
        ('both_richness_ge30_min_richness_weight', '30種以上・少数側の種数重み'),
        ('both_richness_ge100_equal_cells', '100種以上・セル等重み'),
    ]
    colors = {'photo9': '#16608b', 'whitefixed8': '#b35d2e'}
    for outcome, model in [('photo9', 'genus_breadth_sixzones'),
                           ('whitefixed8', 'sixzone_genus_breadth_whitefixed8')]:
        outcome_definitions = definitions + ([('both_nonwhite_richness_ge30_equal_cells',
                                               '8色と同じセル・等重み')]
                                             if outcome == 'photo9' else [])
        for key, label in outcome_definitions:
            row = stats[(stats.outcome == outcome) & (stats.model == model) &
                        (stats.population == 'strict') & (stats.component == 'native_minus_introduced') &
                        (stats.statistic_name == key)]
            if len(row) != 1:
                raise ValueError('Expected one exact statistic: ' + outcome + '/' + key)
            item = row.iloc[0].to_dict()
            item['figure_condition'] = label
            selected.append(item)
    data = pd.DataFrame(selected)
    data.to_csv(out / 'paired_effect_figure_data.csv', index=False)
    fig, ax = plt.subplots(figsize=(12.8, 8.3))
    fig.subplots_adjust(left=.31, right=.75, top=.82, bottom=.20)
    values = data.residual_difference.to_numpy(float)
    padding = max(.012, float(np.ptp(np.r_[values, 0])) * .22)
    ax.set_xlim(min(0, values.min()) - padding, max(0, values.max()) + padding)
    ys = np.array([7.0, 6.0, 5.0, 4.0, 2.5, 1.5, .5])
    labels = []
    for y, (_, row) in zip(ys, data.iterrows()):
        label = '9色' if row.outcome == 'photo9' else '白固定8色'
        labels.append(label + '：' + row.figure_condition + '\n' + f"同一 {int(row.n_cells):,} セル")
        ax.scatter(row.residual_difference, y, s=75, color=colors[row.outcome], zorder=4)
        ax.annotate(f'{row.residual_difference:+.4f}', (row.residual_difference, y),
                    xytext=(0, 10), textcoords='offset points', ha='center', fontsize=10,
                    color=colors[row.outcome])
        if row.statistic_name == 'both_richness_ge30_equal_cells':
            ptext = f"主比較  Holm2 p = {row.p_rank_two_sided_holm2:.3f}"
        else:
            ptext = f"感度分析  両側 p = {row.p_rank_two_sided:.3f}"
        ax.text(1.035, y, ptext, transform=ax.get_yaxis_transform(), va='center', fontsize=10)
    ax.axvline(0, ls='--', lw=1.2, color='#69747c')
    ax.axhline(3.25, lw=.8, color='#d7dce0')
    ax.set_yticks(ys, labels)
    ax.tick_params(axis='y', length=0, pad=14)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(-.1, 7.75)
    ax.set_xlabel('在来−外来の Shannon 偏差差', labelpad=12)
    ax.grid(axis='x', alpha=.15)
    fig.suptitle('在来・外来の色集中：同じ場所での直接比較', x=.055, y=.965,
                 ha='left', fontsize=19, fontweight='bold')
    fig.text(.055, .902, '属・全分布幅・粗い地理親和性を保つ共通帰無モデル。種ごとに同じ色を両群へ適用。',
             ha='left', fontsize=11, color='#46535d')
    fig.text(.055, .105, '正の値：外来群の方が帰無期待に対して相対的に集中。負の値：在来群の方が相対的に集中。', fontsize=10)
    fig.text(.055, .067, '各999反復。点は観測差−帰無平均差。主比較の2検定にHolm補正。感度分析のpは未補正。', fontsize=9.5, color='#46535d')
    fig.text(.055, .043, '9色と白固定8色は適格セルが異なる。4行目では9色の比較場所を8色と統一した。', fontsize=9.5, color='#46535d')
    fig.text(.055, .023, 'Mois 2026-11-22  |  WCVP strict・100 kmセル・写真由来色。滞在期間や進化速度を測った比較ではない。', fontsize=9, color='#66717a')
    for extension in ('png', 'svg'):
        fig.savefig(out / ('MOIS_NATIVE_INTRO_PAIRED_EFFECTS_20261122.' + extension), dpi=180)
    plt.close(fig)

    hues = ['white', 'yellow', 'orange', 'red', 'pink', 'magenta', 'purple', 'blue', 'bronze']
    japanese = ['白', '黄', '橙', '赤', '桃', 'マゼンタ', '紫', '青', 'ブロンズ']
    palette = ['#d8d8d8', '#e3bf3d', '#d98532', '#bd403b', '#df8ca8', '#bd4f9a', '#8960a7', '#4b82b0', '#8d6c4b']
    c = comp[(comp.outcome == 'photo9') & (comp.population == 'strict') &
             (comp.statistic_name == 'both_richness_ge30_equal_cells')].set_index('hue').loc[hues].reset_index()
    if len(c) != 9 or not np.isclose(c.native_mean_cell_species_fraction.sum(), 1) or not np.isclose(c.introduced_mean_cell_species_fraction.sum(), 1):
        raise ValueError('Observed paired hue composition does not sum to one')
    c.to_csv(out / 'paired_composition_figure_data.csv', index=False)
    fig, ax = plt.subplots(figsize=(11.7, 7.4))
    fig.subplots_adjust(left=.16, right=.68, top=.81, bottom=.18)
    y = np.arange(9)
    differences = c.native_minus_introduced_percentage_points.to_numpy(float)
    span = max(1., np.max(np.abs(differences))) * 1.28
    ax.set_xlim(-span, span)
    ax.barh(y, differences, color=palette, edgecolor='#535b62', linewidth=.55, height=.60)
    ax.axvline(0, color='#6f7780', lw=1)
    ax.set_yticks(y, japanese)
    ax.invert_yaxis()
    ax.set_xlabel('平均種割合の在来−外来差（パーセントポイント）', labelpad=12)
    ax.grid(axis='x', alpha=.14)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=12)
    ax.text(1.08, 1.045, '在来       外来', transform=ax.transAxes, fontsize=11, fontweight='bold')
    for i, row in c.iterrows():
        value = row.native_minus_introduced_percentage_points
        ax.annotate(f'{value:+.2f}', (value, i), xytext=(6 if value >= 0 else -6, 0),
                    textcoords='offset points', ha='left' if value >= 0 else 'right', va='center', fontsize=10)
        ax.text(1.08, i, f'{100*row.native_mean_cell_species_fraction:5.2f}%   {100*row.introduced_mean_cell_species_fraction:5.2f}%',
                transform=ax.get_yaxis_transform(), va='center', fontsize=11)
    fig.suptitle('同じ場所で、どの写真色の割合が違うか', x=.055, y=.965,
                 ha='left', fontsize=19, fontweight='bold')
    fig.text(.055, .902, f'在来・外来とも30種以上の同一{int(c.n_cells.iloc[0]):,}セル。セル内の種割合を求め、各セルを等重みで平均。',
             fontsize=11, color='#46535d')
    fig.text(.055, .088, '観測構成比の記述。色別の有意差検定や、導入・選択による因果効果を表す図ではない。', fontsize=10)
    fig.text(.055, .035, 'Mois 2026-11-22  |  WCVP strict・100 kmセル。非白色の8色比較とは適格セル集合が異なる。', fontsize=9, color='#66717a')
    for extension in ('png', 'svg'):
        fig.savefig(out / ('MOIS_NATIVE_INTRO_HUE_COMPOSITION_20261122.' + extension), dpi=180)
    plt.close(fig)
    receipt = {'source_SHA256': {stat_path.name: source_hash(stat_path), comp_path.name: source_hash(comp_path)},
               'figure_code_SHA256': source_hash(__file__),
               'output_SHA256': {p.name: source_hash(p) for p in sorted(out.iterdir())
                                 if p.is_file() and p.name != 'figure_receipt.json'},
               'interpretation': 'point estimates and descriptive composition; no effect confidence intervals drawn'}
    (out / 'figure_receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'figures': [p.name for p in sorted(out.glob('*.png'))]}))


if __name__ == '__main__':
    main()
