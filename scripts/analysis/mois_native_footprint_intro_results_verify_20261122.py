#!/usr/bin/env python3
"""Independent numeric review of MOIS native-footprint and same-cell outputs.

No project analysis module is imported. Source/code/receipt hash validation is
performed separately by the root reviewer; this script checks numerical and
scientific identities from the persisted evidence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


HUES = ('white', 'yellow', 'orange', 'red', 'pink', 'magenta', 'purple', 'blue', 'bronze')
SECTORS = ('Americas_N', 'Americas_S', 'Africa_Eurasia_W_N', 'Africa_Eurasia_W_S',
           'Asia_Oceania_E_N', 'Asia_Oceania_E_S')
PRIMARY = 'both_richness_ge30_equal_cells'
NORTH = 'northern_hemisphere_richness_ge30'
SOUTH = 'southern_hemisphere_richness_ge30'
EQUAL6 = 'six_coordinate_sectors_equal_weight_mean'
TOL = 1e-10
COUNTS = {}
MAXIMUM_ERROR = 0.


def close(label, actual, expected, tolerance=TOL):
    global MAXIMUM_ERROR
    a, b = np.asarray(actual), np.asarray(expected)
    assert a.shape == b.shape, (label, 'shape', a.shape, b.shape)
    assert np.isfinite(a).all() and np.isfinite(b).all(), (label, 'nonfinite')
    error = float(np.max(np.abs(a.astype(float) - b.astype(float)), initial=0))
    assert error <= tolerance, (label, error, tolerance)
    MAXIMUM_ERROR = max(MAXIMUM_ERROR, error)
    COUNTS[label] = COUNTS.get(label, 0) + int(a.size)


def frame(path):
    return pd.read_csv(path, float_precision='round_trip', low_memory=False)


def stat_values(observed, draws):
    draws = np.asarray(draws)
    lower_n = int(np.sum(draws <= observed + 1e-12))
    upper_n = int(np.sum(draws >= observed - 1e-12))
    lower, upper = (1 + lower_n) / (len(draws) + 1), (1 + upper_n) / (len(draws) + 1)
    sd = np.std(draws, ddof=1)
    return dict(n_permutations=len(draws), observed=observed, null_mean=np.mean(draws),
                observed_minus_null=observed - np.mean(draws), null_sd=sd,
                null_mean_monte_carlo_se=sd / np.sqrt(len(draws)),
                p_rank_lower=lower, p_rank_upper=upper,
                p_rank_two_sided=min(1., 2 * min(lower, upper)),
                n_draws_at_or_below_observed=lower_n,
                n_draws_at_or_above_observed=upper_n,
                null_reference_2p5_percentile=np.quantile(draws, .025),
                null_reference_97p5_percentile=np.quantile(draws, .975))


def verify_stats(row, observed, draws):
    values = stat_values(observed, draws)
    for key, value in values.items():
        if key in row and not pd.isna(row[key]):
            close('statistics_' + key, row[key], value)
    return values


def holm(p):
    # Direct rank formulation, independent of implementation vectorization.
    order = sorted(range(len(p)), key=lambda i: p[i])
    adjusted = [None] * len(p)
    running = 0.
    for rank, idx in enumerate(order):
        running = max(running, (len(p) - rank) * p[idx])
        adjusted[idx] = min(1., running)
    return adjusted


def native_masks(table):
    n = table.strict_native_richness.to_numpy()
    lat = table.center_latitude.to_numpy()
    zone = table.coordinate_sector.to_numpy()
    masks = {'richness_ge30': n >= 30, 'richness_ge100': n >= 100}
    for name in SECTORS:
        masks[name + '_richness_ge30'] = (n >= 30) & (zone == name)
    masks[NORTH] = (n >= 30) & np.isin(zone, [s for s in SECTORS if s.endswith('_N')])
    masks[SOUTH] = (n >= 30) & np.isin(zone, [s for s in SECTORS if s.endswith('_S')])
    masks['Arctic_lat_ge66p562_richness_ge30'] = (n >= 30) & (lat >= 66.562)
    masks['north_boreal_lat50to66p562_richness_ge30'] = (n >= 30) & (lat >= 50) & (lat < 66.562)
    return masks


def native_aggregate(values, name, masks):
    if name == EQUAL6:
        return np.mean([values[masks[s + '_richness_ge30']].mean() for s in SECTORS])
    if name == 'south_minus_north':
        return values[masks[SOUTH]].mean() - values[masks[NORTH]].mean()
    if name == 'Arctic_minus_boreal':
        return values[masks['Arctic_lat_ge66p562_richness_ge30']].mean() - values[masks['north_boreal_lat50to66p562_richness_ge30']].mean()
    return values[masks[name]].mean()


def entropy_from_pairs(cell, species, labels, n_cells, n_hues):
    keep = np.ones(len(cell), dtype=bool) if n_hues == 9 else labels[species] > 0
    counts = np.zeros((n_cells, n_hues), dtype=np.int64)
    np.add.at(counts, (cell[keep], labels[species[keep]] - (n_hues == 8)), 1)
    total = counts.sum(axis=1)
    values = np.zeros(n_cells)
    # Summing one colour at a time gives a separate implementation of entropy.
    for j in range(n_hues):
        active = counts[:, j] > 0
        p = counts[active, j] / total[active]
        values[active] -= p * np.log(p)
    return counts, values


def verify_native(directory):
    summary = json.loads((directory / 'native_footprint_summary.json').read_text())
    statistics = frame(directory / 'native_footprint_statistics.csv')
    draws = frame(directory / 'native_footprint_permutation_draws.csv.gz')
    tables = {outcome: frame(directory / (outcome + '_native_footprint_cells.csv.gz'))
              for outcome in ('photo9', 'whitefixed8')}
    series = {key: data.sort_values('draw_id') for key, data in
              draws.groupby(['outcome', 'model', 'statistic_name'], sort=False)}
    assert not draws.duplicated(['outcome', 'model', 'statistic_name', 'draw_id']).any()
    for row in statistics.to_dict('records'):
        table = tables[row['outcome']]
        masks = native_masks(table)
        sample = series[row['outcome'], row['model'], row['statistic_name']]
        close('native_draw_identifiers', sample.draw_id.to_numpy(), np.arange(1, len(sample) + 1), 0)
        observed = native_aggregate(table.strict_native_observed_Shannon.to_numpy(), row['statistic_name'], masks)
        expected = native_aggregate(table['strict_native_expected_Shannon_' + row['model']].to_numpy(), row['statistic_name'], masks)
        close('native_cell_expected_aggregation', row['null_mean'], expected)
        verify_stats(row, observed, sample.draw_value.to_numpy())
        if row['statistic_name'] in ('south_minus_north', 'Arctic_minus_boreal'):
            a, b = (SOUTH, NORTH) if row['statistic_name'] == 'south_minus_north' else (
                'Arctic_lat_ge66p562_richness_ge30', 'north_boreal_lat50to66p562_richness_ge30')
            av = series[row['outcome'], row['model'], a].draw_value.to_numpy()
            bv = series[row['outcome'], row['model'], b].draw_value.to_numpy()
            close('native_paired_regional_draws', sample.draw_value.to_numpy(), av - bv)
    family = summary['multiplicity']['members']
    raw_primary = []
    for member in family:
        outcome, statistic = member['outcome'], member['statistic_name']
        table = tables[outcome]
        observed = native_aggregate(table.strict_native_observed_Shannon.to_numpy(), statistic, native_masks(table))
        values = series[outcome, summary['multiplicity']['primary_model'], statistic].draw_value.to_numpy()
        raw_primary.append(stat_values(observed, values)['p_rank_two_sided'])
    close('native_primary_family_raw_p', [x['p_rank_two_sided'] for x in family], raw_primary)
    close('native_Holm4', [x['p_holm4'] for x in family], holm(raw_primary))

    with np.load(directory / 'strict_native_prepared_arrays.npz') as packed:
        arrays = {key: packed[key] for key in packed.files}
    audit = frame(directory / 'native_species_pool_audit.csv.gz')
    cell, sp, original = arrays['cell'], arrays['species'], arrays['labels']
    n, n_cells = len(original), len(arrays['zone'])
    coverage = np.bincount(sp, minlength=n)
    sectors = np.zeros((n, 6), dtype=np.int64)
    np.add.at(sectors, (sp, arrays['zone'][cell]), 1)
    share = np.divide(sectors.max(axis=1), coverage, out=np.zeros(n), where=coverage > 0)
    affinity = np.full(n, 6, dtype=np.int8)
    affinity[share >= .65] = sectors.argmax(axis=1)[share >= .65]
    affinity[coverage == 0] = 7
    breadth = np.select([coverage == 0, coverage == 1, coverage <= 5, coverage <= 20, coverage <= 100],
                        [-1, 0, 1, 2, 3], default=4)
    close('native_source_species_order', audit.source_species_index.to_numpy(), np.arange(n), 0)
    close('native_independent_coverage', audit.native_occupied_cell_count.to_numpy(), coverage, 0)
    close('native_independent_breadth', audit.native_breadth_bin_index.to_numpy(), breadth, 0)
    close('native_independent_affinity', audit.native_affinity_code.to_numpy(), affinity, 0)
    close('native_independent_affinity_share', audit.native_top_sector_share.to_numpy(), share)
    close('native_independent_donor_ids', audit.native_photo9_donor.to_numpy(), coverage > 0, 0)
    close('native_independent_nonwhite_donor_ids', audit.native_whitefixed8_donor.to_numpy(), (coverage > 0) & (original > 0), 0)
    for outcome, n_hues in (('photo9', 9), ('whitefixed8', 8)):
        table = tables[outcome]
        masks = native_masks(table)
        counts, observed = entropy_from_pairs(cell, sp, original, n_cells, n_hues)
        close('native_independent_counter', counts.sum(axis=1), table.strict_native_richness.to_numpy(), 0)
        close('native_independent_entropy', observed, table.strict_native_observed_Shannon.to_numpy())
        for model in ('native_genus_breadth', 'native_genus_breadth_sixzones'):
            strata = {}
            for idx in np.flatnonzero((coverage > 0) & ((original > 0) if n_hues == 8 else True)):
                key = (int(arrays['genus_codes'][idx]), int(breadth[idx]))
                if model.endswith('sixzones'):
                    key += (int(affinity[idx]),)
                strata.setdefault(key, []).append(idx)
            groups = [np.asarray(strata[k]) for k in sorted(strata)
                      if len(set(original[strata[k]])) > 1]
            mobile = np.zeros(n, dtype=bool)
            if groups:
                mobile[np.concatenate(groups)] = True
            close('native_independent_movable_ids', mobile,
                  audit[outcome + '_' + model + '_eligible_for_label_change'].to_numpy(), 0)
            support = summary['outcomes'][outcome]['models'][model]['exchangeability_support']
            close('native_support_movable_species', support['native_movable_species'], mobile.sum(), 0)
            close('native_support_movable_pairs', support['native_species_cell_pairs_eligible_for_label_change'], coverage[mobile].sum(), 0)
            # Independent first three full-species shuffles establish coherent
            # labels and correctly paired regional statistics from saved arrays.
            rng = np.random.default_rng(20261118 if n_hues == 9 else 20261118 + 803)
            for draw_index in range(3):
                labels = original.copy()
                for group in groups:
                    labels[group] = original[group][rng.permutation(len(group))]
                _, h = entropy_from_pairs(cell, sp, labels, n_cells, n_hues)
                for statistic in ('richness_ge30', 'richness_ge100', 'south_minus_north', EQUAL6):
                    actual = series[outcome, model, statistic].draw_value.iloc[draw_index]
                    close('native_first_three_independent_shuffle_statistics', actual,
                          native_aggregate(h, statistic, masks))
    return {'draw_rows': len(draws), 'statistic_series': len(series), 'statistics_rows': len(statistics),
            'source_species': n, 'native_source_pairs': len(cell), 'independent_full_species_draws': 12,
            'primary_family': family}


def verify_paired(directory):
    summary = json.loads((directory / 'native_intro_samecell_summary.json').read_text())
    statistics = frame(directory / 'native_intro_samecell_statistics.csv')
    draws = frame(directory / 'native_intro_samecell_permutation_draws.csv.gz')
    masks = frame(directory / 'native_intro_samecell_masks_weights.csv.gz')
    composition = frame(directory / 'native_intro_samecell_colour_composition.csv')
    tables = {outcome: frame(directory / ('native_intro_' + outcome + '_samecell_cells.csv.gz'))
              for outcome in ('photo9', 'whitefixed8')}
    for outcome, table in tables.items():
        for population in ('strict', 'moderate'):
            for status in ('native', 'introduced'):
                prefix = population + '_' + status
                hue_names = HUES if outcome == 'photo9' else HUES[1:]
                counts = table[[prefix + '_observed_' + hue + '_species' for hue in hue_names]].to_numpy()
                richness = counts.sum(axis=1)
                observed = np.zeros(len(table))
                for j in range(counts.shape[1]):
                    active = counts[:, j] > 0
                    p = counts[active, j] / richness[active]
                    observed[active] -= p * np.log(p)
                close('paired_entropy_from_saved_colour_counts', observed,
                      table[prefix + '_observed_Shannon'].to_numpy())
                close('paired_richness_from_saved_colour_counts', richness,
                      table[prefix + '_richness'].to_numpy(), 0)
                if outcome == 'whitefixed8':
                    photo9 = tables['photo9']
                    close('paired_original_white_membership_fixed', richness,
                          photo9[prefix + '_richness'].to_numpy() - photo9[prefix + '_observed_white_species'].to_numpy(), 0)
    draw_series = {k: v.sort_values('draw_id') for k, v in draws.groupby(
        ['outcome', 'model', 'population', 'statistic_name'], sort=False)}
    mask_series = {k: v.sort_values('cell_index') for k, v in masks.groupby(
        ['outcome', 'population', 'statistic_name'], sort=False)}
    assert not draws.duplicated(['outcome', 'model', 'population', 'statistic_name', 'draw_id']).any()
    assert not masks.duplicated(['outcome', 'population', 'statistic_name', 'cell_index']).any()
    for key, mask in mask_series.items():
        outcome, population, statistic = key
        table = tables[outcome]
        ids = mask.cell_index.to_numpy()
        w = mask.normalized_common_weight.to_numpy()
        nr = table[population + '_native_richness'].to_numpy()
        ir = table[population + '_introduced_richness'].to_numpy()
        close('paired_weights_sum', w.sum(), 1.)
        close('paired_coordinate_match', mask[['grid_x', 'grid_y']].to_numpy(), table.iloc[ids][['grid_x', 'grid_y']].to_numpy(), 0)
        close('paired_native_richness', mask.native_richness.to_numpy(), nr[ids], 0)
        close('paired_intro_richness', mask.introduced_richness.to_numpy(), ir[ids], 0)
        if statistic not in ('both_richness_ge30_min_richness_weight', 'six_coordinate_sectors_equal_weight'):
            close('paired_all_equal_cell_weights', w, np.full(len(ids), 1 / len(ids)))
        common = (nr >= 30) & (ir >= 30)
        if statistic == PRIMARY:
            close('paired_primary_common_mask', ids, np.flatnonzero(common), 0)
            close('paired_primary_equal_weights', w, np.full(len(ids), 1 / len(ids)))
        elif statistic == 'both_richness_ge30_min_richness_weight':
            close('paired_minrichness_common_mask', ids, np.flatnonzero(common), 0)
            raw = np.minimum(nr[ids], ir[ids])
            close('paired_minrichness_weights', w, raw / raw.sum())
        elif statistic == 'both_richness_ge100_equal_cells':
            close('paired_highrichness_common_mask', ids, np.flatnonzero((nr >= 100) & (ir >= 100)), 0)
        elif statistic == 'both_nonwhite_richness_ge30_equal_cells':
            other = tables['whitefixed8']
            keep = (other[population + '_native_richness'] >= 30) & (other[population + '_introduced_richness'] >= 30)
            close('paired_common8_mask', ids, np.flatnonzero(keep), 0)
        elif statistic == 'six_coordinate_sectors_equal_weight':
            close('paired_equal6_common_mask', ids, np.flatnonzero(common), 0)
            for sector in SECTORS:
                close('paired_equal6_sector_weight', w[table.coordinate_sector.iloc[ids].to_numpy() == sector].sum(), 1 / 6)
        else:
            if statistic.startswith('northern_hemisphere'):
                selected = table.coordinate_sector.str.endswith('_N').to_numpy()
            elif statistic.startswith('southern_hemisphere'):
                selected = table.coordinate_sector.str.endswith('_S').to_numpy()
            else:
                selected = table.coordinate_sector.eq(statistic.replace('_both_richness_ge30', '')).to_numpy()
            close('paired_region_common_mask', ids, np.flatnonzero(common & selected), 0)
    for row in statistics.to_dict('records'):
        outcome, model, population, statistic = (row[k] for k in ('outcome', 'model', 'population', 'statistic_name'))
        sample = draw_series[outcome, model, population, statistic]
        mask = mask_series[outcome, population, statistic]
        table = tables[outcome]
        ids, w = mask.cell_index.to_numpy(), mask.normalized_common_weight.to_numpy()
        native_obs = np.dot(table[population + '_native_observed_Shannon'].to_numpy()[ids], w)
        intro_obs = np.dot(table[population + '_introduced_observed_Shannon'].to_numpy()[ids], w)
        native_exp = np.dot(table[population + '_native_expected_Shannon_' + model].to_numpy()[ids], w)
        intro_exp = np.dot(table[population + '_introduced_expected_Shannon_' + model].to_numpy()[ids], w)
        component = row['component']
        obs, exp = (native_obs, native_exp) if component == 'native' else (
            (intro_obs, intro_exp) if component == 'introduced' else (native_obs - intro_obs, native_exp - intro_exp))
        values = sample[component + '_Shannon'].to_numpy()
        verify_stats(row, obs, values)
        close('paired_cell_expected_aggregation', row['null_mean'], exp)
        close('paired_draw_identifiers', sample.draw_id.to_numpy(), np.arange(1, len(sample) + 1), 0)
        close('paired_synchronous_difference', sample.native_minus_introduced_Shannon.to_numpy(),
              sample.native_Shannon.to_numpy() - sample.introduced_Shannon.to_numpy())
    family = summary['primary_multiplicity']['family']
    raw_primary = []
    for member in family:
        outcome, model, population, statistic = (member[k] for k in ('outcome', 'model', 'population', 'statistic_name'))
        table, mask = tables[outcome], mask_series[outcome, population, statistic]
        ids, w = mask.cell_index.to_numpy(), mask.normalized_common_weight.to_numpy()
        contrast = table[population + '_native_observed_Shannon'].to_numpy() - table[population + '_introduced_observed_Shannon'].to_numpy()
        observed = np.dot(contrast[ids], w)
        values = draw_series[outcome, model, population, statistic].native_minus_introduced_Shannon.to_numpy()
        raw_primary.append(stat_values(observed, values)['p_rank_two_sided'])
    close('paired_primary_family_raw_p', [x['p_rank_two_sided'] for x in family], raw_primary)
    close('paired_Holm2', [x['p_rank_two_sided_holm2'] for x in family], holm(raw_primary))
    for key, data in composition.groupby(['outcome', 'population', 'statistic_name'], sort=False):
        outcome, population, statistic = key
        mask, table = mask_series[key], tables[outcome]
        ids, w = mask.cell_index.to_numpy(), mask.normalized_common_weight.to_numpy()
        close('composition_native_fraction_sum', data.native_mean_cell_species_fraction.sum(), 1.)
        close('composition_intro_fraction_sum', data.introduced_mean_cell_species_fraction.sum(), 1.)
        close('composition_delta_sum', data.native_minus_introduced_percentage_points.sum(), 0.)
        for row in data.to_dict('records'):
            means = {}
            for status in ('native', 'introduced'):
                prefix = population + '_' + status
                fraction = table[prefix + '_observed_' + row['hue'] + '_species'].to_numpy()[ids] / table[prefix + '_richness'].to_numpy()[ids]
                means[status] = np.dot(fraction, w)
                close('composition_independent_cell_weighted_mean', row[status + '_mean_cell_species_fraction'], means[status])
            close('composition_independent_percentage_point_difference', row['native_minus_introduced_percentage_points'], 100 * (means['native'] - means['introduced']))
    return {'draw_rows': len(draws), 'paired_statistic_series': len(draw_series),
            'component_statistic_rows': len(statistics), 'mask_weight_rows': len(masks),
            'composition_rows': len(composition), 'primary_family': family}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--paired', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    native = verify_native(args.native)
    paired = verify_paired(args.paired)
    result = {'status': 'passed', 'project_analysis_modules_imported': False,
              'maximum_absolute_numeric_error': MAXIMUM_ERROR,
              'numeric_comparison_tolerance': TOL,
              'comparison_value_counts': COUNTS, 'native_footprint': native,
              'native_introduced_samecell': paired,
              'receipt_and_code_hash_validation': 'handled separately by root reviewer'}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
