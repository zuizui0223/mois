#!/usr/bin/env python3
"""Strict-native hue nulls whose donors, width and affinity are all native.

This is a bounded extension of the frozen 2026-11-21 global-donor native
analysis. It changes the conditional species-label model, not the observed
native assemblages. It is not a test of selection, residence time or genuine
petal colour. Source geography dependencies are imported only inside run().
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import sys
from importlib import metadata
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

SCHEMA = 'mois_strict_native_footprint_pool_null_20261122_v1'
REFERENCE_SHA = 'c031cf8fe3f80d6cf868ee20477fdb4c550b836435b79ae16d75e771c564f951'
SEED = 20261118
THRESHOLD = .65
TOLERANCE = 1e-12
MODELS = ('native_genus_breadth', 'native_genus_breadth_sixzones')
PRIMARY_MODEL = MODELS[1]
SECTORS = ('Americas_N', 'Americas_S', 'Africa_Eurasia_W_N',
           'Africa_Eurasia_W_S', 'Asia_Oceania_E_N', 'Asia_Oceania_E_S')
SECTOR_GROUPS = tuple(name + '_richness_ge30' for name in SECTORS)
EQUAL_SECTORS = 'six_coordinate_sectors_equal_weight_mean'
NORTH = 'northern_hemisphere_richness_ge30'
SOUTH = 'southern_hemisphere_richness_ge30'
CONTRASTS = {
    'south_minus_north': (SOUTH, NORTH),
    'Arctic_minus_boreal': ('Arctic_lat_ge66p562_richness_ge30',
                           'north_boreal_lat50to66p562_richness_ge30'),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    value = value if isinstance(value, bytes) else Path(value).read_bytes()
    return hashlib.sha256(value).hexdigest()


def build_native_pool(genus_codes, labels, cell, species, zone):
    """Native pair occupancy defines every width and geographic affinity.

    The all-source species catalogue is retained so one label can be applied
    coherently everywhere. Species with no strict-native pair are never donors.
    """
    genus_codes, labels, cell, species, zone = map(
        np.asarray, (genus_codes, labels, cell, species, zone))
    n = len(labels)
    require(genus_codes.shape == labels.shape and labels.ndim == 1 and n > 0,
            'Invalid genus and source label catalogue')
    require(all(np.issubdtype(a.dtype, np.integer) for a in
                (genus_codes, labels, cell, species, zone)), 'Indices must be integers')
    require(cell.shape == species.shape and cell.ndim == 1,
            'Native species and cell vectors differ')
    require(((cell >= 0) & (cell < len(zone))).all() and
            ((species >= 0) & (species < n)).all() and
            ((labels >= 0) & (labels < 9)).all() and
            (genus_codes >= 0).all() and np.isin(zone, np.arange(6)).all(),
            'Invalid native species, colour or coordinate sector index')
    key = cell.astype(np.int64) * n + species
    require(len(np.unique(key)) == len(key), 'Duplicate native species-cell pair')
    coverage = np.bincount(species, minlength=n)
    eligible = coverage > 0
    breadth = np.digitize(coverage, [2, 6, 21, 101]).astype(np.int16)
    breadth[~eligible] = -1
    cross = np.bincount(species.astype(np.int64) * 6 + zone[cell],
                        minlength=n * 6).reshape(n, 6)
    share = cross.max(axis=1) / np.maximum(coverage, 1)
    affinity = np.where(share >= THRESHOLD, cross.argmax(axis=1), 6).astype(np.int8)
    affinity[~eligible] = 7
    primary = genus_codes.astype(np.int64) * 5 + breadth
    regional = primary * 8 + affinity
    primary[~eligible] = regional[~eligible] = -1
    return {'genus_codes': genus_codes, 'coverage': coverage, 'eligible': eligible,
            'breadth_bin': breadth, 'affinity': affinity,
            'top_sector_share': share, 'sector_coverage': cross,
            'codes': {MODELS[0]: primary, MODELS[1]: regional}}


def shuffle_groups(pool, labels, n_hues, geography=True):
    require(n_hues in (8, 9), 'Only photo9 and whitefixed8 are defined')
    labels = np.asarray(labels)
    available = np.flatnonzero(pool['eligible'] & ((labels > 0) if n_hues == 8 else True))
    code = pool['codes'][MODELS[1] if geography else MODELS[0]]
    if not len(available):
        return []
    ordered = available[np.argsort(code[available], kind='stable')]
    splits = np.r_[0, 1 + np.flatnonzero(np.diff(code[ordered])), len(ordered)]
    # Monochromatic strata have a singleton label orbit. Omitting their RNG
    # calls does not change the conditional permutation distribution.
    return [ordered[a:b] for a, b in zip(splits[:-1], splits[1:])
            if b-a >= 2 and len(np.unique(labels[ordered[a:b]])) >= 2]


def draw_labels(original, groups, rng, n_hues):
    original = np.asarray(original)
    ids = np.concatenate(groups) if groups else np.array([], dtype=int)
    require(np.issubdtype(ids.dtype, np.integer) and
            ((ids >= 0) & (ids < len(original))).all() and len(np.unique(ids)) == len(ids),
            'Overlapping or invalid exchange groups')
    require(n_hues in (8, 9), 'Invalid hue dimension')
    if n_hues == 8:
        require((original[ids] > 0).all(), 'White species included in nonwhite exchange')
    labels = original.copy()
    for group in groups:
        labels[group] = original[group][rng.permutation(len(group))]
    return labels


def entropy_counts(cell, species, labels, n_cells, n_hues, original_labels):
    labels = np.asarray(labels)
    require(labels.shape == np.asarray(original_labels).shape and
            np.issubdtype(labels.dtype, np.integer) and
            ((labels >= 0) & (labels < 9)).all(), 'Invalid species hue vector')
    if n_hues == 8:
        require(np.array_equal(labels == 0, np.asarray(original_labels) == 0),
                'White species identities changed')
        require((labels[species] > 0).all(), 'White pair in fixed nonwhite assemblage')
    h = labels[species].astype(np.int64) - (1 if n_hues == 8 else 0)
    counts = np.bincount(np.asarray(cell, dtype=np.int64) * n_hues + h,
                         minlength=n_cells * n_hues).reshape(n_cells, n_hues)
    prob = counts / np.maximum(1, counts.sum(axis=1))[:, None]
    logs = np.zeros_like(prob)
    np.log(prob, out=logs, where=prob > 0)
    return counts, -(prob * logs).sum(axis=1)


def permutation_statistics(observed, draws):
    values = np.asarray(draws, dtype=float)
    require(values.ndim == 1 and len(values) >= 2 and np.isfinite(values).all()
            and np.isfinite(observed), 'Invalid permutation statistics')
    lower = (1 + np.count_nonzero(values <= observed + TOLERANCE)) / (len(values) + 1)
    upper = (1 + np.count_nonzero(values >= observed - TOLERANCE)) / (len(values) + 1)
    sd = float(values.std(ddof=1))
    return {'n_permutations': len(values), 'observed': float(observed),
            'null_mean': float(values.mean()),
            'observed_minus_null': float(observed - values.mean()),
            'null_sd': sd, 'null_mean_monte_carlo_se': sd / np.sqrt(len(values)),
            'null_reference_95pct_interval': np.percentile(values, [2.5, 97.5]).tolist(),
            'interval_interpretation': 'null reference interval, not an effect confidence interval',
            'p_rank_lower': float(lower), 'p_rank_upper': float(upper),
            'p_rank_two_sided': float(min(1, 2 * min(lower, upper))),
            'inference_status': 'conditional_native_species_label_permutation'}


def paired_contrast(observed_a, observed_b, draws_a, draws_b):
    a, b = np.asarray(draws_a), np.asarray(draws_b)
    require(a.ndim == 1 and a.shape == b.shape, 'Paired draws differ in length')
    values = a - b
    stats = permutation_statistics(observed_a - observed_b, values)
    stats.update(observed_contrast=stats['observed'], null_mean_contrast=stats['null_mean'],
                 residual_difference=stats['observed_minus_null'],
                 synchronous_within_model_draws=True)
    return stats, values


def holm_adjust(p):
    p = np.asarray(p, dtype=float)
    require(p.ndim == 1 and len(p) > 0 and np.isfinite(p).all()
            and ((p >= 0) & (p <= 1)).all(), 'Invalid Holm family')
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    result = np.empty(len(p))
    result[order] = adjusted
    return result


def group_means(c, values):
    result = {key: float(values[ids].mean()) for key, ids in c.indices.items()}
    if all(key in result for key in SECTOR_GROUPS):
        result[EQUAL_SECTORS] = float(np.mean([result[key] for key in SECTOR_GROUPS]))
    return result


def simulate_native(c, groups, permutations, seed, model):
    require(permutations >= 2 and all(len(ids) > 0 for ids in c.indices.values()),
            'Insufficient draws or empty summary cells')
    rng = np.random.default_rng(seed)
    sums = np.zeros(c.n_cells)
    draws = {key: np.empty(permutations) for key in group_means(c, c.observed)}
    for i in range(permutations):
        labels = draw_labels(c.original_labels, groups, rng, c.n_hues)
        counts, entropy = entropy_counts(c.cell, c.species, labels, c.n_cells,
                                         c.n_hues, c.original_labels)
        require(np.array_equal(counts.sum(axis=1), c.richness),
                'Fixed native or nonwhite memberships changed')
        sums += entropy
        for key, value in group_means(c, entropy).items():
            draws[key][i] = value
        if (i + 1) % 100 == 0:
            print('MOIS_NATIVE_FOOTPRINT_PROGRESS', c.n_hues, model,
                  i + 1, '/', permutations, flush=True)
    return {'expected': sums / permutations, 'draws': draws}


def support(pool, groups, c):
    mobile = np.concatenate(groups) if groups else np.array([], dtype=int)
    occupied = np.bincount(c.species, minlength=len(c.original_labels))
    eligible = pool['eligible'] & ((c.original_labels > 0) if c.n_hues == 8 else True)
    return {'native_donor_species': int(eligible.sum()), 'variable_donor_strata': len(groups),
            'native_movable_species': len(mobile),
            'native_movable_species_fraction': float(len(mobile) / eligible.sum()),
            'genera_with_variable_donor_strata': int(len(np.unique(pool['genus_codes'][mobile]))),
            'native_species_cell_pairs': len(c.species),
            'native_species_cell_pairs_eligible_for_label_change': int(occupied[mobile].sum()),
            'native_species_cell_pair_fraction_eligible_for_label_change':
                float(occupied[mobile].sum() / len(c.species)),
            'mobility_definition': 'eligible for label change; a draw may retain its original label'}


def summarize(c, result, groups, pool, model, outcome, seed):
    observed = group_means(c, c.observed)
    group_stats, contrasts, frames = {}, {}, []
    for key, value in observed.items():
        stats = permutation_statistics(value, result['draws'][key])
        ids = c.indices['richness_ge30'] if key == EQUAL_SECTORS else c.indices[key]
        require(abs(stats['null_mean'] - group_means(c, result['expected'])[key]) < 1e-10,
                'Per-cell and draw-mean expectations differ')
        stats.update(n_cells=int(len(ids)),
                     mean_native_richness=float(c.richness[ids].mean()),
                     median_native_richness=float(np.median(c.richness[ids])),
                     weighting='equal weight for each of six sector means' if key == EQUAL_SECTORS
                         else 'equal weight for each eligible cell',
                     role='primary family candidate' if model == PRIMARY_MODEL and key == 'richness_ge30'
                         else 'sensitivity or exploratory; no primary multiplicity guarantee')
        group_stats[key] = stats
    for key, (a, b) in CONTRASTS.items():
        stats, values = paired_contrast(observed[a], observed[b],
                                        result['draws'][a], result['draws'][b])
        stats.update(group_A=a, group_B=b, n_cells_A=len(c.indices[a]), n_cells_B=len(c.indices[b]),
                     role='primary family candidate' if model == PRIMARY_MODEL and key == 'south_minus_north'
                         else 'sensitivity or exploratory; no primary multiplicity guarantee')
        require(abs(stats['residual_difference'] - (group_stats[a]['observed_minus_null'] -
                    group_stats[b]['observed_minus_null'])) < 1e-10,
                'Paired regional residual identity failed')
        contrasts[key] = stats
        result['draws'][key] = values
    for key, values in result['draws'].items():
        frames.append(pd.DataFrame({'outcome': outcome, 'model': model, 'cohort': 'strict_native',
                                    'seed': seed, 'statistic_kind': 'contrast' if key in CONTRASTS
                                        else 'group_mean', 'statistic_name': key,
                                    'draw_id': np.arange(1, len(values) + 1), 'draw_value': values}))
    return {'seed': seed, 'exchangeability_support': support(pool, groups, c),
            'group_statistics': group_stats, 'paired_contrasts': contrasts}, frames


def load_reference(path):
    require(digest(path) == REFERENCE_SHA, 'Frozen native 999-reference ZIP SHA differs')
    with ZipFile(path) as z:
        receipt = json.loads(z.read('receipt.json'))
        require(receipt['status'] == 'complete' and receipt['all_reproducibility_gates_passed']
                and receipt['permutations_per_model'] == 999, 'Native reference incomplete')
        for name, sha in receipt['output_SHA256'].items():
            require(digest(z.read(name)) == sha, 'Native reference member differs: ' + name)
        summary = json.loads(z.read('native_georegion_summary.json'))
        tables = {outcome: pd.read_csv(io.BytesIO(z.read(filename)), compression='gzip',
                                      float_precision='round_trip') for outcome, filename in (
            ('photo9', 'native_photo9_georegion_cells.csv.gz'),
            ('whitefixed8', 'native_whitefixed8_georegion_cells.csv.gz'))}
    return summary, tables, {'ZIP_SHA256': REFERENCE_SHA, 'members_SHA256': receipt['output_SHA256']}


def code_provenance():
    root = Path(__file__).resolve().parents[2]
    paths = {Path(__file__).resolve(), root / 'tests/test_native_footprint_pool_null_20261122.py'}
    for name, module in tuple(sys.modules.items()):
        if name.startswith('scripts.analysis.') and getattr(module, '__file__', None):
            path = Path(module.__file__).resolve()
            if path.suffix == '.py' and path.is_relative_to(root):
                paths.add(path)
    return {str(p.relative_to(root)): digest(p) for p in sorted(paths) if p.exists()}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def run(args):
    # Frozen source loaders are only needed for real-data execution; pure
    # estimand and permutation tests do not require geography I/O packages.
    from scripts.analysis import mois_native_hue9_whitefixed8_georegion_null_20261121 as old

    require(args.permutations >= 999, 'At least 999 draws required for the primary analysis')
    out = Path(args.out)
    require(not out.exists(), 'Refuse overwriting a frozen native-footprint output directory')
    previous, reference_cells, reference_integrity = load_reference(args.native999_reference)
    wcvp = Path(args.wcvp)
    require(digest(wcvp / old.MANIFEST) == old.WCVP_RECEIPT_SHA, 'WCVP source receipt changed')
    require(digest(args.tdwg) == old.TDWG_GEO_SHA, 'TDWG3 source polygons changed')
    pairs, species = old.signed_original(args.gbif, args.env)
    x = old.baseline(pairs, species, null_model='genus_breadth')
    coords = pd.DataFrame(np.asarray(x['cell_key'].tolist(), dtype=np.int64),
                          columns=['grid_x', 'grid_y'])
    zone, longitude, latitude = old.geographic_regions(x)
    all_cohorts, source_audit = old.source_guild_rows(pairs, species, coords, wcvp, args.tdwg, coords)
    require(source_audit == previous['native_source_audit'], 'Frozen native status lineage changed')
    source = all_cohorts['strict_native']
    require(source['n_source_species_grid_pairs'] == 2271852, 'Strict-native pair count changed')
    genus_codes, _ = pd.factorize(species.genus.astype(str), sort=True)
    pool = build_native_pool(genus_codes, x['hue'], source['cell_index'], source['species_index'], zone)
    require(int(pool['eligible'].sum()) == 34781, 'Strict-native donor species count changed')
    gates = []
    prepared = {}
    for n, outcome in ((9, 'photo9'), (8, 'whitefixed8')):
        c = old.prepare_assemblage('strict_native', source['cell_index'], source['species_index'],
                                  x['hue'], x['ngrids'], n)
        c.indices = old.summary_indices(c.richness, zone, latitude)
        require(all(len(ids) > 0 for ids in c.indices.values()), 'Empty strict-native summary group')
        ref = reference_cells[outcome]
        gates.extend([
            old.compare_gate(outcome + '_exact_source_coordinate_order', coords.to_numpy(),
                             ref[['grid_x', 'grid_y']].to_numpy(), atol=0),
            old.compare_gate(outcome + '_original_native_richness', c.richness,
                             ref.strict_native_richness.to_numpy(), atol=0),
            old.compare_gate(outcome + '_original_native_observed_entropy', c.observed,
                             ref.strict_native_observed_Shannon.to_numpy()),
        ])
        require(np.array_equal(np.asarray(SECTORS)[zone], ref.coordinate_sector.to_numpy()),
                'Original coordinate sectors changed')
        independent_counts, independent_h = entropy_counts(c.cell, c.species, x['hue'], x['ngrids'], n, x['hue'])
        gates.extend([old.compare_gate(outcome + '_independent_counter', independent_counts, c.counts, atol=0),
                      old.compare_gate(outcome + '_independent_entropy', independent_h, c.observed)])
        for key in ('richness_ge30', 'richness_ge100'):
            expected_n = {(9, 'richness_ge30'): 7231, (9, 'richness_ge100'): 5389,
                          (8, 'richness_ge30'): 5799, (8, 'richness_ge100'): 3051}[n, key]
            require(len(c.indices[key]) == expected_n, 'Frozen native richness threshold changed')
        prepared[outcome] = c
    code_sha = code_provenance()
    print('MOIS_NATIVE_FOOTPRINT_SOURCE_VERIFIED', json.dumps({
        'native_pairs': len(source['cell_index']), 'native_donor_species': int(pool['eligible'].sum()),
        'source_gates_passed': len(gates)}), flush=True)
    tables, outcomes, all_frames, support_members = {}, {}, [], {}
    for n, outcome in ((9, 'photo9'), (8, 'whitefixed8')):
        c = prepared[outcome]
        seed = SEED if n == 9 else SEED + 803
        cells = coords.copy()
        cells['center_longitude'], cells['center_latitude'] = longitude, latitude
        cells['coordinate_sector'] = np.asarray(SECTORS)[zone]
        cells['strict_native_richness'] = c.richness
        cells['strict_native_observed_Shannon'] = c.observed
        for j, hue in enumerate(old.HUES if n == 9 else old.HUES[1:]):
            cells['strict_native_observed_' + hue + '_species'] = c.counts[:, j]
        models = {}
        for model in MODELS:
            groups = shuffle_groups(pool, x['hue'], n, model == PRIMARY_MODEL)
            simulation = simulate_native(c, groups, args.permutations, seed, model)
            models[model], frames = summarize(c, simulation, groups, pool, model, outcome, seed)
            all_frames.extend(frames)
            cells['strict_native_expected_Shannon_' + model] = simulation['expected']
            cells['strict_native_observed_minus_expected_' + model] = c.observed - simulation['expected']
            movable = np.zeros(len(x['hue']), dtype=bool)
            if groups:
                movable[np.concatenate(groups)] = True
            support_members[outcome + '_' + model + '_eligible_for_label_change'] = movable
            # Same source masks and observed values as the frozen global-donor
            # run; this check concerns statistics, not paired cross-model draws.
            legacy_stats = old.permutation_statistics(group_means(c, c.observed)['richness_ge30'],
                                                      simulation['draws']['richness_ge30'])
            current_stats = models[model]['group_statistics']['richness_ge30']
            for key in ('observed', 'null_mean', 'p_rank_lower', 'p_rank_upper', 'p_rank_two_sided'):
                require(abs(legacy_stats[key] - current_stats[key]) <= 1e-12,
                        'Rank statistic differs from frozen verified implementation')
            print('MOIS_NATIVE_FOOTPRINT_MODEL_DONE', outcome, model,
                  json.dumps(current_stats, sort_keys=True), flush=True)
        old_model = 'genus_breadth_sixzones' if n == 9 else 'sixzone_genus_breadth_whitefixed8'
        outcomes[outcome] = {'n_hues': n, 'seed': seed, 'models': models,
                             'previous_global_donor_model': old_model,
                             'previous_global_donor_same_native_evaluation':
                                 previous['outcomes'][outcome]['models'][old_model]['cohorts']['strict_native']}
        tables[outcome] = cells
    family_members = []
    for outcome in ('photo9', 'whitefixed8'):
        for kind, key in (('group_statistics', 'richness_ge30'), ('paired_contrasts', 'south_minus_north')):
            family_members.append((outcome, kind, key,
                outcomes[outcome]['models'][PRIMARY_MODEL][kind][key]))
    adjusted = holm_adjust([stats['p_rank_two_sided'] for _, _, _, stats in family_members])
    family = []
    for i, (outcome, kind, key, stats) in enumerate(family_members):
        stats['p_rank_two_sided_holm4'] = float(adjusted[i])
        stats['role'] = 'primary four-statistic family'
        family.append({'outcome': outcome, 'statistic_kind': kind, 'statistic_name': key,
                       'p_rank_two_sided': stats['p_rank_two_sided'], 'p_holm4': float(adjusted[i])})
    require(code_sha == code_provenance(), 'Analysis code changed during this run')
    source_paths = {'GBIF_zip': args.gbif, 'ENV_zip': args.env, 'native999_reference_zip': args.native999_reference,
                    'WCVP_receipt': wcvp / old.MANIFEST, 'WCVP_exact_names': wcvp / old.NAME_FILE,
                    'WCVP_native_distribution': wcvp / old.DIST_FILE, 'TDWG3_geojson': args.tdwg}
    source_files = {key: {'filename': Path(path).name, 'SHA256': digest(path)}
                    for key, path in source_paths.items()}
    runtime = {'python': platform.python_version(), 'platform': platform.platform()}
    for name in ('numpy', 'pandas', 'pyarrow', 'geopandas', 'shapely', 'pyproj'):
        runtime[name] = metadata.version(name)
    summary = {'schema': SCHEMA,
               'design_status': 'choices recorded before this run; no formal preregistration claim',
               'source_species_unit': 'original accepted GBIF specieskey; WCVP supplies broad regional status only',
               'source_species': len(species), 'source_species_cell_pairs': len(pairs),
               'source_cells': x['ngrids'], 'native_species_cell_pairs': len(source['cell_index']),
               'native_donor_species': int(pool['eligible'].sum()),
               'primary_population': 'strict_native only', 'permutations_per_model': args.permutations,
               'donor_pool': 'species with at least one retained strict-native source cell pair',
               'range_breadth_source': 'number of retained strict-native occupied 100km cells',
               'occupancy_breadth_bins': ['1', '2-5', '6-20', '21-100', '>=101'],
               'affinity_source': 'retained strict-native footprint only', 'affinity_threshold': THRESHOLD,
               'affinity_rule': '>=65% of native cells in one of six sectors; otherwise multizone',
               'sectors': list(SECTORS),
               'sector_definition': 'longitude <-30, -30 to <60, >=60 crossed with signed-northing hemisphere',
               'sector_limit': 'coordinate partitions, not biological realms; TDWG3-specific native colour totals are not held fixed',
               'whitefixed8': 'original white species and nonwhite cell memberships fixed; pool geography comes from all strict-native pairs',
               'shared_seed_does_not_pair_different_models': True,
               'comparison_to_previous': 'same observed native assemblages; donor eligibility, breadth and affinity change together; effect changes cannot be attributed to one component',
               'multiplicity': {'primary_model': PRIMARY_MODEL, 'method': 'Holm4 under arbitrary dependence',
                                'family_size': 4, 'members': family,
                                'other_statistics': 'sensitivity/exploratory; not experiment-wide adjusted'},
               'rank_p_formulas': {'lower': '(1+#draw<=observed+1e-12)/(B+1)',
                                   'upper': '(1+#draw>=observed-1e-12)/(B+1)',
                                   'two_sided': 'min(1,2*min(lower,upper))'},
               'interpretation': 'conditional association of photographic source hue and recorded native distribution; not causal selection, residence time, or validated petal colour',
               'native_status_is_broad_TDWG3_not_individual_wildness': True,
               'source_files': source_files, 'native_source_audit': source_audit,
               'reference_integrity': reference_integrity, 'code_SHA256': code_sha,
               'runtime': runtime, 'outcomes': outcomes}
    # Persist only after source checks and every model have completed.
    out.mkdir(parents=True, exist_ok=False)
    gzip = {'method': 'gzip', 'mtime': 0}
    for outcome, table in tables.items():
        table.to_csv(out / (outcome + '_native_footprint_cells.csv.gz'), index=False, compression=gzip)
    pd.concat(all_frames, ignore_index=True).to_csv(out / 'native_footprint_permutation_draws.csv.gz',
                                                   index=False, compression=gzip)
    rows = []
    for outcome, result in outcomes.items():
        for model, model_result in result['models'].items():
            for kind in ('group_statistics', 'paired_contrasts'):
                for name, stats in model_result[kind].items():
                    row = {'outcome': outcome, 'model': model, 'statistic_kind': kind,
                           'statistic_name': name, **stats}
                    lo, hi = row.pop('null_reference_95pct_interval')
                    row.update(null_reference_2p5_percentile=lo, null_reference_97p5_percentile=hi)
                    rows.append(row)
    pd.DataFrame(rows).to_csv(out / 'native_footprint_statistics.csv', index=False)
    audit = species[['key', 'species', 'genus']].copy().rename(columns={'key': 'specieskey'})
    audit['source_species_index'] = np.arange(len(species))
    audit['original_hue_index'] = x['hue']
    audit['genus_code'] = genus_codes
    audit['native_occupied_cell_count'] = pool['coverage']
    audit['native_breadth_bin_index'] = pool['breadth_bin']
    audit['native_affinity_code'] = pool['affinity']
    audit['native_top_sector_share'] = pool['top_sector_share']
    audit['native_photo9_donor'] = pool['eligible']
    audit['native_whitefixed8_donor'] = pool['eligible'] & (x['hue'] > 0)
    for name, values in support_members.items():
        audit[name] = values
    audit.to_csv(out / 'native_species_pool_audit.csv.gz', index=False, compression=gzip)
    # Compact, receipt-signed derivative for independent pure-numeric replay.
    # Source pairs retain their original catalogue indices and coordinate order.
    np.savez_compressed(out / 'strict_native_prepared_arrays.npz',
                        cell=source['cell_index'].astype(np.int32),
                        species=source['species_index'].astype(np.int32),
                        labels=x['hue'].astype(np.int8), genus_codes=genus_codes.astype(np.int32),
                        grid_x=coords.grid_x.to_numpy(np.int32), grid_y=coords.grid_y.to_numpy(np.int32),
                        zone=zone.astype(np.int8), longitude=longitude, latitude=latitude)
    write_json(out / 'native_footprint_summary.json', summary)
    write_json(out / 'reproducibility_gates.json', {'schema': SCHEMA, 'all_passed': True, 'gates': gates,
                                                  'no_legacy_permutations_rerun': True})
    write_json(out / 'receipt.json', {'schema': SCHEMA, 'status': 'complete', 'source_files': source_files,
                                     'code_SHA256': code_sha, 'runtime': runtime,
                                     'output_SHA256': {p.name: digest(p) for p in sorted(out.iterdir())},
                                     'command_arguments': vars(args), 'permutations_per_model': args.permutations,
                                     'all_reproducibility_gates_passed': True,
                                     'draw_id_interpretation': 'paired within an outcome/model/seed only'})
    print('MOIS_NATIVE_FOOTPRINT_DONE', json.dumps({'output': str(out), 'source_gates_passed': len(gates),
          'saved_draw_rows': sum(len(f) for f in all_frames), 'Holm4': family}, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('gbif', 'env', 'wcvp', 'tdwg', 'native999-reference', 'out'):
        parser.add_argument('--' + key, required=True)
    parser.add_argument('--permutations', type=int, default=999)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
