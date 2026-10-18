"""Acquire unprocessed cached photos and store fixed-model measurements only.

The provider loader, rate limits and ROI-v4 model are unchanged. Finished
256-photo units are reused. Candidate metadata is never counted as an image
or hue; acquisition failures and acquired-but-unclassifiable photos differ.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import pandas as pd

HUES = ('white', 'yellow_orange', 'red_pink', 'blue_purple')
POOL = 'mois_global_angiosperms_broad_open_source_photo_ID_pool.csv'
FIRST = 'mois_broad_first4096_original_open_flower_photo_ID_source.csv'
FIRST_RECEIPT = 'mois_broad_angiosperm_new_taxa_32768_pool_first4096_receipt.json'
SECOND = 'mois_global_angiosperm_second4096_source_only_original_photoids.csv'
SECOND_RECEIPT = 'mois_global_angiosperm_cached_second4096_source_only_receipt.json'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def select_remaining(pool, first, second, count):
    require(count > 0 and count % 256 == 0, 'Use complete 256-photo units')
    require(all(t.photo_id.is_unique and t.measurement_id.is_unique for t in (pool, first, second)),
            'Duplicate original photo IDs')
    excluded = set(first.photo_id) | set(second.photo_id)
    require(set(first.photo_id).isdisjoint(second.photo_id) and excluded.issubset(set(pool.photo_id)),
            'Previous source rosters are overlapping or not in the original pool')
    observations = set(first.observation_id) | set(second.observation_id)
    fresh = pool.loc[~pool.photo_id.isin(excluded) & ~pool.observation_id.isin(observations)]
    require(len(fresh) >= count, 'Not enough unprocessed cached photo records')
    selected = fresh.iloc[:count].copy().reset_index(drop=True)
    require(selected.observation_id.is_unique, 'Repeated original observations in new photo source')
    selected['source_arm'] = 'mois_cached_remaining_after_broad_waves1_2'
    selected['shard'] = (selected.index // 256).astype(str)
    return selected


def prepare(pool_root, second_root, completed_root, output, count=8192):
    pool_root, second_root, completed_root, output = map(Path, (pool_root, second_root, completed_root, output))
    a = json.loads((pool_root / FIRST_RECEIPT).read_text())
    b = json.loads((second_root / SECOND_RECEIPT).read_text())
    done = json.loads((completed_root / 'receipt.json').read_text())
    require(a['schema'] == 'mois_open_global_angiosperm_breadth_first4096_v1'
            and a['original_broad_all_future_source_pool_SHA256'] == digest(pool_root / POOL)
            and a['first_wave_original_4096_csv_sha256'] == digest(pool_root / FIRST),
            'Original cached photo pool fingerprint changed')
    require(b['schema'] == 'mois_angiosperm_cached_second4096_new_species_source_only_v1'
            and b['parent_cached_pool_original_32768_photo_SHA256'] == digest(pool_root / POOL)
            and b['parent_first_wave_original_4096_csv_sha256'] == digest(pool_root / FIRST)
            and b['source_only_second4096_csv_sha256'] == digest(second_root / SECOND)
            and b['original313535_measured_photo_IDs_excluded_via_parent_signed_source'] is True
            and b['other_fifth4096_and_SI4096_reservations_excluded_via_parent'] is True,
            'Previous source exclusions or snapshot identity changed')
    require(done['schema'] == 'mois_cached_wave2_complete_data_only_v1'
            and done['attempted'] == 4096
            and done['source_sha256'] == digest(second_root / SECOND)
            and done['photo_data_sha256'] == digest(completed_root / 'wave2_measured_photos.parquet'),
            'Previous real image wave is not complete and preserved')
    pool, first, second = read_csv(pool_root / POOL), read_csv(pool_root / FIRST), read_csv(second_root / SECOND)
    require((len(pool), len(first), len(second)) == (32768, 4096, 4096), 'Original source counts changed')
    require(first.equals(pool.iloc[:4096].reset_index(drop=True)), 'Original first source rows changed')
    # Verify second-wave metadata came from the same pool before exclusion.
    original = pool.set_index('photo_id')
    immutable = ['measurement_id', 'inat_taxon_id', 'observation_id', 'photo_url_large', 'photo_license']
    require(second.set_index('photo_id')[immutable].equals(original.loc[second.photo_id, immutable]),
            'Previously selected photo metadata does not match the original pool')
    selected = select_remaining(pool, first, second, count)
    require(not output.exists() or not any(output.iterdir()), 'Do not overwrite an existing source plan')
    output.mkdir(parents=True, exist_ok=True)
    selected.to_csv(output / 'candidates.csv', index=False, lineterminator='\n')
    receipt = {'schema': 'mois_cached_remaining_source_v1', 'status': 'SOURCE_ONLY',
               'count': count, 'shards': count // 256, 'data_sha256': digest(output / 'candidates.csv'),
               'no_ecological_analysis': True, 'source_pool_count': len(pool),
               'excluded_previous_photo_ids': len(set(first.photo_id) | set(second.photo_id)),
               'parents': {'pool': digest(pool_root / POOL), 'first': digest(pool_root / FIRST),
                           'second': digest(second_root / SECOND),
                           'second_measurement': done['photo_data_sha256']}}
    (output / 'source.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    for i in range(count // 256):
        source_paths(output, i)
    print('MOIS_NEXT_CACHED_PHOTO_QUEUE_READY', json.dumps(receipt), flush=True)
    return receipt


def source_paths(root, index):
    root = Path(root)
    r = json.loads((root / 'source.json').read_text())
    require(r['schema'] == 'mois_cached_remaining_source_v1' and r['status'] == 'SOURCE_ONLY'
            and r['count'] == r['shards'] * 256 and index in range(r['shards'])
            and r['data_sha256'] == digest(root / 'candidates.csv'), 'Source plan checksum/count mismatch')
    full = read_csv(root / 'candidates.csv')
    require(len(full) == r['count'] and full.photo_id.is_unique and full.measurement_id.is_unique
            and full.observation_id.is_unique and not {'morph', 'hue', 'classified_four_hue'}.intersection(full.columns),
            'Duplicate photo IDs or source phenotype leakage')
    rows = full[full.shard.eq(str(index))].reset_index(drop=True)
    require(len(rows) == 256, 'Incomplete photo unit')
    for x in rows.itertuples():
        require(x.photo_id.isdecimal() and re.fullmatch(
            r'https://inaturalist-open-data\.s3\.amazonaws\.com/photos/' + x.photo_id + r'/large\.(jpg|jpeg|png)',
            x.photo_url_large) and re.fullmatch(r'MOIS-[0-9A-Fa-f]{24}', x.measurement_id)
            and x.image_filename == x.measurement_id + '.jpg', 'Unsafe/mismatched original image URL or path')
    return rows, r


class PhotoBatch:
    def __init__(self, rows, index):
        self.rows, self.index = rows, index

    def get_batch(self, part, batch, size):
        require((part, batch, size) == (self.index, 0, 256), 'Not a single 256-photo unit')
        cols = ['measurement_id', 'image_filename', 'photo_license']
        return self.rows[cols].copy(), self.rows[cols + ['photo_url_large']].copy(), 1


def measure(root, index, runtime, sam, out, runner=None):
    rows, source = source_paths(root, index)
    target = Path(out) / f'shard_{index:02d}'
    report, photos = target / 'receipt.json', target / 'photos.csv'
    production = runner is None
    if report.exists():
        r = json.loads(report.read_text())
        require(r['source_sha256'] == source['data_sha256'] and r['shard'] == index
                and r['attempted'] == 256 and r['photos_sha256'] == digest(photos),
                'Corrupted completed cache; do not overwrite or repeat it')
        require(read_csv(photos).measurement_id.tolist() == rows.measurement_id.tolist(), 'Cached photo IDs changed')
        if production:
            from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
            require(r['method'] == METHOD, 'Cached model settings changed')
        return r
    if production:
        from scripts.acquisition.run_fcp_fresh_image_partition import run_batch
        runner = run_batch
    target.mkdir(parents=True, exist_ok=True)
    sealed = target / 'sealed'; runtime = Path(runtime)
    result = runner(Path('.'), index, 0, 256, sealed,
        source=runtime / 'scripts/analysis/measure_random_photo_first_blind_partition.py',
        detector=runtime / 'data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt',
        sam=Path(sam), indexed_source=PhotoBatch(rows, index))
    require(result['status'] == 'COMPLETE_WORK_UNIT' and result['n_source_photo_ids'] == 256
            and result['n_terminal_photo_ids'] == 256 and result['image_bytes_persisted'] is False
            and result['reused_historical_measurement_csv'] is False
            and result['number_batches_for_this_partition'] == 1, 'Image unit did not complete')
    if production:
        from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
        require(result['method'] == METHOD, 'Original ROI model/source changed')
    terminal = sealed / f'terminal_p{index:03d}_b0000000.csv'
    require(digest(terminal) == result['terminal_csv_sha256'], 'Image measurement fingerprint changed')
    m = read_csv(terminal)
    require(len(m) == 256 and m.measurement_id.is_unique and set(m.measurement_id) == set(rows.measurement_id),
            'Image results missing/duplicating photo identities')
    cols = ['measurement_id', 'morph', 'measurement_status', 'image_sha256']
    cols += [k for k in m.columns if k.startswith('flower_fraction_')]
    data = rows.merge(m[cols], on='measurement_id', validate='one_to_one', sort=False)
    known = data.morph.isin(HUES); failed = data.measurement_status.eq('image_acquisition_failed')
    require(data.morph.isin((*HUES, 'mixed_uncertain')).all()
            and known.eq(data.measurement_status.eq('classified_four_state_morph')).all()
            and int(known.sum()) == result['n_classifiable']
            and int(failed.sum()) == result['n_acquisition_failed'], 'Image/hue status count disagreement')
    data['classified_four_hue'] = known.map({True: 'true', False: 'false'})
    data.to_csv(photos, index=False, lineterminator='\n')
    r = {'schema': 'mois_cached_remaining_image_unit_v1', 'shard': index,
         'attempted': 256, 'downloaded': int((~failed).sum()), 'acquisition_failed': int(failed.sum()),
         'classified': int(known.sum()), 'unclassified': int((~known).sum()),
         'acquired_unclassified': int((~known & ~failed).sum()), 'source_sha256': source['data_sha256'],
         'method': result['method'], 'photos_sha256': digest(photos), 'terminal_sha256': digest(terminal),
         'ecological_models_run': False, 'raw_photo_bytes_retained': False}
    tmp = report.with_suffix('.tmp'); tmp.write_text(json.dumps(r, indent=2, sort_keys=True) + '\n'); tmp.replace(report)
    print('MOIS_CACHED_PHOTO_UNIT_SAVED', json.dumps(r), flush=True)
    return r


def collect(source, workers, output):
    _, plan = source_paths(source, 0)
    frames, receipts = [], []
    for index in range(plan['shards']):
        files = list(Path(workers).rglob(f'shard_{index:02d}/receipt.json'))
        require(len(files) == 1, f'Missing or duplicate completed image unit {index}')
        r = json.loads(files[0].read_text()); f = files[0].parent / 'photos.csv'
        require(r['shard'] == index and r['source_sha256'] == plan['data_sha256']
                and r['attempted'] == 256 and digest(f) == r['photos_sha256'], 'Completed data fingerprint changed')
        data = read_csv(f); known = data.morph.isin(HUES); failed = data.measurement_status.eq('image_acquisition_failed')
        require(len(data) == 256 and int(known.sum()) == r['classified'] and int(failed.sum()) == r['acquisition_failed'],
                'Counts do not match acquired/failed/classified records')
        frames.append(data); receipts.append(r)
    require(len({json.dumps(r['method'], sort_keys=True) for r in receipts}) == 1, 'Mixed original image models')
    data = pd.concat(frames, ignore_index=True)
    original = read_csv(Path(source) / 'candidates.csv')
    require(len(data) == plan['count'] and data.photo_id.is_unique
            and data.measurement_id.tolist() == original.measurement_id.tolist(), 'Combined original photo roster changed')
    output = Path(output); require(not output.exists() or not any(output.iterdir()), 'Never overwrite a saved snapshot')
    output.mkdir(parents=True, exist_ok=True)
    data.to_parquet(output / 'measured_photos.parquet', index=False, compression='zstd')
    totals = {k: sum(r[k] for r in receipts) for k in ('attempted', 'downloaded', 'acquisition_failed',
               'classified', 'unclassified', 'acquired_unclassified')}
    report = {'schema': 'mois_cached_remaining_complete_data_only_v1', **totals,
              'source_sha256': plan['data_sha256'], 'parents': plan['parents'], 'method': receipts[0]['method'],
              'photo_data_sha256': digest(output / 'measured_photos.parquet'), 'ecological_models_run': False,
              'source_arm_preserved': True, 'raw_photo_bytes_retained': False}
    (output / 'receipt.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print('MOIS_REMAINING_POOL_DATA_SAVED', json.dumps(report), flush=True)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); cmd = p.add_subparsers(dest='command', required=True)
    prep = cmd.add_parser('prepare')
    for key in ('pool', 'second', 'completed', 'output'): prep.add_argument('--' + key, type=Path, required=True)
    prep.add_argument('--count', type=int, default=8192)
    work = cmd.add_parser('measure')
    for key in ('source', 'runtime', 'sam', 'output'): work.add_argument('--' + key, type=Path, required=True)
    work.add_argument('--shard', type=int, required=True)
    save = cmd.add_parser('collect')
    for key in ('source', 'workers', 'output'): save.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args()
    if a.command == 'prepare': prepare(a.pool, a.second, a.completed, a.output, a.count)
    elif a.command == 'measure': measure(a.source, a.shard, a.runtime, a.sam, a.output)
    else: collect(a.source, a.workers, a.output)
