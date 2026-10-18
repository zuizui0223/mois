"""Collect cached MOIS wave2 images and fixed-model colours; no ecological analysis.

Uses the existing FCP paced image loader and unchanged ROI-v4 evaluator.
One evaluator invocation handles 256 photos. Completed shards are verified
and reused, and original failed/unknown outcomes are retained.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import pandas as pd

PLAN = 'mois_global_angiosperm_second4096_source_only_original_photoids.csv'
RECEIPT = 'mois_global_angiosperm_cached_second4096_source_only_receipt.json'
HUES = ('white', 'yellow_orange', 'red_pink', 'blue_purple')


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


def verify_source(root, index):
    root = Path(root)
    require(index in range(16), 'Invalid photo shard')
    r = json.loads((root / RECEIPT).read_text())
    require(r['schema'] == 'mois_angiosperm_cached_second4096_new_species_source_only_v1'
            and r['status'] == 'SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES'
            and r['new_photo_pixels_actually_measured'] is False
            and r['second4096_photo_candidates'] == 4096
            and r['first4096_and_new4096_disjoint_photo_IDs'] is True
            and r['original313535_measured_photo_IDs_excluded_via_parent_signed_source'] is True
            and r['other_fifth4096_and_SI4096_reservations_excluded_via_parent'] is True,
            'Unverified or previously used source')
    shard = root / f'mois_broad_secondwave_source_shard_{index:02d}_256.csv'
    require(digest(root / PLAN) == r['source_only_second4096_csv_sha256']
            and digest(shard) == r['source_only_shard_SHA256'][str(index)],
            'Source checksum mismatch')
    full, rows = read_csv(root / PLAN), read_csv(shard)
    require(len(full) == 4096 and full.photo_id.is_unique
            and full.measurement_id.is_unique and full.observation_id.is_unique
            and len(rows) == 256 and rows.photo_id.is_unique
            and rows.measurement_id.is_unique
            and rows.equals(full[full.shard.eq(str(index))].reset_index(drop=True)),
            'Duplicate, incomplete or changed photo roster')
    require(rows.source_snapshot_date.eq(r['parent_original_metadata_snapshot_date']).all()
            and not {'morph', 'hue', 'classified_four_hue'}.intersection(rows.columns),
            'Source dates or unseen colour fields changed')
    for x in rows.itertuples():
        require(re.fullmatch(r'https://inaturalist-open-data\.s3\.amazonaws\.com/photos/'
                             + re.escape(x.photo_id) + r'/large\.(jpg|jpeg|png)', x.photo_url_large)
                and x.image_filename == x.measurement_id + '.jpg',
                'Unsafe or mismatched original photo URL')
    return rows, r


class PhotoBatch:
    def __init__(self, rows, index):
        self.rows, self.index = rows, index

    def get_batch(self, part, batch, size):
        require((part, batch, size) == (self.index, 0, 256), 'Not a single 256-photo unit')
        fields = ['measurement_id', 'image_filename', 'photo_license']
        return self.rows[fields].copy(), self.rows[fields + ['photo_url_large']].copy(), 1


def measure(root, index, runtime, sam, out, runner=None):
    rows, source = verify_source(root, index)
    target = Path(out) / f'shard_{index:02d}'
    receipt_file, photo_file = target / 'receipt.json', target / 'photos.csv'
    if receipt_file.exists():
        r = json.loads(receipt_file.read_text())
        require(r['source_sha256'] == source['source_only_second4096_csv_sha256']
                and r['shard_source_sha256'] == source['source_only_shard_SHA256'][str(index)]
                and r['shard'] == index and r['attempted'] == 256
                and photo_file.is_file() and digest(photo_file) == r['photos_sha256'],
                'Completed cache corrupted; do not replace or repeat it')
        return r
    production = runner is None
    if production:
        from scripts.acquisition.run_fcp_fresh_image_partition import run_batch
        runner = run_batch
    target.mkdir(parents=True, exist_ok=True)
    sealed = target / 'sealed'
    runtime = Path(runtime)
    result = runner(Path('.'), index, 0, 256, sealed,
                    source=runtime / 'scripts/analysis/measure_random_photo_first_blind_partition.py',
                    detector=runtime / 'data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt',
                    sam=Path(sam), indexed_source=PhotoBatch(rows, index))
    require(result['status'] == 'COMPLETE_WORK_UNIT'
            and result['n_source_photo_ids'] == result['n_terminal_photo_ids'] == 256
            and result['image_bytes_persisted'] is False
            and result['reused_historical_measurement_csv'] is False
            and result['number_batches_for_this_partition'] == 1,
            'Image evaluator did not return a complete original photo unit')
    if production:
        from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
        require(result['method'] == METHOD, 'Frozen model/source settings changed')
    terminal = sealed / f'terminal_p{index:03d}_b0000000.csv'
    require(digest(terminal) == result['terminal_csv_sha256'], 'Pixel-result checksum mismatch')
    measured = read_csv(terminal)
    require(len(measured) == 256 and measured.measurement_id.is_unique
            and set(measured.measurement_id) == set(rows.measurement_id), 'Photo outcome IDs changed')
    fields = ['measurement_id', 'morph', 'measurement_status', 'image_sha256']
    fields += [c for c in measured.columns if c.startswith('flower_fraction_')]
    joined = rows.merge(measured[fields], on='measurement_id', validate='one_to_one', sort=False)
    known = joined.morph.isin(HUES)
    require(joined.morph.isin((*HUES, 'mixed_uncertain')).all()
            and known.eq(joined.measurement_status.eq('classified_four_state_morph')).all()
            and int(known.sum()) == result['n_classifiable'], 'Unknown/photo-colour status conflict')
    failed = joined.measurement_status.eq('image_acquisition_failed')
    require(int(failed.sum()) == result['n_acquisition_failed'], 'Photo acquisition count mismatch')
    joined['classified_four_hue'] = known.map({True: 'true', False: 'false'})
    joined.to_csv(photo_file, index=False, lineterminator='\n')
    r = {'schema': 'mois_cached_wave2_image_collection_v1', 'shard': index,
         'attempted': 256, 'downloaded': int((~failed).sum()), 'acquisition_failed': int(failed.sum()),
         'classified': int(known.sum()), 'unclassified': int((~known).sum()),
         'source_sha256': source['source_only_second4096_csv_sha256'],
         'shard_source_sha256': source['source_only_shard_SHA256'][str(index)],
         'photos_sha256': digest(photo_file), 'terminal_sha256': digest(terminal),
         'method': result['method'], 'ecological_models_run': False,
         'raw_photo_bytes_retained': False, 'model_processes_per256': 1}
    temp = receipt_file.with_suffix('.tmp')
    temp.write_text(json.dumps(r, indent=2, sort_keys=True) + '\n'); temp.replace(receipt_file)
    print('MOIS_WAVE2_MEASURED_SHARD', json.dumps(r), flush=True)
    return r


def collect(root, output):
    frames, receipts = [], []
    for i in range(16):
        paths = list(Path(root).rglob(f'shard_{i:02d}/receipt.json'))
        require(len(paths) == 1, f'Missing or repeated completed shard {i}')
        r = json.loads(paths[0].read_text()); f = paths[0].parent / 'photos.csv'
        require(r['shard'] == i and r['attempted'] == 256 and digest(f) == r['photos_sha256'],
                'Corrupt completed data')
        frame = read_csv(f)
        require(len(frame) == 256 and frame.photo_id.is_unique
                and int(frame.morph.isin(HUES).sum()) == r['classified'], 'Outcome count drift')
        receipts.append(r); frames.append(frame)
    require(len({r['source_sha256'] for r in receipts}) == 1
            and len({json.dumps(r['method'], sort_keys=True) for r in receipts}) == 1,
            'Different source or model across workers')
    data = pd.concat(frames, ignore_index=True)
    require(len(data) == 4096 and data.photo_id.is_unique and data.measurement_id.is_unique,
            'Duplicate photo attempts across shards')
    output = Path(output)
    require(not output.exists() or not any(output.iterdir()), 'Do not overwrite data snapshot')
    output.mkdir(parents=True, exist_ok=True)
    f = output / 'wave2_measured_photos.parquet'
    data.to_parquet(f, index=False, compression='zstd')
    r = {'schema': 'mois_cached_wave2_complete_data_only_v1',
         **{k: sum(x[k] for x in receipts) for k in
            ('attempted', 'downloaded', 'acquisition_failed', 'classified', 'unclassified')},
         'source_sha256': receipts[0]['source_sha256'], 'method': receipts[0]['method'],
         'photo_data_sha256': digest(f), 'ecological_models_run': False,
         'source_arm_preserved': True, 'raw_photo_bytes_retained': False}
    (output / 'receipt.json').write_text(json.dumps(r, indent=2, sort_keys=True) + '\n')
    print('MOIS_WAVE2_DATA_SAVED', json.dumps(r), flush=True)
    return r


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    work = commands.add_parser('measure')
    for name in ('source', 'runtime', 'sam', 'output'):
        work.add_argument('--' + name, type=Path, required=True)
    work.add_argument('--shard', type=int, required=True)
    merge = commands.add_parser('collect')
    merge.add_argument('--workers', type=Path, required=True)
    merge.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'measure':
        measure(args.source, args.shard, args.runtime, args.sam, args.output)
    else:
        collect(args.workers, args.output)
