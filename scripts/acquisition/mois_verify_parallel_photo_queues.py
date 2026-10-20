#!/usr/bin/env python3
"""Verify exact pending photo/observation queues are disjoint before parallel RGB.

Reads only source manifests and small worker CSVs. Pending targets are not
measured photos. No image acquisition, colour inference, or source mutation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

SCHEMA = 'mois_full_source_observation_sampling_cohort_v1'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def source_ids(directory):
    directory = Path(directory)
    path = directory / 'manifest.json'
    plan = json.loads(path.read_text(encoding='utf-8'))
    require(plan.get('schema') == SCHEMA and plan.get('colour_measured') is False
            and plan.get('full_frame_photo_ids') == 93209351,
            'Expected full-index pre-colour source queue')
    workers = plan['worker_ids']
    require(workers == list(range(plan['work_units'])), 'Incomplete worker index/count')
    require(set(plan['photo_shard_sha256']) == {str(i) for i in workers},
            'Incomplete worker SHA roster')
    photos, observations, measurements = set(), set(), set()
    for index in workers:
        file = directory / 'shards' / f'shard-{index:03d}.csv'
        require(sha256(file) == plan['photo_shard_sha256'][str(index)],
                'Source shard SHA mismatch')
        with file.open(newline='', encoding='utf-8') as stream:
            reader = csv.DictReader(stream)
            require({'photo_id', 'observation_id', 'measurement_id', 'shard'} <= set(reader.fieldnames or []),
                    'Missing source identity columns')
            count = 0
            for row in reader:
                pid, obs, mid = row['photo_id'], row['observation_id'], row['measurement_id']
                require(pid.isdecimal() and obs.strip() and mid.strip(), 'Missing/invalid source identity')
                require(row['shard'] == str(index), 'Wrong source shard identity')
                require(pid not in photos and mid not in measurements,
                        'Duplicate photo/measurement ID within source queue')
                photos.add(pid); observations.add(obs); measurements.add(mid)
                count += 1
            require(0 < count <= 256, 'Invalid original source work-unit count')
    require(len(photos) == plan['new_photo_attempts_planned'], 'Source photo count mismatch')
    return photos, observations, sha256(path)


def verify_disjoint(new_source, active_source):
    new_photos, new_obs, new_sha = source_ids(new_source)
    active_photos, active_obs, active_sha = source_ids(active_source)
    require(not new_photos & active_photos, 'Concurrent photo ID overlap: refuse repeated download')
    require(not new_obs & active_obs, 'Concurrent observation ID overlap: queues are not independent')
    return {
        'schema': 'mois_parallel_original_photo_queue_disjointness_v1',
        'new_photo_targets': len(new_photos), 'active_photo_targets': len(active_photos),
        'new_observations': len(new_obs), 'active_observations': len(active_obs),
        'overlapping_photo_ids': 0, 'overlapping_observation_ids': 0,
        'new_source_manifest_sha256': new_sha, 'active_source_manifest_sha256': active_sha,
        'active_targets_counted_as_measured_photos': False, 'image_pixels_read': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--new-source', type=Path, required=True)
    parser.add_argument('--active-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'Refuse to overwrite a source disjointness receipt')
    result = verify_disjoint(args.new_source, args.active_source)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print('MOIS_PARALLEL_REAL_PHOTO_QUEUES_DISJOINT', json.dumps(result, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
