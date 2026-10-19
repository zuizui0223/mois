"""Sampling provenance for MOIS, without changing the photo/RGB executor.

Retain the full eligible upstream frame. Prioritize observation clusters with
independent, pre-colour baseline and stratum-top-up draws; every attached
photo inherits its observation's selection. These are acquisition probabilities
within a verified source, NOT probabilities that wild flowers are observed.
Legacy taxonomically selected images do not acquire invented global weights.
"""
from __future__ import annotations
import hashlib
import json
import math
import re
from typing import Any, Mapping

HUES = ('white', 'yellow_orange', 'red_pink', 'blue_purple')
HEX64 = re.compile(r'[0-9a-f]{64}')


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def _probability(value: Any, *, positive: bool = False) -> float:
    _require(not isinstance(value, bool), 'Probability must be numeric, not bool')
    p = float(value)
    _require(math.isfinite(p) and (0 < p <= 1 if positive else 0 <= p <= 1),
             'Probability outside its admissible range')
    return p


def _fingerprint(value: Mapping[str, Any]) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def make_design(frame: Mapping[str, Any], seed: str, baseline: float,
                topups: Mapping[str, float]) -> dict[str, Any]:
    """Freeze rates before outcomes. Count-only documents cannot pass this gate.

    The caller must obtain frame provenance from the independently validated
    full photo-ID index, not fabricate an index-complete flag. Top-up strata
    come from fixed metadata/covariates, never from measured colour. Missing
    environment remains eligible via the strictly positive baseline.
    """
    _require(frame.get('all_photo_ids_indexed') is True,
             'A verified complete photo-ID index is required, not only counts')
    _require(bool(HEX64.fullmatch(str(frame.get('source_sqlite_sha256', ''))))
             and bool(HEX64.fullmatch(str(frame.get('source_counts_sha256', '')))),
             'Missing source-index/count fingerprints')
    n = frame.get('source_photo_candidates')
    _require(isinstance(n, int) and not isinstance(n, bool) and n > 0,
             'Invalid full-frame photo count')
    _require(isinstance(seed, str) and bool(re.fullmatch('[0-9a-f]{32,128}', seed)),
             'Freeze a random hexadecimal seed of at least 128 bits')
    rates = {str(k): _probability(v) for k, v in topups.items()}
    _require(all(k.strip() for k in rates), 'Empty stratum label')
    d = {'schema': 'mois_observation_priority_design_v1',
         'frame_index_sha256': frame['source_sqlite_sha256'],
         'frame_counts_sha256': frame['source_counts_sha256'],
         'frame_photo_count': n, 'seed': seed,
         'baseline_probability': _probability(baseline, positive=True),
         'topup_probabilities': rates, 'sampling_unit': 'observation',
         'photo_cap': None, 'per_species_photo_cap': None,
         'per_observation_photo_cap': None,
         'scope': 'full_upstream_not_world_flora',
         'unselected_records_are_deferred_not_excluded': True,
         'rate_and_stratum_changes_after_inspecting_hues_not_permitted': True}
    d['design_sha256'] = _fingerprint(d)
    return d


def _uniform(seed: str, frame_sha: str, stream: str, observation: str) -> float:
    # Fixed hash ranks permit reproducible increasing-rate milestones. Treat
    # the independently frozen random seed as part of the sampling design.
    payload = json.dumps([seed, frame_sha, stream, observation], separators=(',', ':'))
    bits = int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8], 'big') >> 11
    return bits / (1 << 53)


def priority_record(photo_id: str, observation_id: str, stratum: str,
                    design: Mapping[str, Any]) -> dict[str, Any]:
    """Return provenance for EVERY photo, including deferred ones.

    One frozen stratum per observation is a caller invariant. All photos of
    that observation must be passed together to the executor. If milestones
    retain the seed/frame/strata and only raise p or q, selected sets are nested.
    Never stop halfway through an album and count it as a sampled observation.
    """
    d = dict(design)
    signed = d.pop('design_sha256', None)
    _require(signed == _fingerprint(d), 'Sampling design checksum changed')
    _require(d['schema'] == 'mois_observation_priority_design_v1', 'Unknown design')
    _require(str(photo_id).isdecimal() and int(photo_id) > 0
             and isinstance(observation_id, str) and bool(observation_id.strip())
             and isinstance(stratum, str) and bool(stratum.strip()),
             'Invalid photo, observation or stratum identity')
    p = _probability(d['baseline_probability'], positive=True)
    q = _probability(d['topup_probabilities'].get(stratum, 0.0))
    u = _uniform(d['seed'], d['frame_index_sha256'], 'baseline', observation_id)
    v = _uniform(d['seed'], d['frame_index_sha256'], 'topup', observation_id)
    baseline_draw, topup_draw = u < p, v < q
    inclusion = p + (1.0 - p) * q
    return {'photo_id': str(photo_id), 'observation_id': observation_id,
            'stratum': stratum, 'design_sha256': signed,
            'sampling_unit': 'observation', 'baseline_draw': baseline_draw,
            'topup_draw': topup_draw, 'ready_now': baseline_draw or topup_draw,
            'selection_probability': inclusion,
            'inverse_selection_weight': 1.0 / inclusion,
            'scope': d['scope'], 'all_source_photos_retained': True,
            'source_status': 'metadata_only', 'new_image_measured': False,
            'does_not_correct_observer_or_classifier_selection': True}


def legacy_result_record(row: Mapping[str, Any], cohort: str,
                         counts_sha256: str) -> dict[str, Any]:
    """Sidecar only: keep the original measurements and their source bytes.

    The previous genus/species-first pool was not a probability sample from
    all upstream records. Missing inclusion probability is null, never 1 or
    n_processed/93,209,351. ROI failures do not certify the absence of flowers.
    """
    needed = ('photo_id', 'observation_id', 'observer_id', 'inat_taxon_id',
              'morph', 'measurement_status', 'image_sha256')
    _require(all(k in row for k in needed) and bool(cohort.strip()),
             'Missing collection identity or measurement fields')
    _require(bool(HEX64.fullmatch(counts_sha256)), 'Invalid counts fingerprint')
    _require(str(row['photo_id']).isdecimal() and int(row['photo_id']) > 0
             and bool(str(row['observation_id']).strip()), 'Invalid photo identity')
    known = row['morph'] in HUES
    failed = row['measurement_status'] == 'image_acquisition_failed'
    _require(row['morph'] in (*HUES, 'mixed_uncertain')
             and known == (row['measurement_status'] == 'classified_four_state_morph')
             and not (known and failed), 'Inconsistent hue/measurement status')
    image_sha = str(row['image_sha256'])
    _require(image_sha == '' if failed else bool(HEX64.fullmatch(image_sha)),
             'Acquired-image hash or acquisition failure inconsistent')
    state = ('acquisition_failed' if failed else
             'classified' if known else 'acquired_unclassifiable')
    return {'photo_id': str(row['photo_id']),
            'observation_id': str(row['observation_id']),
            'observer_id': str(row['observer_id']),
            'inat_taxon_id': str(row['inat_taxon_id']), 'source_cohort': cohort,
            'upstream_counts_sha256': counts_sha256,
            'collection_design': 'nonprobability_or_unrecoverable',
            'global_selection_probability': None,
            'global_inverse_probability_weight': None,
            'result_state': state, 'original_measurement_status': row['measurement_status'],
            'botanical_flower_status': 'not_reviewed',
            'flowering_annotation': 'not_recorded',
            'colour_class': row['morph'] if known else None,
            'human_review_status': 'not_reviewed'}
