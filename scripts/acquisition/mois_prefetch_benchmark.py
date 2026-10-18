"""Paired performance check: same original photo bytes, same measured I/O delays.

Download the first 16 source photos once with FCP's original paced loader.
Then replay each measured load duration and identical bytes in AB/BA order:
A=the original sequential run_batch loop; B=the very same loop and evaluator
fed by a one-thread bounded PrefetchLoader. Models, quality gates, source
ordering, temporary image retention and number of model initializations do
not change. This is a small controlled replay, NOT a production throughput
claim. Nothing here adds new photos to the ecological analysis denominator.
"""
from __future__ import annotations
import argparse
import contextlib
import hashlib
import json
from pathlib import Path
import tempfile
import time
import pandas as pd
from scripts.acquisition.mois_prefetch import PrefetchLoader

SOURCE='mois_global_angiosperm_second4096_source_only_original_photoids.csv'
SOURCE_SHA='3483808cc4ae773073e191baac282e8cc87f1dd8996b6369067f89980ab93cf3'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class OriginalPacket:
    def __init__(self, rows): self.rows=rows
    def get_batch(self, part, batch, size):
        if batch!=0 or size!=len(self.rows):
            raise ValueError('Benchmark original photo packet changed')
        cols=['measurement_id','image_filename','photo_license']
        return self.rows[cols].copy(),self.rows[cols+['photo_url_large']].copy(),1


def run(source_root, runtime, sam, output):
    from scripts.acquisition.run_fcp_fresh_image_partition import (
        image_loader, run_batch, PermanentOriginalPhotoError)
    source=source_root/SOURCE
    if digest(source)!=SOURCE_SHA:
        raise ValueError('Pinned completed-wave photo source changed')
    rows=pd.read_csv(source,dtype=str,keep_default_na=False).iloc[:16].copy()
    if len(rows)!=16 or not rows.photo_id.is_unique:
        raise ValueError('Benchmark source requires exactly 16 first original photos')
    if output.exists():
        raise ValueError('Do not overwrite a benchmark')
    output.mkdir(parents=True)
    trials=[];frames=[];captures=[]
    with tempfile.TemporaryDirectory(prefix='mois_paired_pixels_') as tmp:
        tmp=Path(tmp);stored={}
        for i,row in enumerate(rows.itertuples()):
            started=time.perf_counter();error=None;path=tmp/f'{i:02d}.img'
            try: payload=image_loader(row.photo_url_large)
            except PermanentOriginalPhotoError as exc:
                payload=None;error=exc
            elapsed=time.perf_counter()-started
            if payload is not None: path.write_bytes(payload)
            stored[row.photo_url_large]=(path,elapsed,error)
            captures.append(dict(photo_id=row.photo_id,download_seconds=elapsed,
                image_sha256=digest(path) if payload is not None else None,
                permanently_unavailable=error is not None))
        def replay(url):
            path,delay,error=stored[url]
            time.sleep(delay)
            if error is not None: raise error
            return path.read_bytes()
        for repeat,order in enumerate((('sequential','prefetch'),('prefetch','sequential'))):
            for mode in order:
                run_dir=output/f'{repeat}_{mode}';parts=[]
                context=(PrefetchLoader(rows.photo_url_large,replay,tmp,max_items=8,
                         permanent_errors=(PermanentOriginalPhotoError,))
                         if mode=='prefetch' else contextlib.nullcontext(replay))
                started=time.perf_counter()
                with context as loader:
                    for part in range(2):
                        packet=rows.iloc[part*8:(part+1)*8]
                        target=run_dir/f'part{part}'
                        result=run_batch(Path('.'),part,0,len(packet),target,
                            source=runtime/'scripts/analysis/measure_random_photo_first_blind_partition.py',
                            detector=runtime/'data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt',
                            sam=sam,loader=loader,indexed_source=OriginalPacket(packet))
                        terminal=target/f'terminal_p{part:03d}_b0000000.csv'
                        if digest(terminal)!=result['terminal_csv_sha256']:
                            raise ValueError('Original ROI result checksum mismatch')
                        parts.append(pd.read_csv(terminal,dtype=str,keep_default_na=False))
                elapsed=time.perf_counter()-started
                data=pd.concat(parts).sort_values('measurement_id').reset_index(drop=True)
                frames.append(data)
                trial=dict(repeat=repeat,mode=mode,seconds=elapsed,
                           photos_per_minute=len(rows)*60/elapsed,
                           classified=int(data.morph.isin(['white','yellow_orange','red_pink','blue_purple']).sum()))
                if mode=='prefetch':trial['prefetch']=loader.stats
                trials.append(trial)
                print('MOIS_PAIRED_PREFETCH_TRIAL',json.dumps(trial,sort_keys=True),flush=True)
        equal=all(frames[0].equals(f) for f in frames[1:])
        mismatches=[]
        for k,f in enumerate(frames[1:],1):
            for c in frames[0].columns:
                if c not in f or not frames[0][c].equals(f[c]):mismatches.append([k,c])
        elapsed={m:sum(t['seconds'] for t in trials if t['mode']==m)/2
                 for m in ('sequential','prefetch')}
        report=dict(schema='mois_prefetch_paired_real_photo_replay_v1',photos=16,
            source_sha256=SOURCE_SHA,trials=trials,captures=captures,
            mean_seconds=elapsed,speed_ratio=elapsed['sequential']/elapsed['prefetch'],
            time_reduction_fraction=1-elapsed['prefetch']/elapsed['sequential'],
            all_terminal_columns_exactly_equal=equal,mismatched_columns=mismatches,
            both_pairs_faster=all(next(t['seconds'] for t in trials if t['repeat']==i and t['mode']=='prefetch') <
                                 next(t['seconds'] for t in trials if t['repeat']==i and t['mode']=='sequential') for i in range(2)),
            actual_provider_downloads=16,benchmark_replays_do_not_add_new_data=True,
            measured_network_service_times_replayed_not_live_network_speedup=True,
            frozen_evaluator_and_rate_limit_unchanged=True,model_initializations_per_trial=2,
            small_two_unit_benchmark_not_full_256_photo_production_test=True)
        (output/'benchmark.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
        print('MOIS_PREFETCH_BENCHMARK_RESULT',json.dumps(report,sort_keys=True),flush=True)
        if not equal: raise ValueError('Paired output mismatch; do not enable prefetch')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('source','runtime','sam','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();run(a.source,a.runtime,a.sam,a.output)
