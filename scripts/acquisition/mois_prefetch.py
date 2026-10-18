"""Overlap the next photo downloads with the unchanged current ROI evaluation.

One producer calls the original, paced/retrying FCP image loader. No extra
HTTP concurrency or model threads are introduced. A bounded private disk
queue feeds the existing run_batch(loader=...) interface in exactly its
original source order. Images are ephemeral; outputs and quality gates
remain those of the original evaluator. Use one context across several
256-photo units, not one context per image or per unit.
"""
from __future__ import annotations

from collections import deque
from pathlib import Path
import shutil
import tempfile
import threading
import time
from typing import Callable, Iterable


class PrefetchLoader:
    def __init__(self, urls: Iterable[str], loader: Callable[[str], bytes],
                 root: Path, *, max_items: int = 256,
                 max_bytes: int = 512 * 1024**2,
                 max_image_bytes: int = 12 * 1024**2,
                 permanent_errors: tuple[type[Exception], ...] = ()):
        self.urls = tuple(urls)
        if len(set(self.urls)) != len(self.urls):
            raise ValueError('Duplicate original photo URLs')
        if (max_items < 1 or max_image_bytes < 1 or max_bytes < max_image_bytes):
            raise ValueError('Invalid bounded photo queue limits')
        self.loader, self.root = loader, Path(root)
        self.max_items, self.max_bytes = max_items, max_bytes
        self.max_image_bytes = max_image_bytes
        self.permanent_errors = permanent_errors
        self._cv = threading.Condition()
        self._queue = deque()
        self._used_bytes = self._live = self._position = 0
        self._stop = self._done = self._entered = False
        self.stats = dict(download_calls=0, consumed=0, permanent_failures=0,
                          peak_bytes=0, peak_items=0, loader_seconds=0.0,
                          consumer_wait_seconds=0.0)

    def __enter__(self):
        if self._entered:
            raise RuntimeError('Photo prefetch context cannot be reused')
        self._entered = True
        self.root.mkdir(parents=True, exist_ok=True)
        self.directory = Path(tempfile.mkdtemp(prefix='mois_prefetch_', dir=self.root))
        self._thread = threading.Thread(target=self._produce, name='mois-paced-photo-loader')
        self._thread.start()
        return self

    def _produce(self):
        try:
            for index, url in enumerate(self.urls):
                with self._cv:
                    self._cv.wait_for(lambda: self._stop or (
                        self._live < self.max_items and
                        self._used_bytes + self.max_image_bytes <= self.max_bytes))
                    if self._stop:
                        break
                    self._live += 1
                    self._used_bytes += self.max_image_bytes
                path, error, nbytes, fatal = None, None, 0, False
                started = time.perf_counter()
                try:
                    self.stats['download_calls'] += 1
                    payload = self.loader(url)
                    if not isinstance(payload, bytes) or not 0 < len(payload) <= self.max_image_bytes:
                        raise RuntimeError('Original loader returned invalid/unbounded image bytes')
                    nbytes = len(payload)
                    path = self.directory / f'{index:08d}.img'
                    path.write_bytes(payload)
                    del payload
                except Exception as exc:
                    error = exc
                    fatal = not isinstance(exc, self.permanent_errors)
                    if not fatal:
                        self.stats['permanent_failures'] += 1
                    if path is not None:
                        path.unlink(missing_ok=True)
                    path, nbytes = None, 0
                finally:
                    self.stats['loader_seconds'] += time.perf_counter() - started
                with self._cv:
                    self._used_bytes -= self.max_image_bytes - nbytes
                    self.stats['peak_bytes'] = max(self.stats['peak_bytes'], self._used_bytes)
                    self.stats['peak_items'] = max(self.stats['peak_items'], self._live)
                    self._queue.append((url, path, error, nbytes))
                    self._cv.notify_all()
                    if fatal or self._stop:
                        break
        finally:
            with self._cv:
                self._done = True
                self._cv.notify_all()

    def __call__(self, url: str) -> bytes:
        if not self._entered or self._stop:
            raise RuntimeError('Photo prefetch context is not active')
        if self._position >= len(self.urls) or self.urls[self._position] != url:
            raise ValueError('Original photo source order changed')
        started = time.perf_counter()
        with self._cv:
            self._cv.wait_for(lambda: self._queue or self._done)
            self.stats['consumer_wait_seconds'] += time.perf_counter() - started
            if not self._queue:
                raise RuntimeError('Photo producer ended before the expected original ID')
            queued_url, path, error, nbytes = self._queue.popleft()
        try:
            if queued_url != url:
                raise RuntimeError('Photo prefetch identity mismatch')
            self._position += 1
            self.stats['consumed'] += 1
            if error is not None:
                raise error
            return path.read_bytes()
        finally:
            if path is not None:
                path.unlink(missing_ok=True)
            with self._cv:
                self._used_bytes -= nbytes
                self._live -= 1
                self._cv.notify_all()

    def __exit__(self, exc_type, exc, tb):
        with self._cv:
            self._stop = True
            self._cv.notify_all()
        # The original loader has bounded HTTP timeouts/retries. Join it so
        # no producer accesses deleted files or continues after this context.
        self._thread.join()
        shutil.rmtree(self.directory)
        return False


def measure_group(source, indices, runtime, sam, output):
    """Opt-in replacement for the old sequential shell loop over four units."""
    import json
    from scripts.acquisition import mois_collect_remaining_pool as core
    from scripts.acquisition.run_fcp_fresh_image_partition import (
        run_batch, image_loader, PermanentOriginalPhotoError)
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

    indices = list(indices)
    if len(set(indices)) != len(indices):
        raise ValueError('Repeated photo unit in this group')
    pending, reports = [], []
    for index in indices:
        rows, plan = core.source_paths(source, index)
        if (Path(output) / f'shard_{index:02d}' / 'receipt.json').is_file():
            reports.append(core.measure(source, index, runtime, sam, output))
        else:
            pending.append((index, rows))
    if not pending:
        return reports
    urls = [str(u) for _, rows in pending for u in rows.photo_url_large]
    started = time.perf_counter()
    with PrefetchLoader(urls, image_loader, Path(output)/'prefetch_tmp',
                        permanent_errors=(PermanentOriginalPhotoError,)) as prefetched:
        def runner(*args, **kwargs):
            result = run_batch(*args, loader=prefetched, **kwargs)
            if result['method'] != METHOD:
                raise ValueError('Frozen photo-colour model changed')
            return result
        for index, _ in pending:
            reports.append(core.measure(source, index, runtime, sam, output, runner=runner))
        timing = dict(prefetched.stats)
    timing.update(wall_seconds=time.perf_counter()-started,
                  freshly_processed_units=len(pending), reused_units=len(indices)-len(pending),
                  model_and_thresholds_unchanged=True, network_producers=1)
    Path(output).mkdir(parents=True,exist_ok=True)
    (Path(output)/f'prefetch_group_{indices[0]:02d}.json').write_text(
        json.dumps(timing,indent=2,sort_keys=True)+'\n')
    print('MOIS_PREFETCH_GROUP',json.dumps(timing,sort_keys=True),flush=True)
    return reports


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('source','runtime','sam','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--group',type=int,required=True)
    p.add_argument('--units-per-group',type=int,default=4)
    a=p.parse_args()
    if a.group<0 or a.units_per_group<1:
        p.error('Invalid group or units-per-group')
    start=a.group*a.units_per_group
    measure_group(a.source,range(start,start+a.units_per_group),a.runtime,a.sam,a.output)
