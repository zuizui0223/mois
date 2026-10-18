import importlib.util
import threading
import time
from pathlib import Path
import pytest

PATH = Path(__file__).parents[1]/'scripts/acquisition/mois_prefetch.py'

def make(*args, **kwargs):
    assert PATH.is_file(), 'The bounded prefetch implementation is not present'
    spec=importlib.util.spec_from_file_location('prefetch',PATH)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod.PrefetchLoader(*args, **kwargs)


def test_download_next_batch_while_current_batch_is_measured(tmp_path):
    next_downloaded=threading.Event()
    def source(url):
        if url=='b': next_downloaded.set()
        return url.encode()
    with make(['a','b'],source,tmp_path,max_items=2,max_bytes=32,max_image_bytes=16) as p:
        assert p('a')==b'a'
        assert next_downloaded.wait(1), 'Download waits until the next inference ends'
        assert p('b')==b'b'
    assert list(tmp_path.iterdir())==[]


def test_order_and_original_bytes_unchanged_with_one_download_each(tmp_path):
    calls=[]
    def source(url): calls.append(url);return url.encode()
    with make(['a','b','c'],source,tmp_path,max_items=2,max_bytes=32,max_image_bytes=16) as p:
        assert [p(u) for u in ('a','b','c')]==[b'a',b'b',b'c']
        assert p.stats['consumed']==3
    assert calls==['a','b','c']


def test_bounded_prefetch_stops_before_unbounded_network_and_disk(tmp_path):
    calls=[]; two=threading.Event()
    def source(url):
        calls.append(url)
        if len(calls)==2:two.set()
        return b'x'*16
    with make(list('abcd'),source,tmp_path,max_items=2,max_bytes=32,max_image_bytes=16) as p:
        assert two.wait(1)
        time.sleep(.05)
        assert calls==['a','b']
        assert sum(f.stat().st_size for f in p.directory.iterdir())<=32
        assert [p(u) for u in 'abcd']==[b'x'*16]*4
        assert p.stats['peak_bytes']<=32


def test_permanent_missing_photo_is_passed_to_original_handler(tmp_path):
    class Permanent(ValueError):pass
    def source(url):
        if url=='bad':raise Permanent('original object 404')
        return b'actual'
    with make(['bad','good'],source,tmp_path,max_items=2,max_bytes=32,max_image_bytes=16,
              permanent_errors=(Permanent,)) as p:
        with pytest.raises(Permanent,match='404'):p('bad')
        assert p('good')==b'actual'


def test_provider_outage_stops_prefetch_instead_of_fake_missing_photos(tmp_path):
    calls=[]
    class Outage(RuntimeError):pass
    def source(url):
        calls.append(url)
        raise Outage('429 retry budget exhausted')
    with pytest.raises(Outage,match='429'):
        with make(['first','later'],source,tmp_path,max_items=2,max_bytes=32,max_image_bytes=16) as p:
            p('first')
    assert calls==['first']
    assert list(tmp_path.iterdir())==[]


def test_wrong_order_does_not_use_a_different_photo(tmp_path):
    with make(['a','b'],lambda u:u.encode(),tmp_path,max_items=2,max_bytes=32,max_image_bytes=16) as p:
        with pytest.raises(ValueError,match='order'):p('b')
        assert p('a')==b'a'
        assert p('b')==b'b'


def test_cleanup_and_original_inference_exception_is_not_masked(tmp_path):
    with pytest.raises(RuntimeError,match='model failed'):
        with make(list('abcdef'),lambda _:b'x',tmp_path,max_items=2,max_bytes=32,max_image_bytes=16) as p:
            p('a')
            raise RuntimeError('model failed')
    assert list(tmp_path.iterdir())==[]


def test_invalid_or_duplicate_requests_rejected(tmp_path):
    with pytest.raises(ValueError):
        make(['a','a'],lambda _:b'x',tmp_path,max_items=2,max_bytes=32,max_image_bytes=16)
    with pytest.raises(ValueError):
        make(['a'],lambda _:b'x',tmp_path,max_items=2,max_bytes=2,max_image_bytes=16)
