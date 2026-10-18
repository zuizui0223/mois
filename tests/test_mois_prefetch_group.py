from pathlib import Path
import sys
import types
import pandas as pd


def test_completed_units_are_not_downloaded_again(tmp_path,monkeypatch):
    root=Path(__file__).parents[1]
    monkeypatch.syspath_prepend(str(root))
    from scripts.acquisition import mois_prefetch as mod
    loaded=[];measured=[];cached=[]
    method={'source':'unchanged'}
    core=types.ModuleType('scripts.acquisition.mois_collect_remaining_pool')
    def rows(source,i):return pd.DataFrame({'photo_url_large':[f'{i}-a',f'{i}-b']}),{}
    def measure(source,i,runtime,sam,out,runner=None):
        if runner is None:cached.append(i);return {'cached':i}
        measured.append(i);return runner(i)
    core.source_paths=rows;core.measure=measure
    fcp=types.ModuleType('scripts.acquisition.run_fcp_fresh_image_partition')
    def loader(u):loaded.append(u);return u.encode()
    def run_batch(i,loader):
        assert [loader(f'{i}-a'),loader(f'{i}-b')]==[f'{i}-a'.encode(),f'{i}-b'.encode()]
        return {'method':method}
    fcp.image_loader=loader;fcp.run_batch=run_batch;fcp.PermanentOriginalPhotoError=ValueError
    frozen=types.ModuleType('scripts.analysis.mois_reassemble_original_new4096_measurements')
    frozen.METHOD=method
    analysis=types.ModuleType('scripts.analysis');analysis.__path__=[]
    monkeypatch.setitem(sys.modules,core.__name__,core)
    monkeypatch.setitem(sys.modules,fcp.__name__,fcp)
    monkeypatch.setitem(sys.modules,'scripts.analysis',analysis)
    monkeypatch.setitem(sys.modules,frozen.__name__,frozen)
    path=tmp_path/'shard_01';path.mkdir();(path/'receipt.json').write_text('{}')
    mod.measure_group(tmp_path,[0,1,2],tmp_path,tmp_path,tmp_path)
    assert cached==[1]
    assert measured==[0,2]
    assert loaded==['0-a','0-b','2-a','2-b']
    assert not list((tmp_path/'prefetch_tmp').iterdir())
