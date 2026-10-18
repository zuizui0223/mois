"""Second broad 4096 source IDs must not re-scan metadata or repeat first taxa."""
import pandas as pd
import pytest
from scripts.acquisition import mois_cached_angiosperm_second4096_source as m
from scripts.acquisition import mois_global_angiosperm_breadth_source as core


def photo(n,taxon,genus,obs=None,observer=None):
    r={key:"" for key in core.COLUMNS}
    r.update({
        "photo_id":str(n),
        "inat_taxon_id":str(taxon),
        "species":genus+" species"+str(taxon),
        "genus":genus,
        "grid10":f"{n%8:02d}:06",
        "observation_id":str(obs if obs is not None else n+100000),
        "observer_id":str(observer if observer is not None else n+200000),
        "measurement_id":"MOIS-"+str(n).zfill(24),
        "source_arm":"mois_broad_angiosperm_species_first_official_open_S3",
        "source_snapshot_date":"2026-10-18",
        "photo_license":"cc-by",
    })
    return r


def test_next_source_is_disjoint_photo_and_species_and_prefers_new_genus():
    first=pd.DataFrame([photo(100+i,100+i,"Former"+str(i)) for i in range(4)])
    new=[photo(1000+i,1000+i,"Novel"+str(i)) for i in range(14)]
    old_genus_new_species=photo(5000,5000,"Former1")
    pool=pd.DataFrame(first.to_dict("records")+new+[old_genus_new_species])
    selected,status=m.choose_next(pool,first,goal=8)
    assert len(selected)==8
    assert len({r["photo_id"] for r in selected})==8
    assert len({r["genus"] for r in selected})==8
    assert not {r["photo_id"] for r in selected}.intersection(set(first.photo_id))
    assert not {r["inat_taxon_id"] for r in selected}.intersection(
        set(first.inat_taxon_id))
    assert all("novel_genus" in r["source_sampling_priority"] for r in selected)
    assert sum(status.values())==8
    assert sorted(r["shard"] for r in selected)==list(map(str,range(8)))


def test_photographer_reuse_fallback_keeps_species_unique():
    first=pd.DataFrame([photo(1,1,"A"),photo(2,2,"B")])
    fresh=[photo(10+i,10+i,"G"+str(i),observer=9000)
           for i in range(6)]
    chosen,_=m.choose_next(pd.DataFrame(first.to_dict("records")+fresh),first,goal=4)
    assert len(chosen)==4
    assert len({x["inat_taxon_id"] for x in chosen})==4
    assert len({x["observer_id"] for x in chosen})==1
    assert len({x["photo_id"] for x in chosen})==4


def test_rank_is_hue_blind_even_for_untrusted_accidental_extra_columns():
    first=pd.DataFrame([photo(100,100,"Old")])
    new=[photo(1000+i,1000+i,"New"+str(i)) for i in range(8)]
    pool=pd.DataFrame(first.to_dict("records")+new)
    a,_=m.choose_next(pool,first,goal=5)
    pool["morph"]=["white","red_pink"]*4+["white"]
    pool["precipitation"]=1234
    b,_=m.choose_next(pool,first,goal=5)
    assert [x["photo_id"] for x in a]==[x["photo_id"] for x in b]
    assert all("morph" not in x and "precipitation" not in x for x in b)


def test_source_only_not_actual_four_hue_colour():
    assert m.SECOND_N==4096
    assert m.PER_SHARD*m.SHARDS==m.SECOND_N
    assert m.BASE_SOURCE_TAG.startswith("mois-angiosperm-broad-")
    assert "PHOTO IDS" in m.__doc__
    with pytest.raises(ValueError,match="source-locked"):
        m.choose_next(pd.DataFrame([photo(100,1,"A")]),
                      pd.DataFrame([photo(200,2,"B")]),goal=1)
