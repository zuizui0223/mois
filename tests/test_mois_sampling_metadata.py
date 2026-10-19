"""Collection-design contracts, not a new floral ecology model."""
import json
from pathlib import Path
import pytest
from scripts.acquisition import mois_sampling_metadata as m


def frame(indexed=True):
    return {'all_photo_ids_indexed':indexed,'source_photo_candidates':93209351,
            'source_sqlite_sha256':'a'*64,'source_counts_sha256':'b'*64}


def original(status='classified_four_state_morph', hue='white'):
    return {'photo_id':'123','observation_id':'uuid-1','observer_id':'7',
            'inat_taxon_id':'12','morph':hue,'measurement_status':status,
            'image_sha256': '' if status=='image_acquisition_failed' else 'e'*64}


def test_counts_only_never_become_a_sampling_frame():
    with pytest.raises(ValueError, match='photo-ID'):
        m.make_design(frame(False), '1'*32, .1, {'island':.3})


def test_union_probability_and_observation_cluster_decisions():
    d=m.make_design(frame(), '1'*32, .1, {'island':.5})
    a=m.priority_record('1','same-observation','island',d)
    b=m.priority_record('2','same-observation','island',d)
    assert a['selection_probability']==pytest.approx(.55)
    assert a['ready_now']==b['ready_now']
    assert a['baseline_draw']==b['baseline_draw']
    assert a['topup_draw']==b['topup_draw']
    assert a['sampling_unit']=='observation'
    assert a['photo_id']!=b['photo_id']
    assert a['all_source_photos_retained'] is True
    assert a['new_image_measured'] is False


def test_positive_probability_for_every_stratum_and_no_unknown_drop():
    d=m.make_design(frame(), '1'*32, .01, {'island':.2})
    a=m.priority_record('1','obs','unknown_environment',d)
    assert a['selection_probability']==pytest.approx(.01)
    assert a['scope']=='full_upstream_not_world_flora'
    with pytest.raises(ValueError): m.make_design(frame(),'1'*32,0,{})
    with pytest.raises(ValueError): m.make_design(frame(),'1'*32,.1,{'bad':1.1})
    with pytest.raises(ValueError): m.make_design(frame(),'short',.1,{})


def test_progress_to_full_census_is_nested_and_keeps_all_photos():
    early=m.make_design(frame(),'c'*32,.1,{'island':.2})
    later=m.make_design(frame(),'c'*32,.3,{'island':.5})
    full=m.make_design(frame(),'c'*32,1,{})
    for i in range(300):
        a=m.priority_record(str(i+1),'obs'+str(i),'island',early)
        b=m.priority_record(str(i+1),'obs'+str(i),'island',later)
        c=m.priority_record(str(i+1),'obs'+str(i),'island',full)
        assert not a['ready_now'] or b['ready_now']
        assert c['ready_now'] and c['selection_probability']==1
        assert c['all_source_photos_retained']


def test_failure_states_remain_distinct_and_legacy_probability_unknown():
    rows=[original(),original('not_evaluable_roi_or_flip_gate','mixed_uncertain'),
          original('image_acquisition_failed','mixed_uncertain')]
    a,b,c=[m.legacy_result_record(r,'legacy_breadth', 'b'*64) for r in rows]
    assert [r['result_state'] for r in (a,b,c)]==[
        'classified','acquired_unclassifiable','acquisition_failed']
    assert all(r['global_selection_probability'] is None for r in (a,b,c))
    assert all(r['global_inverse_probability_weight'] is None for r in (a,b,c))
    assert b['botanical_flower_status']=='not_reviewed'
    assert c['botanical_flower_status']=='not_reviewed'
    assert a['flowering_annotation']=='not_recorded'


def test_cannot_silently_recode_missing_hue_as_white():
    with pytest.raises(ValueError):
        m.legacy_result_record(original('image_acquisition_failed','white'),'x','b'*64)
    with pytest.raises(ValueError):
        m.legacy_result_record(original('not_evaluable_roi_or_flip_gate','white'),'x','b'*64)


def test_no_taxonomic_or_geographic_inclusion_probability_from_candidate_total():
    r=m.legacy_result_record(original(),'taxon_selected','b'*64)
    assert r['global_selection_probability'] is None
    assert r['collection_design']=='nonprobability_or_unrecoverable'
    d=m.make_design(frame(),'1'*32,.1,{})
    assert d['frame_photo_count']==93209351
    assert d['photo_cap'] is None and d['per_species_photo_cap'] is None


def test_design_tamper_and_incomplete_metadata_rejected():
    d=m.make_design(frame(),'1'*32,.1,{})
    d['baseline_probability']=.2
    with pytest.raises(ValueError,match='checksum'):
        m.priority_record('1','obs','x',d)
    with pytest.raises(ValueError):
        m.legacy_result_record({'photo_id':'1'},'legacy','b'*64)
