"""Contracts for authentic photographic ID union, not synthetic ecology results."""
import pandas as pd
import pytest

from scripts.acquisition import mois_union_FCP594081_first93m18462_real_original_photo_20261020 as m


def previous(pid, tax, hue, obs, arm="core"):
    known=hue in m.HUES
    return {
        "photo_id":str(pid),"inat_taxon_id":str(tax),
        "species":"Campanula punctata" if tax==10 else "Cirsium brevicaule",
        "observation_id":obs,
        "latitude":37.4,"longitude":139.4,"observed_on":"2025-07-01",
        "morph":hue,"classified":known,"white":hue=="white",
        "arm":arm,
        "measurement_status":"classified_four_state_morph" if known else "mixed_uncertain",
        "technical_failure":False,
        "photo_colour_provenance":"original_MOIS_source_ROI_v4_prior350399",
        "source_original_fcp_morph":"",
        "source_fcp_global_classifiable":False,"fcp_source_cohort":""
    }


def new_source(pid, tax, hue, obs, status=None):
    known=hue in m.HUES
    d=previous(pid,tax,hue,obs,arm="global93m_probability_observation")
    d["measurement_status"]=status or (
        "classified_four_state_morph" if known else "mixed_uncertain")
    d["technical_failure"]=d["measurement_status"]=="image_acquisition_failed"
    for col in ["source_original_fcp_morph","source_fcp_global_classifiable",
                "fcp_source_cohort","photo_colour_provenance"]:
        d.pop(col)
    return d


def examples():
    old=pd.DataFrame([
        previous(1,10,"white","oa"),
        previous(2,11,"blue_purple","ob"),
        previous(3,11,"mixed_uncertain","oc","reused_FCP_actual"),
    ])
    new=pd.DataFrame([
        new_source(2,11,"red_pink","ob"),  # photo ID already classified from FCP
        new_source(4,12,"white","od"),  # new classifiable source
        new_source(5,13,"mixed_uncertain","oe","image_acquisition_failed")
    ])
    return old,new


def test_photo_id_join_preserves_old_colors_and_new_real_acquisition_failures():
    old,new=examples()
    combined,report=m.combine(old,new)
    assert report["new_first93m_photo_IDs_already_in_FCP_MOIS"]==1
    assert report["new_photo_IDs_not_previously_measured_in_FCP_MOIS"]==2
    assert report["combined_unique_original_photo_attempts"]==5
    assert report["overlapping_original_photo_ID_both_classified_hue_conflicts_preserved_old"]==1
    assert combined.photo_id.tolist()==["1","2","3","4","5"]
    assert combined.loc[combined.photo_id.eq("2"),"morph"].iloc[0]=="blue_purple"
    assert combined.loc[combined.photo_id.eq("5"),"classified"].iloc[0]==False
    assert combined.loc[combined.photo_id.eq("5"),"technical_failure"].iloc[0]==True
    assert report["combined_unique_original_four_hue_classified"]==3
    assert report["added_nominal_taxa_with_genuine_four_hue_not_in_FCP_MOIS"]==1
    assert report["new_unique_observation_IDs_not_in_old_FCP_MOIS"]==2
    assert combined.iloc[:3].reset_index(drop=True).equals(old.reset_index(drop=True))


def test_no_old_pixel_colour_contradictions_due_to_new_unclassified_photo():
    old,new=examples()
    new.loc[0,"morph"]="mixed_uncertain"
    new.loc[0,"classified"]=False
    new.loc[0,"white"]=False
    result,receipt=m.combine(old,new)
    assert result.loc[result.photo_id.eq("2"),"morph"].iloc[0]=="blue_purple"
    assert receipt["overlapping_original_photo_ID_both_classified_hue_conflicts_preserved_old"]==0
    assert receipt["new_photo_IDs_not_previously_measured_in_FCP_MOIS"]==2


def test_duplicate_p0002_photo_id_cannot_expand_photo_count():
    old,new=examples()
    extra=new.iloc[[1]].copy()
    new=pd.concat([new,extra],ignore_index=True)
    with pytest.raises(ValueError,match="duplications"):
        m.combine(old,new)


def test_old_594081_rows_never_relabelled_in_whole_union():
    old,new=examples()
    original=old.copy(deep=True)
    combined,receipt=m.combine(old,new)
    assert old.equals(original)
    assert receipt["old_genuine_FCP_MOIS_colour_rows_untouched"]
    assert receipt["all_original_photo_classification_failures_remain_missing_not_nonwhite"]
    assert not receipt["no_native_species_or_environmental_modelling_performed"] is False
    assert combined.loc[combined.photo_id.eq("3"),"classified"].iloc[0]==False
