"""Synthetic source identity tests, NOT measured botanical conclusions."""
import pandas as pd
import pytest

from scripts.acquisition.mois_union_FCP447500_and_first93m18462_actual_photos import (
    append,
)


def row(pid,taxon,obs,morph,arm,measured="classified_four_state_morph"):
    classified=morph in {"white","yellow_orange","red_pink","blue_purple"}
    return {
        "photo_id":str(pid), "inat_taxon_id":str(taxon),
        "species":"Campanula punctata" if taxon==10 else "Cirsium brevicaule",
        "observation_id":obs,
        "latitude":38.4,"longitude":137.2,
        "observed_on":"2025-07-12",
        "arm":arm,"morph":morph,"classified":classified,
        "white":morph=="white","measurement_status":measured,
        "technical_failure":measured=="image_acquisition_failed",
    }


def sources():
    fcp=pd.DataFrame([
        row(1,10,"A","white","core"),
        row(2,10,"B","red_pink","core"),
        row(3,10,"C","white","FCP_final_447500_actual_original_ROI"),
        row(4,11,"D","mixed_uncertain","FCP_final_447500_actual_original_ROI",
            "image_acquisition_failed"),
    ])
    fcp["photo_colour_provenance"]="source_synthetic_FCP_only"
    fcp["source_original_fcp_morph"]=fcp.morph
    fcp["fcp_source_cohort"]="wave4"
    fcp["source_fcp_global_classifiable"]=fcp.classified
    first=pd.DataFrame([
        row(1,10,"A","white","core"),
        row(2,10,"B","red_pink","core"),
        # Photo 3 already exists in FCP. The new 93m model says blue.
        # The immutable FCP WHITE label wins and disagreement is reported.
        row(3,10,"C","blue_purple","global93m_probability_observation"),
        row(5,20,"NEW","yellow_orange","global93m_probability_observation"),
    ])
    return fcp,first


def call(fcp,first):
    return append(fcp,first,expected_fcp=4,expected_first=4,
                  expected_base=2,expected_first_new=2)


def test_preserve_existing_pixels_and_only_add_new_ids():
    fcp,first=sources()
    out,report=call(fcp,first)
    assert out.photo_id.tolist()==["1","2","3","4","5"]
    assert len(out)==5 and out.photo_id.is_unique
    assert report["first_full93m_p0002_actual_photo_attempts"]==2
    assert report["first_full93m_p0002_shared_with_FCP_completed_pixel_source"]==1
    assert report["additional_original_photo_IDs_not_previously_in_FCP_or_MOIS"]==1
    assert report["all_genuinely_attempted_unique_original_photo_IDs"]==5
    assert report["all_original_photo_four_hue_classified"]==4
    assert report["additional_first93m_four_hue_classified"]==1
    assert report["first93m_photo_ID_overlap_both_classified_hue_disagreements"]==1
    assert report["new_first93m_classified_nominal_taxa_not_in_FCP_or_old_MOIS"]==1
    assert out.loc[out.photo_id.eq("3"),"morph"].iloc[0]=="white"
    assert out.loc[out.photo_id.eq("4"),"classified"].iloc[0]==False
    assert out.loc[out.photo_id.eq("4"),"white"].iloc[0]==False
    assert out.loc[out.photo_id.eq("5"),"arm"].iloc[0]=="global93m_probability_observation"
    assert out.loc[out.photo_id.eq("5"),"photo_colour_provenance"].iloc[0]==(
        "actual_original_93m_observation_p0002_ROI_v4")
    assert report["original_licensed_photo_images_downloaded_for_this_union"] is False
    assert report["new_native_status_or_environmental_ecology_analyses_run"] is False


def test_original_mois_two_rows_are_immutable_even_when_both_sources_claim_same_id():
    fcp,first=sources()
    first.loc[0,"morph"]="yellow_orange"
    first.loc[0,"white"]=False
    with pytest.raises(ValueError,match="Original 350399 source image pixels"):
        call(fcp,first)


def test_refuse_fake_existing_93m_duplicate_photo_within_new_roi_source():
    fcp,first=sources()
    first.loc[3,"photo_id"]="3"
    with pytest.raises(ValueError,match="invalid or duplicate original photo ID"):
        call(fcp,first)


def test_technical_unknown_is_not_white_or_an_adaptive_phenotype():
    fcp,first=sources()
    first.loc[3,"morph"]="mixed_uncertain"
    first.loc[3,"classified"]=False
    first.loc[3,"white"]=False
    first.loc[3,"measurement_status"]="image_acquisition_failed"
    first.loc[3,"technical_failure"]=True
    out,report=call(fcp,first)
    assert report["all_original_photo_four_hue_classified"]==3
    assert report["unknown_image_ROIs_in_final_union"]==2
    assert report["new_first93m_classified_nominal_taxa_not_in_FCP_or_old_MOIS"]==0
    assert out.loc[out.photo_id.eq("5"),"white"].iloc[0]==False


def test_only_the_same_actual_original_photo_id_is_reconciled_not_taxon_summary():
    fcp,first=sources()
    # An identical name in an independent observation never counts as photo
    # duplication and must remain one actual source photo for each photo ID.
    first.loc[3,"inat_taxon_id"]="10"
    out,report=call(fcp,first)
    assert len(out)==5
    assert report["new_first93m_classified_nominal_taxa_not_in_FCP_or_old_MOIS"]==0
    assert report["first93m_photo_ID_overlap_accepted_taxon_disagreements"]==0


def test_same_overlap_photo_with_different_accepted_taxon_is_logged_not_relabeled():
    fcp,first=sources()
    first.loc[2,"inat_taxon_id"]="99"
    out,report=call(fcp,first)
    assert report["first93m_photo_ID_overlap_accepted_taxon_disagreements"]==1
    assert out.loc[out.photo_id.eq("3"),"inat_taxon_id"].iloc[0]=="10"
