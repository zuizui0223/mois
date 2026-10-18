"""Actual 256 world-wide genus-first ROI pilot: new TAXA, not source photo counts."""
import pandas as pd
import pytest
from scripts.analysis import mois_broad_first256_roi_taxonomic_gain_audit as m


def fixture():
    old=pd.DataFrame({
        "photo_id":["1","2"],"inat_taxon_id":["8","9"],
        "species":["Cirsium sp1","Campanula sp2"]
    })
    p=[];r=[]
    for i in range(256):
        new="Novel"+str(i)
        photo=str(i+10000)
        mid="MOIS-"+str(i).zfill(24)
        p.append({"photo_id":photo,"measurement_id":mid,
                  "inat_taxon_id":str(100+i),"species":new+" alpha",
                  "genus":new,"grid10":f"{i%9:02d}:04",
                  "source_sampling_priority":"new_genus_new_species",
                  "shard":"0"})
        if i<69: hue=("white" if i%2 else "red_pink")
        else: hue="mixed_uncertain"
        status=("classified_four_state_morph" if i<69 else
                "acquisition_failed" if i in range(69,75) else
                "flower_not_detected")
        r.append({"photo_id":photo,"measurement_id":mid,
                  "morph":hue,"measurement_status":status,
                  "classified_four_hue":str(i<69).lower()})
    return old,pd.DataFrame(p),pd.DataFrame(r)


def test_count_only_actual_visible_hues_for_new_species_and_genera():
    old,p,roi=fixture()
    joined,k=m.pilot_metrics(old,p,roi)
    assert len(joined)==256
    assert k["new_photo_hue_classifications"]==69
    assert k["technical_unknown_or_acquisition_failed_attempts"]==187
    assert k["newly_classified_named_taxa_not_in_main_original_source"]==69
    assert k["newly_classified_genera_not_in_main_original_source"]==69
    assert k["source_original_photo_species"]==256
    assert k["source_original_photo_genera"]==256
    assert k["original_ROI_measurement_status_counts"]["acquisition_failed"]==6
    assert k["visible_hue_photo_counts"]["white"]==34
    assert k["visible_hue_photo_counts"]["red_pink"]==35
    assert k["new_taxa_per_1000_original_RGB_photo_attempts"]==69*1000/256


def test_photo_id_mismatch_rejected_without_inventing_extra_colours():
    old,p,roi=fixture()
    roi.loc[10,"photo_id"]="9999999"
    with pytest.raises(ValueError,match="photo IDs"):
        m.pilot_metrics(old,p,roi)


def test_unknown_photo_cannot_be_forced_into_nonwhite_classification():
    old,p,roi=fixture()
    roi.loc[70,"classified_four_hue"]="true"
    with pytest.raises(ValueError,match="unknown original ROI"):
        m.pilot_metrics(old,p,roi)


def test_no_duplicate_pilot_source_ids_or_previous_measured_photos():
    old,p,roi=fixture()
    p.loc[100,"photo_id"]="1"
    with pytest.raises(ValueError,match="do not match"):
        m.pilot_metrics(old,p,roi)
