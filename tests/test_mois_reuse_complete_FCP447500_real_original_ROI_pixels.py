"""Synthetic contracts: only real FCP original measured pixels, no repeated MOIS photos."""
import json
import pandas as pd
import pytest
from scripts.analysis import mois_reuse_complete_FCP447500_real_original_ROI_pixels as mod
from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import load_pixels, sha256


def photo(pid,taxon,morph,valid=True,cohort="discovery",obs=None):
    return {
        "photo_id":str(pid),"inat_taxon_id":str(taxon),
        "species":"Campanula punctata" if taxon==1 else "Cirsium brevifolium",
        "morph":morph,"global_classifiable":valid,
        "fcp_source_cohort":cohort,"observation_id":obs or f"obs{pid}",
        "latitude":32.5,"longitude":140.5,"observed_on":"2025-07-03",
        "measurement_status":"classified_four_state_morph" if valid
                             else "image_acquisition_failed",
    }


def original(path):
    old=[
        photo(1,1,"white"),
        photo(2,2,"red_pink"),
        photo(3,1,"blue_purple")
    ]
    for a in old:a.pop("global_classifiable");a.pop("fcp_source_cohort")
    pd.DataFrame(old).to_parquet(path,index=False)
    return load_pixels(path,"core")


def toy_fcp():
    q=pd.DataFrame([
        photo(1,1,"blue_purple"), # same ID; prior MOIS photo takes precedence
        photo(2,2,"red_pink"),
        photo(4,1,"white"),
        photo(5,2,"yellow_orange"),
        photo(6,2,"white",False), # unclassifiable; still attempted photo pixel
        photo(7,1,"mixed_uncertain",False,obs="obs4"),
    ])
    return q


def test_source_photo_overlap_is_not_repeated_and_original_morph_retained(tmp_path):
    old=original(tmp_path/"old.parquet")
    new=mod.canonical_fcp(toy_fcp())
    all_pixels,reused,report=mod.join(old,new,expected_old=3)
    assert all_pixels.photo_id.tolist()==["1","2","3","4","5","6","7"]
    assert len(reused)==4
    assert report["original_FCP_photo_IDs_already_in_MOIS"]==2
    assert report["new_genuinely_measured_FCP_photo_IDs"]==4
    assert report["new_FCP_four_hue_classified"]==2
    assert report["new_FCP_unclassifiable_original_ROIs"]==2
    assert report["union_total_genuinely_original_photo_IDs"]==7
    assert report["union_total_four_hue_classified"]==5
    assert report["old_source_image_ROI_results_preserved"] is True
    assert report["already_measured_photo_ID_overlap_both_classified_hue_disagreements"]==1
    assert all_pixels.loc[all_pixels.photo_id.eq("1"),"morph"].iloc[0]=="white"
    assert all_pixels.loc[all_pixels.photo_id.eq("6"),"morph"].iloc[0]=="mixed_uncertain"
    assert all_pixels.loc[all_pixels.photo_id.eq("6"),"source_original_fcp_morph"].iloc[0]=="white"
    assert all_pixels.loc[all_pixels.photo_id.eq("6"),"classified"].iloc[0]==False
    assert all_pixels.loc[all_pixels.photo_id.eq("7"),"observation_id"].iloc[0]=="obs4"


def test_original_fcp_full_100_photo_species_and_sha_gate(tmp_path):
    data=[]
    for i in range(5):
        for j in range(100):
            item=photo(1000+100*i+j,i+1,"white" if j%2 else "blue_purple",
                       True,["discovery","validation","prospective","wave4","new_U100"][i])
            item["species"]=f"Campanula alba{i}"
            data.append(item)
    f=tmp_path/"fcp.csv.gz"
    pd.DataFrame(data).to_csv(f,index=False,compression="gzip")
    join_receipt=tmp_path/"source.json"
    join_receipt.write_text(json.dumps({
        "schema":mod.FCP_COMBINE_SCHEMA,
        "total_original_pixel_terminal_photo_rows":500,
        "total_original_100photo_measured_species":5,
        "original_149900_photo_rows_unchanged":True,
        "new_U100_source_reuse_excluded_taxa":3,
        "merged_photo_table_sha256":sha256(f)
    }))
    release=tmp_path/"final.json"
    release.write_text(json.dumps({
        "schema":mod.FCP_SCHEMA,
        "genuinely_measured_photo_rows":500,
        "genuinely_measured_100photo_species":5,
        "genuinely_four_state_classified_photos":500,
        "old_149900_original_photos_unchanged":True,
        "no_new_ecological_statistic_or_classifier":True,
        "source_archive_sha256":sha256(f),
        "source_join_receipt_sha256":sha256(join_receipt)
    }))
    q,r,c=mod.read_fcp_genuine(f,release,join_receipt,
                               expected_rows=500,expected_taxa=5)
    assert len(q)==500
    assert q.inat_taxon_id.nunique()==5
    assert q.fcp_classifiable.all()
    d=pd.read_csv(f);d.loc[0,"morph"]="unknown"
    d.to_csv(f,index=False,compression="gzip")
    with pytest.raises(ValueError,match="altered"):
        mod.read_fcp_genuine(f,release,join_receipt,
                             expected_rows=500,expected_taxa=5)


def test_duplicate_or_missing_original_FCP_photo_IDs_rejected(tmp_path):
    original_old=original(tmp_path/"old.parquet")
    q=toy_fcp()
    q.loc[2,"photo_id"]="1"
    with pytest.raises(ValueError,match="duplicated"):
        mod.join(original_old,mod.canonical_fcp(q),expected_old=3)
