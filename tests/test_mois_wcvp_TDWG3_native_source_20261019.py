"""Source-only tests for Kew WCVP TDWG level-3 native/introduced gating.

WCVP status is not photo point nativeness or island endemism. Tests use
fabricated tiny nomenclatural/distribution rows, not purported source data.
"""
import csv
import zipfile
import pandas as pd

from scripts.acquisition.mois_wcvp_TDWG3_native_source_20261019 import (
    clean_name, collect_wcvp_names, resolve_photo_taxa, source_taxon_table,
    accepted_name_metadata, native_distribution_TDWG3
)


def tiny_archive(tmp_path):
    path=tmp_path/"wcvp.zip"
    fields_names=("plant_name_id","taxon_name","taxon_rank",
                  "taxon_status","accepted_plant_name_id")
    n=[
        {"plant_name_id":"101","taxon_name":"Poa annua","taxon_rank":"Species",
         "taxon_status":"Accepted","accepted_plant_name_id":"101"},
        {"plant_name_id":"102","taxon_name":"Poa minor","taxon_rank":"Species",
         "taxon_status":"Synonym","accepted_plant_name_id":"101"},
        {"plant_name_id":"103","taxon_name":"Poa dubia","taxon_rank":"Species",
         "taxon_status":"Synonym","accepted_plant_name_id":"101"},
        {"plant_name_id":"105","taxon_name":"Poa dubia","taxon_rank":"Species",
         "taxon_status":"Synonym","accepted_plant_name_id":"104"},
        {"plant_name_id":"104","taxon_name":"Poa major","taxon_rank":"Species",
         "taxon_status":"Accepted","accepted_plant_name_id":"104"},
        {"plant_name_id":"106","taxon_name":"Poa annua var. hispida",
         "taxon_rank":"Variety","taxon_status":"Accepted","accepted_plant_name_id":"106"}
    ]
    fields_dist=("plant_name_id","area_code_l3","area","introduced",
                 "extinct","location_doubtful")
    d=[
        {"plant_name_id":"101","area_code_l3":"JAP","area":"Japan",
         "introduced":"0","extinct":"0","location_doubtful":"0"},
        {"plant_name_id":"101","area_code_l3":"NZS","area":"New Zealand South",
         "introduced":"1","extinct":"0","location_doubtful":"0"},
        {"plant_name_id":"101","area_code_l3":"BAD","area":"Mixed",
         "introduced":"1","extinct":"0","location_doubtful":"0"},
        {"plant_name_id":"101","area_code_l3":"BAD","area":"Mixed",
         "introduced":"0","extinct":"0","location_doubtful":"0"},
        {"plant_name_id":"101","area_code_l3":"UNC","area":"uncertain",
         "introduced":"0","extinct":"0","location_doubtful":"1"},
        {"plant_name_id":"101","area_code_l3":"EXT","area":"Extinct",
         "introduced":"0","extinct":"1","location_doubtful":"0"},
        {"plant_name_id":"104","area_code_l3":"JAP","area":"Japan",
         "introduced":"0","extinct":"0","location_doubtful":"0"},
    ]
    from io import StringIO
    def build(rows,fields):
        s=StringIO()
        w=csv.DictWriter(s,fieldnames=fields,lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
        return s.getvalue()
    with zipfile.ZipFile(path,"w",zipfile.ZIP_DEFLATED) as z:
        z.writestr("WCVP_2026/wcvp_names.csv",build(n,fields_names))
        z.writestr("WCVP_2026/wcvp_distribution.csv",build(d,fields_dist))
    return path


def test_source_taxon_table_is_not_photo_majority_colour(tmp_path):
    photo=pd.DataFrame([
        {"photo_id":"1","inat_taxon_id":"A","species":"Poa annua","classified":True,"arm":"core"},
        {"photo_id":"2","inat_taxon_id":"A","species":"Poa annua","classified":False,"arm":"core"},
        {"photo_id":"3","inat_taxon_id":"B","species":"Poa minor","classified":True,"arm":"broad"},
        {"photo_id":"4","inat_taxon_id":"","species":"Poa dubia","classified":True,"arm":"tail"},
        {"photo_id":"5","inat_taxon_id":"D","species":"Poa dubia","classified":False,"arm":"tail"},
    ])
    t=source_taxon_table(photo,expected_photos=5,expected_hues=3)
    assert int(t.original_photos.sum())==5
    assert int(t.classifiable_photos.sum())==3
    assert len(t)==4
    assert t.loc[t.inat_taxon_id.eq("A"),"original_photos"].iloc[0]==2
    assert "Poa annua"==t.loc[t.inat_taxon_id.eq("A"),"match_binomial"].iloc[0]


def test_wcvp_exact_accepted_synonyms_and_ambiguous_protected(tmp_path):
    path=tiny_archive(tmp_path)
    taxa=pd.DataFrame([
        dict(inat_taxon_id="A",source_species="Poa annua",match_binomial="Poa annua",
             source_taxon_inconsistent=False,original_photos=2,classifiable_photos=1),
        dict(inat_taxon_id="B",source_species="Poa minor",match_binomial="Poa minor",
             source_taxon_inconsistent=False,original_photos=1,classifiable_photos=1),
        dict(inat_taxon_id="C",source_species="Poa dubia",match_binomial="Poa dubia",
             source_taxon_inconsistent=False,original_photos=1,classifiable_photos=0),
        dict(inat_taxon_id="",source_species="Poa major",match_binomial="Poa major",
             source_taxon_inconsistent=False,original_photos=1,classifiable_photos=0)
    ])
    with zipfile.ZipFile(path) as z:
        candidates,states,rows,member=collect_wcvp_names(
            z,{"Poa annua","Poa minor","Poa dubia","Poa major"},min_source_rows=2)
        assert rows==6
        assert candidates["Poa annua"]=={"101"}
        assert candidates["Poa minor"]=={"101"}
        assert candidates["Poa dubia"]=={"101","104"}
        matched,selected=resolve_photo_taxa(taxa,candidates,states,expected_photos=5)
        assert selected=={"101"}
        assert list(matched.wcvp_match_state)==[
            "unique_exact_name_to_accepted_id",
            "unique_exact_name_to_accepted_id",
            "multiple_accepted_ids_ambiguous",
            "missing_original_taxon_id"
        ]
        accepted=accepted_name_metadata(z,selected)
        assert accepted=={"101":"Poa annua"}


def test_wcvp_TDWG3_nativeness_flags_without_photo_point_inference(tmp_path):
    path=tiny_archive(tmp_path)
    with zipfile.ZipFile(path) as z:
        accepted={"101":"Poa annua"}
        d,count,matching,member=native_distribution_TDWG3(
            z,set(accepted),accepted,min_source_rows=2)
    assert count==7
    assert matching==6
    statuses=dict(zip(d.tdwg_level3_code,d.native_introduced_state))
    assert statuses["JAP"]=="native_at_TDWG3"
    assert statuses["NZS"]=="introduced_at_TDWG3"
    assert statuses["BAD"]=="contradictory_TDWG3_records"
    assert statuses["UNC"]=="doubtful"
    assert statuses["EXT"]=="extinct"
    assert all(d.native_on_specific_island.eq("not_established_by_TDWG3"))
    assert all(d.native_at_original_photo_lat_lon.eq("not_evaluated_no_TDWG3_polygon"))


def test_taxon_name_uncertainty_unresolved_not_auto_mapped():
    assert clean_name("Poa annua")=="Poa annua"
    assert clean_name("Poa  annua")=="Poa annua"
    assert clean_name("Poa annua var. annua")==""
    assert clean_name("Poa cf. annua")==""
    assert clean_name("Poa sp.")==""
    assert clean_name("Poa")==""
    assert clean_name("poa annua")==""
