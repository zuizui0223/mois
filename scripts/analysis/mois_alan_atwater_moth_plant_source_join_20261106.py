#!/usr/bin/env python3
"""MOIS × independently published Atwater 2025 Table 2 moth–flower relationships.

The seven cited plants have moth flower-visitation and/or moth-body pollen
carriage counts from a Florida sandhill study (Atwater et al., 2025,
DOI 10.1111/een.13446, table 2; data derived from Atwater 2013).
These counts are NOT confirmed stigma deposition or nocturnal dependence.
Join WCVP accepted names, source-original photos and distinct years only.
Never convert technical ROI unknown or direct botanical unknown to nonwhite.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_PINS={
 "canonical":"f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c",
 "taxon":"dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2",
 "names":"9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3",
 "direct":"7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577",
 "viirs2020":"3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3",
 "viirs2024":"1c5f03e5a7c355a0f7045133140c5cf2e39a5802c548f89d1c2f66c8a235db5f",
}
# Public paper Table 2 metadata, not a new GloBI plant interaction.
TABLE2=[
 ("Ageratina aromatica","Asteraceae",54,7),
 ("Eupatorium compositifolium","Asteraceae",29,12),
 ("Liatris tenuifolia","Asteraceae",13,2),
 ("Pityopsis graminifolia","Asteraceae",7,3),
 ("Balduina angustifolia","Asteraceae",9,3),
 ("Dalea pinnata","Fabaceae",20,0),
 ("Eriogonum tomentosum","Polygonaceae",19,1),
]

def require(ok,message):
    if not bool(ok):raise ValueError("MOIS_ATWATER_SOURCE: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def public_paper_ledger():
    d=pd.DataFrame(TABLE2,columns=["plant_name","plant_family","flower_visitor_count_in_paper_table2",
                                   "pollen_transport_count_in_paper_table2"])
    require(d.plant_name.is_unique and len(d)==7,"Study Table 2 species list drifted")
    require(d.flower_visitor_count_in_paper_table2.ge(0).all() and
            d.pollen_transport_count_in_paper_table2.ge(0).all(),"Negative source interactions")
    d["source_doi"]="10.1111/een.13446"
    d["source_study_region"]="Florida sandhill (USA)"
    d["source_evidence_tier"]="reported moth flower visits and body pollen transport (NOT effective pollination)"
    d["source_taxonomy_note"]=np.where(
       d.plant_name.eq("Liatris tenuifolia"),
       "source PDF OCR renders tenuifoia; independently verified accepted spelling tenuifolia",
       "")
    return d

def join_sources(original,names,accepted,direct,viirs2020,viirs2024):
    require(original.photo_id.is_unique and original.observation_id.is_unique,
            "Original photo/observation identity drift")
    mapping=names.loc[names.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
            ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    require(mapping.inat_taxon_id.is_unique,"Duplicated WCVP source taxon ID")
    botanical=accepted[["wcvp_accepted_plant_name_id","wcvp_accepted_name"]].drop_duplicates()
    require(botanical.groupby("wcvp_accepted_plant_name_id").wcvp_accepted_name.nunique().max()==1,
            "Conflicting WCVP accepted names")
    mapping=mapping.merge(botanical.drop_duplicates("wcvp_accepted_plant_name_id"),
           on="wcvp_accepted_plant_name_id",validate="many_to_one",how="inner")
    paper=public_paper_ledger()
    require(paper.plant_name.is_unique,"Duplicate study plant names")
    direct=direct.loc[direct.trait_name.eq("flower_primary_color")&
        direct.quality.isin(("high","medium"))&
        direct.resolution_status.eq("resolved")&
        ~direct.classification.isin(("cultivar_contamination","partial_self_compatibility")),
        ["accepted_species","normalized_value"]].copy()
    require(direct.accepted_species.is_unique,"Duplicate accepted species direct colours")
    colourmap={"white":"white","yellow_orange":"nonwhite","red_pink":"nonwhite",
               "blue_purple":"nonwhite","green_brown_inconspicuous":"nonwhite"}
    direct["known_single_colour"]=direct.normalized_value.map(colourmap)
    x=original.loc[original.source_subarm.eq("historical_main")&
        original.strict9&original.latitude.between(-60,60)&
        original.longitude.between(-180,180)&
        original.observed_on.astype(str).str[:4].between("2019","2024")].copy()
    x=x.merge(mapping[["inat_taxon_id","wcvp_accepted_name"]],on="inat_taxon_id",
              how="left",validate="many_to_one")
    x=x.merge(paper[["plant_name"]],left_on="wcvp_accepted_name",right_on="plant_name",
              how="inner",validate="many_to_one")
    x=x.merge(direct[["accepted_species","known_single_colour"]],
       left_on="wcvp_accepted_name",right_on="accepted_species",
       how="left",validate="many_to_one")
    require(x.photo_id.is_unique,"Study matching multiplied original photo IDs")
    x["year"]=x.observed_on.astype(str).str[:4]
    x["cell005"]=(np.floor((x.longitude+180)/.05).astype(int).astype(str)+"_"+
                  np.floor((x.latitude+90)/.05).astype(int).astype(str))
    require(viirs2020.cell_id.is_unique and viirs2024.cell_id.is_unique,
            "Duplicate source satellite cell")
    require(viirs2020.source_complete.sum()==1044 and len(viirs2024)==302,
            "Original 2020/2024 VIIRS reference cohort drifted")
    years={}
    for year,cellset in [("2020",set(viirs2020.loc[viirs2020.source_complete,"cell_id"])),
                         ("2024",set(viirs2024.cell_id))]:
        t=x.loc[x.year.eq(year)&x.cell005.isin(cellset)]
        known=t.loc[t.known_single_colour.isin(("white","nonwhite"))]
        years[year]={
          "photos_original_strict9":len(t),
          "distinct_accepted_study_species":int(t.wcvp_accepted_name.nunique()),
          "known_direct_colour_species":int(known.wcvp_accepted_name.nunique()),
          "direct_white_species":int(known.loc[known.known_single_colour.eq("white"),
                                              "wcvp_accepted_name"].nunique()),
          "direct_white_photos":int(known.known_single_colour.eq("white").sum()),
          "distinct_satellite_cells_with_atwater_study_species":int(t.cell005.nunique())
        }
    return x,years

def run(paths,out):
    out=Path(out);require(not out.exists(),"Refuse to overwrite original verified results")
    for key,path in paths.items():
        require(sha(path)==SOURCE_PINS[key],"Frozen original-source SHA differs: "+key)
    canonical=pd.read_csv(paths["canonical"],compression="gzip",low_memory=False,
       dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
       usecols=["photo_id","observation_id","inat_taxon_id","source_subarm","strict9",
                "strict9_white","latitude","longitude","observed_on"])
    require(len(canonical)==1032882, "Frozen canonical observation cohort count drift")
    taxa=pd.read_csv(paths["taxon"],compression="gzip",dtype=str,
       usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    accepted=pd.read_csv(paths["names"],compression="gzip",dtype=str,
       usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    direct=pd.read_csv(paths["direct"],compression="gzip",dtype=str,
       usecols=["accepted_species","trait_name","normalized_value","quality",
                "resolution_status","classification"])
    v20=pd.read_csv(paths["viirs2020"],dtype={"cell_id":str})
    v24=pd.read_csv(paths["viirs2024"],dtype={"cell_id":str})
    matched,years=join_sources(canonical,taxa,accepted,direct,v20,v24)
    known=matched.loc[matched.known_single_colour.isin(("white","nonwhite"))]
    paper=public_paper_ledger()
    out.mkdir(parents=True)
    paper.to_csv(out/"public_Atwater_2025_table2_moth_plant_counts.csv",index=False)
    rows=[{"cohort":year,**row} for year,row in years.items()]
    pd.DataFrame(rows).to_csv(out/"MOIS_Atwater_2020_2024_same_year_VIIRS_support.csv",index=False)
    result={
      "schema":"MOIS_20261106_Atwater_public_moth_pollen_evidence_species_support_v1",
      "source_study":"Atwater et al. 2025 Table 2, Florida sandhill; original Atwater 2013",
      "source_doi":"10.1111/een.13446",
      "source_license_note":"study table facts reproduced with scientific citation",
      "original_table_species":len(paper),
      "original_table_flower_visit_rows":int(paper.flower_visitor_count_in_paper_table2.sum()),
      "original_table_pollen_transport_rows":int(paper.pollen_transport_count_in_paper_table2.sum()),
      "source_WCVP_exact_spelling_correction_reviewed_for_Liatris":True,
      "source_WCVP_accepted_source_matching_NOT_global_endemicity":True,
      "2019_2024_strict9_photo_rows":len(matched),
      "2019_2024_distinct_study_species":int(matched.wcvp_accepted_name.nunique()),
      "study_species_found_in_canonical_original_image_photo_IDs":
          int(matched.wcvp_accepted_name.nunique()),
      "2019_2024_direct_white_species":int(
          known.loc[known.known_single_colour.eq("white"),"wcvp_accepted_name"].nunique()),
      "2019_2024_direct_nonwhite_species":int(
          known.loc[known.known_single_colour.eq("nonwhite"),"wcvp_accepted_name"].nunique()),
      "2019_2024_unresolved_direct_colour_species":int(
          matched.loc[~matched.wcvp_accepted_name.isin(set(known.wcvp_accepted_name)),
                      "wcvp_accepted_name"].nunique()),
      "same_year_VIIRS":years,
      "genuinely_measured_nocturnal_pollen_transport_at_the_MOIS_photo_sites":False,
      "original_study_shows_pollen_carrying_NOT_stigma_deposition":True,
      "species_based_global_moth_dependence_not_estimated":True,
      "no_unknown_photo_as_nonwhite_or_white":True,
      "new_images_downloaded_or_ROI_inferred":0,
      "MNAR_tipping_scenarios_run":False,
      "sha256_sources":SOURCE_PINS,
      "output_sha256":{p.name:sha(p) for p in sorted(out.iterdir())}
    }
    (out/"source_receipt.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("MOIS_ATWATER_2025_REAL_SOURCE_JOIN",json.dumps({
      k:result[k] for k in ("original_table_species","2019_2024_strict9_photo_rows",
      "2019_2024_distinct_study_species","2019_2024_direct_white_species",
      "2019_2024_unresolved_direct_colour_species","same_year_VIIRS")}),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in SOURCE_PINS:p.add_argument("--"+key,required=True,type=Path)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run({key:getattr(args,key) for key in SOURCE_PINS},args.out)

if __name__=="__main__":main()
