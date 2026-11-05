#!/usr/bin/env python3
"""Source-backed POSITIVE moth/butterfly flower visit evidence across MOIS flora.

An insect visit is not nocturnal effective pollination. Missing visitor records
are UNKNOWN, never biological absence. Plants with both visitors contribute to
both POSITIVE evidence columns but cannot be treated as independent contrast
groups. Only aggregate counts are published, no licensed species colour rows.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE={
 "canonical":"f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c",
 "names":"dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2",
 "accepted":"9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3",
 "direct":"7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577",
 "globi":"a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994",
 "viirs":"3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3"}
SCALES=(.05,.1,.25,.5,1.,2.)
KNOWN_COLOUR=("white","yellow_orange","red_pink","blue_purple","green_brown_inconspicuous")

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_GLOBI_POSITIVE_SUPPORT: "+msg)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def evidence_flags(source):
    require({"plant_name","guild"} <= set(source),"GloBI interaction source schema drift")
    require(source.guild.isin(("moth_family_visit","butterfly_family_visit")).all(),
            "Not a supported insect family visit channel")
    frame=source[["plant_name","guild"]].drop_duplicates().assign(value=1).pivot(
                index="plant_name",columns="guild",values="value").fillna(0)
    for guild in ("moth_family_visit","butterfly_family_visit"):
        if guild not in frame:frame[guild]=0
    frame=frame.rename(columns={"moth_family_visit":"moth_record",
                "butterfly_family_visit":"butterfly_record"})
    require(frame[["moth_record","butterfly_record"]].ge(0).all().all(),"Invalid positive evidence")
    return frame

def spatial_counts(observed,step):
    """Each accepted species is counted ONCE per cell, even if many photos."""
    require({"longitude","latitude","accepted_species","color_state",
             "moth_record","butterfly_record"}<=set(observed),"Source records absent")
    x=observed.copy()
    x["cell"]=(np.floor((x.longitude+180)/step).astype(int).astype(str)+"_"+
               np.floor((x.latitude+90)/step).astype(int).astype(str))
    x=x.groupby(["cell","accepted_species","color_state"],observed=True,as_index=False).agg(
            moth_record=("moth_record","max"),butterfly_record=("butterfly_record","max"))
    require(not x.duplicated(["cell","accepted_species"]).any(),
            "Original accepted species repeated in same cell")
    known=x.loc[x.color_state.isin(("white","nonwhite"))].copy()
    if not len(known):
        return {"cells_with_known_flower_color":0,"white_moth_and_butterfly_visitor_positive_cells":0,
           "four_evidence_color_classes_min1_cells":0,
           "four_evidence_color_classes_min2_cells":0,
           "four_evidence_color_classes_min3_cells":0,
           "four_classes_no_shared_white_species_cells":0,
           "white_species_shared_moth_butterfly_evidence_cell_rows":0}
    table=known.groupby(["cell","color_state"],observed=True)[
           ["moth_record","butterfly_record"]].sum().unstack("color_state",fill_value=0)
    for role in ("moth_record","butterfly_record"):
        for colour in ("white","nonwhite"):
            if (role,colour) not in table:table[(role,colour)]=0
    c=pd.DataFrame({
       "mw":table[("moth_record","white")],
       "bw":table[("butterfly_record","white")],
       "mn":table[("moth_record","nonwhite")],
       "bn":table[("butterfly_record","nonwhite")]})
    shared=known.loc[known.color_state.eq("white") &
          known.moth_record.eq(1)&known.butterfly_record.eq(1)].groupby("cell").size()
    c["shared_white"]=shared.reindex(c.index,fill_value=0).to_numpy()
    both_white=c.mw.ge(1)&c.bw.ge(1)
    four1=c[["mw","bw","mn","bn"]].ge(1).all(axis=1)
    four2=c[["mw","bw","mn","bn"]].ge(2).all(axis=1)
    four3=c[["mw","bw","mn","bn"]].ge(3).all(axis=1)
    return {"cells_with_known_flower_color":len(c),
            "white_moth_and_butterfly_visitor_positive_cells":int(both_white.sum()),
            "four_evidence_color_classes_min1_cells":int(four1.sum()),
            "four_evidence_color_classes_min2_cells":int(four2.sum()),
            "four_evidence_color_classes_min3_cells":int(four3.sum()),
            "four_classes_no_shared_white_species_cells":int((four1&c.shared_white.eq(0)).sum()),
            "white_species_shared_moth_butterfly_evidence_cell_rows":int(c.shared_white.sum())}

def source_join(canonical,names,accepted,direct,globi):
    require(canonical.photo_id.is_unique and canonical.observation_id.is_unique,
            "Duplicate original photo or observation ID")
    taxa=names.loc[names.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
         ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    require(taxa.inat_taxon_id.is_unique,"Duplicated original taxon key")
    ac=accepted[["wcvp_accepted_plant_name_id","wcvp_accepted_name"]].drop_duplicates()
    require(ac.groupby("wcvp_accepted_plant_name_id").wcvp_accepted_name.nunique().le(1).all(),
            "Conflicting original accepted botanical names")
    taxa=taxa.merge(ac.drop_duplicates("wcvp_accepted_plant_name_id"),
                    on="wcvp_accepted_plant_name_id",how="left",validate="many_to_one")
    colour=direct.loc[direct.trait_name.eq("flower_primary_color")&
        direct.quality.isin(("high","medium"))&
        direct.resolution_status.eq("resolved")&
        ~direct.classification.isin(("cultivar_contamination","partial_self_compatibility")),
        ["accepted_species","normalized_value"]].copy()
    require(colour.accepted_species.is_unique,"Duplicate direct species evidence")
    colour["color_state"]=np.where(colour.normalized_value.eq("white"),"white",
         np.where(colour.normalized_value.isin(KNOWN_COLOUR),"nonwhite","unknown"))
    x=canonical.loc[canonical.source_subarm.eq("historical_main")&
         canonical.strict9 &canonical.latitude.between(-60,60)&
         canonical.longitude.between(-180,180)&
         canonical.observed_on.astype(str).str[:4].between("2019","2024")].copy()
    x=x.merge(taxa[["inat_taxon_id","wcvp_accepted_name"]],on="inat_taxon_id",
         how="left",validate="many_to_one")
    x=x.merge(evidence_flags(globi),left_on="wcvp_accepted_name",right_index=True,
         how="inner",validate="many_to_one")
    x=x.merge(colour[["accepted_species","color_state"]],
         left_on="wcvp_accepted_name",right_on="accepted_species",
         how="left",validate="many_to_one")
    x["accepted_species"]=x.wcvp_accepted_name
    require(x.photo_id.is_unique,"Source join multiplied original photo records")
    return x

def run(input_paths,out):
    out=Path(out);require(not out.exists(),"Refuse to overwrite source results")
    for key,path in input_paths.items():
        require(digest(path)==SOURCE[key],"SHA-locked source differs: "+key)
    c=pd.read_csv(input_paths["canonical"],compression="gzip",low_memory=False,
          dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
          usecols=["photo_id","observation_id","inat_taxon_id","source_subarm",
                   "strict9","latitude","longitude","observed_on"])
    require(len(c)==1032882,"Original canonical population changed")
    names=pd.read_csv(input_paths["names"],compression="gzip",dtype=str,
       usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    accepted=pd.read_csv(input_paths["accepted"],compression="gzip",dtype=str,
       usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    traits=pd.read_csv(input_paths["direct"],compression="gzip",dtype=str,
       usecols=["accepted_species","trait_name","normalized_value",
                "quality","resolution_status","classification"])
    globi=pd.read_csv(input_paths["globi"],compression="gzip",
       usecols=["plant_name","guild"])
    x=source_join(c,names,accepted,traits,globi)
    require(len(x)>=110000 and x.accepted_species.nunique()>3500,
            "Independent source-positive sample suddenly shrank")
    viirs=pd.read_csv(input_paths["viirs"],dtype={"cell_id":str})
    require(len(viirs)==1059 and viirs.source_complete.sum()==1044,
            "2020 existing VIIRS domain changed")
    base=set(viirs.loc[viirs.source_complete,"cell_id"])
    x["site005"]=(np.floor((x.longitude+180)/.05).astype(int).astype(str)+"_"+
                  np.floor((x.latitude+90)/.05).astype(int).astype(str))
    rows=[]
    for degree in SCALES:
        rows.append({"frame":"six_year_all_existing_photo_regions","degrees":degree,
                     **spatial_counts(x,degree)})
    rows.append({"frame":"six_year_photo_at_2020_viirs_1044_sites","degrees":.05,
                 **spatial_counts(x.loc[x.site005.isin(base)],.05)})
    out.mkdir(parents=True)
    result=pd.DataFrame(rows)
    result.to_csv(out/"GloBI_positive_visit_colour_spatial_support.csv",index=False)
    report={
        "schema":"mois_real_GloBI_positive_visit_independent_colour_spatial_support_v1",
        "input_sha256":SOURCE,"photographic_source_2019_2024":"historical_main",
        "original_photo_ids_with_any_positive_guild":len(x),
        "distinct_accepted_species_with_any_guild":int(x.accepted_species.nunique()),
        "distinct_accepted_species_direct_single_colour":
          int(x.loc[x.color_state.isin(("white","nonwhite"))].accepted_species.nunique()),
        "source_GloBI_unique_interaction_rows":61967,
        "GloBI_one_source_family_still_capped":True,
        "two_positive_guild_flags_are_NOT_exclusive":True,
        "same_species_can_contribute_to_both_evidence_columns":True,
        "no_record_is_NOT_negative_pollinator_evidence":True,
        "actual_nocturnality_or_effective_pollination_not_estimated":True,
        "species_extinction_and_causal_ALAN_not_estimated":True,
        "source_2020_light_applied_to_2019_2024_photographs_is_temporally_mismatched":True,
        "no_original_photo_downloads_or_original_RGB_ROI_inference":True,
        "unknown_photographic_flower_color_not_imputed":True,
        "MNAR_tipping_scenarios_not_run":True,
        "result_csv_sha256":digest(out/"GloBI_positive_visit_colour_spatial_support.csv")}
    (out/"source_positive_evidence_receipt.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_ALAN_POSITIVE_EVIDENCE_COHORT_COUNTS",json.dumps(
       {k:report[k] for k in ("original_photo_ids_with_any_positive_guild",
         "distinct_accepted_species_with_any_guild","distinct_accepted_species_direct_single_colour")}),flush=True)
    print("MOIS_ALAN_SOURCE_POSITIVE_SUPPORT",result.to_json(orient="records"),flush=True)
    return result,report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in SOURCE:p.add_argument("--"+k,required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    run({k:getattr(args,k) for k in SOURCE},args.output)
