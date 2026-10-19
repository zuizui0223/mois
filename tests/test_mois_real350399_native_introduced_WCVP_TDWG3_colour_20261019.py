"""Synthetic original-image photographic source/TDWG native-vs-introduced tests.

No simulated value represents observed wild native flower frequency.
"""
import numpy as np
import pandas as pd

from scripts.analysis.mois_real350399_native_introduced_WCVP_TDWG3_colour_20261019 import (
    distinct_original_observations, matched_native_introduced, photo_status_table,
    species_region_table,
)


def row(pid, obs, accepted, genus, region, status, hue, cohort="core"):
    return dict(
        photo_id=str(pid), observation_id=str(obs),
        wcvp_accepted_plant_name_id=str(accepted),
        accepted_genus=genus,
        tdwg_level3_code=region,
        native_introduced_state=status,
        source_cohort=cohort,
        arm=("core" if cohort=="core" else "tail16384"),
        taxon_TDWG3_range_annotation="WCVP_"+status,
        inat_taxon_id=accepted,
        species=genus+" species",
        morph=hue,classified=hue in
            ("white","yellow_orange","red_pink","blue_purple"),
        white=hue=="white"
    )


def source_rows():
    result=[]
    serial=1
    for k in range(12):
        genus=f"G{k}us"
        for nstatus,bit in (("native_at_TDWG3",True),
                            ("introduced_at_TDWG3",False)):
            for tax in range(2):
                pid=serial
                result.append(row(pid,f"obs{pid}",f"k{k}-{nstatus}-{tax}",
                                  genus,"AAA",nstatus,
                                  "white" if bit else "red_pink"))
                serial+=1
    return pd.DataFrame(result)


def test_genus_TDWG3_matched_photo_white_on_accepted_species_equal():
    x=source_rows()
    obs=distinct_original_observations(x)
    assert len(obs)==48
    taxa=species_region_table(obs,"all_source_arms")
    assert len(taxa)==48
    result=matched_native_introduced(taxa,min_species=2,draws=99)
    assert result["matched_genus_region_strata"]==12
    assert result["matched_named_genera"]==12
    assert result["native_minus_introduced_white_pp_equal_genus"]==100
    assert result["native_minus_introduced_white_pp_equal_matched_strata"]==100
    assert result["whole_named_genus_bootstrap95_pp"]==[100,100]
    assert result["matched_accepted_species_by_status"]=={
        "native_at_TDWG3":24,"introduced_at_TDWG3":24}
    assert not result["wild_native_locality_or_island_identity_known"]


def test_original_one_observation_multiple_photos_not_two_photographed_occurrences():
    x=source_rows()
    extra=x.iloc[[0]].copy()
    extra.loc[:,"photo_id"]="other-photo-same-observation"
    extra.loc[:,"morph"]="blue_purple"
    extra.loc[:,"white"]=False
    d=pd.concat((x,extra),ignore_index=True)
    obs=distinct_original_observations(d)
    assert len(obs)==48
    obs1=obs.loc[obs.observation_id.eq("obs1")].iloc[0]
    assert obs1.photos==2
    assert obs1.white_label_fraction==.5
    t=species_region_table(obs,"all_source_arms")
    assert len(t)==48
    one=t.loc[t.wcvp_accepted_plant_name_id.eq("k0-native_at_TDWG3-0")].iloc[0]
    assert one.p_white==.5
    assert one.observation_count==1


def test_two_status_source_photo_classifiability_not_false_nonwhite():
    x=source_rows()
    unknown=x.iloc[[0]].copy()
    unknown.loc[:,"photo_id"]="new-original-image"
    unknown.loc[:,"morph"]="mixed_uncertain"
    unknown.loc[:,"classified"]=False
    unknown.loc[:,"white"]=False
    missing=x.iloc[[0]].copy()
    missing.loc[:,"photo_id"]="unresolved-geographic-point"
    missing.loc[:,"taxon_TDWG3_range_annotation"]="TDWG3_photo_location_unassigned"
    missing.loc[:,"classified"]=False
    missing.loc[:,"morph"]=""
    missing.loc[:,"white"]=False
    d=pd.concat([x,unknown,missing],ignore_index=True)
    qc=photo_status_table(d)
    assert qc.original_source_photo_attempts.sum()==50
    assert qc.classifiable_photo_ROIs.sum()==48
    native=qc[qc.region_species_status.eq("WCVP_native_at_TDWG3")].iloc[0]
    assert native.unknown_or_unclassifiable_photo_ROIs==1
    unassigned=qc[qc.region_species_status.eq("TDWG3_photo_location_unassigned")].iloc[0]
    assert unassigned.classifiable_photo_ROIs==0
    assert np.isnan(unassigned.white_fraction_among_classified)


def test_different_wcvp_species_synonyms_share_one_accepted_species_weight():
    x=source_rows()
    extra=x.iloc[[0]].copy()
    extra.loc[:,"photo_id"]="extra-synonym-photo"
    extra.loc[:,"observation_id"]="extra-synonym-observation"
    extra.loc[:,"inat_taxon_id"]="original-synonym-variant"
    extra.loc[:,"morph"]="red_pink"
    extra.loc[:,"white"]=False
    d=pd.concat((x,extra),ignore_index=True)
    obs=distinct_original_observations(d)
    t=species_region_table(obs,"all_source_arms")
    assert len(t)==48
    a=t.loc[t.wcvp_accepted_plant_name_id.eq("k0-native_at_TDWG3-0")].iloc[0]
    assert a.observation_count==2
    assert a.p_white==.5


def test_no_tdwg3_overlap_does_not_count_as_no_colour_effect():
    x=source_rows()
    x.loc[x.native_introduced_state.eq("introduced_at_TDWG3"),
          "tdwg_level3_code"]="BBB"
    s=species_region_table(distinct_original_observations(x),"all_source_arms")
    result=matched_native_introduced(s,min_species=1,draws=99)
    assert result["matched_genus_region_strata"]==0
    assert result["native_minus_introduced_white_pp_equal_genus"] is None
    assert result["whole_named_genus_bootstrap95_pp"] is None


def test_source_arm_conditioned_comparison_is_not_automatically_pooled():
    x=source_rows()
    x.loc[x.native_introduced_state.eq("introduced_at_TDWG3"),"source_cohort"]="new_species_first"
    obs=distinct_original_observations(x)
    old=species_region_table(obs,"core")
    new=species_region_table(obs,"new_species_first")
    pooled=species_region_table(obs,"all_source_arms")
    assert matched_native_introduced(old,min_species=1,draws=9)["matched_genus_region_strata"]==0
    assert matched_native_introduced(new,min_species=1,draws=9)["matched_genus_region_strata"]==0
    assert matched_native_introduced(pooled,min_species=1,draws=9)["matched_genus_region_strata"]==12
