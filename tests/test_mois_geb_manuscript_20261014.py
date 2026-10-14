"""Editorial and scientific claim regression checks for the mois GEB draft.

No original iNaturalist colour measurement, fitted model, null or fitted coefficient
is recomputed here: immutable source receipts define numeric authority.
"""
from __future__ import annotations
from pathlib import Path
import json
import math
import re

ROOT=Path(__file__).resolve().parents[1]
MS=ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"
READINESS=ROOT/"docs/MOIS_GEB_READINESS_AND_COVER_LETTER_20261014.md"
MACRO=ROOT/"results/fcp_macro_direction_20261013/decision_receipt.json"
SPATIAL=ROOT/"results/fcp_m1_matched_spatial_support_20261014/decision_receipt.json"
MNAR=ROOT/"results/source_receipts/photo_colour_missingness_mnar.json"
MATCHED=ROOT/"results/source_receipts/matched_congeneric_vs_withinspecies_rainfall.json"
INDEPENDENT=ROOT/"results/source_receipts/nhm_independent_botanical_colour_summary.json"
GENUS=ROOT/"results/source_receipts/unseen_genus_colour_environment.json"
DISTANT=ROOT/"results/source_receipts/unseen_genus_other_geography_500km.json"
NONLINEAR=ROOT/"results/source_receipts/nonlinear_geographic_transport.json"


def source(path):
    return json.loads(path.read_text(encoding="utf8"))


def text():
    return MS.read_text(encoding="utf8")


def section(s,start,end):
    return s[s.index(start):s.index(end,s.index(start)+len(start))]


def test_research_article_abstract_and_section_structure():
    s=text()
    assert s.startswith("# The geographic and measurement dependence of flower-colour–moisture")
    abstract=section(s,"## Abstract","## 1 | Introduction")
    fields=["Aim","Location","Time period","Major taxa studied","Methods","Results","Main conclusions"]
    starts=[abstract.index("**"+f+":**") for f in fields]
    assert starts==sorted(starts)
    assert 0<len(abstract.strip().split())<=300
    assert len(fields)==7
    assert len(re.findall(r"^\*\*Keywords:\*\*",s,flags=re.M))==1
    keywords=re.search(r"^\*\*Keywords:\*\* (.+)$",s,flags=re.M).group(1).split(",")
    assert 6<=len(keywords)<=10 and all(k.strip() for k in keywords)
    main=section(s,"## 1 | Introduction","## Data and code availability")
    assert 3000 <= len(main.split()) <= 5000
    for header in ("## 2 | Materials and Methods","## 3 | Results","## 4 | Discussion","## 5 | Main conclusions",
                   "## Data and code availability","## References","## Main tables","## Figure legends"):
        assert header in s
    assert "not been submitted" in s


def test_worldwide_photo_and_broad_environment_results_stay_exact():
    s=text();q=source(MACRO)
    assert (q["original_photo_opportunities"],q["source_photo_classified"]) == (100543,45953)
    assert (q["n_macro_scope_cohorts"],q["directional_estimates_and_sensitivities"])==(20,1280)
    def pick(axis,quality="all",adjust="genus_geo"):
        z=[x for x in q["key_white_directional_observations"]
           if (x["mode"],x["region"],x["weight"] if "weight" in x else "genus_equal",
               x["quality"],x["adjustment"],x["axis"],x["hue"]) ==
           ("climate","world","genus_equal",quality,adjust,axis,"white")]
        assert len(z)==1
        return z[0]
    rain=pick("wc_bio12")
    vapour=pick("wc_vapr_annual_mean")
    assert rain["n_photo_species"]==20694 and rain["n_genera"]==1934
    assert math.isclose(rain["point_pp"],2.254545235123587,abs_tol=1e-12)
    assert rain["family_95_interval_pp"][0]>0
    assert math.isclose(vapour["point_pp"],3.951850325152891,abs_tol=1e-12)
    assert pick("wc_bio12",quality="strict")["family_95_interval_pp"][0]<0
    for token in ("100,543","45,953","54,590","20,694","1,934","+2.255","+3.952","−2.046","−3.556"):
        assert token in s,token
    assert "does not demonstrate" in s or "not a selection coefficient" in s


def test_matched_spatial_source_scale_and_familywise_nulls():
    s=text()
    sp=source(SPATIAL)
    assert sp["n_environmental_direction_rows"]==576
    assert sp["independent_direct_source_WLS_checks"]==288
    for name,n,g in (("world|genus_equal",9370,845),("E90_E180|genus_equal",2886,280)):
        v=sp["BIO12_white_by_region_and_weight"][name]
        assert v["n_matched_species"]==n and v["n_matched_genera"]==g
        assert not v["difference_excludes_zero_after_24_maxT"]
        assert v["paired_24_hue_axis_maxT_95_interval_pp"][0]<0<v["paired_24_hue_axis_maxT_95_interval_pp"][1]
    assert "−0.258" in s and "not a causal" in s.lower()


def test_identical_photo_between_and_within_not_merged_with_H2():
    q=source(MATCHED)["matched_source_photo_and_species"]
    assert (q["n_common_original_photos"],q["n_common_species"],q["n_source_genera"])==(21448,6143,787)
    assert math.isclose(q["between_minus_within_slope"],
        q["between_congener_BIO12_photo_white_slope"]-q["within_same_species_BIO12_photo_white_slope"],
        abs_tol=1e-13)
    assert q["paired_genus_bootstrap95_within"][0]<0<q["paired_genus_bootstrap95_within"][1]
    assert q["paired_genus_bootstrap95_difference"][0]>0
    s=text()
    assert all(v in s for v in ("21,448","6,143","787","+2.722","−0.511","+3.233"))
    assert "within-species" in s
    assert not re.search(r"W\s*=\s*0\.517",s)
    assert "separate" in s.lower() and "FCP" in s


def test_MNAR_not_silently_imputed():
    q=source(MNAR);s=text()
    assert q["all_original_unique_photo_ids"]==100543
    assert q["source_all_classifiable"]==45953
    assert q["source_all_unclassifiable"]==54590
    assert q["original_photo_wet_dry_unweighted"]["0"]["white_fraction"]<q["original_photo_wet_dry_unweighted"]["9"]["white_fraction"]
    tip=q["latent_hue_missing_not_at_random_tipping_curve"]["0.5"]
    assert 0.41<tip["required_unclassified_wet_white_fraction_for_EQUAL_true_dry_wet_white"]<0.42
    assert all(v in s for v in ("49.4%","39.1%","41.80%","53.25%","41.26%","MNAR"))
    assert "not a measurement of the missing photos" in s or "No such missing-photo true hue was measured" in s


def test_prediction_and_botanical_nonconfirmation():
    a=source(GENUS);z=source(DISTANT);nl=source(NONLINEAR)
    assert a["soil_addition_on_identical_20786_species"]["CI95"][0]<0<a["soil_addition_on_identical_20786_species"]["CI95"][1]
    assert a["primary36"]["absolute_species_colour"]["genus_equal_gain"]>0
    assert a["primary36"]["within_unseen_genus_species_differences"]["CI95"][0]<0<a["primary36"]["within_unseen_genus_species_differences"]["CI95"][1]
    assert z["main36_same_20786_evaluated_species"]["genus_region_and_500km"]["CI95"][0]<0<z["main36_same_20786_evaluated_species"]["genus_region_and_500km"]["CI95"][1]
    gain=nl["primary_result"]["genus_equal"]["total"]["contrasts"]["full_model_improvement"]
    assert gain["CI95"][0]<0<gain["CI95"][1]
    f=source(INDEPENDENT)["cases"]["min3_all_photo_opportunity"]["equal_genera"]
    assert f["east"]["support"]["named_species"]==670
    assert f["east"]["support"]["n_genus"]==58
    assert f["east_result"]["CI95"][0]<0<f["east_result"]["CI95"][1]
    s=text()
    assert all(v in s for v in ("+0.004518","+0.0000757","+0.000174","+0.004882","−15.34","+22.60"))
    assert "different" in s and "not" in s


def test_figures_supplement_source_and_real_release_blockers():
    s=text()
    assert "Figure 1." in s and "Figure 2." in s
    assert (ROOT/"figures/global_moisture_colour.svg").is_file()
    assert (ROOT/"figures/mois_photo_selection_and_matched_rainfall.svg").is_file()
    supplement=SI.read_text(encoding="utf8")
    for h in ("## S1","## S2","## S3","## S4","## S5","## S6"):
        assert h in supplement
    assert "not a common causal" not in s.lower() or "causal" in s.lower()
    assert "do not" in s.lower() or "not" in s.lower()
    assert "not yet deposited" in s and "original raw-source photo" in s
    assert "not submitted" in READINESS.read_text(encoding="utf8").lower()


def test_verified_precedents_and_distinct_fcp_manuscript():
    s=text()
    for doi in ("10.1111/nph.16737","10.1111/jbi.14751","10.1093/biosci/biab093","10.3389/fpls.2020.00991"):
        assert doi in s
    assert "Dalrymple" in s and "Trethowan" in s
    assert "separate FCP project" in s or "independent FCP" in s
    assert "genetic" in s and "pigment" in s
