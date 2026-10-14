"""Reproducible aggregate scientific checks of between-/within-genus moisture scales.

Original licensed source photo/environment ZIPs are not in public GitHub CI.
The raw-data script has been run with SHA-verified source ZIPs separately.
"""
from __future__ import annotations
import csv,importlib.util,json,math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/"scripts/analysis/replay_mois_genus_between_within_moisture_20261015.py"
RECEIPT=ROOT/"results/mois_genus_between_within_moisture_20261015.json"
LEDGER=ROOT/"results/fcp_macro_direction_20261013/direction_ledger.csv"


def source():
    return json.loads(RECEIPT.read_text(encoding="utf8"))


def replay_module():
    spec=importlib.util.spec_from_file_location("genus_moisture",CODE)
    a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
    return a


def test_synthetic_group_moments_do_not_confuse_levels():
    m=replay_module()
    rng=np.random.default_rng(571)
    g=np.repeat(np.arange(10),4)
    x_g=np.linspace(-1,1,10)
    data=pd.DataFrame({
       "genus":["g"+str(i) for i in g],
       "abs_lat":rng.normal(size=len(g)),
       "lon_sin":rng.normal(size=len(g)),
       "lon_cos":rng.normal(size=len(g)),
       "wc_bio12":x_g[g]+rng.normal(scale=.4,size=len(g)),
       "hue_white":(g>=5).astype(float),
    })
    norms={key:{"mean":float(data[key].mean()),"sd":float(data[key].std(ddof=0))}
           for key in (*m.GEO,"wc_bio12")}
    aw,bw,ab,bb,_=m.genus_moments(data,norms,"wc_bio12")
    assert aw.shape==(10,4,4) and ab.shape==(10,5,5)
    # A genus-fixed phenotype has no within-genus species response.
    betas=m.regress(aw.sum(axis=0),bw.sum(axis=0),ab.sum(axis=0),bb.sum(axis=0))
    assert abs(betas[0])<1e-10
    assert betas[1]>0.1
    assert abs(betas[2]-betas[1]+betas[0])<1e-10


def test_source_support_and_within_slope_equal_original_24_colour_ledger():
    q=source()
    assert q["schema"]=="mois_same_source_genus_between_within_moisture_v1"
    assert q["source_cohort"]["n_primary_taxa"]==20694
    assert q["source_cohort"]["n_primary_genera"]==1934
    assert q["source_cohort"]["n_primary_photos"]==38498
    assert q["estimates"]["n_bootstrap"]==1999
    assert q["estimates"]["bootstrap_genus_unit"]==1934
    with LEDGER.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    for axis,result in q["estimates"]["axes"].items():
        subset=[row for row in rows if
                row["region"]=="world" and row["mode"]=="climate" and
                row["quality"]=="all" and row["adjustment"]=="genus_geo" and
                row["weight"]=="genus_equal" and row["anchor"]==axis and
                row["hue"]=="white"]
        assert len(subset)==1
        old=float(subset[0]["direction_pp_per_env_global_sd"])
        within=result["within_nominal_genus"]["point_pp_per_original_photo_sd"]
        between=result["between_genus_means"]["point_pp_per_original_photo_sd"]
        difference=result["between_minus_within"]["point_pp_per_original_photo_sd"]
        assert math.isclose(old,within,abs_tol=1e-9)
        assert math.isclose(difference,between-within,abs_tol=1e-9)
        assert 0<within<between
        assert result["within_nominal_genus"]["whole_genus_percentile_CI95_pp"][0]>0
        assert result["between_genus_means"]["whole_genus_percentile_CI95_pp"][0]>0
        assert result["between_minus_within"]["whole_genus_percentile_CI95_pp"][0]>0
        # The original simultaneous interval remains distinct from new ordinary intervals.
        assert float(subset[0]["maxT_ci_lo_pp"])>0
    assert q["estimates"]["axes"]["wc_bio12"]["within_nominal_genus"]["point_pp_per_original_photo_sd"]>2.25
    assert q["estimates"]["axes"]["wc_bio12"]["between_genus_means"]["point_pp_per_original_photo_sd"]>8.55


def test_genus_depth_subsets_and_claim_limits():
    q=source()
    thresholds=q["descriptive_genus_support_thresholds_without_new_bootstrap_intervals"]
    assert [row["min_classified_climate_species_per_genus"] for row in thresholds]==[3,5,10,20]
    assert [(row["n_species"],row["n_genera"]) for row in thresholds]==[(20694,1934),(18164,1184),(14130,555),(9675,225)]
    for x in thresholds:
        for axis in ("wc_bio12","wc_vapr_annual_mean"):
            d=x["point_coefficients_pp"][axis]
            assert 0<d["within"]<d["between"]
            assert math.isclose(d["between_minus_within"],d["between"]-d["within"],abs_tol=1e-10)
    limits=" ".join(q["limitations"]).lower()
    assert "not a causal decomposition" in limits
    assert "genus labels are not a phylogeny" in limits
    assert "fitness" in limits and "pigment costs" in limits
    assert '"photo_id"' not in RECEIPT.read_text(encoding="utf8")
    assert '"inat_taxon_id"' not in RECEIPT.read_text(encoding="utf8")


def test_geb_claims_not_fitness_and_si_included():
    ms=(ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md").read_text(encoding="utf8")
    si=(ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md").read_text(encoding="utf8")
    memo=(ROOT/"docs/MOIS_GENUS_TURNOVER_AND_CONGENERIC_MOISTURE_20261015.md").read_text(encoding="utf8")
    assert "### S2.6 |" in si
    assert all(token in ms for token in ("+8.555","+2.255","+6.300","genus-mean","within-genus"))
    for token in ("20,694","1,934","38,498","+8.555","+9.370","+6.300","+5.419","1,999"):
        assert token in si and token in memo
    assert "not a causal attribution" in ms.lower()
    assert "not" in memo.lower() and "genetic" in memo.lower()
