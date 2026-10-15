"""Figure 5 guards: reproducible ecological synthesis from existing source aggregates.

No new fitting or classifying. Main between- and within-genus slope inference,
white-photo frequency, and genus signed contributions remain different targets.
"""
from __future__ import annotations

import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WHITE = ROOT/"results/mois_original_photo_species_genus_balance_20261014.json"
SCALES = ROOT/"results/mois_genus_between_within_moisture_20261015.json"
SIGNS = ROOT/"results/mois_congeneric_moisture_genus_sign_mosaic_20261015.json"
SVG = ROOT/"figures/mois_white_moisture_ecological_scales.svg"
SCRIPT = ROOT/"scripts/render_mois_figure5_ecological_scales.mjs"
MANUSCRIPT = ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md"


def test_same_source_receipts_and_distinct_denominators():
    w = json.loads(WHITE.read_text())
    scales = json.loads(SCALES.read_text())
    signs = json.loads(SIGNS.read_text())
    assert (w["n_classified_photos"], w["n_classified_white_photos"]) == (45953, 21574)
    assert (w["n_species_with_classified_photo"],w["n_genera_with_classified_photo"]) == (25064, 5369)
    assert (scales["source_cohort"]["n_primary_taxa"],scales["source_cohort"]["n_primary_genera"])==(20694,1934)
    assert (signs["source_support"]["n_primary_species"],signs["source_support"]["n_primary_genera"])==(20694,1934)
    assert math.isclose(w["point_white_proportion"]["photo_equal"],21574/45953,abs_tol=1e-14)
    assert .46<w["point_white_proportion"]["species_equal"]<.47
    assert .49<w["point_white_proportion"]["genus_equal"]<.50
    assert w["min_classified_species_thresholds"][3]["min_classified_species_per_genus"]==5
    assert .46<w["min_classified_species_thresholds"][3]["genus_equal_white"]<.47


def test_scale_ci_and_sign_mosaic_match_frozen_source():
    scales=json.loads(SCALES.read_text())
    signs=json.loads(SIGNS.read_text())
    axes=scales["estimates"]["axes"]
    for axis, expected in (("wc_bio12",2.254545235),("wc_vapr_annual_mean",3.951850325)):
        a=axes[axis]
        white=a["within_nominal_genus"]
        mean=a["between_genus_means"]
        assert math.isclose(white["point_pp_per_original_photo_sd"],expected,abs_tol=1e-8)
        assert 0<white["whole_genus_percentile_CI95_pp"][0]<white["point_pp_per_original_photo_sd"]
        assert 0<mean["whole_genus_percentile_CI95_pp"][0]<mean["point_pp_per_original_photo_sd"]
        assert mean["point_pp_per_original_photo_sd"]>white["point_pp_per_original_photo_sd"]
        sig=signs["cases"]["3"]["axes"][axis]
        assert sig["n_positive_contributing_genera"]+sig["n_negative_contributing_genera"]+sig["n_no_within_genus_photo_white_variation"]==1934
        assert math.isclose(sig["source_wls_slope_pp_per_environment_sd"],expected,abs_tol=1e-6)


def test_svg_well_formed_readable_and_no_sensitive_source_rows():
    raw=SVG.read_text(encoding="utf-8")
    doc=ET.fromstring(raw)
    assert doc.tag.endswith("svg")
    assert doc.attrib["width"]=="1560"
    assert doc.attrib["height"]=="870"
    assert doc.attrib["role"]=="img"
    ns="{http://www.w3.org/2000/svg}"
    assert doc.find(ns+"title") is not None
    assert doc.find(ns+"desc") is not None
    assert len(doc.findall(ns+"text"))>=35
    assert len(doc.findall(ns+"rect"))>=20
    assert len(doc.findall(ns+"line"))>=20
    for literal in ("GLOBAL NEAR-PARITY","20,694","1,934","830","719","385",
                    "Rain","Vapour","46.95%","46.82%","49.89%","46.63%"):
        assert literal in raw, literal
    for forbidden in ('inat_taxon_id','observation_id','photo_id','gps_coordinates','href="http','<image'):
        assert forbidden not in raw
    assert "no inferred pigment fitness" in raw
    assert "not original maxT" in raw
    source=SCRIPT.read_text(encoding="utf-8")
    for p in (WHITE,SCALES,SIGNS):
        assert p.name in source
    assert "Figure 5 disagrees with source receipts" in source


def test_article_figure_five_exists_and_keeps_reproductive_fitness_out_of_claims():
    ms=MANUSCRIPT.read_text(encoding="utf-8")
    assert ms.count("**Figure 5.")==1
    assert "(Fig. 5a)" in ms
    assert "Fig. 5b;" in ms and "Fig. 5c;" in ms
    assert "not global biological morph frequencies" in ms
    assert "not genus-specific statistically significant selection effects" in ms
    assert "figures/mois_white_moisture_ecological_scales.svg" in ms
    figleg=ms[ms.index("**Figure 5."):ms.index("## Supporting information")]
    assert "ordinary paired whole-genus bootstrap" in figleg
    assert "not the original 24-contrast maxT" in figleg
    assert "not global biological morph frequencies" in figleg
    assert "not genus-specific statistically significant selection effects" in figleg
    main=ms[ms.index("## 1 | Introduction"):ms.index("## Data and code availability")]
    assert len(main.split())<=5000
