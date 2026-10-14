"""Aggregate scientific and synthetic invariants for genus colour clustering.

The original attributed photo/site ZIP is not part of public CI: exact original
source replay is performed locally with SHA256-verified owner-only data.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/mois_original_genus_clustering_20261015.json"
SOURCE = ROOT / "scripts/analysis/mois_replay_genus_clustering_20261015.py"


def module():
    spec = importlib.util.spec_from_file_location("genus_cluster_source", SOURCE)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def test_illustrative_genus_clustering_math_and_boundaries():
    m = module()
    p = np.array([1,1,1,1,1,0,0,0,0,0],dtype=float)
    codes = np.array([0]*5+[1]*5,dtype=int)
    counts = np.array([5,5],dtype=int)
    np.testing.assert_allclose(m.clustering_metrics(p,codes,counts),[1,1,1])
    # Within each genus exactly the same white/coloured fraction -> no genus R2.
    p2 = np.array([1,1,0,0,0,1,1,0,0,0],dtype=float)
    assert abs(m.clustering_metrics(p2,codes,counts)[0]) < 1e-12


def test_real_source_aggregate_is_genus_clustered_not_balancing_selection():
    a = json.loads(RESULT.read_text(encoding="utf-8"))
    assert a["schema"] == "mois_photo_genus_clustering_source_perm_v1"
    assert a["original_photo_sha256"] == m_sha()
    assert a["data_n"] == {
        "unique_original_photos":100543,
        "source_classifiable_photos":45953,
        "taxa_with_classifiable_photos":25064,
        "genera_with_classifiable_photos":5369
    }
    comparisons=a["comparisons"]
    assert set(comparisons) == {"union_global","union_conditioned","breadth_global","breadth_conditioned"}
    for c in comparisons.values():
        assert c["permutations"] == 1999
        assert c["n_species"] >= 14000 and c["n_genera"] >= 1500
        for v in c["metrics"].values():
            assert v["upper_tail_permutation_p"] == 0.0005
            assert v["observed"] > v["null_95"][1] > v["null_95"][0]
    u = comparisons["union_global"]
    assert (u["n_species"],u["n_genera"],u["n_genera_ge5"]) == (20725,1935,1185)
    assert abs(u["metrics"]["genus_r2_species_white_fraction"]["observed"] - .343762761) < 1e-9
    assert abs(u["metrics"]["genus_r2_species_white_fraction"]["null_mean"] - .093329281) < 1e-9
    assert abs(u["metrics"]["genera_nearly_white_or_nonwhite_ge5_fraction"]["observed"] - 277/1185) < 1e-9
    assert comparisons["union_conditioned"]["conditioning_detail"]["n_blocks"] == 46
    assert comparisons["breadth_global"]["n_species"] == 14639
    assert "neither identifies phylogenetic signal" in a["interpretation"]
    for forbidden in ("photo_id","inat_taxon_id","longitude","original_coordinate"):
        assert '"'+forbidden+'"' not in RESULT.read_text(encoding="utf-8")


def m_sha():
    return "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"


def test_manuscript_marks_exploratory_nominal_genus_scope():
    ms=(ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md").read_text(encoding="utf-8")
    si=(ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md").read_text(encoding="utf-8")
    memo=(ROOT/"docs/MOIS_GENUS_COLOUR_CLUSTERING_AND_BALANCE_20261015.md").read_text(encoding="utf-8")
    assert "34.4%" in ms and "9.3%" in ms
    assert "### S2.5 |" in si
    for v in ("1,935","20,725","34.4%","9.3%","23.38%","2.19%","P*=0.0005"):
        assert v in si
    assert "no fitness" in memo.lower()
    assert "phylogenetic" in memo.lower()
