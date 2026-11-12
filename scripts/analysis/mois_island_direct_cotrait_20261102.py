#!/usr/bin/env python3
"""Species-direct floral-colour association with other direct traits, by rank.

Read SHA-verified source artifacts; export only aggregate counts/differences.
Underlying Island direct-trait rows have third-party rights. No image work,
hypothetical MNAR, imputation, or inference of pollinator guild.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ISLAND_ZIP = "5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd"
MOIS_ZIP = "3eb62f10607c13893e410dc1b7abdf9f90d4f1e730b08edecb245669823d1c63"
TRAIT_SHA = "7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577"
REGION_SHA = "a138ab56c0233b769f41c56c7eee7df82f09365d052983c1460f74f8cbf881f0"
COLOR_STATES = {
    "green_brown_inconspicuous": 1,
    "white": 0, "yellow_orange": 0, "red_pink": 0, "blue_purple": 0,
}
RESPONSE_STATES = {
    "small_flower": ("flower_size_class", {
        "very_small": 1, "small": 1, "medium": 0, "large": 0, "very_large": 0,
    }),
    "self_compatible": ("self_incompatibility", {"SC": 1, "SI": 0}),
    "open_radial_form": ("floral_form", {
        "open_radial": 1,
        **{x: 0 for x in ("reduced_wind", "bell_campanulate", "tubular", "salverform",
             "funnel_trumpet", "urn_urceolate", "brush_puff", "composite_head",
             "papilionaceous", "bilabiate", "spurred")},
    }),
    "predominantly_selfing": ("mating_system", {
        "predominantly_selfing": 1, "obligate_selfing": 1,
        "predominantly_outcrossing": 0,
    }),
    "cleistogamy_present": ("cleistogamy", {
        "obligate": 1, "facultative": 1, "absent": 0,
    }),
}
SECTORS = ("arctic", "north_subpolar", "north_temperate", "south_temperate",
           "island_cached", "mainland_cached")


def ensure(ok: bool, msg: str) -> None:
    if not bool(ok):
        raise ValueError("MOIS_DIRECT_COTRAIT: " + msg)


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for blob in iter(lambda: f.read(4 << 20), b""):
            h.update(blob)
    return h.hexdigest()


def read_sources(island_zip: Path, mois_zip: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ensure(sha_file(island_zip) == ISLAND_ZIP, "Island artifact ZIP SHA changed")
    ensure(sha_file(mois_zip) == MOIS_ZIP, "MOIS artifact ZIP SHA changed")
    with zipfile.ZipFile(island_zip) as z:
        blob = z.read("direct_species_trait_ledger.csv.gz")
        receipt = json.loads(z.read("integration_summary.json"))
        ensure(hashlib.sha256(blob).hexdigest() == TRAIT_SHA and
               receipt["output_sha256"]["direct_species_trait_ledger.csv.gz"] == TRAIT_SHA,
               "Island direct trait SHA disagrees with frozen scientific source")
        ensure(receipt["denominator_species"] == 106295 and receipt["filled"] == 222688,
               "Original Island scientific frame changed")
        source = pd.read_csv(io.BytesIO(blob), compression="gzip", dtype=str, low_memory=False)
    with zipfile.ZipFile(mois_zip) as z:
        blob = z.read("accepted_family_species_region.csv.gz")
        receipt = json.loads(z.read("family_holdout_receipt.json"))
        ensure(hashlib.sha256(blob).hexdigest() == REGION_SHA and
               receipt["out_SHA256"]["accepted_family_species_region.csv.gz"] == REGION_SHA,
               "Original MOIS regional species source SHA changed")
        taxa = pd.read_csv(io.BytesIO(blob), compression="gzip",
                           dtype={"accepted_id": str}, low_memory=False)
    ensure(len(taxa) == 30054 and taxa.accepted_id.nunique() == 22712,
           "MOIS original accepted species frame drifted")
    return source, taxa, {"Island": ISLAND_ZIP, "MOIS": MOIS_ZIP,
                          "Island_direct_trait": TRAIT_SHA, "MOIS_species_region": REGION_SHA}


def build_taxon_cotrait(source: pd.DataFrame, taxa: pd.DataFrame) -> pd.DataFrame:
    required = {"accepted_species", "trait_name", "normalized_value", "quality",
                "resolution_status", "classification"}
    ensure(required <= set(source.columns), "Island trait source schema changed")
    needed = {"accepted_id", "accepted_name", "genus", "family", "sector"}
    ensure(needed <= set(taxa.columns), "MOIS accepted species schema changed")
    ensure(not taxa.duplicated(["sector", "accepted_id"]).any(), "MOIS sector species duplicates")
    ensure(taxa.accepted_id.notna().all() and taxa.accepted_name.notna().all(),
           "Missing accepted identity")
    gr = taxa.groupby("accepted_id").agg(names=("accepted_name", "nunique"),
                                          families=("family", "nunique"),
                                          genera=("genus", "nunique"))
    ensure((gr == 1).all().all(), "Inconsistent taxon identity across overlapping sectors")
    unique = taxa.drop_duplicates("accepted_id")[["accepted_id", "accepted_name", "genus", "family"]].copy()
    ensure(unique.genus.eq(unique.accepted_name.str.split().str[0]).all(),
           "Genus and accepted species conflict")
    relevant = {"flower_primary_color", *(key for key, _ in RESPONSE_STATES.values())}
    source = source.loc[source.trait_name.isin(relevant)].copy()
    ensure(source.accepted_species.notna().all() and
           not source.duplicated(["accepted_species", "trait_name"]).any(),
           "Invalid/duplicate species-direct trait evidence")
    valid = (source.quality.isin(("high", "medium")) &
             source.resolution_status.eq("resolved") &
             ~source.classification.isin(("cultivar_contamination", "partial_self_compatibility")))
    source = source.loc[valid, ["accepted_species", "trait_name", "normalized_value"]]
    species = source.pivot(index="accepted_species", columns="trait_name", values="normalized_value")
    out = unique.join(species, on="accepted_name")
    ensure(out.accepted_id.is_unique and len(out) == 22712, "Source to accepted ID join multiplied rows")
    out["green_brown"] = out.flower_primary_color.map(COLOR_STATES).astype("float64")
    for label, (field, mapping) in RESPONSE_STATES.items():
        out[label] = out[field].map(mapping).astype("float64")
    return out


def paired_summary(frame: pd.DataFrame, label: str, scope: str) -> dict:
    """One accepted species per scope; no missing phenotype is made negative."""
    ensure(frame.accepted_id.is_unique, "A species counted twice within a scope")
    colour = frame.green_brown.notna()
    doc = frame[label].notna()
    known = frame.loc[colour & doc].copy()
    c1 = known.loc[known.green_brown.eq(1)]
    c0 = known.loc[known.green_brown.eq(0)]
    z = {"scope": scope, "trait": label,
         "total_accepted_species": len(frame),
         "direct_colour_evaluable": int(colour.sum()),
         "green_brown_single_state": int(frame.green_brown.eq(1).sum()),
         "other_single_colour": int(frame.green_brown.eq(0).sum()),
         "response_covered_among_green_brown": int((frame.green_brown.eq(1) & doc).sum()),
         "response_covered_among_other_colour": int((frame.green_brown.eq(0) & doc).sum()),
         "response_coverage_rate_green_brown":
             (len(c1) / int(frame.green_brown.eq(1).sum()) if frame.green_brown.eq(1).any() else None),
         "response_coverage_rate_other_colour":
             (len(c0) / int(frame.green_brown.eq(0).sum()) if frame.green_brown.eq(0).any() else None),
         "jointly_evaluable_accepted_species": len(known),
         "green_brown_evaluable": len(c1), "other_colour_evaluable": len(c0),
         "green_brown_and_positive": int(c1[label].sum()),
         "other_colour_and_positive": int(c0[label].sum()),
         "green_brown_positive_rate": float(c1[label].mean()) if len(c1) else None,
         "other_colour_positive_rate": float(c0[label].mean()) if len(c0) else None,
         "species_equal_diff_pp":
             (100 * (float(c1[label].mean()) - float(c0[label].mean()))
              if len(c1) and len(c0) else None),
         "taxon_inference_or_test_significance": False}
    for rank in ("genus", "family"):
        means = known.groupby([rank, "green_brown"], observed=True)[label].agg(["mean", "size"])
        m = means["mean"].unstack("green_brown")
        n = means["size"].unstack("green_brown")
        matched = m.dropna(subset=[0.0, 1.0]) if {0.0, 1.0}.issubset(m.columns) else m.iloc[:0]
        z["both_colour_" + rank + "_count"] = len(matched)
        z["within_" + rank + "_equal_diff_pp"] = (
            100 * float((matched[1.0] - matched[0.0]).mean()) if len(matched) else None)
        if len(matched):
            mm = n.loc[matched.index]
            z["within_" + rank + "_species_support"] = int((mm[1.0] + mm[0.0]).sum())
            z["within_" + rank + "_green_brown_support"] = int(mm[1.0].sum())
        else:
            z["within_" + rank + "_species_support"] = 0
            z["within_" + rank + "_green_brown_support"] = 0
    ensure(z["response_covered_among_green_brown"] == len(c1) and
           z["response_covered_among_other_colour"] == len(c0),
           "Response documentation and pair counts conflict")
    return z


def run(island_zip: Path, mois_zip: Path, output: Path):
    output = Path(output)
    ensure(not output.exists(), "Refuse to overwrite previous outputs")
    source, regions, inputs = read_sources(Path(island_zip), Path(mois_zip))
    frame = build_taxon_cotrait(source, regions)
    rows = []
    for label in RESPONSE_STATES:
        rows.append(paired_summary(frame, label, "global_unique_accepted_species"))
    for sector in SECTORS:
        subset = frame.loc[frame.accepted_id.isin(
            set(regions.loc[regions.sector.eq(sector), "accepted_id"]))]
        for label in RESPONSE_STATES:
            rows.append(paired_summary(subset, label, sector))
    final = pd.DataFrame(rows)
    ensure(len(final) == len(RESPONSE_STATES) * (len(SECTORS) + 1),
           "Expected cross-trait count changed")
    output.mkdir(parents=True)
    name = "floral_colour_cotrait_direct_summary.csv"
    final.to_csv(output / name, index=False, float_format="%.12g")
    report = {"schema": "mois_WCVP_Island_direct_species_cotrait_source_audit_v1",
              "source_sha256": inputs, "original_MOIS_accepted_species": len(frame),
              "sector_names": list(SECTORS), "response_names": list(RESPONSE_STATES),
              "region_memberships_overlap_and_are_not_independent": True,
              "restricted_species_row_evidence_not_redistributed": True,
              "outcomes_are_distinct_not_composite_dullness_score": True,
              "green_brown_is_species_description_not_optical_contrast": True,
              "photo_colour_ROI_does_not_validate_anatomical_flower": True,
              "missing_and_multistate_traits_not_converted_to_zero": True,
              "within_genus_comparisons_require_both_colour_groups": True,
              "no_confidence_intervals_or_multiple_trait_significance_claims": True,
              "SC_does_not_imply_autonomous_selfing": True,
              "no_new_photo_download_inference_or_MNAR_scenarios": True,
              "output_sha256": {name: sha_file(output / name)}}
    (output / "floral_cotrait_receipt.json").write_text(
        json.dumps(report, sort_keys=True, indent=2) + "\n")
    print("MOIS_DIRECT_COTRAIT_SUCCESS", json.dumps({
         "accepted_species": len(frame), "rows": len(final),
         "global": final.loc[final.scope.eq("global_unique_accepted_species"),
             ["trait", "jointly_evaluable_accepted_species", "green_brown_evaluable",
              "species_equal_diff_pp", "both_colour_genus_count",
              "within_genus_equal_diff_pp"]].to_dict("records")}, sort_keys=True))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--island-zip", type=Path, required=True)
    p.add_argument("--mois-zip", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    run(a.island_zip, a.mois_zip, a.output)


if __name__ == "__main__":
    main()
