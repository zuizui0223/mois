#!/usr/bin/env python3
"""Recheck a deterministic, non-inferential GBIF specimen coverage pilot.

Requires the original 2960-row Dryad flower-colour CSV (sha256 locked).
Uses only the public GBIF occurrence search endpoint, unauthenticated, up to
30 queries (one per taxon, limit=1). Outputs aggregate counts with no taxon
names, GBIF IDs or specimen coordinates. NOT a GBIF occurrence download/DOI,
raster extraction, climate test or sample-representative ecological estimate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

SOURCE_SHA256 = "a253308785e4cbd0e361b3ca04cdfdae29c523eefa843375cebf8030ee0874af"
PILOT_NAMES_SHA256 = "d025092677a9377265d48d1484fa8ab9555788a76b3c3b08c04fdb2729f6b0be"
GBIF_ENDPOINT = "https://api.gbif.org/v1/occurrence/search"
FLOWER_HUES = {"black", "green", "orange", "pink", "purple", "red", "white", "yellow"}


def sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def pilot_sample(path: Path) -> tuple[list[dict], list[str]]:
    raw = path.read_bytes()
    if sha256(raw) != SOURCE_SHA256:
        raise ValueError("Original Dryad source SHA256 differs: refuse alternative sample")
    rows = list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    if len(rows) != 2960:
        raise ValueError("Expected exactly 2,960 published trait rows")
    by_clade = defaultdict(list)
    for row in rows:
        colour = row["flower_color"].strip().lower()
        if colour not in FLOWER_HUES:
            raise ValueError("Unrecognized original flower hue: do not silently label nonwhite")
        by_clade[row["clade"]].append(row)
    pairs = []
    for _, group in sorted(by_clade.items(), key=lambda item: (-len(item[1]), item[0])):
        white = sorted(x["species"].replace("_", " ") for x in group if x["flower_color"].strip().lower() == "white")
        coloured = sorted(x["species"].replace("_", " ") for x in group if x["flower_color"].strip().lower() != "white")
        if len(group) < 20 or not white or not coloured:
            continue
        pairs.append({"clade_size": len(group), "white": white[0], "nonwhite": coloured[0]})
        if len(pairs) == 15:
            break
    names = [value for group in pairs for value in (group["white"], group["nonwhite"])]
    if len(names) != 30 or sha256(("\n".join(names) + "\n").encode("utf-8")) != PILOT_NAMES_SHA256:
        raise ValueError("Deterministic source-name selection not identical to locked pilot")
    return pairs, names


def species_search_url(name: str) -> str:
    params = urllib.parse.urlencode({
        "scientificName": name,
        "basisOfRecord": "PRESERVED_SPECIMEN",
        "hasCoordinate": "true",
        "year": "2010,2026",
        "limit": 1,
    })
    return GBIF_ENDPOINT + "?" + params


def fetch_one(name: str) -> dict:
    req = urllib.request.Request(species_search_url(name), headers={
        "Accept": "application/json",
        "User-Agent": "mois-source-coverage-pilot/1.0 (public ecological research)"
    })
    # Avoid aggressive retries; fail clearly rather than manufacturing a zero.
    with urllib.request.urlopen(req, timeout=25) as response:
        body = json.load(response)
    if "count" not in body or "results" not in body:
        raise ValueError("GBIF API output schema unexpected")
    return body


def analyze(pairs: list[dict], fetcher=fetch_one, *, delay: float = 0.45) -> dict:
    result = []
    for g, pair in enumerate(pairs, 1):
        pairout = {"rank": g, "source_clade_species_n": pair["clade_size"]}
        for classlabel in ("white", "nonwhite"):
            name = pair[classlabel]
            body = fetcher(name)
            count = body["count"]
            if not isinstance(count, int) or count < 0:
                raise ValueError("Negative or non-integer GBIF search count")
            first = body["results"][0] if body["results"] else None
            if count > 0 and not first:
                raise ValueError("Positive count but no first indexed specimen record")
            if count == 0 and first:
                raise ValueError("Zero count but indexed specimen record returned")
            # This is a first-result audit, not a rigorous taxonomic or identity crosswalk.
            exact = None if first is None else (
                (first.get("species") or "").strip().casefold() == name.strip().casefold()
            )
            pairout[classlabel] = {"indexed_specimens": count, "first_result_species_exact": exact}
            if delay > 0:
                time.sleep(delay)
        result.append(pairout)

    arr = [p[t] for p in result for t in ("white", "nonwhite")]
    paired_both5 = [p for p in result if min(p["white"]["indexed_specimens"], p["nonwhite"]["indexed_specimens"]) >= 5]
    paired_exact = [p for p in result if p["white"]["first_result_species_exact"] is True and p["nonwhite"]["first_result_species_exact"] is True]
    paired_both = [p for p in paired_exact if min(p["white"]["indexed_specimens"], p["nonwhite"]["indexed_specimens"]) >= 5]
    counts = Counter()
    counts["nonzero_named_searches"] = sum(z["indexed_specimens"] > 0 for z in arr)
    counts["searches_with_at_least_5_raw_indexed_occurrences"] = sum(z["indexed_specimens"] >= 5 for z in arr)
    counts["searches_with_zero"] = sum(z["indexed_specimens"] == 0 for z in arr)
    counts["first_result_exact_binomial_matching_source"] = sum(z["first_result_species_exact"] is True for z in arr)
    counts["first_result_accepted_species_mismatch"] = sum(z["first_result_species_exact"] is False for z in arr)
    counts["n_hits_over_all_30_search_queries"] = sum(z["indexed_specimens"] for z in arr)
    for k, label in (("white", "white"), ("nonwhite", "nonwhite")):
        counts[label + "_nonzero"] = sum(p[k]["indexed_specimens"] > 0 for p in result)
        counts[label + "_first_exact"] = sum(p[k]["first_result_species_exact"] is True for p in result)
    counts["paired_clades_both_names_raw_indexed_count_at_least_5"] = len(paired_both5)
    counts["paired_clades_both_first_accepted_species_exact"] = len(paired_exact)
    counts["paired_clades_both_first_exact_and_raw_count_at_least_5"] = len(paired_both)
    return {"source_clade_sizes_in_sample_order": [p["clade_size"] for p in pairs],
            "n_clades_sampled": len(pairs), "n_species_sampled": len(arr),
            "totals": dict(counts), "anonymous_pairs": result,
            "not_a_climate_result": True,
            "limits": ("GBIF name-filter counts can include synonymous accepted taxa and do not "
                       "establish 5 spatially distinct specimen sites. First-result identity "
                       "is merely a warning flag; citable occurrence downloads remain required.")}


def self_test():
    pairs = [{"clade_size": 20, "white": "Genus albiflora", "nonwhite": "Genus purpurea"}]
    def mock(name):
        return {"count": 8, "results": [{"species": name}]} if "albiflora" in name else {
            "count": 4, "results": [{"species": "Other purpurea"}]}
    report = analyze(pairs, fetcher=mock, delay=0)
    assert report["totals"]["nonzero_named_searches"] == 2
    assert report["totals"]["searches_with_at_least_5_raw_indexed_occurrences"] == 1
    assert report["totals"]["first_result_exact_binomial_matching_source"] == 1
    assert report["totals"]["first_result_accepted_species_mismatch"] == 1
    assert report["totals"]["paired_clades_both_first_exact_and_raw_count_at_least_5"] == 0
    assert "scientificName=" in species_search_url("Genus albiflora")
    print("PASS synthetic GBIF search-count and synonym guard; NO live source queried")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--traits", type=Path, help="Locally authorized official final_dataset.csv")
    parser.add_argument("--output", type=Path, help="Aggregate-only JSON, not taxon/site rows")
    parser.add_argument("--delay", type=float, default=0.45)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.traits is None:
        parser.error("--traits official public Dryad CSV required")
    pairs, _ = pilot_sample(args.traits)
    summary = analyze(pairs, delay=args.delay)
    summary.update({"source_traits_sha256": SOURCE_SHA256,
                    "ordered_selected_names_sha256": PILOT_NAMES_SHA256,
                    "query_endpoint": GBIF_ENDPOINT,
                    "filter": "PRESERVED_SPECIMEN;hasCoordinate=true;year=2010,2026;limit=1",
                    "not_a_citable_gbif_bulk_download": True})
    payload = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    print(payload)


if __name__ == "__main__":
    main()
