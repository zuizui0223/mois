#!/usr/bin/env python3
"""Source-verified taxonomic clustering of photographed flower colour.

This is a retrospective permutation test of NOMINAL-genus photocolour
clustering in the historically frozen FCP photo table. No genetic, lifetime
fitness, balancing-selection, or species-phylogeny parameter is estimated.
Original photos, species IDs, names, sites and photo locations are never
written to the report. The original ZIP must be held and checked locally.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd

SOURCE_SHA256 = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
BRE = "fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz"
CEL = "fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz"
USE_COLS = ["inat_taxon_id", "photo_id", "species", "morph", "measurement_status", "wc_bio12", "latitude"]
VALID_HUES = {"white", "yellow_orange", "red_pink", "blue_purple"}
N_PERM = 1999
SEEDS = {"union_global":2026101501, "union_conditioned":2026101502,
         "breadth_global":2026101503, "breadth_conditioned":2026101504}


def file_sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(part)
    return digest.hexdigest()


def load_source(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if file_sha(path) != SOURCE_SHA256:
        raise ValueError("Wrong original photo archive SHA256")
    with zipfile.ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise ValueError("Original ZIP integrity check failed")
        def load(member: str) -> pd.DataFrame:
            with archive.open(member) as compressed, gzip.GzipFile(fileobj=compressed) as stream:
                return pd.read_csv(stream, usecols=USE_COLS, low_memory=False)
        breadth, cell = load(BRE), load(CEL)
    if (len(breadth),len(cell)) != (42111,85337):
        raise AssertionError("Source photographic opportunity count drift")
    combined = pd.concat([breadth,cell],ignore_index=True)
    if combined.photo_id.duplicated().sum() != 26905:
        raise AssertionError("Historical photo ID overlap drift")
    if combined.groupby("photo_id").inat_taxon_id.nunique().max() != 1:
        raise AssertionError("Conflicting source photo taxon identity")
    unique = combined.drop_duplicates("photo_id").copy()
    if len(unique) != 100543 or unique.inat_taxon_id.nunique() != 42111:
        raise AssertionError("Deduplicated source population drift")
    return unique, breadth


def species_photographic_fractions(photo: pd.DataFrame) -> pd.DataFrame:
    classified = photo.loc[
        photo.measurement_status.eq("classified_four_state_morph") & photo.morph.isin(VALID_HUES)
    ].copy()
    classified["genus"] = classified.species.str.split().str[0]
    classified["white"] = classified.morph.eq("white").astype(float)
    return (classified.groupby(["inat_taxon_id","genus"],sort=True)
             .agg(p=("white","mean"), n=("white","size"),rain=("wc_bio12","mean"),lat=("latitude","mean"))
             .reset_index())


def prepare_groups(species: pd.DataFrame, min_species: int = 3) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    counts = species.groupby("genus").size()
    eligible = species.loc[species.genus.isin(counts.index[counts>=min_species])].copy().reset_index(drop=True)
    code,_ = pd.factorize(eligible.genus,sort=True)
    count = np.bincount(code)
    assert count.min() >= min_species
    return eligible,code,count


def clustering_metrics(p: np.ndarray, group_idx: np.ndarray, group_n: np.ndarray) -> np.ndarray:
    mu = p.mean()
    avg = np.bincount(group_idx,weights=p,minlength=len(group_n))/group_n
    total = np.sum((p-mu)**2)
    if total <= 0:
        raise ValueError("No species-level visible-colour variation")
    between = np.sum(group_n*(avg-mu)**2)
    near = (avg<=.10)|(avg>=.90)
    return np.asarray([between/total, near.mean(), near[group_n>=5].mean()])


def groups_for_permutation(species: pd.DataFrame, conditioning: str) -> tuple[list[np.ndarray],dict]:
    n = len(species)
    if conditioning == "global":
        return [np.arange(n)],{"n_blocks":1,"photo_depth_categories":"not conditioned", "rainfall_quintile_categories":"not conditioned", "absolute_latitude_categories":"not conditioned"}
    if conditioning != "rainfall_latitude_photo_depth":
        raise ValueError("Unexpected fixed source conditioning")
    rain = species.rain.to_numpy(float)
    okay = np.isfinite(rain)
    edges = np.quantile(rain[okay],[.2,.4,.6,.8])
    rainbin = np.searchsorted(edges,rain,side="right")
    rainbin[~okay]=5
    lat=species.lat.to_numpy(float)
    latbin=np.digitize(np.abs(lat),[30,60]);latbin[~np.isfinite(lat)]=3
    depthbin=np.digitize(species.n.to_numpy(),[2,3])
    block=((rainbin*4+latbin)*3+depthbin).astype(int)
    return [np.flatnonzero(block==b) for b in np.unique(block)],{
        "n_blocks":len(np.unique(block)),
        "photo_depth_categories":"1, 2, >=3 classifiable photos/taxon",
        "rainfall_quintile_categories":"five original species-mean rainfall quintiles plus unknown",
        "absolute_latitude_categories":"0-30, 30-60, >=60 degrees plus unknown"
    }


def permutation_test(species: pd.DataFrame, conditioning: str, seed: int, reps: int = N_PERM) -> dict:
    selected,idx,n=prepare_groups(species,min_species=3)
    x=selected.p.to_numpy(float)
    obs=clustering_metrics(x,idx,n)
    blocks, detail=groups_for_permutation(selected,conditioning)
    rng=np.random.default_rng(seed)
    null=np.empty((reps,len(obs)),float)
    for j in range(reps):
        y=x.copy()
        for block in blocks:
            if len(block)>1:
                y[block]=x[rng.permutation(block)]
        null[j]=clustering_metrics(y,idx,n)
    qs=np.quantile(null,[.025,.975],axis=0)
    p_upper=(np.count_nonzero(null>=obs,axis=0)+1)/(reps+1)
    metric_names=("genus_r2_species_white_fraction", "genera_nearly_white_or_nonwhite_ge3_fraction",
                  "genera_nearly_white_or_nonwhite_ge5_fraction")
    return {
        "n_species":len(selected),"n_genera":len(n),"n_genera_ge5":int(np.sum(n>=5)),
        "species_equal_white_fraction":round(float(x.mean()),9),
        "conditioning":conditioning,"conditioning_detail":detail,
        "permutations":reps,"seed":seed,
        "metrics":{
            key:{"observed":round(float(obs[i]),9),
                 "null_mean":round(float(null[:,i].mean()),9),
                 "null_95":np.round(qs[:,i],9).tolist(),
                 "upper_tail_permutation_p":round(float(p_upper[i]),7)}
            for i,key in enumerate(metric_names)}
    }


def build_report(zip_path: Path) -> dict:
    allphoto,breadth=load_source(zip_path)
    union_species=species_photographic_fractions(allphoto)
    breadth_species=species_photographic_fractions(breadth)
    if (len(union_species),union_species.genus.nunique(),len(breadth_species),breadth_species.genus.nunique())!=(25064,5369,18457,4559):
        raise ValueError("Original source species/genus four-hue cohort drift")
    comparisons={}
    for key,species,conditioning in (
        ("union_global",union_species,"global"),
        ("union_conditioned",union_species,"rainfall_latitude_photo_depth"),
        ("breadth_global",breadth_species,"global"),
        ("breadth_conditioned",breadth_species,"rainfall_latitude_photo_depth"),
    ):
        comparisons[key]=permutation_test(species,conditioning,SEEDS[key])
    primary=comparisons["union_global"]["metrics"]
    if not (primary["genus_r2_species_white_fraction"]["observed"] > .30):
        raise AssertionError("Unexpected original-source taxonomic-colour clustering shift")
    return {
        "schema":"mois_photo_genus_clustering_source_perm_v1",
        "date":"2026-10-15",
        "original_photo_sha256":SOURCE_SHA256,
        "photographic_source":"original 42111 breadth and 85337 taxon-cell rows; 26905 photo-ID overlaps removed; source four-hue class labels unchanged",
        "data_n":{"unique_original_photos":len(allphoto),
                  "source_classifiable_photos":int((allphoto.measurement_status.eq("classified_four_state_morph") & allphoto.morph.isin(VALID_HUES)).sum()),
                  "taxa_with_classifiable_photos":len(union_species),
                  "genera_with_classifiable_photos":union_species.genus.nunique()},
        "analysis":"Taxon-level mean of classifiable photo-white(1)/coloured(0). Restrict to genera with >=3 classifiable taxa; weighted genus-mean explained fraction of across-taxon colour variance, and fractions of eligible genera with genus mean <=0.1 or >=0.9. Source-locked species-value permutation keeps genus sample sizes and p_i intact. Four exploratory tests: photo-union and original one-photo breadth, with global and rainfall/latitude/photo-depth conditional permutations. These are post hoc, not programme-wide error controlled.",
        "near_pure_threshold":"Taxon-mean classified-photo white <=10% or >=90%; not genetically fixed morphological states",
        "comparisons":comparisons,
        "interpretation":"Visible photographed flower-colour values cluster much more strongly within nominal genera than expected from reallocating the same observed species colour fractions among size-matched genera, including in a coarse environment/photo-depth stratification. This alone neither identifies phylogenetic signal beyond genus labels, nor local sympatric balance, nor a pigment synthesis cost or rare-morph fitness advantage. The global near-half genus-equal mean masks heterogeneous genus-level observational composition.",
        "publication_boundary":"No raw/source photo IDs, names, sites, photograph attributions or licensed original tables exported. Original photos and classification mechanism not remeasured; null permutations concern retrospective species fractions and observed genus memberships only. Raw original ZIP is needed for source replay; public CI tests aggregate-only and synthetic properties."
    }


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--photo-zip",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=ap.parse_args()
    report=build_report(args.photo_zip)
    if args.write:
        args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        print("Wrote source-only aggregate genus clustering",args.output)
    else:
        other=json.loads(args.output.read_text(encoding="utf-8"))
        assert other==report, "Frozen genus-clustering source receipt differs from source ZIP re-analysis"
        print("PASS original genus clustering photo source and all four null checks")


if __name__=="__main__":
    main()
