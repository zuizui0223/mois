#!/usr/bin/env python3
"""Source-only MOIS audit of 268 season-pairs and a blinded five-taxon photo reservation.

No pixels, flower stages, climatic predictors or flower colour are read.
Public iNaturalist source photographs retain their original licence and attribution;
individual photo manifest is uploaded as an Actions artifact, not Git-tracked.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

SOURCE_SHA = "58e037a210e23d3ff5d6f04d392665ada5ee9ac5149d0676672c239cf092df01"
SITE_SHA = "7691c7cf975da4e73494ad3fa0f2cf8df52302e84206017dd525051cd672b432"
PRIORITY_SHA = "ef95638da60b010c8671877854163132be304f5b7245edf9c5e2f9536265102a"
MOIS_SHA = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
PILOT_SHA = "d5c4bdcd632d0de8f2d20a663674e87f35e323a42b896291efefb93c46023053"
SEED = "MOIS-5taxa-independent-observer-siteyear-source-20261017-v1"
MORPH_STATUS = "source-only pre-colour photo opportunity"

AUDIT_COLUMNS = (
    "inat_taxon_id", "species", "site_cell_lat_0p25", "site_cell_lon_0p25",
    "season_year", "local_season_a", "local_season_b",
    "source_independent_obs_a", "source_independent_obs_b",
    "source_photos_original_a", "source_photos_original_b",
    "excluded_prior_mois_and_pilot_a", "excluded_prior_mois_and_pilot_b",
    "distinct_photo_source_observations_after_exclusion_a",
    "distinct_photo_source_observations_after_exclusion_b",
    "distinct_observers_a", "distinct_observers_b",
    "cross_season_overlapping_observer_count",
    "max_balanced_independent_observer_disjoint_per_season",
    "min_original_observations_per_season",
    "no_biological_hue_or_flowering_status_used",
)
PHOTO_COLUMNS = (
    "inat_taxon_id", "species", "photo_id", "observation_id", "observer_id",
    "observed_on", "photo_license", "photo_creator_id", "attribution",
    "photo_url_large", "source_snapshot_date", "site_latcell", "site_loncell",
    "season_year", "selected_season", "source_rank", "wave_photo_stage",
)


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()


def assert_sha(path: Path, expected: str) -> None:
    observed = file_sha(path)
    if observed != expected:
        raise ValueError(f"Input/source checksum mismatch for {path.name}: {observed}")


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def season(month: int, lat: float) -> str:
    if abs(lat) < 23.5:
        return "tropical_no_four_season_label"
    s = ("winter" if month in (12,1,2) else "spring" if month in (3,4,5)
         else "summer" if month in (6,7,8) else "autumn")
    if lat < 0:
        s = {"winter":"summer", "spring":"autumn", "summer":"winter",
             "autumn":"spring"}[s]
    return s


def balanced_observer_capacity(a: set[str], b: set[str]) -> int:
    """Max same number per season while forbidding cross-season observer reuse."""
    return min(len(a), len(b), (len(a | b)) // 2)


def rank_photo(pid: str) -> str:
    return hashlib.sha256((SEED + "|" + pid).encode("utf-8")).hexdigest()


def old_mois_photo_ids(path: Path) -> set[str]:
    assert_sha(path, MOIS_SHA)
    ids = set()
    members = [
        "fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz",
        "fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz",
    ]
    with zipfile.ZipFile(path) as z:
        for member in members:
            with gzip.open(io.BytesIO(z.read(member)), "rt", newline="") as f:
                for r in csv.DictReader(f):
                    ids.add(r["photo_id"])
    if len(ids) != 100543:
        raise ValueError("Historical MOIS source IDs do not total 100543")
    return ids


def write_rows(path: Path, fields, rows) -> None:
    with path.open("x", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k:r.get(k,"") for k in fields})


def run(source: Path, site: Path, site_receipt: Path, priority: Path,
        old_mois_zip: Path, pilot: Path, pilot_receipt: Path, outdir: Path) -> dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Never overwrite an original frozen source/photo reservation")
    assert_sha(source, SOURCE_SHA)
    assert_sha(site, SITE_SHA)
    assert_sha(priority, PRIORITY_SHA)
    assert_sha(pilot, PILOT_SHA)
    document = json.loads(site_receipt.read_text(encoding="utf-8"))
    if (document.get("status") != "COMPLETE_ORIGINAL_SOURCE_CALENDAR_CENSUS_NO_FLOWER_PIXELS"
        or document.get("original_source_site_month_file_sha256") != SITE_SHA
        or document.get("n_original_photo_ids") != 3430502
        or document.get("source_photo_opportunity_is_not_flowering_detection") is not True):
        raise ValueError("Not authentic 2026-09-27 iNaturalist photo calendar")
    pilot_data = json.loads(pilot_receipt.read_text(encoding="utf-8"))
    if (pilot_data.get("original_photo_attempts_sha256") != PILOT_SHA
        or pilot_data.get("n_original_photo_attempts") != 64
        or pilot_data.get("open_flower_botanically_validated") is not False):
        raise ValueError("Pilot original observed-on/source identities not frozen")
    old_ids = old_mois_photo_ids(old_mois_zip)
    prior = {r["photo_id"] for r in read_csv(pilot)}
    if len(prior) != 64:
        raise ValueError("Pilot photo-ID exclusion count must be 64")

    rank = read_csv(priority)
    if len(rank) != 268 or len({r["inat_taxon_id"] for r in rank}) != 39:
        raise ValueError("Original independently generated 268/39 source priority has changed")
    by_key = {}
    for q in rank:
        key = (q["inat_taxon_id"],int(q["site_cell_lat_0p25"]),
               int(q["site_cell_lon_0p25"]),int(q["season_year"]))
        if key in by_key:
            raise ValueError("Duplicate ranked source site-year")
        by_key[key] = q
    allowed = set()
    with site.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            ym = r["year_month"]
            try:
                y,m = int(ym[:4]),int(ym[5:7])
                ilat,ilon = int(r["lat_0p25_cell"]),int(r["lon_0p25_cell"])
            except (ValueError,TypeError):
                continue
            sy = y + (m == 12)
            key = (r["inat_taxon_id"], ilat, ilon, sy)
            q = by_key.get(key)
            if q is None:
                continue
            ss = season(m, -90 + (ilat + .5)/4)
            if (ss in (q["local_season_a"],q["local_season_b"])
                and int(r["distinct_source_observations"]) >= 2
                and int(r["distinct_source_observers"]) >= 2):
                allowed.add((r["inat_taxon_id"],ilat,ilon,ym))

    best = defaultdict(lambda: defaultdict(dict))
    all_observations = defaultdict(lambda: defaultdict(set))
    raw_season = Counter()
    exclude_season = Counter()
    n_rows = 0
    with source.open(newline="",encoding="utf-8") as f:
        for p in csv.DictReader(f):
            n_rows += 1
            tax=p["inat_taxon_id"]
            try:
                acc=float(p["positional_accuracy_m"])
                lat,lon=float(p["latitude"]),float(p["longitude"])
                ym=p["observed_on"][:7]
                y,m=int(ym[:4]),int(ym[5:7])
            except (ValueError,TypeError):
                continue
            if not (math.isfinite(acc) and 0 < acc <= 5000 and
                    math.isfinite(lat) and math.isfinite(lon)):
                continue
            ilat,ilon=int(math.floor((lat+90)*4)),int(math.floor((lon+180)*4))
            if not 0<=ilat<=719 or not 0<=ilon<=1439:
                continue
            sy = y + (m == 12)
            key = (tax,ilat,ilon,sy)
            q = by_key.get(key)
            if q is None or (tax,ilat,ilon,ym) not in allowed:
                continue
            ss=season(m,lat)
            if ss not in (q["local_season_a"],q["local_season_b"]):
                continue
            raw_season[(key,ss)] += 1
            if p["photo_id"] in old_ids or p["photo_id"] in prior:
                exclude_season[(key,ss)] += 1
                continue
            user,obs,pid = p["observer_id"],p["observation_id"],p["photo_id"]
            if not user or not obs or not pid:
                continue
            all_observations[key][ss].add(obs)
            existing=best[key][ss].get(user)
            pid_rank=rank_photo(pid)
            if existing is None or pid_rank < existing["source_rank"]:
                best[key][ss][user]={**p, "source_rank":pid_rank,
                                     "site_latcell":ilat,"site_loncell":ilon,
                                     "season_year":sy,"selected_season":ss,
                                     "wave_photo_stage":MORPH_STATUS}
    if n_rows != 3430502:
        raise ValueError("Exact 3.43M real photo-source census incomplete")

    audits=[]
    for key,q in by_key.items():
        a,b = q["local_season_a"],q["local_season_b"]
        A,B = set(best[key][a]),set(best[key][b])
        audits.append({
            "inat_taxon_id":q["inat_taxon_id"],"species":q["species"],
            "site_cell_lat_0p25":key[1],"site_cell_lon_0p25":key[2],
            "season_year":key[3],"local_season_a":a,"local_season_b":b,
            "source_independent_obs_a":q["independent_source_observations_a"],
            "source_independent_obs_b":q["independent_source_observations_b"],
            "source_photos_original_a":raw_season[(key,a)],
            "source_photos_original_b":raw_season[(key,b)],
            "excluded_prior_mois_and_pilot_a":exclude_season[(key,a)],
            "excluded_prior_mois_and_pilot_b":exclude_season[(key,b)],
            "distinct_photo_source_observations_after_exclusion_a":len(all_observations[key][a]),
            "distinct_photo_source_observations_after_exclusion_b":len(all_observations[key][b]),
            "distinct_observers_a":len(A),"distinct_observers_b":len(B),
            "cross_season_overlapping_observer_count":len(A&B),
            "max_balanced_independent_observer_disjoint_per_season":
                balanced_observer_capacity(A,B),
            "min_original_observations_per_season":
                q["min_obs_per_season"],
            "no_biological_hue_or_flowering_status_used":True,
        })
    audits.sort(key=lambda r:(-r["max_balanced_independent_observer_disjoint_per_season"],
              -min(r["distinct_observers_a"],r["distinct_observers_b"]),
              -int(r["min_original_observations_per_season"]),
              int(r["inat_taxon_id"]),r["site_cell_lat_0p25"],r["site_cell_lon_0p25"]))
    target=[]
    chosen_taxa=set()
    for q in audits:
        if q["max_balanced_independent_observer_disjoint_per_season"] < 40:
            break
        if q["inat_taxon_id"] not in chosen_taxa:
            chosen_taxa.add(q["inat_taxon_id"])
            target.append(q)
    if len(target) != 5:
        raise ValueError("Five observer-balanced taxa not available")

    selected=[]
    global_users=set()
    elig=[]
    for q in target:
        key=(q["inat_taxon_id"],int(q["site_cell_lat_0p25"]),
             int(q["site_cell_lon_0p25"]),int(q["season_year"]))
        a,b=q["local_season_a"],q["local_season_b"]
        users={s:{u:r for u,r in best[key][s].items() if u not in global_users}
               for s in (a,b)}
        shared=sorted(set(users[a])&set(users[b]),
            key=lambda u:hashlib.sha256((SEED+"|"+str(key)+"|"+u).encode()).hexdigest())
        solo={a:{u:r for u,r in users[a].items() if u not in users[b]},
              b:{u:r for u,r in users[b].items() if u not in users[a]}}
        choice={s:sorted(solo[s].values(),key=lambda r:r["source_rank"])[:40]
                for s in (a,b)}
        for u in shared:
            wants=[s for s in (a,b) if len(choice[s]) < 40]
            if not wants:
                break
            use=min(wants,key=lambda s:(len(choice[s]),s))
            choice[use].append(users[use][u])
        if not all(len(choice[s])==40 for s in (a,b)):
            raise ValueError("Independent observer-balanced colour-blind target unachievable")
        for s in (a,b):
            selected.extend(choice[s])
        source_users={r["observer_id"] for s in (a,b) for r in choice[s]}
        if len(source_users)!=80 or global_users & source_users:
            raise ValueError("Observer appears twice within/across taxa or seasons")
        global_users.update(source_users)
        elig.append({
            "inat_taxon_id":q["inat_taxon_id"],"species":q["species"],
            "latcell":key[1],"loncell":key[2],"year":key[3],
            "season_a":a,"season_b":b,
            "obs_a":q["source_independent_obs_a"],"obs_b":q["source_independent_obs_b"],
            "photos_a":raw_season[(key,a)],"photos_b":raw_season[(key,b)],
            "min_obs":int(q["min_original_observations_per_season"]),
            "min_photos":min(raw_season[(key,a)],raw_season[(key,b)]),
            "original_photo_candidates_a":raw_season[(key,a)],
            "original_photo_candidates_b":raw_season[(key,b)],
            "excluded_historical_mois_photo_ids_a":exclude_season[(key,a)],
            "excluded_historical_mois_photo_ids_b":exclude_season[(key,b)],
            "excluded_previous_2025_pilot_a":0,
            "excluded_previous_2025_pilot_b":0,
            "unique_source_observers_a":len(users[a]),
            "unique_source_observers_b":len(users[b]),
            "shared_observers_between_source_seasons":len(shared),
            "reserved_observer_disjoint_photo_ids_a":len(choice[a]),
            "reserved_observer_disjoint_photo_ids_b":len(choice[b]),
            "reserved_source_is_not_flowering_photo":True,
        })
    if len(selected)!=400 or len(global_users)!=400:
        raise ValueError("400 original photo IDs/observers not distinct")
    ids={r["photo_id"] for r in selected}
    if len(ids)!=400 or ids & old_ids or ids & prior:
        raise ValueError("Selection reuses historical MOIS or known pilot photo")

    outdir.mkdir(parents=True,exist_ok=True)
    audit_path=outdir/"all_268_actual_observer_disjoint_source_only_audit.csv"
    write_rows(audit_path, AUDIT_COLUMNS, audits)
    photo_path=outdir/"first_five_taxa_temporal_colour_blind_original_photo_manifest.csv"
    selected.sort(key=lambda r:(r["inat_taxon_id"],r["season_year"],
                   r["site_latcell"],r["selected_season"],r["source_rank"]))
    write_rows(photo_path,PHOTO_COLUMNS,selected)
    elig_path=outdir/"first_five_taxa_temporal_source_only_eligibility.csv"
    write_rows(elig_path,tuple(elig[0]),elig)
    summary={
       "status":"ORIGINAL_PHOTO_SOURCE_ONLY_NOT_A_BLOOM_OR_COLOUR",
       "n_source_original_photos":n_rows,"n_ranked_source_siteyear_pairs":len(audits),
       "n_ranked_taxa":len({r["inat_taxon_id"] for r in audits}),
       "actual_capacity_by_minimum_observer_disjoint_per_season":{
          str(t):{"siteyears":sum(r["max_balanced_independent_observer_disjoint_per_season"]>=t for r in audits),
                  "taxa":len({r["inat_taxon_id"] for r in audits if
                             r["max_balanced_independent_observer_disjoint_per_season"]>=t})}
          for t in (5,10,15,20,24,30,40,50,64,80)},
       "selected_taxa":len(target),"selected_original_photo_ids":len(selected),
       "selected_distinct_source_observers":len(global_users),
       "selected_photo_manifest_sha256":file_sha(photo_path),
       "observer_capacity_audit_sha256":file_sha(audit_path),
       "selected_siteyear_eligibility_sha256":file_sha(elig_path),
       "did_not_open_or_reclassify_photo_pixels":True,
       "did_not_infer_flower_opening_or_white_fraction":True,
    }
    (outdir/"mois_observer_disjoint_source_photo_selection_receipt.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return summary


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    for arg in ("source","site","site-receipt","priority","mois-zip","pilot","pilot-receipt","outdir"):
        ap.add_argument("--"+arg,required=True,type=Path)
    a=ap.parse_args()
    print("MOIS_ORIGINAL_SOURCE_OBSERVER_DISJOINT_WAVE",
          json.dumps(run(a.source,a.site,a.site_receipt,a.priority,
                         a.mois_zip,a.pilot,a.pilot_receipt,a.outdir),sort_keys=True))
