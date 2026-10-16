#!/usr/bin/env python3
"""Post-ROI FCP flowering-season audit on an EXACT original monthly photo-ID source.

Source metadata are read BEFORE looking up new ROI-v4 outcomes. The original
observed_on, observation ID, taxon and coordinates are authoritative; old
149,900 photographic colour outcomes never enter this script.

Crucial: an unmeasured or ROI-unclassifiable image is NOT a non-flowering
observation. The measured, valid-date, confidently classifiable flower-ROI
subset is the only denominator for descriptive white proportions. A white+
yellow same-ROI organ-mixed image is excluded from the primary denominator.
Independent source observations (not photo counts) are the white-rate unit.
This produces descriptive month/site QC, never a causal climate estimate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sqlite3
from collections import defaultdict
from datetime import date
from pathlib import Path

SNAPSHOT = "2026-09-27"
MORPHS = {"white", "yellow_orange", "red_pink", "blue_purple"}
OK = "classified_four_state_morph"
SOURCE_STATUS = "COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS"
ROI_STATES = {
    "SEALED_PHOTO_SCREEN_UNBIASED_NOT_FINAL": "new_photo_table_sha256",
    "FIRST_GENUINE_MONTHLY_ORIGINAL_GEO_MATCHED_PIXELS_MEASURED_NOT_1000_SPECIES":
        "original_geographic_source_join_terminal_sha256",
}
SOURCE_COLS = {"inat_taxon_id","species","photo_id","observation_id","observer_id",
               "latitude","longitude","observed_on","source_snapshot_date"}
ROI_COLS = {"inat_taxon_id","species","photo_id","observation_id","observer_id",
            "morph","measurement_status","flower_fraction_white","flower_fraction_yellow"}
OUTCOLS = ("inat_taxon_id","species","year_month","calendar_month","hemisphere",
           "season","lat5","lon5","dated_original_photo_opportunities",
           "new_ROI_photos_with_known_date","classifiable_flower_ROIs",
           "same_ROI_white_yellow_organ_mixed","unclassifiable_ROIs",
           "usable_flower_photo_labels","independent_usable_observations",
           "white_observations","nonwhite_observations",
           "contradictory_multiphoto_observations","distinct_usable_observers",
           "white_rate_per_observation_descriptive")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def source_day(raw: str) -> date | None:
    try:
        # Do not turn unknown day precision, creation/upload date or 00 dates
        # into an invented day/month of observation.
        value=str(raw or "").strip()
        if len(value)!=10 or value[4]!="-" or value[7]!="-":
            return None
        day=date.fromisoformat(value)
        return day if day<=date.fromisoformat(SNAPSHOT) else None
    except ValueError:
        return None


def geo(raw_lat: str, raw_lon: str) -> tuple[str,str,int,int]:
    try:
        lat,lon=float(raw_lat),float(raw_lon)
        if not (math.isfinite(lat) and math.isfinite(lon) and
                -90<=lat<=90 and -180<=lon<=180):
            raise ValueError("Coordinate outside globe")
    except (ValueError,TypeError) as e:
        raise ValueError("Original photo source has invalid geographic coordinates") from e
    hemisphere = "tropical" if abs(lat)<23.5 else ("north" if lat>=0 else "south")
    lat5=min(35,int(math.floor((lat+90)/5)))
    lon5=int(math.floor(((lon+180)%360)/5))
    return hemisphere, f"{lat:.7f},{lon:.7f}", lat5, lon5


def season(month: int, hemisphere: str) -> str:
    """Local meteorological seasons; tropical sites retain month only.

    DJF/MAM/JJA/SON refer to months rather than the local season, so label
    summer/winter by hemisphere explicitly instead of miscalling SH July DJF.
    """
    if hemisphere=="tropical":
        return "tropical_month_only"
    if month in (12,1,2): northern="winter"
    elif month in (3,4,5): northern="spring"
    elif month in (6,7,8): northern="summer"
    elif month in (9,10,11): northern="autumn"
    else:
        raise ValueError("Impossible calendar month")
    if hemisphere=="north":
        return "local_"+northern
    southern={"winter":"summer","spring":"autumn",
              "summer":"winter","autumn":"spring"}
    return "local_"+southern[northern]


def verified_input(source: Path, source_receipt: Path, measured: Path,
                   roi_receipt: Path, *, synthetic: bool=False) -> tuple[dict,dict]:
    a=json.loads(source_receipt.read_text(encoding="utf-8"))
    b=json.loads(roi_receipt.read_text(encoding="utf-8"))
    if (a.get("status")!=SOURCE_STATUS or
        a.get("source_snapshot_date")!=SNAPSHOT or
        a.get("selected_all_photo_IDs_csv_sha256")!=sha(source) or
        a.get("historical_flower_colour_CSV_imported") is not False or
        a.get("original_image_pixels_opened") is not False):
        raise ValueError("Missing complete SHA-verified original dated photo-ID source")
    if synthetic:
        if b.get("status")!="SYNTHETIC_ORIGINAL_GEOGRAPHIC_PAIR_TEST_NOT_BIOLOGY":
            raise ValueError("Synthetic test must carry a synthetic ROI receipt")
        key="original_geographic_source_join_terminal_sha256"
    else:
        if b.get("status") not in ROI_STATES:
            raise ValueError("Only genuinely remeasured, sealed original ROI-v4 photos are eligible")
        key=ROI_STATES[b["status"]]
        if (b.get("old_measured_colours_imported") is not False and
            b.get("historical_measured_flower_colours_reused") is not False):
            raise ValueError("Historical flower colour measurements forbidden")
        if b["status"]=="SEALED_PHOTO_SCREEN_UNBIASED_NOT_FINAL" and (
                b.get("original_ROI_v4_method_attested") is not True):
            raise ValueError("Unverified original image measurement method")
    if b.get(key)!=sha(measured):
        raise ValueError("Original ROI terminal file SHA differs from receipt")
    return a,b


def analyze(source:Path,source_receipt:Path,measured:Path,roi_receipt:Path,
            output:Path,*,_synthetic_only:bool=False)->dict:
    if output.exists() and any(output.iterdir()):
        raise ValueError("Never overwrite signed original season/flower QC")
    src,roi=verified_input(source,source_receipt,measured,roi_receipt,
                           synthetic=_synthetic_only)
    output.mkdir(parents=True,exist_ok=True)
    db_path=output/"_monthly_join.sqlite"
    completed=False
    try:
        with sqlite3.connect(db_path) as db:
            db.execute("""CREATE TABLE roi(
                photo_id TEXT PRIMARY KEY, taxon TEXT, species TEXT,
                observation TEXT, observer TEXT, morph TEXT, status TEXT,
                white REAL, yellow REAL, observed_on TEXT, lat TEXT, lon TEXT,
                matched INT NOT NULL DEFAULT 0)""")
            db.execute("""CREATE TABLE source_ids(photo_id TEXT PRIMARY KEY)""")
            db.execute("""CREATE TABLE usable_observations(
                key TEXT, observation TEXT, observer TEXT, morph TEXT)""")
            n_measured=0
            with measured.open(newline="",encoding="utf-8") as f:
                reader=csv.DictReader(f)
                if not ROI_COLS.issubset(reader.fieldnames or []):
                    raise ValueError("Measured terminal missing frozen ROI-v4 or source join fields")
                for row in reader:
                    pid=row["photo_id"]
                    if not pid.isascii() or not pid.isdecimal():
                        raise ValueError("ROI terminal photo ID invalid")
                    good=row["measurement_status"]==OK
                    if good!=(row["morph"] in MORPHS):
                        raise ValueError("Invalid classifiable morphology/status pairing")
                    if not good and row["morph"]!="mixed_uncertain":
                        raise ValueError("Unclassifiable ROI cannot be a biological colour")
                    white=yellow=None
                    if good:
                        white=float(row["flower_fraction_white"])
                        yellow=float(row["flower_fraction_yellow"])
                        if not (math.isfinite(white) and math.isfinite(yellow) and
                                0<=white<=1 and 0<=yellow<=1):
                            raise ValueError("Invalid white/yellow ROI colour fractions")
                    if "global_classifiable" in row and row["global_classifiable"].strip().lower()!=str(good).lower():
                        raise ValueError("Global classifiability conflicts with original ROI status")
                    try:
                        db.execute("INSERT INTO roi VALUES (?,?,?,?,?,?,?,?,?,?,?,?,0)",
                            (pid,row["inat_taxon_id"],row["species"],
                             row["observation_id"],row["observer_id"],row["morph"],
                             row["measurement_status"],white,yellow,
                             row.get("observed_on",""),row.get("latitude",""),
                             row.get("longitude","")))
                    except sqlite3.IntegrityError as e:
                        raise ValueError("Duplicated measured ROI photo ID") from e
                    n_measured+=1
            db.commit()
            monthly=defaultdict(lambda:[0,0,0,0,0,0])
            names={}
            n_source=n_dated=n_missing=0
            n_classifiable=0
            n_organ_mixed=0
            n_undated_measured=0
            with source.open(newline="",encoding="utf-8") as f:
                reader=csv.DictReader(f)
                if not SOURCE_COLS.issubset(reader.fieldnames or []):
                    raise ValueError("Original source lacks observation date/ID/position")
                for row in reader:
                    n_source+=1
                    pid=row["photo_id"];tax=row["inat_taxon_id"]
                    if (row["source_snapshot_date"]!=SNAPSHOT or not pid.isascii()
                        or not pid.isdecimal() or not tax.isdecimal()
                        or not row["observation_id"] or not row["observer_id"]):
                        raise ValueError("Original dated source taxon/photo/observer identity invalid")
                    try:
                        db.execute("INSERT INTO source_ids VALUES (?)",(pid,))
                    except sqlite3.IntegrityError as e:
                        raise ValueError("Original source repeats a photo ID") from e
                    if tax in names and names[tax]!=row["species"]:
                        raise ValueError("Source taxon species name drift")
                    names[tax]=row["species"]
                    day=source_day(row["observed_on"])
                    hemisphere,coord,lat5,lon5=geo(row["latitude"],row["longitude"])
                    old=db.execute("SELECT taxon,species,observation,observer,morph,status,white,yellow,observed_on,lat,lon FROM roi WHERE photo_id=?",
                                   (pid,)).fetchone()
                    if old:
                        if (old[0:4]!=(tax,row["species"],row["observation_id"],row["observer_id"]) or
                            (old[8] and old[8]!=row["observed_on"])):
                            raise ValueError("ROI outcome source identity/date differs from original monthly record")
                        if old[9] and old[10]:
                            if (abs(float(old[9])-float(row["latitude"]))>1e-6 or
                                abs(float(old[10])-float(row["longitude"]))>1e-6):
                                raise ValueError("ROI source coordinates do not match monthly original")
                        db.execute("UPDATE roi SET matched=1 WHERE photo_id=?",(pid,))
                    if day is None:
                        n_missing+=1
                        if old:n_undated_measured+=1
                        continue
                    n_dated+=1
                    key=(tax,day.strftime("%Y-%m"),hemisphere,lat5,lon5)
                    stat=monthly[key]
                    stat[0]+=1
                    if not old:continue
                    stat[1]+=1
                    if old[5]!=OK:
                        stat[4]+=1
                        continue
                    stat[2]+=1
                    n_classifiable+=1
                    mixed=old[6]>=0.10 and old[7]>=0.10
                    if mixed:
                        stat[3]+=1
                        n_organ_mixed+=1
                        continue
                    stat[5]+=1
                    obs_key=json.dumps(key,separators=(",",":"))
                    db.execute("INSERT INTO usable_observations VALUES (?,?,?,?)",
                               (obs_key,row["observation_id"],row["observer_id"],old[4]))
            db.commit()
            actual_measured=db.execute("SELECT COUNT(*) FROM roi WHERE matched=1").fetchone()[0]
            if actual_measured!=n_measured:
                raise ValueError("New ROI terminal includes photo IDs outside the exact original source")
            if n_source!=src.get("n_selected_original_photo_IDs"):
                raise ValueError("Complete authentic source photo-ID count changed")
            if len(names)!=src.get("enrolled_outcome_blind_taxa"):
                raise ValueError("Enrolled original taxon frame was truncated")
            # Deduplicate multiple photographs from the SAME observation.
            # Conflicting photo labels or multiple observer IDs cannot support
            # an independent biological white/non-white observation.
            observation_counts=defaultdict(lambda:[0,0,0,set()])
            query="""SELECT key,observation,COUNT(DISTINCT observer),
                           COUNT(DISTINCT morph),MIN(morph),MIN(observer)
                     FROM usable_observations
                     GROUP BY key,observation"""
            for key,obs,n_owner,n_morph,morph,owner in db.execute(query):
                k=tuple(json.loads(key))
                rec=observation_counts[k]
                if n_morph!=1 or n_owner!=1:
                    rec[2]+=1
                else:
                    rec[0]+=int(morph=="white")
                    rec[1]+=int(morph!="white")
                    rec[3].add(owner)
            path=output/"original_dated_month_site_white_qc.csv"
            with path.open("x",newline="",encoding="utf-8") as f:
                writer=csv.DictWriter(f,fieldnames=OUTCOLS,lineterminator="\n")
                writer.writeheader()
                for k,stat in sorted(monthly.items()):
                    tax,ym,hemisphere,lat5,lon5=k
                    w,other,contradict,observers=observation_counts[k]
                    count=w+other
                    writer.writerow(dict(zip(OUTCOLS,(
                        tax,names[tax],ym,int(ym[-2:]),hemisphere,
                        season(int(ym[-2:]),hemisphere),lat5,lon5,
                        stat[0],stat[1],stat[2],stat[3],stat[4],stat[5],
                        count,w,other,contradict,len(observers),
                        f"{w/count:.8f}" if count else ""))))
        report={
            "schema":"fcp_original_dated_flowering_month_site_white_qc_v1",
            "status":("SYNTHETIC_TEST_NO_BIOLOGICAL_FINDING" if _synthetic_only else
                      "DESCRIPTIVE_ORIGINAL_ROI_MONTHLY_PHENOTYPE_NOT_GLOBAL_COMPLETION"),
            "source_snapshot_date":SNAPSHOT,
            "source_manifest_sha256":sha(source),
            "original_source_receipt_sha256":sha(source_receipt),
            "new_roi_terminal_sha256":sha(measured),
            "new_roi_receipt_sha256":sha(roi_receipt),
            "monthly_site_QC_csv_sha256":sha(path),
            "source_taxa":len(names),
            "all_original_photos_in_supplied_source":n_source,
            "original_photos_with_valid_observation_day":n_dated,
            "original_photos_with_invalid_or_unknown_observation_day":n_missing,
            "new_roi_terminal_photos":n_measured,
            "new_roi_photos_without_valid_observation_day":n_undated_measured,
            "new_roi_classifiable_flower_photos_with_valid_day":n_classifiable,
            "same_roi_white_yellow_organ_mixed_excluded_from_primary":n_organ_mixed,
            "original_source_is_not_confirmed_flowers":True,
            "roi_failed_is_not_called_nonflowering":True,
            "white_fraction_uses_independent_original_observations_not_photos":True,
            "calendar_date_was_never_used_to_select_or_classify_original_pixels":True,
            "not_genetic_polymorphism_or_environmental_causation":True,
            "full_42111_taxon_source_photo_census_complete":False,
            "all_93209351_photo_pixels_classified":False,
        }
        (output/"original_dated_month_site_receipt.json").write_text(
            json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
            encoding="utf-8")
        completed=True
        return report
    finally:
        db_path.unlink(missing_ok=True)
        if not completed:
            (output/"original_dated_month_site_white_qc.csv").unlink(missing_ok=True)
            (output/"original_dated_month_site_receipt.json").unlink(missing_ok=True)


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    for flag in ("original-source","original-source-receipt",
                 "new-roi-terminal","new-roi-receipt","outdir"):
        parser.add_argument("--"+flag,type=Path,required=True)
    a=parser.parse_args(argv)
    print(json.dumps(analyze(
        a.original_source,a.original_source_receipt,
        a.new_roi_terminal,a.new_roi_receipt,a.outdir),sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
