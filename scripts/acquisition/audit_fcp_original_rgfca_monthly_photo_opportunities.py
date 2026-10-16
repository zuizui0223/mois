#!/usr/bin/env python3
"""Census ALL observed-on dates in a validated RGFCA original photo-ID source.

Original source opportunity, NOT a flowering-image detector or a colour
phenotype. Each candidate photo remains eligible, including missing dates.
The site-month overlap audit uses only independent original observations,
original observers and conservative <=5 km recorded position accuracy;
none of those data enter the frozen ROI-v4 image evaluator.

The 1,000 RGFCA taxa are a first source screen, not a cap on all 42,111 taxa.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sqlite3
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

SOURCE_DATE = "2026-09-27"
ORIGINAL_1000 = 1000
ORIGINAL_PHOTOS = 3_430_502
SOURCE_STATE = "COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS"
REQUIRED = (
    "inat_taxon_id","species","photo_id","observation_id","observer_id",
    "latitude","longitude","observed_on","positional_accuracy_m",
    "source_snapshot_date",
)
MONTHS = ("inat_taxon_id","species","year_month","hemisphere",
          "all_original_date_known_photo_ids","distinct_original_observations",
          "distinct_original_observers","unmeasured_source_photo_opportunities_only")
SITES = ("inat_taxon_id","species","year_month","lat_0p25_cell",
         "lon_0p25_cell","all_source_photo_ids_with_5km_accuracy",
         "distinct_source_observations","distinct_source_observers",
         "date_and_site_are_original_not_pixel_inferred")
TAXA = ("inat_taxon_id","species","all_original_photo_ids",
        "photos_with_valid_original_observed_day",
        "photos_without_valid_original_observed_day",
        "distinct_original_months_with_dated_candidates",
        "dated_photos_at_positions_accuracy_0_to_5km_exclusive_zero",
        "independent_site_years_with_two_supported_observation_months",
        "do_not_call_any_candidate_a_blooming_flower")


def digest(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def precise_day(value:str)->str:
    value=str(value or "").strip()
    if len(value)!=10 or value[4]!="-" or value[7]!="-":
        return ""
    try:
        day=date.fromisoformat(value)
    except ValueError:
        return ""
    return day.strftime("%Y-%m") if day<=date.fromisoformat(SOURCE_DATE) else ""


def source_coordinates(lat,lon,acc)->tuple[str,int,int,bool]:
    try:
        x,y=float(lat),float(lon)
    except (ValueError,TypeError) as err:
        raise ValueError("Original monthly source location invalid") from err
    if not (math.isfinite(x) and math.isfinite(y) and
            -90<=x<=90 and -180<=y<=180):
        raise ValueError("Original monthly coordinates outside globe")
    hemisphere=("tropical" if abs(x)<23.5 else "north" if x>=0 else "south")
    latcell=min(719,int(math.floor((x+90)*4)))
    loncell=int(math.floor(((y+180)%360)*4))
    try:
        accuracy=float(acc)
        precise=math.isfinite(accuracy) and 0<accuracy<=5000
    except (ValueError,TypeError):
        precise=False
    return hemisphere,latcell,loncell,precise


def execute(manifest:Path,source_receipt:Path,out:Path,*,_small_test:bool=False)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Do not overwrite complete real dated-source opportunity audit")
    doc=json.loads(source_receipt.read_text(encoding="utf-8"))
    if (doc.get("status")!=SOURCE_STATE or
        doc.get("source_snapshot_date")!=SOURCE_DATE or
        doc.get("selected_all_photo_IDs_csv_sha256")!=digest(manifest) or
        doc.get("original_image_pixels_opened") is not False or
        doc.get("historical_flower_colour_CSV_imported") is not False or
        doc.get("all_enrolled_species_match_independent_full_source_photo_ID_counts") is not True):
        raise ValueError("Not the exact complete dated original 1000 photo-ID source")
    if not _small_test and (
            doc.get("n_selected_original_photo_IDs")!=ORIGINAL_PHOTOS or
            doc.get("enrolled_outcome_blind_taxa")!=ORIGINAL_1000):
        raise ValueError("Full original 1000 species and photo count gate not satisfied")

    out.mkdir(parents=True,exist_ok=True)
    db_path=out/"_complete_dated_source_audit.sqlite"
    complete=False
    try:
        with sqlite3.connect(db_path) as db:
            db.execute("PRAGMA journal_mode=DELETE")
            db.execute("PRAGMA synchronous=NORMAL")
            db.execute("PRAGMA cache_size=-96000")
            db.execute("CREATE TABLE ids(photo_id TEXT PRIMARY KEY)")
            db.execute("""CREATE TABLE observations (
                observation TEXT PRIMARY KEY, taxon TEXT, species TEXT,
                observer TEXT, ym TEXT, hemisphere TEXT,
                ilat INTEGER, ilon INTEGER, precise INTEGER)""")
            db.execute("""CREATE TABLE months(
                taxon TEXT, ym TEXT, hemisphere TEXT, n_photos INTEGER,
                PRIMARY KEY(taxon,ym,hemisphere))""")
            db.execute("""CREATE TABLE sites(
                taxon TEXT, ym TEXT, ilat INTEGER, ilon INTEGER,
                n_photos INTEGER, PRIMARY KEY(taxon,ym,ilat,ilon))""")
            source_counts=Counter()
            dated_counts=Counter()
            site_counts=Counter()
            names={}
            n_source=n_dated=n_unknown=n_precise=0
            with manifest.open(newline="",encoding="utf-8") as file:
                reader=csv.DictReader(file)
                if not set(REQUIRED).issubset(reader.fieldnames or []):
                    raise ValueError("Original monthly photo-ID source no longer includes dates")
                for row in reader:
                    pid=row["photo_id"]
                    taxon=row["inat_taxon_id"]
                    ident=row["observation_id"]
                    observer=row["observer_id"]
                    if (row["source_snapshot_date"]!=SOURCE_DATE or
                        not pid.isascii() or not pid.isdecimal() or
                        not taxon.isascii() or not taxon.isdecimal() or
                        not ident or not observer or not row["species"]):
                        raise ValueError("Original monthly observation/photo source identity invalid")
                    if taxon in names and names[taxon]!=row["species"]:
                        raise ValueError("Original taxon name drift")
                    names[taxon]=row["species"]
                    ym=precise_day(row["observed_on"])
                    hem,ilat,ilon,precise=source_coordinates(
                        row["latitude"],row["longitude"],
                        row["positional_accuracy_m"])
                    try:
                        db.execute("INSERT INTO ids VALUES (?)",(pid,))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError("Original source photo ID repeated globally") from exc
                    db.execute("""INSERT OR IGNORE INTO observations VALUES(?,?,?,?,?,?,?,?,?)""",
                        (ident,taxon,row["species"],observer,ym,hem,ilat,ilon,int(precise)))
                    n_source+=1
                    source_counts[taxon]+=1
                    if not ym:
                        n_unknown+=1
                        continue
                    n_dated+=1
                    dated_counts[taxon]+=1
                    db.execute("""INSERT INTO months VALUES (?,?,?,1)
                        ON CONFLICT(taxon,ym,hemisphere) DO UPDATE
                        SET n_photos=n_photos+1""",(taxon,ym,hem))
                    if precise:
                        n_precise+=1
                        site_counts[taxon]+=1
                        db.execute("""INSERT INTO sites VALUES(?,?,?,?,1)
                            ON CONFLICT(taxon,ym,ilat,ilon) DO UPDATE
                            SET n_photos=n_photos+1""",(taxon,ym,ilat,ilon))
                    if n_source%20000==0:
                        db.commit()
            db.commit()
            if n_source!=doc.get("n_selected_original_photo_IDs") or (
                  len(names)!=doc.get("enrolled_outcome_blind_taxa")):
                raise ValueError("Dated photo manifest and original source sample counts drifted")
            if n_source!=n_dated+n_unknown:
                raise ValueError("Original date eligibility accounting incomplete")
            if n_precise>n_dated:
                raise ValueError("Original precise georeferences exceed dated photos")

            monthfile=out/"all_original_1000_source_photo_observation_months.csv"
            month_counts=defaultdict(set)
            with monthfile.open("x",encoding="utf-8",newline="") as file:
                writer=csv.DictWriter(file,fieldnames=MONTHS,lineterminator="\n")
                writer.writeheader()
                sql="""SELECT p.taxon,p.ym,p.hemisphere,p.n_photos,
                       COUNT(o.observation),COUNT(DISTINCT o.observer)
                       FROM months p LEFT JOIN observations o
                       ON p.taxon=o.taxon AND p.ym=o.ym
                       AND p.hemisphere=o.hemisphere
                       GROUP BY p.taxon,p.ym,p.hemisphere
                       ORDER BY p.taxon,p.ym,p.hemisphere"""
                for tax,ym,hem,photos,obs,observers in db.execute(sql):
                    if not (0<obs<=photos and 0<observers<=obs):
                        raise ValueError("Original month independent observation counts failed")
                    month_counts[tax].add(ym)
                    writer.writerow(dict(zip(MONTHS,
                        (tax,names[tax],ym,hem,photos,obs,observers,True))))

            sitefile=out/"all_original_1000_source_month_0p25degree_5km_accuracy.csv"
            supported=defaultdict(set)
            with sitefile.open("x",encoding="utf-8",newline="") as file:
                writer=csv.DictWriter(file,fieldnames=SITES,lineterminator="\n")
                writer.writeheader()
                sql="""SELECT p.taxon,p.ym,p.ilat,p.ilon,p.n_photos,
                       COUNT(o.observation),COUNT(DISTINCT o.observer)
                       FROM sites p LEFT JOIN observations o
                       ON p.taxon=o.taxon AND p.ym=o.ym AND p.ilat=o.ilat
                       AND p.ilon=o.ilon AND o.precise=1
                       GROUP BY p.taxon,p.ym,p.ilat,p.ilon
                       ORDER BY p.taxon,p.ym,p.ilat,p.ilon"""
                for tax,ym,ilat,ilon,photos,obs,observers in db.execute(sql):
                    if not (0<obs<=photos and 0<observers<=obs):
                        raise ValueError("Original local source photo/observation depth mismatch")
                    if obs>=2 and observers>=2:
                        supported[(tax,ym[:4],ilat,ilon)].add(ym)
                    writer.writerow(dict(zip(SITES,
                        (tax,names[tax],ym,ilat,ilon,photos,obs,observers,True))))
            supported_siteyears=Counter()
            for (tax,year,ilat,ilon),months in supported.items():
                if len(months)>=2:
                    supported_siteyears[tax]+=1
            taxfile=out/"all_original_1000_taxa_date_and_site_sampling_coverage.csv"
            with taxfile.open("x",newline="",encoding="utf-8") as file:
                writer=csv.DictWriter(file,fieldnames=TAXA,lineterminator="\n")
                writer.writeheader()
                for tax in sorted(names,key=int):
                    writer.writerow(dict(zip(TAXA,(
                        tax,names[tax],source_counts[tax],dated_counts[tax],
                        source_counts[tax]-dated_counts[tax],len(month_counts[tax]),
                        site_counts[tax],supported_siteyears[tax],True))))
            if sum(source_counts.values())!=n_source or sum(dated_counts.values())!=n_dated:
                raise ValueError("Whole-source month coverage totals fail")

        result={
            "schema":"fcp_exact_original_rgfca_1000_all_source_observed_months_v1",
            "status":("SYNTHETIC_ONLY_NOT_BIOLOGICAL_PHOTOS" if _small_test else
                      "COMPLETE_ORIGINAL_SOURCE_CALENDAR_CENSUS_NO_FLOWER_PIXELS"),
            "snapshot_date":SOURCE_DATE,
            "original_monthly_manifest_sha256":digest(manifest),
            "original_monthly_source_receipt_sha256":digest(source_receipt),
            "n_taxa":len(names),
            "n_original_photo_ids":n_source,
            "n_photo_ids_with_full_observed_on_date":n_dated,
            "n_photo_ids_without_valid_observed_on_date":n_unknown,
            "n_dated_source_photos_5km_location_accuracy":n_precise,
            "n_taxa_with_original_dated_photo":sum(v>0 for v in dated_counts.values()),
            "n_taxa_with_repeated_source_calendar_months":sum(len(v)>=2 for v in month_counts.values()),
            "n_taxa_with_supported_0p25_degree_multiple_month_same_year":sum(v>0 for v in supported_siteyears.values()),
            "n_supported_geocell_years_two_months_two_obs_two_observers_each":sum(supported_siteyears.values()),
            "original_source_months_file_sha256":digest(monthfile),
            "original_source_site_month_file_sha256":digest(sitefile),
            "original_taxa_date_coverage_sha256":digest(taxfile),
            "source_photo_opportunity_is_not_flowering_detection":True,
            "source_date_not_image_file_timestamp_or_upload_date":True,
            "geo_precision_without_source_accuracy_is_not_fabricated":True,
            "all_original_image_pixels_classified":False,
            "new_flower_colour_or_white_frequency_estimated":False,
            "not_42111_source_photo_census_or_genetic_polymorphism":True,
        }
        (out/"original_dated_source_opportunity_receipt.json").write_text(
            json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        complete=True
        return result
    finally:
        db_path.unlink(missing_ok=True)
        if not complete:
            for filename in (
                "all_original_1000_source_photo_observation_months.csv",
                "all_original_1000_source_month_0p25degree_5km_accuracy.csv",
                "all_original_1000_taxa_date_and_site_sampling_coverage.csv",
                "original_dated_source_opportunity_receipt.json",
            ):
                (out/filename).unlink(missing_ok=True)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-monthly-photo-ids",type=Path,required=True)
    p.add_argument("--original-monthly-source-receipt",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args(argv)
    print(json.dumps(execute(args.original_monthly_photo_ids,
                             args.original_monthly_source_receipt,args.outdir),
                     sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
