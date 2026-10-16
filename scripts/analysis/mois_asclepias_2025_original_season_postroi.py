#!/usr/bin/env python3
"""MOIS post-blind 2025 seasonal ROI audit and exact TerraClimate month-cell keys.

FCP is sole authority for 2026-09-27 original photo IDs/measurements. This
module only verifies the source terminal and coordinates after original blind
ROI-v4; it does not establish anthesis or estimate a seasonal white effect.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,zipfile
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path

ORIG_RUN=38908308446
PIXEL_RUN=38909164068
ROI_SHA="1b16cdabd2936567dee87b8344019f40f315600c6510dfc1f25a03835024b6a0"
SNAPSHOT="2026-09-27"
STATUS="FIRST_ORIGINAL_SEASONAL_PIXELS_MEASURED_BOTANICAL_ANTHESIS_NOT_VALIDATED"
MORPHS={"white","yellow_orange","red_pink","blue_purple"}
CLIMATE_FIELDS=("year_month","tc_lat_idx","tc_lon_idx","tc_lat_center","tc_lon_center",
                "n_photo_observations","n_classifiable_roi","n_white_roi_only",
                "n_chromatic_roi_only","n_mixed_floral_organs","n_unclassifiable_roi",
                "anthesis_is_unverified")
MONTH_FIELDS=("observed_year_month","local_season","n_original_photo_attempts","n_independent_observers",
              "n_automated_four_state","n_mixed_white_yellow_roi","n_unusable_roi",
              "n_photo_white_unverified","n_photo_chromatic_unverified",
              "n_botanically_confirmed_open_flower")
PHOTO_FIELDS=("inat_taxon_id","species","photo_id","measurement_id","observation_id",
              "observer_id","observed_on","observed_year_month","local_season",
              "latitude","longitude","positional_accuracy_m","tc_lat_idx","tc_lon_idx",
              "roi_measurement_status","roi_morph_unverified","same_roi_organ_mix",
              "botanical_open_flower_review","no_white_rate_until_open_flower_confirmation")


def sha(raw:bytes)->str:return hashlib.sha256(raw).hexdigest()


def grid(lat:float,lon:float)->tuple[int,int,str,str]:
    if not math.isfinite(lat) or not math.isfinite(lon) or not (-90<=lat<=90 and -180<=lon<=180):
        raise ValueError("Original date source coordinates outside sphere")
    latid=min(4319,int(math.floor((90-lat)*24)))
    lonid=int(math.floor(((lon+180)%360)*24))
    return latid,lonid,f"{90-(latid+.5)/24:.8f}",f"{-180+(lonid+.5)/24:.8f}"


def read_csv(raw:bytes)->list[dict]:
    with io.StringIO(raw.decode("utf-8")) as stream:
        reader=csv.DictReader(stream)
        if not reader.fieldnames:raise ValueError("Original CSV contains no header")
        return list(reader)


def write_csv(out:Path,fields:tuple[str,...],rows:list[dict])->None:
    with out.open("x",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        writer.writeheader();writer.writerows(rows)


def run(artifact:Path,out:Path,*,_synthetic=False)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Mois seasonal analysis cannot overwrite sealed original inputs")
    with zipfile.ZipFile(artifact) as z:
        source_doc=json.loads(z.read("signed_season_selection_origin_run.json"))
        sel=json.loads(z.read("source/frozen/original_season_photo_freeze_receipt.json"))
        m=json.loads(z.read("measure/terminal_p000_b0000000.json"))
        receipt=json.loads(z.read("postblind/first_original_source_matched_season_64_roi_receipt.json"))
        source_raw=z.read("postblind/first_original_source_matched_season_64_photo_roi_QC.csv")
        original_raw=z.read("measure/terminal_p000_b0000000.csv")
        blind=z.read("source/frozen/season_first_original_blind_worker.csv")
    if not _synthetic and (source_doc.get("id")!=ORIG_RUN or
          source_doc.get("head_sha")!="55a43469450d86145d2a19e3f70c67e633f7fddf" or
          source_doc.get("status")!="completed" or source_doc.get("conclusion")!="success"):
        raise ValueError("Wrong original FCP 128-photo season selection job")
    if (sel.get("source_photo_queue_sha256") is None or
        sel.get("status")!="REAL_ORIGINAL_SEASONAL_SOURCE_PHOTOS_FROZEN_NO_PIXELS" or
        sel.get("no_climate_or_colour_outcome_used_to_select_photos") is not True or
        sel.get("observer_ids_unique_across_seasons") is not True or
        sel.get("file_sha256",{}).get("season_first_original_blind_worker.csv")!=sha(blind) or
        m.get("status")!="COMPLETE_WORK_UNIT" or
        m.get("terminal_csv_sha256")!=sha(original_raw) or
        receipt.get("status")!=STATUS or
        receipt.get("terminal_sha256")!=sha(original_raw) or
        receipt.get("postblind_join_sha256")!=sha(source_raw) or
        receipt.get("original_photos_measured")!=len(read_csv(source_raw)) or
        receipt.get("no_genetic_polymorphism_or_flowering_colour_ecological_effect_estimated") is not True):
        raise ValueError("Original seasonal FCP source/pixel/method SHA or ecological boundary drift")
    if not _synthetic and (len(read_csv(source_raw))!=64 or sha(source_raw)!=ROI_SHA):
        raise ValueError("Not the exact original 64 September27 Asclepias photo ROI")
    original=read_csv(source_raw)
    original_receipt=receipt.get("first_work_unit_by_season",{})
    attempted=Counter();by_month=defaultdict(lambda:Counter());by_climate=defaultdict(lambda:Counter())
    photo_out=[];seen_photos=set();seen_obs=set();seen_observers=set()
    unique_grids=set();season_observers=defaultdict(set)
    for r in original:
        tid,species=r["inat_taxon_id"],r["species"]
        if not _synthetic and (tid!="47911" or species!="Asclepias syriaca"):
            raise ValueError("Original paired season source species differs")
        pid,measurement,obs,observer=r["photo_id"],r["measurement_id"],r["observation_id"],r["observer_id"]
        if (not pid.isascii() or not pid.isdecimal() or not measurement or not obs or not observer
            or pid in seen_photos or obs in seen_obs or observer in seen_observers):
            raise ValueError("Duplicate photo, observation or source observer across original matched seasons")
        seen_photos.add(pid);seen_obs.add(obs);seen_observers.add(observer)
        if r["source_snapshot_date"]!=SNAPSHOT:raise ValueError("Changed frozen FCP archive date")
        day=date.fromisoformat(r["observed_on"])
        ym=day.strftime("%Y-%m")
        if day.year!=2025 or ym!=r["observed_year_month"] or r["season_year"]!="2025":
            raise ValueError("Original summer/autumn season-year mismatched observation date")
        season=("summer" if day.month in (6,7,8) else "autumn" if day.month in (9,10,11) else None)
        if season!=r["local_season"] or season is None:
            raise ValueError("Non-2025 summer/autumn photo admitted")
        lat,lon=float(r["latitude"]),float(r["longitude"])
        acc=float(r["positional_accuracy_m"])
        if not math.isfinite(acc) or not 0<acc<=5000:raise ValueError("Mois climate cannot invent precise source geolocation")
        ilat,ilon,clat,clon=grid(lat,lon)
        if str(int(math.floor((lat+90)*4)))!=r["site_cell_lat_0p25"] or str(int(math.floor(((lon+180)%360)*4)))!=r["site_cell_lon_0p25"]:
            raise ValueError("Original temporal photo not same frozen 0.25-degree site")
        good=r["measurement_status"]=="classified_four_state_morph"
        morph=r["morph"]
        if good and morph not in MORPHS:raise ValueError("Four-state ROI classified without valid morph")
        if not good and morph!="mixed_uncertain":raise ValueError("ROI failure relabelled as flower colour")
        mix=r["same_roi_white_yellow_organ_mixed"]=="True"
        if r["same_roi_white_yellow_organ_mixed"] not in ("True","False") or mix and not good:
            raise ValueError("Mixed-organ flag inconsistent with original ROI stage")
        if good:
            white=float(r["flower_fraction_white"]);yellow=float(r["flower_fraction_yellow"])
            if not (0<=white<=1 and 0<=yellow<=1 and mix==(white>=.10 and yellow>=.10)):
                raise ValueError("Original palette composition is not preserved")
        if r["flower_stage_verified_by_independent_botanists"]!="False" or r["flower_roi_organ_human_review"]!="UNREVIEWED":
            raise ValueError("Mois cannot promote unreviewed original ROI to flowering-stage evidence")
        attempted[season]+=1;season_observers[season].add(observer)
        group=by_month[(ym,season)];climate=by_climate[(ym,ilat,ilon,clat,clon)]
        for g in (group,climate):
            g["photos"]+=1
            g["classifiable"]+=int(good)
            g["mixed"]+=int(mix)
            g["unusable"]+=int(not good)
            g["white"]+=int(good and not mix and morph=="white")
            g["chromatic"]+=int(good and not mix and morph!="white")
        unique_grids.add((ilat,ilon))
        photo_out.append(dict(zip(PHOTO_FIELDS,(tid,species,pid,measurement,obs,observer,r["observed_on"],
             ym,season,r["latitude"],r["longitude"],r["positional_accuracy_m"],ilat,ilon,
             r["measurement_status"],morph,mix,"UNREVIEWED",True))))
    if len(season_observers["summer"]&season_observers["autumn"]):
        raise ValueError("Original matched-seasons observers not disjoint")
    if not _synthetic and (attempted!={"summer":32,"autumn":32} or
                         sum(g["classifiable"] for g in by_month.values())!=18):
        raise ValueError("Original real 32/32 measured seasonal ROI receipt not reproduced")
    for season,count in attempted.items():
        expected=original_receipt[season]
        sample={k:v for k,v in by_month.items() if k[1]==season}
        countfour=sum(c["classifiable"] for c in sample.values())
        white=sum(c["white"] for c in sample.values())
        chrom=sum(c["chromatic"] for c in sample.values())
        mixed=sum(c["mixed"] for c in sample.values())
        if (count!=expected["attempted_original_images"] or
            countfour!=expected["automated_four_state_classifiable"] or
            white!=expected["classifiable_white_photo_labels"] or
            chrom!=expected["classifiable_chromatic_photo_labels"] or
            mixed!=expected["same_flower_white_yellow_organ_mix"] or
            expected["true_open_flower_photo_count_verified"]!=0):
            raise ValueError("Real original FCP by-season automated photo status drift")
    out.mkdir(parents=True,exist_ok=True)
    photos=out/"original_photo_observed_date_and_automated_colour_UNREVIEWED.csv"
    months=out/"observed_month_botanical_and_roi_missingness.csv"
    request=out/"source_month_terraclimate_v11_point_grid_keys.csv"
    write_csv(photos,PHOTO_FIELDS,photo_out)
    monthrows=[]
    for (ym,season),c in sorted(by_month.items()):
        monthrows.append(dict(zip(MONTH_FIELDS,(ym,season,c["photos"],c["photos"],c["classifiable"],
                           c["mixed"],c["unusable"],c["white"],c["chromatic"],0))))
    write_csv(months,MONTH_FIELDS,monthrows)
    reqrows=[]
    for (ym,ilat,ilon,clat,clon),c in sorted(by_climate.items()):
        reqrows.append(dict(zip(CLIMATE_FIELDS,(ym,ilat,ilon,clat,clon,c["photos"],
                     c["classifiable"],c["white"],c["chromatic"],c["mixed"],c["unusable"],True))))
    write_csv(request,CLIMATE_FIELDS,reqrows)
    result={
        "schema":"mois_actual_original_2025_asclepias_postblind_seasonal_roi_to_climate_grid_v1",
        "status":("SYNTHETIC_NOT_REAL_FLOWER_COLOUR" if _synthetic else
                   "REAL_DATE_AND_ROI_ATTEMPTS_NO_BOTANICAL_WHITE_CLIMATE_EFFECT"),
        "source_snapshot_date":SNAPSHOT,"source_fcp_run":ORIG_RUN,"first_pixels_fcp_run":PIXEL_RUN,
        "source_original_postblind_sha256":sha(source_raw),
        "source_method_terminal_sha256":sha(original_raw),
        "n_original_photo_attempts":len(original),
        "n_original_unique_observers":len(seen_observers),
        "n_original_source_grid_cells":len(unique_grids),
        "n_exact_2025_grid_and_month_climate_requests":len(reqrows),
        "n_2025_original_observed_months":len(monthrows),
        "by_season":{
            season:{"original_photos":attempted[season],
                    "automated_roi_four_state":sum(c["classifiable"] for (ym,s),c in by_month.items() if s==season),
                    "original_photo_white_labels_unreviewed":sum(c["white"] for (ym,s),c in by_month.items() if s==season),
                    "original_photo_chromatic_labels_unreviewed":sum(c["chromatic"] for (ym,s),c in by_month.items() if s==season),
                    "mixed_white_yellow_organs":sum(c["mixed"] for (ym,s),c in by_month.items() if s==season)}
            for season in ("summer","autumn")},
        "source_photographs_are_not_all_flowers":True,
        "open_flower_botanically_validated":False,
        "no_original_source_pixels_measured_by_mois":True,
        "no_climate_values_extracted_or_effect_estimated":True,
        "original_month_climate_grid_csv_sha256":sha(request.read_bytes()),
        "month_roi_opportunity_qc_sha256":sha(months.read_bytes()),
        "original_photo_attempts_sha256":sha(photos.read_bytes()),
    }
    (out/"mois_original_asclepias_source_rois_seasons_receipt.json").write_text(
       json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fcp-season64-artifact",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    print(json.dumps(run(args.fcp_season64_artifact,args.outdir),sort_keys=True))


if __name__=="__main__":main()
