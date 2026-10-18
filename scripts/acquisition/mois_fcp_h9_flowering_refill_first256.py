#!/usr/bin/env python3
"""MOIS photo acquisition speed: FCP H9 one random flowering page per FAILED taxon.

The FCP random-photo-first method uses iNaturalist phenology term_id=12,
term_value_id=13, research-grade georeferenced unobscured observations,
position <=5000m, eligible photo licences, page 1/per_page <=200, API
pacing and bounded retry. DO NOT copy FCP's ~100 photos/species design.
Here query ONLY broad global MOIS pilot taxa whose earlier genuine RGB
ROI-v4 result was unclassifiable. Keep newly classifiable 69 taxa untouched.
No new photo is assigned a hue until a SECOND actual image ROI run succeeds.
"""
from __future__ import annotations

import argparse,csv,hashlib,json,math,re
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

import pandas as pd

from scripts.acquisition import mois_global_angiosperm_breadth_source as broad
from scripts.analysis import mois_broad_first256_roi_taxonomic_gain_audit as pilot

FCP_CLIENT_BLOB="a46341c98aa7ada422ed0e2799ccb863f01a3626"
FCP_QUERY_RESEARCH="research"
FCP_FLOWER_TERM_ID=12
FCP_FLOWER_TERM_VALUE=13
FCP_MAX_POSITIONAL_ACCURACY=5000
FCP_PAGE_SIZE=200
FCP_ALLOWED_PHOTO_LICENSES=("cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa")
SALT="MOIS_FIRST256_UNKNOWN_H9_FLOWERING_FCP_API_FIRST_PAGE_20261018_V1"
PILOT_ROI="mois_broad_first256_real_photo_ROI_hue_taxon_quality.csv"
PILOT_RECEIPT="mois_broad_first256_photo_RBG_novel_taxa_genera_QC_receipt.json"
FROZEN_SOURCE="mois_fcp_H9_original_flowering_photo_refills_first256_unknown.csv"
FROZEN_RECEIPT="mois_fcp_H9_flowering_refills_photo_source_provenance.json"
AUDIT="mois_fcp_H9_original_187_unclassifiable_photo_species_source_outcomes.csv"
MIN_FAILURES=187
OUT_COLS=(
    "measurement_id","photo_id","inat_taxon_id","species","observation_id",
    "observer_id","observed_on","latitude","longitude","photo_license",
    "photo_url_large","image_filename","original_failed_photo_id",
    "original_failed_measurement_status","source_sampling_arm","source_annotation",
    "source_protocol","api_source_taxon_id","api_observation_id",
    "fcp_client_blob_sha","source_order","source_original_observation_uuid",
    "source_snapshot_date"
)


def require(ok,msg):
    if not ok:raise ValueError(msg)


def sha(path:Path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def gitsource_blob(path:Path):
    body=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(body)).encode()+b"\0"+body).hexdigest()


def fcp_h9_params(taxon_id):
    """Identical FCP H9 flower-grade photographic eligibility, taxon-specific."""
    require(str(taxon_id).isdecimal(),"H9 rescue taxon must be original numeric iNat ID")
    return {
        "taxon_id":int(taxon_id),
        "quality_grade":FCP_QUERY_RESEARCH,
        "photos":"true",
        "geo":"true",
        "rank":"species",
        "term_id":FCP_FLOWER_TERM_ID,
        "term_value_id":FCP_FLOWER_TERM_VALUE,
        "acc_below":FCP_MAX_POSITIONAL_ACCURACY,
        "obscuration":"none",
        "photo_license":",".join(FCP_ALLOWED_PHOTO_LICENSES),
        "order_by":"random",
        "per_page":FCP_PAGE_SIZE,
        "page":1,
    }


def geographic_position(item):
    geo=item.get("geojson") or {}
    p=geo.get("coordinates") or []
    if not isinstance(p,(list,tuple)) or len(p)!=2:return None
    try:
        lng=float(p[0]);lat=float(p[1])
        if math.isfinite(lng) and math.isfinite(lat) and -180<=lng<=180 and -90<=lat<=90:
            return lat,lng
    except (ValueError,TypeError): pass
    return None


def candidate_from_response(item,source,known_photo_ids):
    """Local independent guards for FCP API preselection, without hue lookup."""
    taxon=item.get("taxon") or {}
    user=item.get("user") or {}
    if (str(taxon.get("id") or "")!=str(source["inat_taxon_id"])
        or str(taxon.get("rank") or "")!="species"
        or str(taxon.get("name") or "").strip().casefold()!=
            str(source["species"]).strip().casefold()
        or str(item.get("quality_grade") or "")!="research"
        or bool(item.get("obscured"))
        or item.get("geoprivacy") not in (None,"","open")
        or not user.get("id") or not item.get("id")):
        return None
    place=geographic_position(item)
    if place is None:return None
    try:
        positional=float(item["positional_accuracy"])
        observed=str(item["observed_on"])
        actualday=date.fromisoformat(observed)
        if not (0<=positional<=FCP_MAX_POSITIONAL_ACCURACY
                and date(1800,1,1)<=actualday<=date.today()):
            return None
    except (KeyError,TypeError,ValueError,OverflowError):return None
    if str(item.get("uuid") or "")==str(source.get("observation_id") or ""):
        return None
    # The observation's flowering annotation applies to the observation
    # and may not guarantee that each attached PHOTO shows an open petal.
    options=[]
    for position,p in enumerate(item.get("photos") or []):
        if not isinstance(p,dict):continue
        pid=str(p.get("id") or "")
        license=str(p.get("license_code") or "").lower()
        u=str(p.get("url") or "")
        if not pid.isdecimal() or pid in known_photo_ids or (
            license not in FCP_ALLOWED_PHOTO_LICENSES):
            continue
        extension=Path(urlsplit(u).path).suffix.lower().removeprefix(".")
        if extension not in ("jpg","jpeg","png"):
            continue
        options.append((position,pid,license,extension))
    if not options:return None
    cover=min(options,key=lambda o:(o[0],int(o[1])))
    _,pid,license,ext=cover
    return {
        "photo_id":pid,
        "inat_taxon_id":str(source["inat_taxon_id"]),
        "species":str(source["species"]),
        "observation_id":str(item["id"]),
        "source_original_observation_uuid":str(item.get("uuid") or ""),
        "observer_id":str(user["id"]),
        "observed_on":observed,
        "latitude":f"{place[0]:.7f}",
        "longitude":f"{place[1]:.7f}",
        "photo_license":license,
        "photo_url_large":
           f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}",
        "original_failed_photo_id":str(source["photo_id_source"]),
        "original_failed_measurement_status":str(source["measurement_status"]),
        "source_sampling_arm":"fcp_h9_flowering_annotated_rescue_after_unknown",
        "source_annotation":"flowers_and_fruits:flowers_(12:13)",
        "source_protocol":"FCP_H9_one_random_API_page_per_failed_angiosperm_species",
        "api_source_taxon_id":str(source["inat_taxon_id"]),
        "api_observation_id":str(item["id"]),
        "fcp_client_blob_sha":FCP_CLIENT_BLOB,
    }


def choose_one_h9_photo(rows,source,excluded_photo_ids):
    """API one-page 200 observations; deterministic outcome-blind pick."""
    eligible=[]
    for o in rows:
        x=candidate_from_response(o,source,excluded_photo_ids)
        if x is not None:
            rank=hashlib.sha256(
                (SALT+"|"+str(source["inat_taxon_id"])+"|"+
                 x["photo_id"]).encode()).hexdigest()
            eligible.append((rank,x))
    if not eligible:return None
    return min(eligible,key=lambda a:a[0])[1]


def load_pilot(prior,pilot_release,wide):
    existing,_,_,_,_,_ = broad.old_photo_source(prior)
    p=pd.read_csv(wide/broad.FIRST,dtype=str,keep_default_na=False)
    rec=json.loads((wide/broad.RECEIPT).read_text())
    require(len(p)==4096 and p.photo_id.is_unique
            and sha(wide/broad.FIRST)==rec["first_wave_original_4096_csv_sha256"]
            and rec["first_wave_named_taxa"]==4096
            and len(existing)==313535,
            "Global species-first 4096 photo source not pinned")
    quality=json.loads((pilot_release/PILOT_RECEIPT).read_text())
    real=pilot_release/PILOT_ROI
    require(real.is_file() and quality["source_outcomes_csv_SHA256"]==sha(real)
            and quality["actual_pilot_metrics"]["new_photo_hue_classifications"]==69,
            "First 256 original RGB ROI pilot source SHA/colour count changed")
    tested=pd.read_csv(real,dtype=str,keep_default_na=False)
    require(len(tested)==256 and tested.measurement_id.is_unique
            and tested.photo_id_source.is_unique
            and set(tested.measurement_id)==set(p.loc[p.shard.eq("0"),"measurement_id"]),
            "Pilot colour results do not map to fixed 256 photo IDs")
    original=pd.DataFrame(tested)
    nclasses=int(tested.classified_four_hue.eq("true").sum())
    require(nclasses==69 and
            int(tested.morph.isin(("white","yellow_orange","red_pink","blue_purple")).sum())==69,
            "Unknown vs photo hue reconciliation changed")
    unknown=tested.loc[~tested.classified_four_hue.eq("true")].copy()
    require(len(unknown)==MIN_FAILURES
            and unknown.inat_taxon_id.nunique()==MIN_FAILURES,
            "Pilot failed ROI named angiosperm species denominator changed")
    excluded=set(existing)|set(p.photo_id)
    return unknown,excluded,rec,quality


def run(previous,oldpilot,broadsrc,fifth,season,output,client=None,
        *,maximum_species=None):
    require(not output.exists() or not any(output.iterdir()),
            "Do not overwrite an original photo-source rescue snapshot")
    unknown,excluded,freezerec,oldrec=load_pilot(previous,oldpilot,broadsrc)
    fids,frec=broad.reservation(fifth,broad.EXCLUDED_OLD,
        broad.EXCLUDED_OLD_RECEIPT,"first_wave_original_4096_csv_sha256",4096)
    sids,srec=broad.reservation(season,broad.EXCLUDED_SI,
        broad.EXCLUDED_SI_RECEIPT,"first_wave_original_4096_csv_sha256",4096)
    require(fids.isdisjoint(sids) and
            not (fids & excluded) and not (sids & excluded),
            "Fifth/season photo reservation overlaps original broad main source")
    excluded.update(fids)
    excluded.update(sids)
    if maximum_species is not None:
        require(int(maximum_species)>=1,
                "Pilot maximum rescue species must be positive")
        unknown=unknown.iloc[:int(maximum_species)].copy()
    # The original 32768 candidate pool is also reserved for new taxa
    # acquisition and should not be silently removed from that future pool.
    pool=broadsrc/broad.POOL
    require(sha(pool)==freezerec["original_broad_all_future_source_pool_SHA256"],
            "Large future source pool changed since original-source freeze")
    candidates=pd.read_csv(pool,usecols=["photo_id"],dtype=str,
                           keep_default_na=False)
    excluded.update(candidates.photo_id)
    statuses=Counter();selected=[];audits=[];seen=set()
    if client is None:
        from fcp_pipeline.random_photo_pool import InaturalistObservationClient
        client=InaturalistObservationClient(
            user_agent="mois-h9-flowering-original-photo-refill/1.0 "
                       "(github.com/zuizui0223/mois)",
            request_interval_seconds=1.05,max_retries=4,timeout_seconds=45)
    for i,row in enumerate(unknown.to_dict("records")):
        tid=str(row["inat_taxon_id"])
        params=fcp_h9_params(tid)
        try:
            response=client.observations(params)
            results=response.get("results")
            require(isinstance(results,list),
                    "iNaturalist response has no list of original observation results")
            candidate=choose_one_h9_photo(results,row,excluded|seen)
            status="new_flowering_annotated_photo_source" if candidate else (
                "no_extra_photo_in_annotated_first_page")
            results_n=len(results)
        except Exception as exc:
            candidate=None
            status="provider_error"
            results_n=0
        statuses[status]+=1
        audits.append({
            "source_taxon_id":tid,
            "source_failed_photo_id":row["photo_id_source"],
            "source_failed_ROI_status":row["measurement_status"],
            "fcp_H9_api_random_page":1,
            "photo_candidates_returned":results_n,
            "source_rescue_state":status,
            "new_original_photo_id":candidate["photo_id"] if candidate else "",
            "phenology_term":"12:13",
            "new_visible_hue_measured":False,
        })
        if candidate:
            pid=candidate["photo_id"];seen.add(pid)
            candidate["source_order"]=str(i)
            mid=hashlib.sha256((SALT+"|"+pid).encode()).hexdigest()[:24].upper()
            candidate["measurement_id"]="MOIS-"+mid
            candidate["image_filename"]="MOIS-"+mid+".jpg"
            candidate["source_snapshot_date"]=freezerec["source_snapshot_date"]
            selected.append(candidate)
        if (i+1)%20==0:
            print("MOIS_FCP_H9_FLOWER_REFILL_SOURCE_PROGRESS",json.dumps(
                {"attempted_failed_source_species":i+1,
                 "new_original_annotated_photo_ids":len(selected),
                 "no_annotated_alternative":statuses["no_extra_photo_in_annotated_first_page"],
                 "provider_errors":statuses["provider_error"]}),flush=True)
    attempts=len(unknown)
    require(attempts==sum(statuses.values())
            and statuses["provider_error"]<=max(2,int(.05*attempts)),
            "iNaturalist API unavailable too often; do not interpret provider failures as missing flowers")
    require(len(selected)==len({x["photo_id"] for x in selected})
            and len(selected)==len({x["inat_taxon_id"] for x in selected}),
            "New source flower-annotated photographic result duplicates IDs/taxa")
    output.mkdir(parents=True,exist_ok=True)
    source=output/FROZEN_SOURCE
    with source.open("x",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=OUT_COLS,lineterminator="\n")
        writer.writeheader();writer.writerows(selected)
    audits_path=output/AUDIT
    pd.DataFrame(audits).to_csv(audits_path,index=False,lineterminator="\n")
    report={
        "schema":"mois_fcp_H9_flowering_metadata_rescue_original_unclassifiable256_pilot_v1",
        "source_method":"FCP_H9_frozen_random_photo_first_api_one_page_per_species",
        "FCP_H9_query_UTC":datetime.now(timezone.utc).isoformat(),
        "source_fifth4096_photo_reservation_SHA256":
            frec["first_wave_original_4096_csv_sha256"],
        "source_season4096_photo_reservation_SHA256":
            srec["first_wave_original_4096_csv_sha256"],
        "FCP_client_code_git_blob_sha":FCP_CLIENT_BLOB,
        "request_params_template":fcp_h9_params("123"),
        "original_256_pilot_completed_RGB_classifications":69,
        "original_256_pilot_unknown_hue_photos":MIN_FAILURES,
        "rescued_species_API_attempts":attempts,
        "rescued_flowering_annotated_new_photo_SOURCE_ONLY_IDs":len(selected),
        "source_rescue_status_counts":dict(statuses),
        "source_original_313535_and_broad_pool_photo_IDs_excluded":True,
        "source_original_256_successful_taxa_not_requeried":True,
        "FCP_100_photo_per_species_requirement_not_applied":True,
        "iNaturalist_flowering_annotation_is_OBSERVATION_not_per_photo_certainty":True,
        "new_photo_pixels_actually_measured":False,
        "new_original_photo_four_hues_classified":0,
        "original_photo_pre_colour_source_selection_blind_to_hue":True,
        "rescue_sampling_stratum_not_unbiased_global_white_fraction":True,
        "photo_colour_unknown_is_not_nonwhite":True,
        "separate_future_and_season_photo_IDs_reserved":True,
        "source_new_photo_candidates_SHA256":sha(source),
        "source_query_audit_SHA256":sha(audits_path),
        "old_first256_RGB_source_SHA256":oldrec["source_outcomes_csv_SHA256"],
        "old_broad_original_4096_photo_source_SHA256":
            freezerec["first_wave_original_4096_csv_sha256"],
        "no_FCP_within_species_colour_polymorphism_50km_estimator":True,
        "no_biochemical_UV_anthocyanin_pigment_cost_or_flower_fitness":True,
    }
    (output/FROZEN_RECEIPT).write_text(
        json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("MOIS_FCP_H9_FIRST256_UNKNOWN_SPECIES_FLOWERING_ANNOTATED_NEW_SOURCE_ONLY",
          json.dumps({"queried_failed_original_ROI_species":attempts,
                      "new_photo_ID_candidates":len(selected),
                      "provider_errors":statuses["provider_error"]},
                     sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original313535-photo-release",type=Path,required=True)
    p.add_argument("--original256-hue-result-release",type=Path,required=True)
    p.add_argument("--broad4096-source-release",type=Path,required=True)
    p.add_argument("--fifth4096-photo-reservation",type=Path,required=True)
    p.add_argument("--season4096-photo-reservation",type=Path,required=True)
    p.add_argument("--max-query-species",type=int)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.original313535_photo_release,a.original256_hue_result_release,
        a.broad4096_source_release,a.fifth4096_photo_reservation,
        a.season4096_photo_reservation,a.output,
        maximum_species=a.max_query_species)


if __name__=="__main__":
    main()
