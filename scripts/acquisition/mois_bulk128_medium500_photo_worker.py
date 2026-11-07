#!/usr/bin/env python3
"""Bulk128 real 500px original-photo image acquisition and flower-region colour.

One signed source compute shard (part=0..127) contains four frozen FCP
photo-ID work packets, each <=256, usually 1024 original observations/photos.
Convert their *metadata* to one signed original 93m row sequence; do NOT
rescan the 93m index or infer any 9-palette/petal-level colour. Fetch each
photo's 500px S3 variant once; run ONE warm source-pinned YOLO11 flower box
and preserve signed per-photo status, Lab box descriptors and failed results.

A selected P1/P2/P3/P4 cohort is biased by the previous 1m-photo ROI results.
Separate its sampling arm from the purely diagnostic source-reference cohort.
A box may include leaves and backgrounds, and no green/brown/unknown is
silently promoted to botanical petal colour.
"""
from __future__ import annotations

import argparse,csv,gzip,hashlib,io,json,math
from collections import Counter,defaultdict
from pathlib import Path

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS,obs_draw,usable_date,usable_point,
)
from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    COLUMNS_LAB,LAB_SCHEMA,WarmFlowerBoxDetector,Pacer,download_medium,
    scan_source,sha_file,
)
from scripts.acquisition.mois_flower_box_lab_colour import PROTOCOL

SCHEMA="mois_bulk128_real_medium500_photo_flower_box_lab_from_signed_photo_roster_v1"
SOURCE_SCHEMA="mois_original93m_tier10_bulk128_probability_reference_plus_gap_v1"
TIER=10
EXPECTED_PHOTOS=131_072
PARTS=128
PACKETS=4
REFERENCE="mois_bulk128_tier10_photo_reference_v1"
TARGETED="mois_bulk128_tier10_taxon_region_season_gap_v1"
SHARD_MANIFEST_SHA="7e975e701281ba812451acd48be7e91a7d17f2bbe6f9f53ad15ef4f8db9e557e"
SOURCE_ARCHIVE_SHA="61c0eee377f6edead471436d3f75b92700e54687e6b9dd0a417819a472c23dcb"
WEIGHT_SHA="f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5"
SELECTED_ROLE="bulk128_selected"
ALLOWED_STATUS={"candidate","negative","uncertain","error"}

def require(ok,message):
    if not bool(ok):
        raise ValueError("MOIS_BULK128_MEDIUM: "+message)

def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()

def source_photo_ext(row):
    pid=str(row["photo_id"])
    url=str(row["photo_url_large"])
    prefix="https://inaturalist-open-data.s3.amazonaws.com/photos/"+pid+"/large."
    require(url.startswith(prefix),"Original image URL must match frozen numeric photo ID")
    ext=url.removeprefix(prefix).lower()
    require(ext in ("jpg","jpeg","png"),"Original image file extension unsupported")
    return ext

def build_signed_input(source,part,out):
    """Return exact source 93m field rows plus reference/targeted provenance."""
    source,out=Path(source),Path(out)
    m=json.loads((source/"manifest.json").read_text())
    require(m.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and m.get("sampling_schema")==SOURCE_SCHEMA
            and m.get("source_128_compute_shards")==PARTS
            and m.get("new_photo_attempts_planned")==EXPECTED_PHOTOS
            and m.get("P5_targeted_photo_attempts")==0
            and m.get("selected_reference_photo_attempts")==32768
            and m.get("selected_targeted_photo_attempts")==98304
            and m.get("source_tier")==TIER
            and m.get("all_upstream_128_source_taxon_shards_verified"),
            "This is not the exact 128-part image-source cohort")
    require(isinstance(part,int) and 0<=part<PARTS,"Invalid bulk128 part")
    rp=source/f"part-{part:03d}.receipt.json"
    require(sha_file(rp)==m["per_128_shard_receipts_sha256"][str(part)],
            "Original photo metadata part receipt SHA changed")
    proof=json.loads(rp.read_text())
    provenance=source/"provenance"/f"part-{part:03d}.csv"
    require(sha_file(provenance)==proof["source_photo_provenance_SHA256"],
            "Original per-photo stage/stratum source provenance SHA changed")
    with provenance.open(newline="",encoding="utf-8") as fh:
        prev=list(csv.DictReader(fh))
    prov={r["photo_id"]:r for r in prev}
    require(len(prev)==len(prov)==proof["selected_new_original_photo_IDs"],
            "Original photo source provenance has absent or repeated photo ID")
    require(all(r["acquisition_design_cohort"] in ("reference","targeted")
                and r["previous_photo_ID_and_observation_overlap"]=="False"
                for r in prev),
            "Original reference and targeted cohorts have unknown identity/previous overlap")
    ordered=[]
    source_roles={}
    source_taxa={}
    source_observations=set()
    for unit in range(part*PACKETS,(part+1)*PACKETS):
        file=source/"shards"/f"shard-{unit:03d}.csv"
        require(sha_file(file)==m["photo_shard_sha256"][str(unit)]
                ==proof["source_photo_worker_unit_SHA256"][str(unit)],
                "Original selected image worker unit is not SHA-pinned")
        with file.open(newline="",encoding="utf-8") as f:
            partrows=list(csv.DictReader(f))
        require(0<len(partrows)<=256
                and all(x["shard"]==str(unit) for x in partrows),
                "Original 128-part work unit size or group identity invalid")
        for r in partrows:
            pid,oid,tax=r["photo_id"],r["observation_id"],r["inat_taxon_id"]
            require(pid.isdecimal() and tax.isdecimal() and oid
                    and pid in prov and
                    prov[pid]["observation_id"]==oid and
                    prov[pid]["inat_taxon_id"]==tax,
                    "Before-photo source taxon or observation mismatched original roster")
            require(pid not in source_roles and oid not in source_observations,
                    "Original photo or original observation duplicated in a source part")
            source_observations.add(oid)
            source_roles[pid]=prov[pid]["acquisition_design_cohort"]
            source_taxa[pid]=tax
            ext=source_photo_ext(r)
            z={x:"" for x in COLUMNS}
            z.update(
                photo_id=pid,observation_uuid=oid,inat_taxon_id=tax,
                observer_id=r["observer_id"],photographer_id="",
                observed_on=r["observed_on"],latitude=r["latitude"],
                longitude=r["longitude"],positional_accuracy=r["positional_accuracy"],
                photo_license=r["photo_license"],image_extension=ext,
                frozen_loader_supported="true",
                valid_observed_date=str(usable_date(r["observed_on"])).lower(),
                valid_geopoint=str(usable_point(r["latitude"],r["longitude"])).lower(),
                observation_inclusion_score=format(obs_draw(oid),".17g"),
                tier=str(TIER),
            )
            ordered.append(z)
    n=len(ordered)
    require(n==len(prev)==len(source_roles)==len(source_observations)
            and n==proof["selected_new_original_photo_IDs"]
            and Counter(source_roles.values())==Counter({
                "reference":proof["selected_photo_attempts_by_design"].get("reference",0),
                "targeted":proof["selected_photo_attempts_by_design"].get("targeted",0)}),
            "Source photo roster, observation UUID and design arm counts not conserved")
    dest=out/"selected-original-photo-source.csv.gz"
    out.mkdir(parents=True,exist_ok=True)
    temporary=out/"selected-original-photo-source.csv.gz.partial"
    require(not temporary.exists(),"Interrupted previous metadata selection; diagnose before continuing")
    with temporary.open("wb") as binary:
        with gzip.GzipFile(filename="",fileobj=binary,mode="wb",mtime=0,compresslevel=3) as zipped:
            with io.TextIOWrapper(zipped,encoding="utf-8",newline="") as writer:
                w=csv.DictWriter(writer,fieldnames=COLUMNS,lineterminator="\n")
                w.writeheader();w.writerows(ordered)
    if dest.exists():
        require(sha_file(dest)==sha_file(temporary),
                "Previously scanned P1-P4/diagnostic source original photo bytes changed on resume")
        temporary.unlink()
    else:
        temporary.replace(dest)
    return {
        "source_file":str(dest),
        "source_input_SHA256":sha_file(dest),
        "provenance_SHA256":sha_file(provenance),
        "part_receipt_SHA256":sha_file(rp),
        "original_photo_ID_count":n,
        "source_roles":source_roles,
        "source_taxa":source_taxa,
        "source_provenance_rows":prov,
    }

def reconcile_pixels(measured,roster,scan):
    measured=Path(measured)
    packets=sorted(measured.glob("packet-*.csv.gz"))
    require(len(packets)==scan["packets_sealed"],
            "Some signed new S3 photo packet or ROI proxy rows missing")
    roles=roster["source_roles"]
    source_taxa=roster["source_taxa"]
    base={name:Counter() for name in ("reference","targeted")}
    identifiers=set()
    joined=[]
    for packet in packets:
        rp=packet.with_suffix("").with_suffix(".receipt.json")
        require(rp.is_file(),"Original packet SHA receipt missing")
        receipt=json.loads(rp.read_text())
        require(receipt["schema"]==LAB_SCHEMA
                and receipt["source_role"]==SELECTED_ROLE
                and receipt["source_gzip_sha256"]==roster["source_input_SHA256"]
                and receipt["output_sha256"]==sha_file(packet),
                "Original 500px image packet was reclassified or changed")
        rows=0
        with gzip.open(packet,"rt",encoding="utf-8",newline="") as f:
            reader=csv.DictReader(f)
            require(tuple(reader.fieldnames or ())==COLUMNS_LAB,
                    "Detected photo-box Lab measurements changed column contract")
            for r in reader:
                pid=r["photo_id"]
                require(pid in roles and pid not in identifiers
                        and r["inat_taxon_id"]==source_taxa[pid]
                        and r["screen_status"] in ALLOWED_STATUS
                        and r["source_photo_inclusion_probability"]==""
                        and r["source_arm"]==
                        "bulk128_tier10_source_selected_reference_and_gap_SEPARATE",
                        "Original photo source identity/detection/estimator contradicted")
                identifiers.add(pid)
                role=roles[pid]
                base[role]["photo_attempts"]+=1
                base[role]["medium_decoded"]+=bool(r["image_pixel_sha256"])
                base[role]["flower_boxes"]+=r["screen_status"]=="candidate"
                base[role]["screen_negative"]+=r["screen_status"]=="negative"
                base[role]["technical_unknown"]+=r["screen_status"] in ("uncertain","error")
                base[role]["Lab_modal_descriptor"]+=bool(r["floral_box_main_L_star"])
                p=roster["source_provenance_rows"][pid]
                joined.append({
                    "photo_id":pid,"inat_taxon_id":r["inat_taxon_id"],
                    "observation_id":r["observation_id"],
                    "design_cohort":role,
                    "original_priority":p["original_source_priority"],
                    "latitude_region":p["latitude_region"],
                    "local_season":p["recorded_local_season"],
                    "screen_status":r["screen_status"],
                    "medium_decoded":bool(r["image_pixel_sha256"]),
                    "floral_box_main_L_star":r["floral_box_main_L_star"],
                    "floral_box_main_C_star":r["floral_box_main_C_star"],
                    "floral_box_main_rgb_hex":r["floral_box_main_rgb_hex"],
                    "unverified_floral_box_not_petal_ROI":True,
                })
                rows+=1
        require(rows==receipt["photo_ids"],
                "Original packet photo row count not conserved")
    require(identifiers==set(roles)
            and sum(v["photo_attempts"] for v in base.values())==
            roster["original_photo_ID_count"],
            "Full SHA-checked source-selected original photos not all attempted")
    return {k:dict(v) for k,v in base.items()},joined

def run(source,output,part,weights,*,model_sha=WEIGHT_SHA,rps=.5,
        workers=8,batch=256,detector_override=None,download_override=None):
    source,output,weights=Path(source),Path(output),Path(weights)
    require(not output.is_file(),
            "Original scanned pixel output must be a directory, not a file")
    require(0<rps<=1 and 1<=workers<=12 and batch==256,
            "Provider rate, original S3 decode workers or fixed source packet altered")
    output.mkdir(parents=True,exist_ok=True)
    signed=build_signed_input(source,part,output)
    measured=output/"medium500"
    detector=(
        WarmFlowerBoxDetector(weights,model_sha,0,"item",device="cpu",conf=.02)
        if detector_override is None else detector_override
    )
    pace=Pacer(rps)
    download=download_override or (lambda url:download_medium(url,pace))
    summary=scan_source(
        signed["source_file"],measured,
        expected_source_sha256=signed["source_input_SHA256"],
        expected_count=signed["original_photo_ID_count"],
        part=part,tier=TIER,source_role=SELECTED_ROLE,
        detector=detector,download=download,workers=workers,
        batch_size=batch,feature_mode="floral_box_lab_v2",
    )
    require(summary["complete_preselected_bulk128_photo_roster_not_full93m_census"]
            and summary["actual_unique_source_photo_ids_attempted"]==
            signed["original_photo_ID_count"]
            and summary["source_role"]==SELECTED_ROLE
            and summary["screening_only_not_verified_flower_colour"]
            and summary["outcome_conditioned_supplement_not_population_unbiased"]
            and not summary["all_original_input_photos_attempted"],
            "Genuine 500px preselected bulk originals were confused with full world census")
    roles,rows=reconcile_pixels(measured,signed,summary)
    roll=Counter(r["original_priority"] for r in rows if r["design_cohort"]=="targeted")
    require(set(roll).isdisjoint({"P5_existing_10plus_strict9","unknown"}),
            "Old well-measured P5 taxa entered missing-colour targeted original cohort")
    output_table=output/"photo_ID_matched_reference_vs_targeted_boxLab.csv.gz"
    temp_output=output/"photo_ID_matched_reference_vs_targeted_boxLab.csv.gz.partial"
    require(not temp_output.exists(),
            "Interrupted old photo-status aggregation must be diagnosed")
    with temp_output.open("wb") as binary:
        with gzip.GzipFile(filename="",fileobj=binary,mode="wb",mtime=0,
                           compresslevel=3) as zipped:
            with io.TextIOWrapper(zipped,encoding="utf-8",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
                w.writeheader();w.writerows(rows)
    if output_table.exists():
        require(sha_file(output_table)==sha_file(temp_output),
                "Existing original photo classifications would change on restart")
        temp_output.unlink()
    else:
        temp_output.replace(output_table)
    report={
        "schema":SCHEMA,"source_tier":TIER,"original_source_part":part,
        "expected_original_photos":signed["original_photo_ID_count"],
        "actually_medium_photo_IDs_attempted":summary["actual_unique_source_photo_ids_attempted"],
        "real_medium_photo_pixels_decoded":summary["medium_image_pixels_decoded"],
        "real_medium_photo_bytes_transferred":summary["medium_image_payload_bytes_downloaded"],
        "real_srgb_floral_box_not_verified_petals":True,
        "original_FCP_four_hue_full_image_ROI_inferred":False,
        "real_strict9_flower_pixel_vectors_generated":0,
        "photo_download_errors_preserved_as_unknown":True,
        "reference_and_adaptive_targeted_are_distinct":True,
        "by_design":roles,
        "new_targeted_original_priority":dict(roll),
        "source_128_photo_part_receipt_sha256":signed["part_receipt_SHA256"],
        "selected_original_photo_roster_SHA256":signed["source_input_SHA256"],
        "original_source_provenance_sha256":signed["provenance_SHA256"],
        "original_model_detector_weight_sha256":model_sha,
        "original_source_metadata_snapshot_SHA256":SHARD_MANIFEST_SHA,
        "actual_pixel_terminal_photo_and_cohort_SHA256":sha_file(output_table),
        "raw_licensed_image_bytes_archived":False,
        "stage1_selected_photo_roster_fully_attempted":True,
        "model_only_box_sRGB_CIELAB_not_petals":PROTOCOL,
        "full93m_source_original_photo_census_complete":False,
        "prior_93m_original_upstream_frame_intact":True,
    }
    final_receipt=output/"bulk128_medium_photo_completion_receipt.json"
    encoded=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n"
    if final_receipt.exists():
        require(final_receipt.read_text()==encoded,
                "Original 500px photo outcomes or source provenance changed after completion")
    else:
        final_receipt.write_text(encoded)
    print("MOIS_BULK128_TRUE_500PX_ORIGINAL_PHOTOS_STAGE1_FINISHED",
          json.dumps({k:report[k] for k in ("original_source_part",
          "expected_original_photos","actually_medium_photo_IDs_attempted",
          "real_medium_photo_pixels_decoded","by_design","new_targeted_original_priority")},
          sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--part",required=True,type=int)
    p.add_argument("--weights",required=True,type=Path)
    p.add_argument("--weights-sha256",default=WEIGHT_SHA)
    p.add_argument("--requests-per-second",type=float,default=.5)
    p.add_argument("--workers",type=int,default=8)
    a=p.parse_args()
    run(a.source,a.output,a.part,a.weights,model_sha=a.weights_sha256,
        rps=a.requests_per_second,workers=a.workers)

if __name__=="__main__":
    main()
