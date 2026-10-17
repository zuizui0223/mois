#!/usr/bin/env python3
"""MOIS reuse FCP wave4's 49,800 genuinely original ROI-v4 image terminal rows.

Do not import FCP within-species D, W, 50km spatial/repeated-polymorphism
statistics; import original photo ID, old ROI hue label, site/date, and status
only. All source overlap handled by photo-ID dedup, and the prior MOIS 256504
data are never overwritten. This file deliberately refuses source-only FCP
metadata without all 256 finished pixel partition receipts.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json
from collections import Counter
from pathlib import Path

HUES=("white","yellow_orange","red_pink","blue_purple")
EXPECTED_FCP_PHOTOS=49800
EXPECTED_FCP_TAXA=498
EXPECTED_SOURCE_SHA="8d98b14856dd2539b54e1db287aa08adea9c191e112223ff380cb49ac0f8a554"
SOURCE_SCHEMA="fcp_same_149900_method_wave4_real_photo_terminals_v1"
BASE_COLS=("photo_id","inat_taxon_id","species","morph","source_original_measured",
 "old_mois_or_new_fcp","measurement_status","latitude","longitude",
 "photo_license","observation_id","observed_on","time_observed_at",
 "wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay",
 "soil_available_water_proxy","wc_bio1","wc_bio5","wc_bio12","wc_bio15")


def require(ok,msg):
    if not ok:raise ValueError(msg)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def read_fcp_wave4(result_root:Path):
    candidate=list(result_root.rglob("wave4_original_method_measured_photos.csv"))
    receipts=list(result_root.rglob("terminal_result.json"))
    require(len(candidate)==len(receipts)==1,
            "Original FCP complete wave4 photo pixels/results not available")
    original,receipt=candidate[0],json.loads(receipts[0].read_text())
    require(receipt.get("schema")==SOURCE_SCHEMA and
            receipt.get("wave4_new_photo_ids")==EXPECTED_FCP_PHOTOS and
            receipt.get("wave4_100photo_species")==EXPECTED_FCP_TAXA and
            receipt.get("terminal_partitions")==256 and
            receipt.get("source_metadata_sha256")==EXPECTED_SOURCE_SHA and
            receipt.get("measurement_original_fixed_colour_worker") is True and
            receipt.get("no_floral_colour_imputed") is True and
            receipt.get("original_results_changed") is False,
            "FCP wave4 actual image pixel/immutable source receipt invalid")
    require(receipt.get("measured_sha256")==sha(original),
            "FCP original 49,800 photo pixel table checksum invalid")
    with original.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        required={"photo_id","inat_taxon_id","species","observation_id",
                  "latitude","longitude","morph","measurement_status",
                  "global_classifiable","photo_license","observed_on"}
        require(required.issubset(set(reader.fieldnames or [])),
                "FCP old original photo metadata/hue fields lost")
        records=list(reader)
    photo_ids=[r["photo_id"] for r in records]
    require(len(records)==len(set(photo_ids))==EXPECTED_FCP_PHOTOS and
            len({r["inat_taxon_id"] for r in records})==EXPECTED_FCP_TAXA,
            "Original 49800 photos/498 taxa have missing/duplicated identifiers")
    classes=Counter()
    for r in records:
        status=r["measurement_status"]
        cls=r["morph"] in HUES
        require((status=="classified_four_state_morph")==cls and
                str(r["global_classifiable"]).strip().lower()==str(cls).lower(),
                "FCP actual original image classifier status/morph mismatch")
        require(r["photo_id"].isdigit() and r["inat_taxon_id"].isdigit() and
                r["photo_license"] and r["species"].strip(),
                "Original FCP picture/source identity missing")
        classes["four_hue_classified" if cls else "uncertain"]+=1
    require(classes["four_hue_classified"]==receipt["wave4_classifiable_photo_rows"] and
            classes["uncertain"]==receipt["wave4_nonclassifiable_photo_rows"],
            "Full original FCP photo-level four-hue totals disagree with source")
    return records,receipt


def merge(prior:Path,prior_receipt:Path,source:Path,outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Never overwrite published/released original photohue source")
    old=json.loads(prior_receipt.read_text())
    expected=old.get("all_original_unique_photos")
    classified=old.get("all_four_hue_classified_photos")
    require(expected==256504 and isinstance(classified,int) and
            sha(prior)==old.get("all_real_source_photo_hue_csv_sha256") and
            old.get("original_model_SHA_verified_for_every_photo") is True,
            "Prior 2×4096 original photo pixel source not verified")
    new,receipt=read_fcp_wave4(source)
    prior_ids=set()
    overlap=Counter()
    taxonomy_mismatch=Counter()
    outdir.mkdir(parents=True,exist_ok=True)
    output=outdir/"mois_original_photo_hues_plus_real_FCP_wave4_49800.csv.gz"
    count=0
    classified_count=0
    with output.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw,mode="wb",mtime=0,filename="") as gz:
            with io.TextIOWrapper(gz,encoding="utf-8",newline="") as f:
                w=csv.DictWriter(f,fieldnames=BASE_COLS,lineterminator="\n")
                w.writeheader()
                old_dict={}
                with gzip.open(prior,"rt",encoding="utf-8",newline="") as src:
                    rd=csv.DictReader(src)
                    require(tuple(rd.fieldnames or ())==BASE_COLS,
                            "Prior 256504 original source changed photographic fields")
                    for row in rd:
                        pid=row["photo_id"]
                        require(pid and pid not in prior_ids,
                                "Existing MOIS source has repeated original photographic IDs")
                        prior_ids.add(pid)
                        old_dict[pid]=(row["morph"],row["inat_taxon_id"])
                        w.writerow(row)
                        count+=1
                        classified_count+=row["morph"] in HUES
                require(count==expected and classified_count==classified,
                        "Prior already measured original photo/hue denominators changed")
                for row in new:
                    pid=row["photo_id"]
                    if pid in prior_ids:
                        overlap["source_id_already_measured"]+=1
                        old_hue,old_tax=old_dict[pid]
                        overlap["both_four_hue"]=overlap["both_four_hue"]+int(
                            old_hue in HUES and row["morph"] in HUES)
                        overlap["colour_mismatch"]=overlap["colour_mismatch"]+int(
                            old_hue in HUES and row["morph"] in HUES and old_hue!=row["morph"])
                        taxonomy_mismatch["taxon_id"]=taxonomy_mismatch["taxon_id"]+int(
                            old_tax!=row["inat_taxon_id"])
                        continue
                    output_row={k:"" for k in BASE_COLS}
                    for k in BASE_COLS:
                        if k in row and row[k] is not None:
                            output_row[k]=row[k]
                    output_row["source_original_measured"]="fcp_wave4_49800_original_photo_ROI_v4"
                    output_row["old_mois_or_new_fcp"]="fcp_additional_ROI"
                    w.writerow(output_row)
                    prior_ids.add(pid)
                    count+=1
                    classified_count+=row["morph"] in HUES
                    overlap["actually_new_original_photo_ids"]+=1
                    overlap["actually_new_classified_four_hue"]+=row["morph"] in HUES
    require(count==256504+overlap["actually_new_original_photo_ids"] and
            overlap["source_id_already_measured"]+
            overlap["actually_new_original_photo_ids"]==EXPECTED_FCP_PHOTOS,
            "Original FCP wave4 source photo ID accounting not reconciled")
    rec={
       "schema":"mois_reuse_complete_FCP_wave4_49800_original_photo_pixels_v1",
       "status":"ORIGINAL_IMAGE_PIXELS_CLASSIFIED_SAME_ROI_ALGORITHM_NOT_FCP_POLYMORPHISM",
       "original_mois_unique_photo_ids":expected,
       "original_mois_four_hue_classified":classified,
       "original_FCP_wave4_photo_ids":EXPECTED_FCP_PHOTOS,
       "original_FCP_wave4_named_taxa":EXPECTED_FCP_TAXA,
       "original_FCP_wave4_four_hue_classified":receipt["wave4_classifiable_photo_rows"],
       "original_FCP_wave4_photo_ids_overlapping_MOIS":overlap["source_id_already_measured"],
       "new_original_photo_ids":overlap["actually_new_original_photo_ids"],
       "new_four_hue_classified_original_photo_ids":overlap["actually_new_classified_four_hue"],
       "source_shared_photo_original_colour_disagreements":overlap["colour_mismatch"],
       "source_shared_photo_original_taxon_id_disagreements":taxonomy_mismatch["taxon_id"],
       "all_photo_ids_after_photo_id_deduplication":count,
       "all_photo_hues_classified_after_original_id_deduplication":classified_count,
       "source_FCP_wave4_terminal_result_SHA256":sha(next(iter(source.rglob("terminal_result.json")))),
       "source_FCP_wave4_full_measured_original_csv_sha256":receipt["measured_sha256"],
       "original_MOIS_256504_source_SHA256":sha(prior),
       "source_MOIS_photo_extension_SHA256":sha(output),
       "no_FCP_intraspecific_D_W_or_50km_polymorphism_result_imported":True,
       "data_are_photographic_colour_labels_not_genotypes_or_adaptation":True,
       "all_256_original_FCP_wave4_pixel_partitions_required":True,
       "same_iNaturalist_source_and_classifier_not_independent_replication":True
    }
    target=outdir/"mois_original_FCP_wave4_real_photo_reuse_receipt.json"
    target.write_text(json.dumps(rec,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_ACTUAL_REUSED_FCP_49800_ORIGINAL_MEASURED_PHOTO_LABELS",
          json.dumps({k:rec[k] for k in (
              "new_original_photo_ids",
              "new_four_hue_classified_original_photo_ids",
              "all_photo_ids_after_photo_id_deduplication",
              "all_photo_hues_classified_after_original_id_deduplication")},sort_keys=True))
    return rec


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-MOIS256504-photo-gzip",type=Path,required=True)
    p.add_argument("--source-MOIS256504-receipt",type=Path,required=True)
    p.add_argument("--FCP-wave4-complete-original-ROI-artifact",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    merge(a.source_MOIS256504_photo_gzip,a.source_MOIS256504_receipt,
          a.FCP_wave4_complete_original_ROI_artifact,a.output)
