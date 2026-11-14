#!/usr/bin/env python3
"""Original-image RGB/ROI-v4 pilot on licensed GBIF NON-iNaturalist media.

A *distinct* observational source from MOIS' frozen 93.2M iNaturalist frame.
No pixels are downloaded by 'source'. 'work' loads ONLY approved HTTPS
public image hosts, uses original pinned FCP YOLO11/EfficientSAM ROI-v4 and
seals per-image SHA; the old MOIS 42,111 atlas is not modified. 'audit'
checks actual terminals, hue fractions, URL/image SHA collisions.

Non-iNat URLs are NOT evidence of non-reposted images: cross-frame image-SHA
comparison and botanical flower visibility QA are REQUIRED BEFORE an
external candidate can count as a new strictly covered MOIS source taxon.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, HTTPRedirectHandler, build_opener
from urllib.error import HTTPError, URLError

import numpy as np
import pandas as pd

from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import (
    safe_image, BASELINE, require, source_sha
)
from scripts.acquisition.mois_zero9_gbif_other_cc_media_union import (
    SCHEMA as UNION_SCHEMA
)
from scripts.analysis.mois_species_first_palette9 import PALETTE, class_labels

SCHEMA="mois_zero9_gbif_other_CC_original_photo_RGB9_pilot_v1"
COHORT="mois_zero9_gbif_other_CC_photo_source_v1"
MAX_BYTES=12*1024*1024


def validate_https_public_host(url,allowed_hosts):
    parsed=urlsplit(str(url))
    host=(parsed.hostname or "").lower()
    require(parsed.scheme=="https" and parsed.port in (None,443)
            and parsed.username is None and parsed.password is None and
            host in allowed_hosts and host not in ("localhost","127.0.0.1")
            and not parsed.fragment and parsed.path and
            "inaturalist" not in host and not host.endswith(".internal"),
            "Refuse external image URL except source-audited exact public HTTPS host")
    return host


def approved_media_rows(union_root, allowed_hosts, quota=32):
    """Conservative one photo / one independently sourced GBIF occurrence / taxon."""
    union_root=Path(union_root)
    require(1<=quota<=96 and bool(allowed_hosts),
            "Bound original external-photo cohort and require manual host allowlist")
    receipt=json.loads((union_root/"receipt.json").read_text())
    require(receipt.get("schema")==UNION_SCHEMA and
            receipt.get("source_original_42111_taxon_photo_atlas_sha256")==BASELINE and
            receipt.get("stratified_hash_draw_taxa")==304 and
            receipt.get("photo_pixels_downloaded")==0 and
            receipt.get("new_strict9_numerical_photo_measurements")==0,
            "External image source has no signed 304-taxon source-only certificate")
    file=union_root/"mois_CC_noninat_GBIF_unique_image_metadata_304taxa.csv.gz"
    require(source_sha(file)==receipt["outputs_sha256"][file.name],
            "External original image URL source digest changed")
    a=pd.read_csv(file,compression="gzip",dtype=str,keep_default_na=False)
    needed={"inat_taxon_id","source_binomial","gbif_COL_XR_accepted_species_key",
            "gbif_occurrence_key","gbif_dataset_key","image_url","image_license"}
    require(needed<=set(a.columns),"Missing GBIF accepted taxon or image source identity")
    require(a.image_url.is_unique,"Duplicate external media URL was not quarantined")
    selection=[]
    for row in a.to_dict("records"):
        if not safe_image(row["image_url"],row["image_license"]):
            continue
        host=(urlsplit(row["image_url"]).hostname or "").lower()
        if host not in allowed_hosts:continue
        validate_https_public_host(row["image_url"],allowed_hosts)
        if not row["gbif_occurrence_key"] or not row["gbif_dataset_key"]:
            continue
        row["media_host"]=host
        row["fixed_draw"]=hashlib.sha256((
            "MOIS-EXTERNAL-GBIF-ROI9-FROZEN-1|"+
            row["inat_taxon_id"]+"|"+row["image_url"]).encode()).hexdigest()
        selection.append(row)
    a=pd.DataFrame(selection)
    require(len(a)>0,"No GBIF licensed photo URL on an explicitly approved host")
    a=a.sort_values(["fixed_draw"],kind="stable")
    a=a.drop_duplicates(subset=["gbif_occurrence_key"],keep="first")
    a=a.drop_duplicates(subset=["inat_taxon_id"],keep="first").head(quota)
    require(len(a)>0 and a.inat_taxon_id.is_unique and
            a.gbif_occurrence_key.is_unique and a.image_url.is_unique,
            "Original externally licensed taxon/photo/observation source collisions")
    a["measurement_id"]=a.image_url.map(lambda v:
        "MOIS-CCGBIF-ROI9-"+hashlib.sha256(v.encode()).hexdigest()[:24])
    a["image_filename"]=a.measurement_id+".jpg"
    a["photo_url_large"]=a.image_url
    a["photo_license"]=a.image_license
    a["observation_id"]="GBIF-occurrence:"+a.gbif_occurrence_key
    a["photo_id"]="GBIF-CC-media:"+a.image_url.map(lambda v:
        hashlib.sha256(v.encode()).hexdigest()[:24])
    a["shard"]="0"
    return a.drop(columns=["fixed_draw"])


def source(union_root,output,*,quota=32,allowed_hosts=()):
    output=Path(output)
    require(not output.exists(),"Do not overwrite a signed external source-photo manifest")
    allowed_hosts=set(x.lower() for x in allowed_hosts)
    rows=approved_media_rows(union_root,allowed_hosts,quota=quota)
    (output/"shards").mkdir(parents=True)
    p=output/"shards"/"shard-000.csv"
    rows.to_csv(p,index=False,lineterminator="\n")
    m={
        "schema":COHORT,
        "upstream_304_taxon_union_receipt_sha256":source_sha(Path(union_root)/"receipt.json"),
        "source_real_original_photo_urls_SHA256":source_sha(p),
        "source_historical_original_MOIS_photo_taxon_atlas_SHA256":BASELINE,
        "source_taxon_count":len(rows),
        "original_external_image_count_planned":len(rows),
        "source_manual_approved_public_image_hosts":sorted(allowed_hosts),
        "source_unique_GBIF_occurrence_per_species":True,
        "photo_colour_measured":False,
        "botanical_flower_visibility_checked":False,
        "original_image_sha_overlap_with_MOIS_1231244_prior_pixels_checked":False,
        "downstream_appending_to_iNaturalist_original_frame_allowed":False,
        "raw_image_pixels_persisted":False,
    }
    (output/"manifest.json").write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
    print("MOIS_OTHER_CCGBIF_SOURCE_FROZEN_NOT_COLOUR_MEASURED",
          json.dumps({k:m[k] for k in (
            "source_taxon_count","original_external_image_count_planned",
            "source_manual_approved_public_image_hosts")},sort_keys=True),flush=True)
    return m


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):
        return None


def load_photo(url,allowed_hosts):
    """Bounded, no-redirect, explicit-public-host network loader for one licensed photo."""
    from PIL import Image,UnidentifiedImageError
    from io import BytesIO
    from scripts.acquisition.run_fcp_fresh_image_partition import (
        PermanentOriginalPhotoError, ProviderTransientAccessError)
    validate_https_public_host(url,allowed_hosts)
    opener=build_opener(NoRedirect())
    req=Request(url,headers={"User-Agent":
        "MOIS-ecology-CCGBIF-ROI9-research/1.0 (source-only no public image reuse)",
        "Accept":"image/jpeg,image/png,image/webp"})
    last=None
    for i in range(3):
        try:
            with opener.open(req,timeout=32) as response:
                typ=str(response.headers.get("Content-Type","")).split(";")[0].strip().lower()
                require(typ in ("image/jpeg","image/png","image/webp","image/jpg"),
                        "Remote host returned non-image/unlicensed bytes")
                buf=response.read(MAX_BYTES+1)
            if not buf or len(buf)>MAX_BYTES:
                raise PermanentOriginalPhotoError("CC source photo absent/exceeds 12 MiB hard cap")
            try:
                with Image.open(BytesIO(buf)) as im:
                    if im.format not in ("JPEG","PNG","WEBP"):
                        raise PermanentOriginalPhotoError("Unsupported original format")
                    if im.width*im.height>60_000_000:
                        raise PermanentOriginalPhotoError("Decode-bomb pixel dimensions")
                    im.verify()
            except (OSError,UnidentifiedImageError,Image.DecompressionBombError) as exc:
                raise PermanentOriginalPhotoError("Undecodable CC source original image") from exc
            time.sleep(0.75)  # polite per-provider original still-image acquisition
            return buf
        except HTTPError as exc:
            last=exc
            if exc.code in (404,410,301,302,303,307,308):
                raise PermanentOriginalPhotoError("Remote original photo no longer a direct image") from exc
            if exc.code not in (408,425,429,500,502,503,504):
                raise ProviderTransientAccessError("Public licensed media HTTP denied "+str(exc.code)) from exc
        except (URLError,TimeoutError) as exc:
            last=exc
        if i<2:time.sleep(1+i*2)
    raise ProviderTransientAccessError("Public licensed image provider temporary failure: "+str(last))


def work(source_dir,runtime,sam,output):
    from scripts.acquisition.run_fcp_fresh_image_partition import run_batch
    from scripts.acquisition.mois_full_index_real_original_photo_measure import PhotoBatch,HUES
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
    source_dir=Path(source_dir);output=Path(output)
    require(not output.exists(),"Original CCGBIF photo cohort must never be silently repeated")
    manifest=json.loads((source_dir/"manifest.json").read_text())
    require(manifest.get("schema")==COHORT
            and manifest.get("photo_colour_measured") is False
            and manifest.get("downstream_appending_to_iNaturalist_original_frame_allowed") is False,
            "Not an authenticated original CCGBIF image cohort")
    p=source_dir/"shards"/"shard-000.csv"
    require(source_sha(p)==manifest["source_real_original_photo_urls_SHA256"],
            "CCGBIF photo URL changed after frozen ID cohort")
    rows=pd.read_csv(p,dtype=str,keep_default_na=False)
    allowed=set(manifest["source_manual_approved_public_image_hosts"])
    require(1<=len(rows)<=96 and rows.photo_id.is_unique and
            rows.observation_id.is_unique and rows.inat_taxon_id.is_unique and
            rows.image_url.is_unique,
            "Original external source image/taxon/GBIF occurrence duplicated")
    for row in rows.itertuples(index=False):
        require(safe_image(row.image_url,row.image_license),"CC licence or external media host invalid")
        validate_https_public_host(row.image_url,allowed)
    output.mkdir()
    sealed=output/"sealed"
    runtime=Path(runtime);sam=Path(sam)
    fingerprint=run_batch(
        Path("."),0,0,len(rows),sealed,
        source=runtime/"scripts/analysis/measure_random_photo_first_blind_partition.py",
        detector=runtime/"data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt",
        sam=sam,
        loader=lambda u:load_photo(u,allowed),
        indexed_source=PhotoBatch(rows,0))
    require(fingerprint["status"]=="COMPLETE_WORK_UNIT"
            and fingerprint["n_source_photo_ids"]==len(rows)
            and fingerprint["n_terminal_photo_ids"]==len(rows)
            and fingerprint["method"]==METHOD
            and fingerprint["image_bytes_persisted"] is False,
            "Original frozen FCP ROI-v4 model or external image cohort changed")
    terminal=sealed/"terminal_p000_b0000000.csv"
    require(source_sha(terminal)==fingerprint["terminal_csv_sha256"],
            "Authenticated real external image terminal SHA mismatch")
    measured=pd.read_csv(terminal,dtype=str,keep_default_na=False)
    require(measured.measurement_id.is_unique and len(measured)==len(rows),
            "External source image terminal missing one or more outcomes")
    frame=rows.merge(measured,on="measurement_id",how="inner",validate="one_to_one",sort=False)
    required=["morph","measurement_status","image_sha256"]
    require(all(k in frame.columns for k in required),"Original image SHA missing from ROI4 terminal")
    frame.to_csv(output/"external_original_photo9_and_unknowns.csv.gz",index=False,compression="gzip")
    report={
        "schema":SCHEMA,"source_CCGBIF_manifest_SHA256":source_sha(source_dir/"manifest.json"),
        "source_original_media_shard_SHA256":source_sha(p),
        "source_frozen_FCP_work_unit_receipt_SHA256":source_sha(sealed/"terminal_p000_b0000000.json"),
        "actual_external_original_photo_attempts":len(rows),
        "FCP_frozen_method_verified":True,
        "actual_external_original_image_acquisition_failures":fingerprint["n_acquisition_failed"],
        "actual_FCP_fourstate_classifiable":fingerprint["n_classifiable"],
        "actual_external_image_sha256_present":int(frame.image_sha256.ne("").sum()),
        "actual_external_duplicate_image_SHA256":int(frame.loc[
            frame.image_sha256.ne(""),"image_sha256"].duplicated(keep=False).sum()),
        "original_external_pixels_persisted":False,
        "source_iNaturalist_93m_image_sha_overlap_checked":False,
        "botanical_flower_organ_true_visible_checked":False,
        "external_hue_ready_to_append_as_new_source_taxa":False,
        "GBIF_location_specific_flower_colour_proven":False,
        "output_numeric_photo_terminal_SHA256":source_sha(output/"external_original_photo9_and_unknowns.csv.gz"),
    }
    (output/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_NONINAT_CCGBIF_ORIGINAL_FCP_ROI9_WORKER",
          json.dumps(report,sort_keys=True),flush=True)
    return report


def audit(source_dir,measured_dir,output):
    output=Path(output);source_dir=Path(source_dir);measured_dir=Path(measured_dir)
    require(not output.exists(),"Cannot overwrite independent original external ROI9 audit")
    manifest=json.loads((source_dir/"manifest.json").read_text())
    worker=json.loads((measured_dir/"receipt.json").read_text())
    require(worker.get("schema")==SCHEMA and worker.get("FCP_frozen_method_verified")
            and worker.get("source_CCGBIF_manifest_SHA256")==source_sha(source_dir/"manifest.json"),
            "Actual photo worker and signed source identity changed")
    file=measured_dir/"external_original_photo9_and_unknowns.csv.gz"
    require(source_sha(file)==worker["output_numeric_photo_terminal_SHA256"],
            "Actual original CCGBIF photo colour terminal altered after worker")
    d=pd.read_csv(file,dtype=str,keep_default_na=False)
    assert len(d)==worker["actual_external_original_photo_attempts"]
    require(d.inat_taxon_id.is_unique and d.observation_id.is_unique and
            d.photo_id.is_unique and d.measurement_id.is_unique,
            "External media original photo/GBIF occurrence/taxon overlap")
    cols=["flower_fraction_"+x for x in PALETTE]
    x=d[cols].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    strict=d.morph.isin(("white","yellow_orange","red_pink","blue_purple")).to_numpy(bool)
    strict=(strict & np.isfinite(x).all(axis=1) &
            (x>=0).all(axis=1) & (x<=1.00001).all(axis=1) &
            (np.abs(np.nan_to_num(x).sum(axis=1)-1)<=0.001))
    labels=np.full(len(d),"",dtype=object)
    labels[strict]=class_labels(x[strict])
    require(int(strict.sum())==worker["actual_FCP_fourstate_classifiable"],
            "True numerical external original FCP9 ROI does not match model receipt")
    raw=d.loc[strict & d.image_sha256.ne(""),"image_sha256"]
    duplicates=int(raw.duplicated(keep=False).sum())
    result={
        "schema":"mois_zero9_noninat_CCGBIF_external_RGB9_independent_audit_v1",
        "actual_external_original_photo_attempts":len(d),
        "actual_valid_numeric_nine_palette_vectors":int(strict.sum()),
        "actual_clear_photo_nine_palette_labels":int((labels!="").sum()),
        "actual_technical_unresolved_or_failed":int((~strict).sum()),
        "unique_external_source_taxa_with_numeric_photo9":int(d.loc[strict,"inat_taxon_id"].nunique()),
        "external_original_image_sha_duplicates":duplicates,
        "external_photo_SHA_overlap_prior_iNaturalist_frame_certified":False,
        "botanical_flower_only_pixel_true_verified":False,
        "allowed_to_add_external_hues_to_MOIS_iNaturalist_frame":False,
        "original_worker_terminal_SHA256":source_sha(file),
        "GBIF_local_occurrence_hue_measured":False,
    }
    d["strict9_external_original_photo_numeric"]=strict
    d["dominant_hue9_photo_label"]=labels
    for c in cols:
        d.loc[~strict,c]=""
    output.mkdir()
    out=output/"external_CCGBIF_raw_photo9_audited.csv.gz"
    d.to_csv(out,index=False,compression="gzip")
    result["audited_nine_fraction_photo_outcome_SHA256"]=source_sha(out)
    (output/"receipt.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("MOIS_EXTERNAL_NONINAT_CCGBIF_ORIGINAL_RGB9_AUDIT",
          json.dumps(result,sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="mode",required=True)
    src=sub.add_parser("source")
    src.add_argument("--union",required=True)
    src.add_argument("--output",required=True)
    src.add_argument("--quota",type=int,default=32)
    src.add_argument("--allow-host",action="append",required=True)
    workp=sub.add_parser("work")
    workp.add_argument("--source",required=True)
    workp.add_argument("--runtime",required=True)
    workp.add_argument("--sam",required=True)
    workp.add_argument("--output",required=True)
    auditp=sub.add_parser("audit")
    auditp.add_argument("--source",required=True)
    auditp.add_argument("--measured",required=True)
    auditp.add_argument("--output",required=True)
    a=p.parse_args()
    if a.mode=="source": source(a.union,a.output,quota=a.quota,allowed_hosts=a.allow_host)
    elif a.mode=="work": work(a.source,a.runtime,a.sam,a.output)
    else: audit(a.source,a.measured,a.output)
