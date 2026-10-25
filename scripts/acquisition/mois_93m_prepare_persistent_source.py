#!/usr/bin/env python3
"""Materialize actual signed full93m photo-ID TAR groups into bounded scan roots.

Reads each of the 224 immutable tier×8part metadata archives ONE time.
Never fetches or decodes original photograph pixels. All original photo IDs,
observations and taxa remain in exactly one 4,096-photo source chunk. Result
is directly consumable by mois_93m_persistent_worker.py across many groups.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile
import tempfile

from scripts.acquisition.mois_full93m_source_part_chunks import run as split_one, sha
from scripts.acquisition.mois_onepass_full93m_target_tiers import RESULT_SCHEMA

SCHEMA="mois_all93m_source_group_once_to_persistent_bounded_chunk_roots_v1"
REGISTRY_SCHEMA="mois_complete_full93m_original_photo_target_registry_v1"

def require(ok, message):
    if not ok:
        raise ValueError("MOIS_PERSISTENT_SOURCE_PREP: "+message)

def prepare_group(registry, group_receipt, archive, output, *,
                  tier,group,chunk_size=4096,
                  expected_world_photos=93_209_351,parts_per_group=8):
    registry,group_receipt,archive,output=map(
        Path,(registry,group_receipt,archive,output))
    require(not output.exists(),"Existing authenticated source group must not be rewritten")
    require(0<=tier<14 and 0<=group<128//parts_per_group and
            1<=parts_per_group<=128 and 1<=chunk_size<=8192,
            "Malformed source tier or group split")
    reg=json.loads(registry.read_text(encoding="utf-8"))
    proof=json.loads(group_receipt.read_text(encoding="utf-8"))
    require(reg.get("schema")==REGISTRY_SCHEMA
            and reg.get("exact_original_photo_IDs_in_all_disjoint_tiers")
                ==expected_world_photos
            and reg.get("native_vs_introduced_not_hard_filtered") is True
            and reg.get("source_only_NO_photo_pixels_or_hue_measured") is True,
            "Not the complete signed global photo-ID registry")
    require(proof.get("schema")==RESULT_SCHEMA
            and proof.get("source_snapshot_SHA256")
                ==reg["full_source_manifest_sha256"],
            "Original 93m source group and complete registry differ")
    require(archive.name==f"tier-{tier:02d}-group-{group:02d}.tar"
            and reg["published_tier_archive_SHA256"][archive.name]==sha(archive),
            "Published source tier archive SHA is absent or mismatched")
    parts={int(p["part"]):p for p in proof["parts"]}
    target=set(range(group*parts_per_group,(group+1)*parts_per_group))
    require(set(parts)==target,
            "Missing original photo-ID source partition in archival group")
    expected=sum(int(p["counts_by_disjoint_tier"][str(tier)])
                 for p in parts.values())
    require(expected==int(proof["counts_by_disjoint_tier"][str(tier)]),
            "Original photo count not identical to group receipt")
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-93m-source-chunks-",dir=output.parent))
    staging_tar=stage/"one_source_part.csv.gz"
    roots=[]
    try:
        with tarfile.open(archive,mode="r:") as pack:
            members={m.name:m for m in pack.getmembers() if m.isfile()}
            names={f"tier-{tier:02d}/part-{part:03d}.csv.gz" for part in target}
            require(set(members)==names,
                    "Source archive includes unknown or missing metadata photo-ID parts")
            total=0
            for part_no in sorted(parts):
                one=parts[part_no]
                expected_part=int(one["counts_by_disjoint_tier"][str(tier)])
                with pack.extractfile(members[f"tier-{tier:02d}/part-{part_no:03d}.csv.gz"]) as inp:
                    with staging_tar.open("wb") as out:
                        shutil.copyfileobj(inp,out,length=4<<20)
                require(sha(staging_tar)==one["ranked_target_csv_sha256"][str(tier)],
                        "Original tier/part metadata gzip SHA mismatch")
                root=stage/f"tier-{tier:02d}"/f"part-{part_no:03d}"
                if expected_part:
                    result=split_one(
                        staging_tar,root,
                        expected_source_sha256=one["ranked_target_csv_sha256"][str(tier)],
                        expected_photo_count=expected_part,tier=tier,part=part_no,
                        chunk_size=chunk_size)
                    require(result["actual_original_photo_IDs"]==expected_part,
                            "Photo IDs lost inside source partition chunking")
                    roots.append((part_no,expected_part,str(
                        (output/f"tier-{tier:02d}"/f"part-{part_no:03d}").resolve())))
                else:
                    require(not root.exists(),
                            "A source part with no IDs must not produce fabricated pixels")
                total+=expected_part
            staging_tar.unlink()
        require(total==expected,"Original source SHA count drift within group")
        roots_file=stage/"source_roots.txt"
        roots_file.write_text("".join(name+"\n" for _,_,name in roots),
                              encoding="utf-8")
        result={
            "schema":SCHEMA,
            "source_world_original_photo_IDs":expected_world_photos,
            "signed_registry_sha256":sha(registry),
            "signed_group_receipt_sha256":sha(group_receipt),
            "signed_archive_sha256":sha(archive),
            "original_source_tier":tier,
            "source_archive_group":group,
            "signed_original_photo_parts":parts_per_group,
            "nonempty_photo_parts":len(roots),
            "actual_signed_source_photo_IDs_accounted":total,
            "chunk_max_original_photo_IDs":chunk_size,
            "photo_root_paths_sha256":sha(roots_file),
            "original_source_image_downloads":0,
            "actual_image_pixels_inferred":0,
            "native_taxon_season_photo_exclusions":0,
            "source_93m_full_pixels_completed":False,
        }
        (stage/"source_prep_receipt.json").write_text(
            json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(stage,output)
        return result
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise

def collect_roots(group_directories, output, *,
                  world_total=93_209_351, require_world_complete=True):
    output=Path(output)
    require(not output.exists(),"Cannot overwrite old complete source-root roster")
    seen=set()
    total=0
    names=[]
    for group_root in sorted(map(Path,group_directories)):
        proof=json.loads((group_root/"source_prep_receipt.json").read_text())
        require(proof.get("schema")==SCHEMA and
                proof.get("source_world_original_photo_IDs")==world_total
                and proof.get("actual_image_pixels_inferred")==0,
                "Unverified source metadata group")
        root_file=group_root/"source_roots.txt"
        require(sha(root_file)==proof["photo_root_paths_sha256"],
                "Source chunk root paths changed")
        for name in root_file.read_text(encoding="utf-8").splitlines():
            folder=Path(name)
            receipt=json.loads((folder/"chunk_receipt.json").read_text())
            identity=(receipt["original_tier"],receipt["original_part"])
            require(identity not in seen,"Duplicate photo-ID source tier/part")
            seen.add(identity)
            names.append(str(folder.resolve()))
            total+=receipt["actual_original_photo_IDs"]
        # Include zero-photo partitions in the 1792 source coverage check.
    if require_world_complete:
        # The prepared-source group receipts prove all 8 original part sources,
        # including technical zero-photo partitions; each tier has 16 groups.
        groups={
            (json.loads((Path(root)/"source_prep_receipt.json").read_text())[
                "original_source_tier"],
             json.loads((Path(root)/"source_prep_receipt.json").read_text())[
                "source_archive_group"])
            for root in group_directories
        }
        require(len(groups)==14*16 and total==world_total,
                "Complete worldwide scan roots require all 224 source archives")
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text("".join(n+"\n" for n in sorted(names)),encoding="utf-8")
    return {
        "world_source_photo_IDs":world_total,
        "prepared_nonempty_signed_source_part_roots":len(seen),
        "signed_photo_IDs_in_root_roster":total,
        "all_world_original_photo_metadata_available":
            require_world_complete and total==world_total,
        "all_actual_medium_pixel_scans_completed":False,
        "roots_file_sha256":sha(output),
    }

def main():
    a=argparse.ArgumentParser(description=__doc__)
    sp=a.add_subparsers(dest="command",required=True)
    p=sp.add_parser("prepare-group")
    p.add_argument("--registry",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--archive",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--tier",type=int,required=True)
    p.add_argument("--group",type=int,required=True)
    p.add_argument("--chunk-size",type=int,default=4096)
    c=sp.add_parser("collect-roots")
    c.add_argument("--group-roots-file",type=Path,required=True)
    c.add_argument("--output",type=Path,required=True)
    c.add_argument("--partial",action="store_true")
    args=a.parse_args()
    if args.command=="prepare-group":
        res=prepare_group(args.registry,args.receipt,args.archive,args.output,
                          tier=args.tier,group=args.group,chunk_size=args.chunk_size)
    else:
        roots=[s.strip() for s in args.group_roots_file.read_text().splitlines()
               if s.strip()]
        res=collect_roots(roots,args.output,
                          require_world_complete=not args.partial)
    print("MOIS_93M_SIGNED_SOURCE_ONCE_ONLY_PREPARED",json.dumps(res,sort_keys=True))

if __name__=="__main__":
    main()
