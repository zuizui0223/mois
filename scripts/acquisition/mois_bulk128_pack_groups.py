#!/usr/bin/env python3
"""Pack already signed 128 photo groups into 16 tiny download packets.

This is DATA ONLY and never downloads original images. Each of the 128 warm
workers retrieves just one of 16 per-eight-part group tarballs instead of
re-downloading the whole 128-shard job's tens of MB provenance package.
"""
from __future__ import annotations

import argparse,gzip,hashlib,io,json,tarfile
from pathlib import Path

SCHEMA="mois_bulk128_compact16_read_only_group_bundles_v1"

def require(ok,msg):
    if not ok:raise ValueError("BULK128_PACK: "+msg)

def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()

def sha(path):
    return sha_bytes(Path(path).read_bytes())

def pack(root,output):
    root,output=Path(root),Path(output)
    require(not output.exists(),"Refuse overwrite of previous source-bundle export")
    manifest_path=root/"manifest.json"
    m=json.loads(manifest_path.read_text())
    require(m.get("sampling_schema")==
            "mois_original93m_tier10_bulk128_probability_reference_plus_gap_v1"
            and m.get("source_128_compute_shards")==128
            and m.get("work_units")==512
            and m.get("group_ids")==list(range(128))
            and m.get("P5_targeted_photo_attempts")==0
            and m.get("original_source_photo_IDs_and_observation_UUIDs_all_unique"),
            "Not the SHA-verified true 128-photo-source manifest")
    output.mkdir(parents=True)
    bundles={}
    for group in range(16):
        files=[]
        for part in range(group*8,(group+1)*8):
            rp=root/f"part-{part:03d}.receipt.json"
            require(sha(rp)==m["per_128_shard_receipts_sha256"][str(part)],
                    "Per-original-photo-shard source receipt changed")
            receipt=json.loads(rp.read_text())
            require(receipt["part"]==part
                    and receipt["schema"]==m["sampling_schema"]
                    and receipt["P5_targeted_original_photo_count"]==0,
                    "P5 or unrecognized source shard identity found")
            files.append((f"part-{part:03d}.receipt.json",rp))
            p=root/"provenance"/f"part-{part:03d}.csv"
            require(sha(p)==receipt["source_photo_provenance_SHA256"],
                    "Original photo/taxon/season provenance SHA mismatch")
            files.append((f"provenance/part-{part:03d}.csv",p))
            for unit in range(part*4,(part+1)*4):
                f=root/"shards"/f"shard-{unit:03d}.csv"
                require(sha(f)==m["photo_shard_sha256"][str(unit)]
                        ==receipt["source_photo_worker_unit_SHA256"][str(unit)],
                        "Original signed source worker photo-ID roster altered")
                files.append((f"shards/shard-{unit:03d}.csv",f))
        output_path=output/f"bulk128-source-group-{group:02d}.tar.gz"
        with output_path.open("wb") as raw:
            with gzip.GzipFile(filename="",mode="wb",fileobj=raw,
                               compresslevel=4,mtime=0) as zipped:
                with tarfile.open(fileobj=zipped,mode="w|") as tar:
                    for name,p in sorted(files):
                        data=p.read_bytes()
                        info=tarfile.TarInfo(name)
                        info.size=len(data)
                        info.mtime=0
                        info.uid=info.gid=0
                        info.uname=info.gname=""
                        info.mode=0o644
                        tar.addfile(info,io.BytesIO(data))
        bundles[output_path.name]=sha(output_path)
    report={
      "schema":SCHEMA,
      "original_source_manifest_SHA256":sha(manifest_path),
      "original_128_shards_packed_exactly_once":True,
      "group_count":len(bundles),
      "original_work_units":512,
      "bundle_SHA256":bundles,
      "new_original_photo_attempt_targets":m["new_photo_attempts_planned"],
      "no_original_raw_image_bytes_acquired":True,
    }
    (output/"bundle_receipt.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_BULK128_SOURCE_PACK16_FOR_EFFICIENT_PARALLEL_IMAGE_READS",
          json.dumps({"groups":len(bundles),
                      "new_photo_IDs":m["new_photo_attempts_planned"],
                      "sum_bytes":sum((output/x).stat().st_size for x in bundles)},
                     sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    opts=p.parse_args()
    pack(opts.source,opts.output)

if __name__=="__main__":
    main()
