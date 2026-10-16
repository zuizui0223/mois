#!/usr/bin/env python3
"""Offline SHA-pinned extraction of EXACT WorldClim 2.1 10m precipitation ZIP.

Only operates on locally acquired original WorldClim bytes. It does not
download files, open original plant photos or measure floral colour.
"""
from __future__ import annotations
import argparse
import json
import re
import zipfile
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    from scripts.analysis.mois_monthly_rainfall_flower_colour_20261016 import (
        MANIFEST_SCHEMA, NORMALS, MonthlyNormals, sha256_file
    )
except ModuleNotFoundError:
    from mois_monthly_rainfall_flower_colour_20261016 import (
        MANIFEST_SCHEMA, NORMALS, MonthlyNormals, sha256_file
    )

ARCHIVE_BASENAME='wc2.1_10m_prec.zip'
ARCHIVE_SOURCE='https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_prec.zip'
PATH_PATTERN=re.compile(r'(?:wc2\.1_10m_prec/)?(wc2\.1_10m_prec_(\d\d)\.tif)$')
STATUS='TWELVE_WORLDCLIM_MONTHLY_RAINFALL_NORMALS_BYTE_PINNED_NO_BIOLOGICAL_RESULTS'
MAX_UNCOMPRESSED=2_000_000_000
MAX_PER_FILE=256_000_000


def seal(archive:Path, archive_sha256:str, output:Path, *, _synthetic_small=False)->dict:
    stage=output.with_name(output.name+'.INCOMPLETE')
    if output.exists() or stage.exists():
        raise ValueError('Refuse overwrite of a sealed/interrupted monthly raster source')
    if not _synthetic_small and archive.name!=ARCHIVE_BASENAME:
        raise ValueError('Official WorldClim monthly 10m precip ZIP basename required')
    if len(archive_sha256)!=64 or sha256_file(archive)!=archive_sha256:
        raise ValueError('Original monthly WorldClim provider ZIP SHA256 mismatch')
    with zipfile.ZipFile(archive) as z:
        matching={}
        total=0
        for item in z.infolist():
            if item.is_dir():
                continue
            pos=PATH_PATTERN.fullmatch(item.filename)
            if (pos is None or item.file_size<=0 or item.file_size>MAX_PER_FILE
                or item.flag_bits&1 or item.compress_size<=0):
                raise ValueError('Unexpected, potentially unsafe or non-precip original WorldClim ZIP member')
            month=int(pos.group(2))
            if month<1 or month>12 or month in matching:
                raise ValueError('Duplicate/out-of-range WorldClim monthly GeoTIFF member')
            matching[month]=(pos.group(1),item)
            total+=item.file_size
            if total>MAX_UNCOMPRESSED:
                raise ValueError('Original zip uncompressed raster budget exceeded')
        if set(matching)!=set(range(1,13)):
            raise ValueError('Provider monthly archive must contain exactly 12 original 10m precipitation rasters')
        bad=z.testzip()
        if bad is not None:
            raise ValueError('Original provider ZIP CRC failed on '+bad)
        stage.mkdir(parents=True)
        months={}
        try:
            for month in range(1,13):
                basename,member=matching[month]
                target=stage/basename
                with z.open(member) as inp,target.open('xb') as out:
                    while True:
                        block=inp.read(4*1024*1024)
                        if not block:break
                        out.write(block)
                months[f'{month:02}']={'path':basename,'sha256':sha256_file(target)}
            manifest={
                'schema':MANIFEST_SCHEMA,'source':NORMALS,
                'provider_original_url':ARCHIVE_SOURCE,
                'original_provider_zip_sha256':archive_sha256,
                'units':'mm/month','resolution':'10m','months':months,
                'source_photo_colour_outcomes_read':False
            }
            manifest_path=stage/'monthly_precipitation_manifest.json'
            manifest_path.write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n',encoding='utf8')
            verified=MonthlyNormals(manifest_path,sha256_file(manifest_path),
                                    test_small=_synthetic_small)
            result={
                'schema':'mois_byte_sealed_original_worldclim21_precipitation_zip_v1',
                'status':STATUS,
                'input_archive_sha256':archive_sha256,
                'climatology':'1970–2000 (not photo-year weather)',
                'monthly_data_manifest_sha256':sha256_file(manifest_path),
                'n_byte_verified_monthly_rasters':len(months),
                'n_gis_aligned_monthly_rasters':12,
                'source_raster_width':verified.width,
                'source_raster_height':verified.height,
                'source_original_geotiff_sha256':verified.raster_hashes,
                'original_plant_photos_opened':False,
                'worldwide_flowering_season_or_colour_result':False,
                'requires_independent_original_site_bio12_reconciliation':True,
                'synthetic_test_only':bool(_synthetic_small)
            }
            (stage/'monthly_source_receipt.json').write_text(
                json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf8')
            stage.rename(output)
            return result
        except Exception:
            # Failed staging remains quarantined; no completed source receipt.
            raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-worldclim-zip',required=True,type=Path)
    parser.add_argument('--original-zip-sha256',required=True)
    parser.add_argument('--output-private-raster-source',required=True,type=Path)
    args=parser.parse_args()
    result=seal(args.original_worldclim_zip,args.original_zip_sha256,
                args.output_private_raster_source)
    print(json.dumps({k:result[k] for k in (
        'status','input_archive_sha256','monthly_data_manifest_sha256',
        'n_byte_verified_monthly_rasters')},sort_keys=True))

if __name__=='__main__':
    main()
