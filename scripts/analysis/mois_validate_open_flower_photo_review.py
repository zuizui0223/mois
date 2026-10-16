#!/usr/bin/env python3
"""Independent two-reviewer anthesis/flower-organ QC for original FCP pixels.

Never infer open flowers from ROI-v4 colour, observation-month or an unreturned
iNaturalist annotation. Open anthesis of a *visible* flower plus an ROI on a
petal/tepal is primary, while other floral structures are a separate stratum.
A photograph of buds/fruit or an unresolved ROI is not scored as white/nonwhite.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,io,zipfile
from pathlib import Path

REVIEW_COLUMNS=("measurement_id","photo_id","species","observation_url","photo_url_large",
    "photo_license","r1_reviewer_id","r1_open_flower_stage","r1_roi_tissue",
    "r2_reviewer_id","r2_open_flower_stage","r2_roi_tissue",
    "adjudicator_id","adjudicated_stage","adjudicated_roi_tissue","review_notes")
STAGE={"open_flower","buds_only","fruit_only","no_flower_visible","uncertain"}
TISSUE={"petal_tepal","other_open_flower_organ","spikelet_glume_bract","fruit_seed",
        "leaf_stem","background","uncertain"}


def digest(b:bytes)->str:return hashlib.sha256(b).hexdigest()


def rows(buffer:bytes)->list[dict]:
    with io.StringIO(buffer.decode('utf-8')) as f:
        return list(csv.DictReader(f))


def read_receipt(zip_file:zipfile.ZipFile)->tuple[list[dict],list[dict],dict]:
    source_file='measured/first_original_geographically_matched_photo_terminal.csv'
    acq_file='prepared/first_matched_pair_original_photo_source_acquisition.csv'
    report=json.loads(zip_file.read('measured/first_geo_pair_real_pixel_receipt.json'))
    if report['status']!='FIRST_GENUINE_MONTHLY_ORIGINAL_GEO_MATCHED_PIXELS_MEASURED_NOT_1000_SPECIES' or report['old_measured_colours_imported'] is not False:
        raise ValueError('No genuine frozen original photo pixel method attested')
    raw=zip_file.read(source_file)
    if digest(raw)!=report['original_geographic_source_join_terminal_sha256']:
        raise ValueError('Changed original photo pixels CSV or source identities')
    measured=rows(raw)
    if len(measured)!=report['new_original_photo_pixels_processed']:
        raise ValueError('Original real ROI terminal photos incomplete')
    acq=rows(zip_file.read(acq_file))
    if len(acq)!=len(measured) or len({q['measurement_id'] for q in acq})!=len(acq):
        raise ValueError('Original acquisition source photo IDs missing or repeated')
    if {r['measurement_id'] for r in acq}!={r['measurement_id'] for r in measured}:
        raise ValueError('Source photo IDs do not match original pixel results')
    return measured,acq,report


def make_blank(zip_path:Path,output:Path)->dict:
    if output.exists():raise ValueError('Never overwrite original botanist review form')
    with zipfile.ZipFile(zip_path) as z:measured,acq,report=read_receipt(z)
    urls={x['measurement_id']:x for x in acq}
    with output.open('x',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=REVIEW_COLUMNS,lineterminator='\n');writer.writeheader()
        for r in measured:
            a=urls[r['measurement_id']]
            assert a['photo_license']==r.get('photo_license',a['photo_license']) or 'photo_license' not in r
            writer.writerow({
                'measurement_id':r['measurement_id'],'photo_id':r['photo_id'],
                'species':r['species'],
                'observation_url':'https://www.inaturalist.org/observations/'+r['observation_id'],
                'photo_url_large':a['photo_url_large'],'photo_license':a['photo_license'],
            })
    return {'original_photos':len(measured),'original_source_photo_terminal_sha256':report['original_geographic_source_join_terminal_sha256'],
            'review_form_sha256':digest(output.read_bytes()),'original_unreviewed_flower_stages':len(measured),
            'not_true_nonflower_or_white_rate':True}


def decide(stage1,tissue1,stage2,tissue2):
    if not all((stage1,tissue1,stage2,tissue2)):return 'unreviewed'
    if stage1 not in STAGE or stage2 not in STAGE or tissue1 not in TISSUE or tissue2 not in TISSUE:
        raise ValueError('Unknown predeclared botanical flowering review code')
    if (stage1,tissue1)!=(stage2,tissue2):return 'needs_independent_adjudication'
    if stage1=='uncertain' or tissue1=='uncertain':return 'unresolved_flower_stage_or_roi'
    if stage1!='open_flower':return 'non_open_flower_reviewed_not_colour'
    if tissue1=='petal_tepal':return 'petal_tepal_open_flower_reviewed'
    if tissue1=='other_open_flower_organ':return 'other_flower_organ_secondary_not_primary_colour'
    return 'noncomparable_bract_fruit_foliage_or_background'


def audit(review:Path,zip_path:Path,outdir:Path)->dict:
    if outdir.exists() and any(outdir.iterdir()):raise ValueError('Never overwrite original phenotype human review')
    with zipfile.ZipFile(zip_path) as z:measured,acq,orig=read_receipt(z)
    expected={r['measurement_id']:r for r in measured}
    urls={r['measurement_id']:r for r in acq}
    out=[];seen=set();counts={}
    with review.open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=REVIEW_COLUMNS:
            raise ValueError('Botanical review form exposes forbidden outcome columns or wrong schema')
        for r in reader:
            mid=r['measurement_id']
            if mid in seen or mid not in expected:
                raise ValueError('Unknown or duplicated source original photo in botanical review')
            seen.add(mid)
            b=expected[mid];a=urls[mid]
            if (r['photo_id']!=b['photo_id'] or r['species']!=b['species'] or
                r['photo_url_large']!=a['photo_url_large'] or
                r['observation_url']!='https://www.inaturalist.org/observations/'+b['observation_id'] or
                r['photo_license']!=a['photo_license']):
                raise ValueError('Original image/source identity changed since selection')
            r1=r['r1_reviewer_id'].strip();r2=r['r2_reviewer_id'].strip()
            stage1=r['r1_open_flower_stage'].strip();stage2=r['r2_open_flower_stage'].strip()
            t1=r['r1_roi_tissue'].strip();t2=r['r2_roi_tissue'].strip()
            if any((stage1,t1)) and not r1 or any((stage2,t2)) and not r2 or r1 and r2 and r1==r2:
                raise ValueError('Two independent botanical reviewers required')
            status=decide(stage1,t1,stage2,t2)
            adjudicator=r['adjudicator_id'].strip();astage=r['adjudicated_stage'].strip();atissue=r['adjudicated_roi_tissue'].strip()
            if any((adjudicator,astage,atissue)):
                if not all((adjudicator,astage,atissue)) or adjudicator in (r1,r2) or status!='needs_independent_adjudication':
                    raise ValueError('Only third independent reviewer may resolve conflicting full reviews')
                status=decide(astage,atissue,astage,atissue)
            classifiable=b['measurement_status']=='classified_four_state_morph'
            w=float(b['flower_fraction_white']) if classifiable else 0
            y=float(b['flower_fraction_yellow']) if classifiable else 0
            organ_mix=classifiable and w>=.1 and y>=.1
            if status=='petal_tepal_open_flower_reviewed' and (not classifiable or organ_mix):
                status='anthesis_reviewed_but_original_photo_colour_not_usable'
            counts[status]=counts.get(status,0)+1
            out.append({'measurement_id':mid,'photo_id':b['photo_id'], 'species':b['species'],
                'botanical_review_status':status,'method_four_state_classifiable':classifiable,
                'same_roi_white_and_yellow_organs_mixed':organ_mix,
                'eligible_original_open_petaloid_flower_photo_colour':status=='petal_tepal_open_flower_reviewed',
                'is_flowers_not_fruits_inferred_by_automated_colour':False})
    if seen!=set(expected):raise ValueError('Missing 64 source photos in the independent review')
    outdir.mkdir(parents=True,exist_ok=True)
    dest=outdir/'audited_first_original64_botanical_flower_stage.csv'
    with dest.open('x',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(out[0]),lineterminator='\n');writer.writeheader();writer.writerows(out)
    status='BOTANICAL_REVIEW_INCOMPLETE_NO_FLOWER_WHITE_RATE' if counts.get('unreviewed',0)+counts.get('needs_independent_adjudication',0)+counts.get('unresolved_flower_stage_or_roi',0) else 'BOTANICAL_STAGE_REVIEW_COMPLETED_NOT_CAUSAL'
    report={'schema':'fcp_first_64_original_photo_botanical_stage_adjudication_v1','status':status,
        'n_original_pixels':len(expected),'botanical_statuses':counts,
        'n_reviewed_primary_flower_organs':sum(x['eligible_original_open_petaloid_flower_photo_colour'] for x in out),
        'unreviewed_not_treated_as_nonflower':True,'original_image_pixels_not_reopened':True,
        'no_white_rate_without_complete_bloom_and_organ_validation':True,
        'prior_colour_model_not_used_to_label_flowering':True,
        'input_review_sheet_sha256':digest(review.read_bytes()),
        'source_original_pixel_sha256':orig['original_geographic_source_join_terminal_sha256'],
        'botanical_review_qc_sha256':digest(dest.read_bytes())}
    (outdir/'audited_first_original64_botanical_flower_stage_receipt.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--real-first64-archive',type=Path,required=True)
    p.add_argument('--review-sheet',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(audit(a.review_sheet,a.real_first64_archive,a.outdir),sort_keys=True))
if __name__=='__main__':main()
