#!/usr/bin/env python3
"""Offline, two-reviewer botanical petal/background visible-colour instrument.

This creates a review tool only, never a botanical label or image-derived
contrast. Reviewers must load images and mark petal and nearby-background
pixels independently. Source IDs remain strict9 UNKNOWN throughout.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

SCHEMA = 'mois_existing128_two_blind_botanical_petal_background_srgb_v1'
ROSTER_SHA = '5fef860e4d2d0ff4b2dfec850253e4a08c8dc87878647d6712bc61e4b75dcb56'
ROSTER_N = 128
ROSTER_COLUMNS = ('photo_id', 'observation_id', 'inat_taxon_id', 'source_design',
                  'latitude_zone', 'local_season', 'postoutcome_source_gap',
                  'screen_status', 'technical_diagnostic_class', 'review_priority')
REVIEW_COLUMNS = (ROSTER_COLUMNS + ('source_roster_sha256', 'reviewer_slot',
                  'reviewer_id', 'flower_visible', 'petal_visible', 'background_visible',
                  'photo_file_sha256', 'image_width', 'image_height', 'image_file_name',
                  'fg_points_json', 'bg_points_json', 'fg_r', 'fg_g', 'fg_b',
                  'bg_r', 'bg_g', 'bg_b', 'deltaE76', 'deltaLstar',
                  'sRGB_luminance_ratio', 'notes'))
CHOICES = {'yes', 'no', 'uncertain', ''}


def require(ok, message):
    if not bool(ok):
        raise ValueError('MOIS_OPTICAL_REVIEW: ' + message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda: f.read(4 << 20), b''):
            h.update(part)
    return h.hexdigest()


def _srgb_linear(value):
    v = float(value) / 255
    require(0 <= v <= 1 and math.isfinite(v), 'RGB channel out of range')
    return v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4


def _colour_to_lab(rgb):
    require(len(rgb) == 3, 'RGB requires 3 channels')
    r, g, b = [_srgb_linear(q) for q in rgb]
    # IEC sRGB / D65 / CIE 1931 2 degree white; pixel colours are uncalibrated.
    x = (0.4124564*r + 0.3575761*g + 0.1804375*b)/0.95047
    y = (0.2126729*r + 0.7151522*g + 0.0721750*b)/1.
    z = (0.0193339*r + 0.1191920*g + 0.9503041*b)/1.08883
    delta = 6/29
    def f(t):
        return t**(1/3) if t > delta**3 else t/(3*delta**2)+4/29
    fx, fy, fz = f(x), f(y), f(z)
    return (116*fy-16, 500*(fx-fy), 200*(fy-fz))


def colour_contrast(fg, bg):
    p, q = _colour_to_lab(fg), _colour_to_lab(bg)
    diff = math.dist(p, q)
    yl = [sum(w*_srgb_linear(c) for w, c in zip((.2126,.7152,.0722), rgb))
          for rgb in (fg,bg)]
    ratio = (max(yl)+.05)/(min(yl)+.05)
    return dict(deltaE76=float(diff),deltaLstar=float(p[0]-q[0]),
                sRGB_luminance_ratio=float(ratio))


def _load_roster(path, expected_sha=None, expected_n=ROSTER_N):
    path = Path(path)
    if expected_sha is not None:
        require(sha(path)==expected_sha,'Frozen source roster SHA changed')
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        reader=csv.DictReader(f)
        require(tuple(reader.fieldnames or ())==ROSTER_COLUMNS,
                'Unexpected source roster schema')
        rows=list(reader)
    require(len(rows)==expected_n, 'Frozen review roster count changed')
    for field in ('photo_id','observation_id','inat_taxon_id'):
        values=[row[field] for row in rows]
        require(all(values) and len(set(values))==len(values),
                'Source roster duplicate or missing '+field)
    require(all(x['photo_id'].isdecimal() for x in rows), 'Invalid real photo IDs')
    return rows


def build(roster, template, out, *, expected_sha=ROSTER_SHA, expected_n=ROSTER_N):
    roster,template,out=map(Path,(roster,template,out))
    require(not out.exists(), 'Refuse to overwrite a botanical review packet')
    rows=_load_roster(roster, expected_sha, expected_n)
    markup=template.read_text(encoding='utf-8')
    require('__SOURCE_ROWS__' in markup and '__SOURCE_SHA__' in markup and
            '__REVIEWER_SLOT__' in markup, 'Review template is incomplete')
    source=sha(roster)
    out.mkdir(parents=True)
    files={}
    for slot in 'AB':
        ordered=sorted(rows,key=lambda x:hashlib.sha256(
            ('MOIS_128_CONTRAST_REVIEW_BLIND_'+slot+'|'+x['photo_id']).encode()).hexdigest())
        payload=json.dumps(ordered,ensure_ascii=False,separators=(',',':'))
        payload=payload.replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
        page=markup.replace('__SOURCE_ROWS__',payload).replace('__SOURCE_SHA__',source).replace('__REVIEWER_SLOT__',slot)
        target=out/('mois_manual_contrast_reviewer_'+slot+'.html')
        target.write_text(page,encoding='utf-8')
        files[target.name]=sha(target)
    receipt={'schema':SCHEMA,'source_roster_sha256':source,
             'source_ROI9_original_status':'strict9_unknown_unchanged',
             'n_original_unknown_ROI9_photo_ids':len(rows),
             'original_photo_ids_unique':True,'original_observation_UUIDs_unique':True,
             'original_taxon_IDs_unique':True,'botanical_confirmed_count':0,
             'optical_contrast_valid_count':0,
             'independent_reviewers_needed':2,
             'only_manual_photo_file_loading_not_automatic_fetch':True,
             'no_new_image_acquisition_by_generator':True,
             'no_botanical_labels_or_colours_fabricated':True,
             'original_strict9_unchanged':True,
             'image_bytes_archived':False,'derived_image_lab_is_not_calibrated_spectral_reflectance':True,
             'deltaE76_requires_manually_verified_petal_and_local_background':True,
             'outputs_sha256':files}
    (out/'mois128_optical_review_blank_packet_receipt.json').write_text(
        json.dumps(receipt,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    return receipt


def _point_valid(row, which):
    raw=row.get(which+'_points_json','')
    try:
        points=json.loads(raw)
    except (TypeError,ValueError) as exc:
        raise ValueError('MOIS_OPTICAL_REVIEW: invalid '+which+' point array') from exc
    require(isinstance(points,list) and len(points)==3,'Three '+which+' points required')
    require(len({(p[0],p[1]) for p in points})==3,'Duplicate '+which+' point')
    width,height=int(row['image_width']),int(row['image_height'])
    require(10<=width<=12000 and 10<=height<=12000,'Image dimensions out of range')
    for p in points:
        require(isinstance(p,list) and len(p)==5,'Missing image-pixel RGB record')
        x,y,*rgb=p
        require(all(isinstance(v,int) for v in p),'Noninteger pixel or RGB record')
        require(2<=x<width-2 and 2<=y<height-2,'point outside image margins')
        require(all(0<=v<=255 for v in rgb),'RGB point channel outside range')
    return points


def _median(xs):
    return sorted(xs)[len(xs)//2]


def _validated_vote(row):
    for col in ('flower_visible','petal_visible','background_visible'):
        require(row[col] in CHOICES,'Unexpected botanical review choice')
    valid = all(row[x]=='yes' for x in ('flower_visible','petal_visible','background_visible'))
    keys = ('fg_points_json','bg_points_json','fg_r','fg_g','fg_b',
            'bg_r','bg_g','bg_b','deltaE76','deltaLstar','sRGB_luminance_ratio')
    if not valid:
        require(not any(str(row.get(x,'')).strip() for x in keys),
                'Unverified plant organ cannot have contrast metrics')
        return None
    digest=row.get('photo_file_sha256','').lower()
    require(len(digest)==64 and all(c in '0123456789abcdef' for c in digest),
            'Human-selected photo bytes are not SHA-recorded')
    fg,bg=_point_valid(row,'fg'),_point_valid(row,'bg')
    require(not (set((x[0],x[1]) for x in fg) & set((x[0],x[1]) for x in bg)),
            'shared foreground/background pixel')
    for name,points in [('fg',fg),('bg',bg)]:
        values=tuple(_median([p[j] for p in points]) for j in (2,3,4))
        declared=tuple(int(row[name+'_'+x]) for x in ('r','g','b'))
        require(values==declared,'Original point median changed: '+name)
    derived=colour_contrast(tuple(int(row['fg_'+k]) for k in ('r','g','b')),
                            tuple(int(row['bg_'+k]) for k in ('r','g','b')))
    for key,val in derived.items():
        try:
            measured=float(row[key])
        except (TypeError,ValueError) as exc:
            raise ValueError('MOIS_OPTICAL_REVIEW: derived colour metric missing') from exc
        require(math.isfinite(measured) and abs(measured-val)<.002,
                'derived colour metric cannot be reproduced: '+key)
    return derived


def reconcile(roster, first, second, *, expected_roster_sha=None):
    """Private local validator; never relabel original ROI9 or source-release ID.

    Returns photo ID rows only to the coordinator locally, never for public upload.
    Blank/uncertain/disagreeing cases cannot become verified contrast numbers.
    """
    keys={x['photo_id']:x for x in roster}
    require(len(keys)==len(roster) and len(roster)>0,'Duplicate original source IDs')
    require(len(first)==len(second)==len(roster),'Reviewer files lost original photos')
    def check(records,slot):
        require(all(r.get('reviewer_slot')==slot for r in records),
                'Reviewer slot changed')
        names={r.get('reviewer_id','').strip() for r in records}
        require(len(names)==1 and '' not in names,'Reviewer identity not complete')
        index={r['photo_id']:r for r in records}
        require(len(index)==len(records) and index.keys()==keys.keys(),
                'Reviewer source photo IDs altered')
        for photo_id,row in index.items():
            for col in ROSTER_COLUMNS:
                require(row.get(col)==keys[photo_id].get(col),
                        'Reviewer source identity changed: '+col)
            if expected_roster_sha is not None:
                require(row.get('source_roster_sha256')==expected_roster_sha,
                        'Reviewer source roster SHA changed')
            _validated_vote(row)
        return index,next(iter(names))
    a,aid=check(first,'A');b,bid=check(second,'B')
    require(aid!=bid,'Two independent reviewer identities required')
    output=[];totals=dict(confirmed=0, disagreement=0,pending=0, nonpetal=0, different_bytes=0)
    for photo_id in sorted(keys,key=int):
        x,y=a[photo_id],b[photo_id]
        state_x=tuple(x[t] for t in ('flower_visible','petal_visible','background_visible'))
        state_y=tuple(y[t] for t in ('flower_visible','petal_visible','background_visible'))
        complete_x=all(s in ('yes','no') for s in state_x)
        complete_y=all(s in ('yes','no') for s in state_y)
        if not complete_x or not complete_y:
            totals['pending']+=1;continue
        if state_x!=state_y:
            totals['disagreement']+=1;continue
        if state_x!=('yes','yes','yes'):
            totals['nonpetal']+=1;continue
        if x['photo_file_sha256']!=y['photo_file_sha256']:
            totals['different_bytes']+=1;continue
        ax=_validated_vote(x);by=_validated_vote(y)
        ensure=ax is not None and by is not None
        require(ensure,'Positive voting requires numeric image samples')
        totals['confirmed']+=1
        output.append({'photo_id':photo_id,'observation_id':keys[photo_id]['observation_id'],
                       'inat_taxon_id':keys[photo_id]['inat_taxon_id'],
                       'source_design':keys[photo_id]['source_design'],
                       'latitude_zone':keys[photo_id]['latitude_zone'],
                       'local_season':keys[photo_id]['local_season'],
                       'deltaE76_mean':(ax['deltaE76']+by['deltaE76'])/2,
                       'deltaE76_reviewer_gap':abs(ax['deltaE76']-by['deltaE76']),
                       'deltaLstar_mean':(ax['deltaLstar']+by['deltaLstar'])/2,
                       'luminance_ratio_mean':(ax['sRGB_luminance_ratio']+by['sRGB_luminance_ratio'])/2,
                       'source_photo_SHA256':x['photo_file_sha256']})
    return {'schema':SCHEMA,'original_photo_ids':len(roster),
            'technical_ROI9_unknown_still_unknown':len(roster),
            'confirmed_optical_photo_ids':totals['confirmed'],
            'unresolved_or_disagreement':totals['pending']+totals['disagreement'],
            'different_image_bytes':totals['different_bytes'],
            'confirmed_nonpetal_cases':totals['nonpetal'],
            'confirmed_records':output,
            'claim_botanically_verified_colour_distribution':False,
            'claim_world_native_flora_conspicuousness':False,
            'requires_independent_source_photo_file_authentication_beyond_human_identity_check':True}


def _read_rows(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='task',required=True)
    build_p=sub.add_parser('build')
    build_p.add_argument('--roster',type=Path,required=True)
    build_p.add_argument('--template',type=Path,required=True)
    build_p.add_argument('--output',type=Path,required=True)
    review_p=sub.add_parser('validate-private')
    review_p.add_argument('--roster',type=Path,required=True)
    review_p.add_argument('--reviewer-a',type=Path,required=True)
    review_p.add_argument('--reviewer-b',type=Path,required=True)
    review_p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.task=='build':
        print('MOIS_OFFLINE_BLIND_OPTICAL_REVIEW_PACKET_READY',json.dumps(build(
            args.roster,args.template,args.output),sort_keys=True))
    else:
        original=_load_roster(args.roster,ROSTER_SHA,ROSTER_N)
        verified=reconcile(original,_read_rows(args.reviewer_a),_read_rows(args.reviewer_b),
                           expected_roster_sha=ROSTER_SHA)
        require(not args.output.exists(),'Refuse overwrite local private human review')
        args.output.mkdir(parents=True)
        outcome={k:v for k,v in verified.items() if k!='confirmed_records'}
        (args.output/'PRIVATE_review_aggregate.json').write_text(json.dumps(outcome,indent=2)+'\n')
        if verified['confirmed_records']:
            with (args.output/'PRIVATE_photo_background_contrast_rows.csv').open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(verified['confirmed_records'][0]))
                writer.writeheader();writer.writerows(verified['confirmed_records'])
        print('MOIS_PRIVATE_REVIEW_COMPLETENESS',json.dumps(outcome,sort_keys=True))

if __name__=='__main__':main()
