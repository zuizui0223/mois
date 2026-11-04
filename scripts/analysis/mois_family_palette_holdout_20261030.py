#!/usr/bin/env python3
"""Signed Kew WCVP accepted-family versus accepted-genus geographic transfer.

Uses existing measured ROI9 numeric traits. Original Kew species names are
read as metadata; full Kew archive and original photos are not redistributed.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from pathlib import Path
import pandas as pd
import numpy as np

METRICS=('strict9_white','strict9_blue','strict9_purple','strict9_bronze','palette_entropy')
PAIRS=(('arctic','north_subpolar'),('south_temperate','north_temperate'),('island_cached','mainland_cached'))
PINS={'wcvp_zip':'d32ea2b3a85e489b14e83bcc9eae7274532e1d113753f7be290d4b2dfde573fa',
      'species':'8fef269fba22e084e4ff859ca43629ce5cb1a87e3fa781ce39e9e060d342c150',
      'genus_predictions':'3c13b98b9d9c589da0ecfad3a14c52b173cbc7d0ed94adab62114627c7f4850f'}

def require(ok,msg):
    if not bool(ok):raise ValueError('MOIS_FAMILY_HOLDOUT: '+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def read_accepted_families(archive,accepted_ids):
    """Stream accepted species names with reliable families, never photo native status."""
    accepted_ids=set(map(str,accepted_ids));out={};rows=0
    with zipfile.ZipFile(archive) as z:
        candidates=[x for x in z.namelist() if Path(x).name.lower()=='wcvp_names.csv']
        require(len(candidates)==1,'Exactly one official Kew name table required')
        with z.open(candidates[0]) as raw:
            txt=io.TextIOWrapper(raw,encoding='utf-8-sig',newline='')
            rd=csv.DictReader(txt,delimiter='|')
            require({'plant_name_id','taxon_name','taxon_status','taxon_rank','family','genus'}<=set(rd.fieldnames or ()),
                    'Unexpected official WCVP family header or delimiter')
            for row in rd:
                rows+=1;pid=str(row['plant_name_id']).strip()
                if pid not in accepted_ids:continue
                if row['taxon_status'].strip().casefold()!='accepted' or row['taxon_rank'].strip().casefold()!='species':
                    continue
                name=str(row['taxon_name']).strip()
                family=str(row['family']).strip()
                genus=str(row['genus']).strip()
                if not name or not family or not genus or not name.startswith(genus+' '):continue
                result=(pid,name,genus,family)
                if pid in out:
                    require(out[pid]==result,'Conflicting official family or accepted name for same plant ID')
                else:out[pid]=result
    table=pd.DataFrame(sorted(out.values()),columns=['accepted_id','accepted_name','genus','family'])
    stats={'accepted_species_ids_requested':len(accepted_ids),
           'accepted_species_with_verified_family':len(table),
           'unresolved_accepted_ids':len(accepted_ids)-len(table),
           'official_WCVP_name_rows_scanned':rows,'official_WCVP_delimiter':'pipe',
           'source_region_or_native_status_assigned':False}
    return table,stats

def join_family(species,families):
    require({'sector','accepted_id','genus'}<=set(species),'Original accepted species source invalid')
    require({'accepted_id','accepted_name','genus','family'}<=set(families),'Original Kew family mapping missing')
    require(not species.duplicated(['sector','accepted_id']).any(),'Repeated species sector')
    require(families.accepted_id.is_unique,'Duplicate accepted family ID')
    j=species.merge(families,on='accepted_id',how='left',validate='many_to_one',suffixes=('_source','_wcvp'))
    found=j.family.notna()
    require(j.loc[found,'genus_source'].eq(j.loc[found,'genus_wcvp']).all(),
            'WCVP genus drift from signed prior taxon mapping')
    mapped=j.loc[found].copy().rename(columns={'genus_source':'genus'})
    mapped=mapped.drop(columns=['genus_wcvp'])
    conflicts=mapped.groupby('genus').family.nunique()
    require(conflicts.le(1).all(),'Cross-family genus identity disagreement')
    stats={'accepted_species_sector_rows':len(j),
           'family_resolved_species_sector_rows':len(mapped),
           'family_unresolved_species_sector_rows':len(j)-len(mapped),
           'photo_records_with_family':int(mapped.n.sum()) if 'n' in mapped else None,
           'distinct_accepted_families':int(mapped.family.nunique())}
    return mapped,stats

def family_composition(species,focal,reference,metric):
    a=species.loc[species.sector.eq(focal)].groupby('family')[metric].mean()
    b=species.loc[species.sector.eq(reference)].groupby('family')[metric].mean()
    require(len(a)>0 and len(b)>0,'No accepted families in regional comparison')
    common=a.index.intersection(b.index);one=a.index.difference(b.index);two=b.index.difference(a.index)
    wa=1/len(a);wb=1/len(b)
    terms={'common_family_trait_pp':100*float(((a.loc[common]-b.loc[common])*(wa+wb)/2).sum()),
           'common_family_weight_pp':100*float(((a.loc[common]+b.loc[common])*(wa-wb)/2).sum()),
           'focal_only_families_pp':100*float(a.loc[one].sum()*wa),
           'reference_only_families_pp':-100*float(b.loc[two].sum()*wb)}
    delta=100*float(a.mean()-b.mean())
    require(abs(sum(terms.values())-delta)<1e-9,'Accepted family composition does not reconcile')
    return {'focal':focal,'reference':reference,'metric':metric,
            'focal_families':len(a),'reference_families':len(b),
            'shared_families':len(common),'family_equal_difference_pp':delta,
            'shared_family_only_difference_pp':100*float((a.loc[common]-b.loc[common]).mean()) if len(common) else None,
            **terms}

def score_family_holdout(species,focal,reference,metrics=METRICS,min_species=5):
    require(min_species>=2,'Family must have independent training species')
    a=species.loc[species.sector.eq(focal)].copy()
    b=species.loc[species.sector.eq(reference)].copy()
    require(a.accepted_id.is_unique and b.accepted_id.is_unique,
            'Family source duplicated accepted species')
    known=set(b.accepted_id);new=a.loc[~a.accepted_id.isin(known)]
    means=b.groupby('family',as_index=False).agg({'accepted_id':'nunique',**{k:'mean' for k in metrics}})
    means=means.rename(columns={'accepted_id':'training_species_in_family',
                                **{k:'pred_'+k for k in metrics}})
    supported=means.loc[means.training_species_in_family.ge(min_species)]
    scored=new.merge(supported,on='family',how='inner',validate='many_to_one').sort_values('accepted_id').reset_index(drop=True)
    require(not scored.accepted_id.isin(known).any() and scored.accepted_id.is_unique,
            'Family target-train leakage')
    scored.insert(0,'focal',focal);scored.insert(1,'reference',reference)
    result={'focal':focal,'reference':reference,'novel_target_species':len(new),
            'supported_novel_species':len(scored),
            'family_supported_fraction_novel':len(scored)/len(new) if len(new) else None,
            'family_supported_fraction_all_target':len(scored)/len(a) if len(a) else None,
            'family_training_species_minimum':min_species,
            'supported_families':int(scored.family.nunique())}
    return result,scored

def compare_family_genus(family_pred,genus_pred,metrics=METRICS):
    keys=['focal','reference','accepted_id']
    require(not family_pred.duplicated(keys).any() and not genus_pred.duplicated(keys).any(),
            'Duplicate heldout species')
    selected=genus_pred[keys+['pred_'+m for m in metrics]].rename(
        columns={'pred_'+m:'genus_pred_'+m for m in metrics})
    matched=family_pred.merge(selected,on=keys,how='inner',validate='one_to_one')
    results={'compared_novel_species':len(matched),'metrics':{}}
    for metric in metrics:
        if len(matched):
            obs=matched[metric].to_numpy(float)
            pf=matched['pred_'+metric].to_numpy(float)
            pg=matched['genus_pred_'+metric].to_numpy(float)
            results['metrics'][metric]={
                'family_mse':float(np.mean((obs-pf)**2)),
                'genus_mse':float(np.mean((obs-pg)**2)),
                'family_mae':float(np.mean(np.abs(obs-pf))),
                'genus_mae':float(np.mean(np.abs(obs-pg)))}
        else:
            results['metrics'][metric]={'family_mse':None,'genus_mse':None,
                                        'family_mae':None,'genus_mae':None}
    return results,matched

def run(accepted_species,genus_predictions,wcvp_zip,out):
    out=Path(out);require(not out.exists(),'Immutable family output exists')
    inputs={'species':Path(accepted_species),'genus_predictions':Path(genus_predictions),
            'wcvp_zip':Path(wcvp_zip)}
    checks={k:sha(p) for k,p in inputs.items()}
    require(checks==PINS,'SHA-pinned original WCVP or genus stage source changed')
    source=pd.read_csv(inputs['species'],compression='gzip',dtype={'accepted_id':str},low_memory=False)
    genus=pd.read_csv(inputs['genus_predictions'],compression='gzip',dtype={'accepted_id':str},low_memory=False)
    require(len(source)==30054 and source.sector.nunique()==6 and len(genus)==6508,
            'Unexpected source genus stage input')
    official,wd=read_accepted_families(inputs['wcvp_zip'],source.accepted_id.unique())
    mapped,coverage=join_family(source,official)
    require(len(mapped)>26000 and mapped.family.nunique()>100,
            'Incomplete official family source coverage')
    per_family=[];family_scoring=[];matched=[];heldout=[]
    for a,b in PAIRS:
        for metric in METRICS:
            per_family.append(family_composition(mapped,a,b,metric))
        for focal,ref in ((a,b),(b,a)):
            d,fr=score_family_holdout(mapped,focal,ref,METRICS,5)
            gg=genus.loc[genus.focal.eq(focal)&genus.reference.eq(ref)]
            comparison,pairs=compare_family_genus(fr,gg,METRICS)
            d['common_family_genus_novel_species']=comparison['compared_novel_species']
            for metric in METRICS:
                for k,v in comparison['metrics'][metric].items():
                    d[metric+'_'+k]=v
            family_scoring.append(d);heldout.append(fr);matched.append(pairs)
    out.mkdir(parents=True)
    datasets={'accepted_family_species_region.csv.gz':mapped,
              'family_composition_contrasts.csv':pd.DataFrame(per_family),
              'family_vs_genus_holdout.csv':pd.DataFrame(family_scoring),
              'family_heldout_predictions.csv.gz':pd.concat(heldout,ignore_index=True),
              'paired_family_genus_species_predictions.csv.gz':pd.concat(matched,ignore_index=True)}
    for name,frame in datasets.items():
        if name.endswith('.gz'):
            frame.to_csv(out/name,index=False,compression={'method':'gzip','mtime':0})
        else:frame.to_csv(out/name,index=False)
    report={
        'schema':'mois_original_WCVP_accepted_family_vs_genus_region_heldout_v1',
        'input_SHA256':checks,'WCVP_official_accepted_family_source':wd,
        'family_region_coverage':coverage,'family_composition_terms':len(per_family),
        'six_directional_heldout':family_scoring,
        'family_family_source_only_not_functional_pollen_guild':True,
        'green_and_flower_background_colour_not_measured':True,
        'no_photo_download_or_roi_reinference':True,'MNAR_scenario_not_run':True,
        'no_native_island_photo_status_assigned':True,
        'out_SHA256':{p.name:sha(p) for p in sorted(out.iterdir())}}
    (out/'family_holdout_receipt.json').write_text(
        json.dumps(report,sort_keys=True,indent=2,allow_nan=False)+'\n')
    print('MOIS_ORIGINAL_WCVP_FAMILY_GENUS_REGION_HELDOUT_COMPLETE',json.dumps({
        'mapped_region_species':len(mapped),
        'families':coverage['distinct_accepted_families'],
        'holdout_directions':len(family_scoring)},sort_keys=True),flush=True)
    for row in family_scoring:
        print('FAMILY_HELDOUT',row['focal'],'from',row['reference'],
            'family_supported',row['supported_novel_species'],
            'paired_genus_family',row['common_family_genus_novel_species'],
            'family_mse',row['strict9_white_family_mse'],
            'genus_mse',row['strict9_white_genus_mse'],flush=True)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for x in ('accepted_species','genus_predictions','wcvp_zip','out'):
        p.add_argument('--'+x.replace('_','-'),dest=x,required=True)
    a=p.parse_args()
    run(a.accepted_species,a.genus_predictions,a.wcvp_zip,a.out)
