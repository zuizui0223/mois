"""Accumulate global environment batches without erasing prior measurements.

2020-only collector output is not interchangeable with another observation year.
Terminal quality-incomplete results are retained, not automatically retried.
Repeated identical terminal results are idempotent; conflicts fail closed.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

STATUS_COLUMNS=[f'viirs_{m}_status' for m in (1,4,7,10)]+['population_global1_0_status']
TERMINAL={'OK','NO_DATA','LOW_CLOUD_COUNT'}
TERMINAL_STATES={'new_environment_complete','queried_environment_quality_incomplete'}
DEFAULTS={'new_environment_status':'not_queried_this_batch','new_viirs_valid_months':0,
          'new_viirs_fourmonth_median':np.nan,'new_population_density_2020':np.nan}
MUTABLE={'collection_state',*DEFAULTS}

def _equal(a,b):
    if pd.isna(a) and pd.isna(b):return True
    if pd.isna(a) or pd.isna(b):return False
    if isinstance(a,(int,float,np.number)) and isinstance(b,(int,float,np.number)):
        # CSV median roundtrips may differ by a few IEEE-754 ULPs. Preserve
        # the stored value; this tolerance is never used to alter source data.
        if not(np.isfinite(a) and np.isfinite(b)):return False
        return abs(float(a)-float(b))<=4*np.spacing(max(abs(float(a)),abs(float(b)),1.0))
    return bool(a==b)

def _derived(row):
    statuses=[row.get(c,'NOT_ATTEMPTED') for c in STATUS_COLUMNS]
    vals=[]
    for m in (1,4,7,10):
        if row.get(f'viirs_{m}_status')=='OK':
            v=pd.to_numeric(row.get(f'viirs_{m}_value'),errors='coerce')
            if not np.isfinite(v):raise ValueError('Successful radiance is missing')
            vals.append(float(v))
    pop=pd.to_numeric(row.get('population_global1_0_value'),errors='coerce')
    pop_ok=row.get('population_global1_0_status')=='OK' and np.isfinite(pop) and pop>=0
    if row.get('population_global1_0_status')=='OK' and not pop_ok:
        raise ValueError('Successful population invalid')
    if any(s not in TERMINAL for s in statuses):state='environment_retryable_error'
    elif len(vals)>=3 and pop_ok:state='new_environment_complete'
    else:state='queried_environment_quality_incomplete'
    return {'collection_state':state,'new_environment_status':state,
            'new_viirs_valid_months':len(vals),
            'new_viirs_fourmonth_median':float(np.median(vals)) if len(vals)>=3 else np.nan,
            'new_population_density_2020':float(pop) if pop_ok else np.nan}

def merge_results(queue,values):
    if not queue.target_id.is_unique or not values.target_id.is_unique:
        raise ValueError('Duplicate source or result target ID')
    if not values.target_id.isin(queue.target_id).all():
        raise ValueError('Result not in source target frame')
    if not values.target_id.astype(str).str.startswith('2020|').all():
        raise ValueError('This population collector is restricted to same-year 2020')
    a=queue.set_index('target_id');b=values.set_index('target_id')
    for col in ('query_longitude','query_latitude'):
        if not np.allclose(a.loc[b.index,col].astype(float),b[col].astype(float),rtol=0,atol=1e-9):
            raise ValueError('Coordinates disagree')
    if not a.loc[b.index,'land_class'].eq(b.land_class).all():
        raise ValueError('Landmass category changed')
    out=queue.copy()
    for col,default in DEFAULTS.items():
        if col not in out:out[col]=default
    positions=pd.Series(out.index,index=out.target_id).to_dict()
    newly_queried=complete=quality=retry=replayed=updated_retry=0
    changed=[]
    for target,row in b.iterrows():
        i=positions[target];previous=out.at[i,'collection_state'];result=_derived(row)
        if previous in TERMINAL_STATES:
            if not all(_equal(out.at[i,col],v) for col,v in result.items()):
                raise ValueError('Conflicting terminal measurement; preserve previous source: '+target)
            replayed+=1
            continue
        if previous not in {'pending_new_VIIRS','environment_retryable_error'}:
            raise ValueError('Cannot replace a legacy/unverified source state: '+str(previous))
        if previous=='pending_new_VIIRS':newly_queried+=1
        else:updated_retry+=1
        for col,v in result.items():out.at[i,col]=v
        changed.append(i)
        state=result['collection_state']
        complete+=state=='new_environment_complete'
        quality+=state=='queried_environment_quality_incomplete'
        retry+=state=='environment_retryable_error'
    for col in queue.columns:
        if col not in MUTABLE and not out[col].equals(queue[col]):
            raise AssertionError('Modified source field '+col)
    untouched=queue.index.difference(changed)
    pd.testing.assert_frame_equal(queue.loc[untouched],out.loc[untouched,queue.columns])
    summary={'queue_targets':len(queue),'input_result_targets':len(values),
       'new_queried_targets':newly_queried,'updated_retry_targets':updated_retry,
       'new_complete_targets':int(complete),'quality_incomplete_targets':int(quality),
       'retryable_targets':int(retry),'replayed_targets':replayed,
       'unchanged_target_count':len(untouched),'source_fields_preserved':True,
       'prior_environment_measurements_preserved':True,
       'new_photos_or_ROI_successes':0,'population_response_year_not_echoed':True,
       'total_complete_targets':int(out.collection_state.eq('new_environment_complete').sum()),
       'total_quality_incomplete_targets':int(out.collection_state.eq('queried_environment_quality_incomplete').sum())}
    return out,summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--queue',type=Path,required=True);p.add_argument('--values',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise ValueError('Refuse to overwrite an existing batch output')
    a.out.mkdir(parents=True)
    q=pd.read_csv(a.queue,dtype={'target_id':str,'cell_id':str,'land_id':str},low_memory=False)
    v=pd.read_csv(a.values,dtype={'target_id':str},low_memory=False)
    result,receipt=merge_results(q,v)
    target=a.out/'GLOBAL_environment_queue_updated.csv.gz'
    result.to_csv(target,index=False,compression={'method':'gzip','compresslevel':1,'mtime':0})
    receipt['queue_input_sha256']=hashlib.sha256(a.queue.read_bytes()).hexdigest()
    receipt['values_input_sha256']=hashlib.sha256(a.values.read_bytes()).hexdigest()
    receipt['queue_output_sha256']=hashlib.sha256(target.read_bytes()).hexdigest()
    (a.out/'queue_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,sort_keys=True))
