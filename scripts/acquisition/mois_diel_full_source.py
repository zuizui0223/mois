"""Select uncapped experimental-species observations from all 53.65m source rows.

No flower-stage assertion or original photograph download. Literature species
traits label photographed plant distributions, not the organ in each photograph.
"""
from pathlib import Path
import argparse,datetime,gzip,hashlib,io,json,subprocess,tarfile
import numpy as np
import pandas as pd

TAG='mois-full-93209351-indexed-photo-observation-source-v1'
MANIFEST_SHA='22f7b6ddd9770e5ef51a093e547e99a25cf704b98de7c4e102d6b1bd99d6ebd9'
SUMS_SHA='6edcf3b807884ddab0e467dd66e981fdc464ef92c8263ebd9015faa888bcbcb4'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

def select_chunk(d,targets):
    x=d.copy();day=pd.to_datetime(x.observed_on,errors='coerce',format='%Y-%m-%d')
    x['year']=day.dt.year;x['month']=day.dt.month
    for c in ['latitude','longitude','positional_accuracy']:x[c]=pd.to_numeric(x[c],errors='coerce')
    x=x.loc[x.year.between(2014,2024)&x.latitude.between(-60,60)&x.longitude.between(-180,180,inclusive='left')].copy()
    x['year']=x.year.astype(int);x['month']=x.month.astype(int)
    x['coordinate_precision_5km']=x.positional_accuracy.between(0,5000)
    x['gx025']=np.floor((x.longitude+180)/.25).astype(int)
    x['gy025']=np.floor((x.latitude+90)/.25).astype(int)
    x['cell025']=x.gx025.astype(str)+'_'+x.gy025.astype(str)
    bg=x.groupby(['cell025','year'],as_index=False).agg(n_observations=('observation_uuid','size'),n_precise_observations=('coordinate_precision_5km','sum'))
    return x.loc[x.inat_taxon_id.isin(targets)].copy(),bg

def run(group,target_file,out,cache):
    out=Path(out);cache=Path(cache);out.mkdir(parents=True,exist_ok=True);cache.mkdir(parents=True,exist_ok=True)
    target=pd.read_csv(target_file,dtype={'inat_taxon_id':str})
    assert target.inat_taxon_id.is_unique
    targets=set(target.inat_taxon_id)
    for name in ['manifest.json','SHA256SUMS.txt',f'source-{group:02d}.tar']:
        if not (cache/name).exists():subprocess.run(['gh','release','download',TAG,'--repo','zuizui0223/mois','--pattern',name,'--dir',str(cache)],check=True)
    assert sha(cache/'manifest.json')==MANIFEST_SHA and sha(cache/'SHA256SUMS.txt')==SUMS_SHA
    manifest=json.loads((cache/'manifest.json').read_text());shards={v['part']:v for v in manifest['shards']}
    sums={line.split()[1]:line.split()[0] for line in (cache/'SHA256SUMS.txt').read_text().splitlines()}
    archive=cache/f'source-{group:02d}.tar';assert sha(archive)==sums[archive.name]
    selected=[];background=[];counts={};seen_parts=set()
    with tarfile.open(archive) as tf:
        for member in tf:
            name=member.name.removeprefix('./')
            if not name.startswith('observations/part-') or not name.endswith('.csv.gz'):continue
            part=int(Path(name).name[5:8]);assert part//8==group and part not in seen_parts;seen_parts.add(part)
            raw=tf.extractfile(member).read();assert hashlib.sha256(raw).hexdigest()==shards[part]['observations_sha256']
            n=0
            for d in pd.read_csv(io.BytesIO(raw),compression='gzip',dtype={'observation_uuid':str,'inat_taxon_id':str,'observer_id':str},chunksize=200000):
                assert (pd.to_numeric(d.inat_taxon_id)%128==part).all()
                n+=len(d);s,b=select_chunk(d,targets);selected.append(s);background.append(b)
            assert n==shards[part]['observation_count'];counts[str(part)]=n
    assert seen_parts==set(range(group*8,(group+1)*8))
    d=pd.concat(selected,ignore_index=True);assert d.observation_uuid.is_unique
    b=pd.concat(background,ignore_index=True).groupby(['cell025','year'],as_index=False)[['n_observations','n_precise_observations']].sum()
    d.to_csv(out/'selected_observations.csv.gz',index=False);b.to_csv(out/'all_plant_cell_year_effort.csv.gz',index=False)
    receipt={'schema':'mois_full_uncapped_diel_species_source_v1','group':group,'manifest_sha256':MANIFEST_SHA,'archive_sha256':sums[archive.name],'target_sha256':sha(target_file),'source_observation_rows_scanned':sum(counts.values()),'part_counts':counts,'selected_observations':len(d),'selected_taxa':d.inat_taxon_id.nunique(),'selected_precise':int(d.coordinate_precision_5km.sum()),'source_years':[2014,2024],'per_species_cap':None,'ROI_success_gate':False,'photo_stage_verified':False,'new_original_pixels':0,'outputs':{p.name:sha(p) for p in out.glob('*.csv.gz')},'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--group',type=int,required=True,choices=range(16));p.add_argument('--targets',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--cache',type=Path,required=True);a=p.parse_args();run(a.group,a.targets,a.out,a.cache)
