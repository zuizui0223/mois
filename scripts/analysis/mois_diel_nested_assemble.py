"""Build experimental-species subset in the broad primary environmental frame."""
from pathlib import Path
import argparse,json,hashlib
import pandas as pd

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run(source,broad,out):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    receipts=sorted(source.glob('*/receipt.json'));assert len(receipts)==16
    data=[];n=0
    for p in receipts:
        r=json.loads(p.read_text());n+=r['source_observation_rows_scanned']
        target=Path('config/diel_experimental_target_species_20261105.csv').read_bytes().replace(b'\r\n',b'\n')
        assert r['target_sha256']==hashlib.sha256(target).hexdigest()
        f=p.parent/'selected_observations.csv.gz';assert sha(f)==r['outputs'][f.name]
        x=pd.read_csv(f,dtype={'inat_taxon_id':str})
        x=x.loc[x.coordinate_precision_5km & x.year.between(2015,2024)]
        data.append(x[['inat_taxon_id','cell025','year','observation_uuid']])
    assert n==53652241
    x=pd.concat(data);assert x.observation_uuid.is_unique
    traits=pd.read_csv('config/diel_experimental_target_species_20261105.csv',dtype={'inat_taxon_id':str})
    assert traits.inat_taxon_id.is_unique
    x=x.merge(traits,on='inat_taxon_id',validate='many_to_one')
    keys=['cell025','year'];frame=pd.read_csv(Path(broad)/'primary_cell_year_source.csv.gz')
    counts=x.groupby(['species']+keys,as_index=False).size().rename(columns={'size':'n_observations'})
    counts=counts.merge(frame[keys],on=keys,validate='many_to_one')
    support=counts.groupby('species').agg(n_cell_years=('cell025','size'),n_cells=('cell025','nunique'),n_years=('year','nunique'),n_observations=('n_observations','sum')).reset_index()
    support['eligible']=support.n_cells.ge(20)&support.n_years.ge(3)
    support=support.merge(traits.drop(columns='inat_taxon_id').drop_duplicates('species'),on='species',validate='one_to_one')
    support.to_csv(out/'experimental_species_support_before_environment.csv',index=False)
    counts=counts.loc[counts.species.isin(support.loc[support.eligible,'species'])]
    counts=counts.merge(traits.drop(columns='inat_taxon_id').drop_duplicates('species'),on='species',validate='many_to_one')
    counts.to_csv(out/'experimental_species_cell_year_source.csv.gz',index=False)
    r={'source_observations_scanned':n,'experimental_precise_2015_2024_observations':len(x),'experimental_species_in_broad_frame':len(support),'eligible_species':int(support.eligible.sum()),'eligible_species_cell_years':len(counts),'eligibility':'at least 20 distinct broad-primary cells and 3 years before inspecting ALAN','night_share_interpretation':'descriptive night/(night+day) treatment outcome; not corrected for complete-exclusion baseline and not proof of effective nocturnal pollinator identity'}
    (out/'nested_source_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--broad',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.source,a.broad,a.out)
