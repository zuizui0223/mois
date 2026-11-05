"""Assemble all signed source groups before environmental selection/modeling."""
from pathlib import Path
import json,argparse,hashlib
import pandas as pd
import numpy as np
from scripts.acquisition.mois_diel_full_source import sha,MANIFEST_SHA

def assemble(source,out):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    paths=sorted(source.glob('*/receipt.json'));assert len(paths)==16
    receipts=[json.loads(p.read_text()) for p in paths]
    assert {r['group'] for r in receipts}==set(range(16))
    assert sum(r['source_observations_scanned'] for r in receipts)==53652241
    assert all(r['source_manifest_sha256']==MANIFEST_SHA for r in receipts)
    cells=[];genera=[];effort=[]
    for p,r in zip(paths,receipts):
        for name,digest in r['outputs'].items():assert sha(p.parent/name)==digest
        cells.append(pd.read_csv(p.parent/'documented_flower_colour_cell_year_counts.csv.gz',dtype={'cell025':str}))
        genera.append(pd.read_csv(p.parent/'documented_flower_colour_genus_cell_year_counts.csv.gz',dtype={'cell025':str,'genus':str}))
        effort.append(pd.read_csv(p.parent/'all_source_plants_cell_year_effort.csv.gz',dtype={'cell025':str}))
    c=pd.concat(cells).groupby(['cell025','year'],as_index=False)[['n_species','n_white_species']].sum()
    g=pd.concat(genera).groupby(['cell025','year','genus'],as_index=False)[['n_species','n_white_species']].sum()
    b=pd.concat(effort).groupby(['cell025','year'],as_index=False)[['n_observations','n_precise_observations']].sum()
    n=g.groupby(['cell025','year']).genus.nunique().rename('n_genera').reset_index()
    c=c.merge(n,on=['cell025','year'],validate='one_to_one').merge(b,on=['cell025','year'],validate='one_to_one')
    assert c.n_white_species.between(0,c.n_species).all()
    c['eligible']=c.n_species.ge(10)&c.n_genera.ge(5)&c.year.between(2015,2024)
    c['white_species_fraction']=c.n_white_species/c.n_species
    c.to_csv(out/'all_documented_species_cell_years.csv.gz',index=False)
    g.to_csv(out/'all_documented_genus_cell_years.csv.gz',index=False)
    primary=c.loc[c.eligible].copy()
    xy=primary.cell025.str.split('_',expand=True).astype(int)
    primary['longitude']=(xy[0]+.5)*.25-180;primary['latitude']=(xy[1]+.5)*.25-90
    primary['block5']=np.floor(primary.longitude/5).astype(int).astype(str)+'_'+np.floor(primary.latitude/5).astype(int).astype(str)
    primary['region10']=np.floor(primary.longitude/10).astype(int).astype(str)+'_'+np.floor(primary.latitude/10).astype(int).astype(str)
    primary.to_csv(out/'primary_cell_year_source.csv.gz',index=False)
    primary[['cell025','year']].sort_values(['year','cell025']).to_csv(out/'environment_queue.csv',index=False)
    receipt={'schema':'mois_full_broad_flora_source_assembled_v1','source_run':41982925538,'source_observations_scanned':53652241,'selected_observations_2014_2024':sum(r['selected_observations_2014_2024'] for r in receipts),'precise_selected_observations_2015_2024':sum(r['precise_selected_observations_2015_2024'] for r in receipts),'target_accepted_species':7730,'species_with_precise_observations':sum(r['species_with_precise_observations'] for r in receipts),'all_cell_years':len(c),'primary_cell_years':len(primary),'primary_distinct_cells':primary.cell025.nunique(),'primary_geographic_blocks5':primary.block5.nunique(),'primary_year_counts':primary.year.value_counts().sort_index().to_dict(),'primary_north_south_rows':{'north':int(primary.latitude.ge(0).sum()),'south':int(primary.latitude.lt(0).sum())},'exposure_not_yet_joined':True,'unknown_occurrence_not_absence':True,'photo_specific_flowering_not_verified':True,'species_colour_not_ROI_colour':True,'outputs':{p.name:sha(p) for p in out.glob('*.csv*')}}
    (out/'assembly_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True);a=p.parse_args();assemble(a.source,a.out)
