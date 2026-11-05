"""Broad documented single-colour flora from the full observation frame.

Publish aggregate cell/year and genus/cell/year counts only, never original
third-party species trait rows. One source taxon per accepted species avoids
duplicating species across taxon-modulo partitions.
"""
from pathlib import Path
import argparse,gzip,hashlib,io,json,subprocess,tarfile
import pandas as pd
from scripts.acquisition.mois_diel_full_source import TAG,MANIFEST_SHA,SUMS_SHA,sha,select_chunk

TRAIT_SHA='7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577'
NAME_SHA='dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2'
ACCEPTED_SHA='9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3'
COUNT_SHA='ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9'
KNOWN=['white','yellow_orange','red_pink','blue_purple','green_brown_inconspicuous']

def direct_colours(d):
    return d.loc[d.trait_name.eq('flower_primary_color')&d.quality.isin(['high','medium'])&d.resolution_status.eq('resolved')&~d.classification.isin(['cultivar_contamination','partial_self_compatibility'])&d.normalized_value.isin(KNOWN)].copy()

def aggregate_occurrence_species(d):
    if d.groupby('wcvp_id').white.nunique().gt(1).any():raise ValueError('Conflicting species colour')
    d=d.drop_duplicates(['cell025','year','wcvp_id'])
    c=d.groupby(['cell025','year'],as_index=False).agg(n_species=('wcvp_id','size'),n_white_species=('white','sum'))
    g=d.groupby(['cell025','year','genus'],as_index=False).agg(n_species=('wcvp_id','size'),n_white_species=('white','sum'))
    return c,g

def targets_from_sources(traits,names,accepted,counts):
    for p,s in [(traits,TRAIT_SHA),(names,NAME_SHA),(accepted,ACCEPTED_SHA),(counts,COUNT_SHA)]:assert sha(p)==s
    d=direct_colours(pd.read_csv(traits,dtype=str));assert d.accepted_species.is_unique
    n=pd.read_csv(names,dtype=str);n=n.loc[n.wcvp_match_state.eq('unique_exact_name_to_accepted_id')]
    a=pd.read_csv(accepted,dtype=str)[['wcvp_accepted_plant_name_id','wcvp_accepted_name']].drop_duplicates()
    assert a.wcvp_accepted_plant_name_id.is_unique
    k=pd.read_csv(counts,dtype={'inat_taxon_id':str})
    x=n.merge(a,on='wcvp_accepted_plant_name_id',validate='many_to_one').merge(d,left_on='wcvp_accepted_name',right_on='accepted_species',validate='many_to_one').merge(k[['inat_taxon_id','n_unique_source_photo_ids']],on='inat_taxon_id',validate='one_to_one')
    # Choose the most supported nominal source taxon, without exposure inspection.
    x=x.sort_values(['n_unique_source_photo_ids','inat_taxon_id'],ascending=[False,True]).drop_duplicates('wcvp_accepted_plant_name_id')
    x=x.loc[x.n_unique_source_photo_ids>=100].copy()
    x['wcvp_id']=x.wcvp_accepted_plant_name_id;x['white']=x.normalized_value.eq('white').astype(int)
    x['genus']=x.accepted_species.str.split().str[0]
    assert x.inat_taxon_id.is_unique and x.wcvp_id.is_unique
    return x[['inat_taxon_id','wcvp_id','white','genus','n_unique_source_photo_ids']]

def run(group,traits,names,accepted,counts,out,cache):
    out=Path(out);cache=Path(cache);out.mkdir(parents=True,exist_ok=True);cache.mkdir(parents=True,exist_ok=True)
    targets=targets_from_sources(traits,names,accepted,counts);local=targets.loc[pd.to_numeric(targets.inat_taxon_id).mod(128).floordiv(8).eq(group)]
    for name in ['manifest.json','SHA256SUMS.txt',f'source-{group:02d}.tar']:
        if not (cache/name).exists():subprocess.run(['gh','release','download',TAG,'--repo','zuizui0223/mois','--pattern',name,'--dir',str(cache)],check=True)
    assert sha(cache/'manifest.json')==MANIFEST_SHA and sha(cache/'SHA256SUMS.txt')==SUMS_SHA
    manifest=json.loads((cache/'manifest.json').read_text());shards={v['part']:v for v in manifest['shards']}
    sums={line.split()[1]:line.split()[0] for line in (cache/'SHA256SUMS.txt').read_text().splitlines()}
    archive=cache/f'source-{group:02d}.tar';assert sha(archive)==sums[archive.name]
    selected=[];background=[];total=0;observations=0;precise=0;parts=set();support=[]
    with tarfile.open(archive) as tf:
        for member in tf:
            name=member.name.removeprefix('./')
            if not name.startswith('observations/part-') or not name.endswith('.csv.gz'):continue
            part=int(Path(name).name[5:8]);assert part//8==group and part not in parts;parts.add(part)
            raw=tf.extractfile(member).read();assert hashlib.sha256(raw).hexdigest()==shards[part]['observations_sha256'];n=0
            for d in pd.read_csv(io.BytesIO(raw),compression='gzip',dtype={'observation_uuid':str,'inat_taxon_id':str,'observer_id':str},chunksize=200000):
                assert (pd.to_numeric(d.inat_taxon_id)%128==part).all();n+=len(d)
                s,b=select_chunk(d,set(local.inat_taxon_id));background.append(b);observations+=len(s)
                # Primary period uses one consistent WorldPop Global2 product.
                s=s.loc[s.coordinate_precision_5km & s.year.between(2015,2024)]
                precise+=len(s);s=s.merge(local,on='inat_taxon_id',validate='many_to_one')
                selected.append(s[['cell025','year','wcvp_id','white','genus']].drop_duplicates())
                support.append(s[['wcvp_id','cell025','year']].drop_duplicates())
            assert n==shards[part]['observation_count'];total+=n
    assert parts==set(range(group*8,(group+1)*8))
    d=pd.concat(selected,ignore_index=True).drop_duplicates(['cell025','year','wcvp_id'])
    c,g=aggregate_occurrence_species(d)
    c.to_csv(out/'documented_flower_colour_cell_year_counts.csv.gz',index=False)
    g.to_csv(out/'documented_flower_colour_genus_cell_year_counts.csv.gz',index=False)
    bg=pd.concat(background,ignore_index=True).groupby(['cell025','year'],as_index=False)[['n_observations','n_precise_observations']].sum()
    bg.to_csv(out/'all_source_plants_cell_year_effort.csv.gz',index=False)
    sp=pd.concat(support).drop_duplicates().groupby('wcvp_id').agg(n_site_years=('cell025','size'),n_sites=('cell025','nunique'),n_years=('year','nunique'))
    # No labels or original species-level trait records in this coverage histogram.
    histogram=sp.groupby(['n_sites','n_years']).size().reset_index(name='species_count')
    histogram.to_csv(out/'species_geographic_support_histogram.csv',index=False)
    receipt={'schema':'mois_broad_documented_colour_full_source_v1','group':group,'source_observations_scanned':total,'eligible_accepted_species_global':len(targets),'eligible_white_species_global':int(targets.white.sum()),'source_photo_opportunities_global':int(targets.n_unique_source_photo_ids.sum()),'eligible_species_in_group':len(local),'selected_observations_2014_2024':observations,'precise_selected_observations_2015_2024':precise,'species_with_precise_observations':int(d.wcvp_id.nunique()),'species_cell_year_rows':len(d),'source_trait_sha256':TRAIT_SHA,'source_manifest_sha256':MANIFEST_SHA,'only_aggregate_trait_outputs':True,'no_new_photo_pixels':True,'photo_specific_flowering_verified':False,'missing_species_records_not_absences':True,'outputs':{p.name:sha(p) for p in out.glob('*.csv*')}}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--group',type=int,choices=range(16),required=True)
    for k in ['traits','names','accepted','counts','out','cache']:p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();run(a.group,a.traits,a.names,a.accepted,a.counts,a.out,a.cache)
