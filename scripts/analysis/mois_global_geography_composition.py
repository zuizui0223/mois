"""Exploratory MOIS flower-photo geography; species turnover is part of the response.
GSHHG source matches island's corrected coastline. Unknown hues stay missing.
"""
from pathlib import Path
import argparse, hashlib, json, zipfile
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import statsmodels.formula.api as smf
hues=('white','yellow_orange','red_pink','blue_purple')
SHA_GSHHG='8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41'
P='primary313535_fullfourth4096_original_photo_hues.parquet'
E='primary313535_fullfourth4096_annual_worldclim.parquet'
S='primary313535_original_partial_soil_terrain.parquet'
M='primary313535_fullfourth4096_featurelake_manifest.json'
SEEDS={'africa':(20,0),'eurasia':(80,50),'north_america':(-100,40),'south_america':(-60,-15),'australia':(135,-25)}
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def require(x,message):
    if not x:raise ValueError(message)
def read_core(root):
    rec=json.loads((root/M).read_text());frames=[]
    for f in (P,E,S):
        require(digest(root/f)==rec['files_sha256'][f],f'Source checksum differs: {f}')
        frames.append(pd.read_parquet(root/f))
    p,e,s=frames
    require(len(p)==313535 and p.photo_id.is_unique,'Unexpected measured core')
    require(p.photo_id.tolist()==e.photo_id.tolist()==s.photo_id.tolist(),'Photo environment IDs differ')
    require(p.morph.isin(hues).sum()==151239 and p.morph.eq('white').sum()==71427,'Frozen colour outcomes changed')
    d=p.copy()
    for dst,src in [('heat','wc2_bio05_10m'),('rain','wc2_bio12_10m'),('solar','annual_solar_Wm2')]:d[dst]=pd.to_numeric(e[src],errors='coerce')
    d['elevation']=pd.to_numeric(s['wc2_elev_2p5m'],errors='coerce')
    for c in ('latitude','longitude'):d[c]=pd.to_numeric(d[c],errors='coerce')
    return d,rec
def load_land(archive,out):
    require(digest(archive)==SHA_GSHHG,'Not the coastline source used by island')
    with zipfile.ZipFile(archive) as z:
        for n in z.namelist():
            if Path(n).name.startswith(('GSHHS_h_L1.','GSHHS_h_L5.')):(out/Path(n).name).write_bytes(z.read(n))
    land=gpd.read_file(out/'GSHHS_h_L1.shp')
    require({'sibling_id','area','geometry'}<=set(land.columns),'Unexpected GSHHG schema')
    ids=set()
    for name,(lon,lat) in SEEDS.items():
        hit=land[land.geometry.covers(Point(lon,lat))]
        require(len(hit)>0 and hit.sibling_id.nunique()==1,f'Continental seed not unique: {name}')
        ids.add(hit.sibling_id.iloc[0])
    return land,ids
def classify_land(d,land,continental_ids):
    valid=d.latitude.between(-90,90)&d.longitude.between(-180,180)
    points=gpd.GeoDataFrame(d.loc[valid,['latitude','longitude']].copy(),geometry=gpd.points_from_xy(d.loc[valid,'longitude'],d.loc[valid,'latitude']),crs=4326)
    # GeoPandas within uses inverse contains and indexes points. Unlike point-first
    # intersects, this avoids repeatedly testing huge coastline vertices per point.
    # Exact boundary locations are deliberately unassigned, rather than snapped.
    hits=gpd.sjoin(points,land[['sibling_id','area','geometry']],how='left',predicate='within')
    nu=hits.groupby(level=0).sibling_id.nunique();ambiguous=nu[nu>1].index
    hits=hits[~hits.index.isin(ambiguous)];hits=hits[~hits.index.duplicated(keep='first')]
    out=pd.DataFrame(index=d.index);out['land_type']='unassigned';out['landmass']='';out['area_km2']=np.nan
    known=hits[hits.sibling_id.notna()]
    out.loc[known.index,'land_type']=np.where(known.sibling_id.isin(continental_ids),'mainland','island')
    out.loc[known.index,'landmass']=known.sibling_id.astype(str)
    areas=land.groupby('sibling_id')['area'].sum()
    out.loc[known.index,'area_km2']=known.sibling_id.map(areas).astype(float)
    return out
def community_table(d):
    x=d[d.morph.isin(hues)&d.land_type.isin(['island','mainland'])].copy()
    x=x[x.inat_taxon_id.astype(str).ne('')]
    x['ilat']=np.floor(x.latitude).astype(int);x['ilon']=np.floor(x.longitude).astype(int)
    for h in hues:x['p_'+h]=x.morph.eq(h).astype(float)
    key=['landmass','land_type','ilat','ilon']
    fields=['p_'+h for h in hues]+['heat','rain','solar','elevation','area_km2','latitude','longitude']
    tax=x.groupby(key+['inat_taxon_id'],observed=True).agg({**{f:'mean' for f in fields},'photo_id':'size'})
    cells=tax.groupby(level=key).agg({**{f:'mean' for f in fields},'photo_id':'sum'})
    cells['n_taxa']=tax.groupby(level=key).size();cells=cells.rename(columns={'photo_id':'n_photos'}).reset_index()
    cells['region']=(np.floor((cells.latitude+90)/10).astype(int).astype(str)+'_'+np.floor((cells.longitude+180)/30).astype(int).astype(str))
    cells['block']=(np.floor((cells.latitude+90)/10).astype(int).astype(str)+'_'+np.floor((cells.longitude+180)/10).astype(int).astype(str))
    cells['island']=cells.land_type.eq('island').astype(int)
    return cells
def summarize(d,label):
    known=d[d.morph.isin(hues)]
    return {'group':label,'attempts':len(d),'classified':len(known),'unknown':len(d)-len(known),'taxa_classified':int(known.inat_taxon_id.nunique()),**{h+'_count':int(known.morph.eq(h).sum()) for h in hues},**{h+'_fraction':float(known.morph.eq(h).mean()) if len(known) else None for h in hues}}
def fit_contrast(cells):
    c=cells[(cells.n_taxa>=5)&(cells.n_photos>=10)].dropna(subset=['heat','rain','solar','elevation']).copy()
    both=c.groupby('region').island.nunique();c=c[c.region.isin(both[both==2].index)].copy()
    support={'cells':len(c),'island_cells':int(c.island.sum()),'mainland_cells':int((1-c.island).sum()),'geographic_strata':c.region.nunique(),'spatial_clusters':c.block.nunique(),'photos':int(c.n_photos.sum())}
    if len(c)<50 or c.block.nunique()<10:return [],support,c
    for name in ('heat','solar','elevation'):
        sd=c[name].std();require(sd>0,f'Constant predictor {name}');c['z_'+name]=(c[name]-c[name].mean())/sd
    lograin=np.log1p(c.rain.clip(lower=0));c['z_rain']=(lograin-lograin.mean())/lograin.std()
    result=[]
    for label,f in [('geography_only','p_white ~ island + C(region)'),('geography_climate_altitude','p_white ~ island + z_heat + z_rain + z_solar + z_elevation + C(region)')]:
        model=smf.ols(f,c).fit(cov_type='cluster',cov_kwds={'groups':c.block},use_t=True);ci=model.conf_int().loc['island']
        result.append({'model':label,'island_minus_mainland_pp':float(100*model.params['island']),'nominal_spatial_cluster_95_low_pp':float(100*ci.iloc[0]),'nominal_spatial_cluster_95_high_pp':float(100*ci.iloc[1]),'nominal_p':float(model.pvalues['island']),'R2_in_sample':float(model.rsquared),**support})
    return result,support,c
def run(core,coast,out):
    require(not out.exists() or not any(out.iterdir()),'Output already exists');out.mkdir(parents=True,exist_ok=True)
    d,source=read_core(core);print('MOIS_GEO_CORE_LOADED',len(d),flush=True)
    land,continents=load_land(coast,out);print('MOIS_GEO_LAND_LOADED',len(land),flush=True)
    geo=classify_land(d,land,continents);d=pd.concat([d,geo],axis=1)
    print('MOIS_GEO_ASSIGNMENTS',d.land_type.value_counts().to_dict(),flush=True)
    raw=[summarize(d,'all_original')]+[summarize(g,k) for k,g in d.groupby('land_type')]
    for label,mask in [('islands_under_1000km2',d.land_type.eq('island')&d.area_km2.between(5,1000)),('islands_1000_100000km2',d.land_type.eq('island')&d.area_km2.between(1000,100000,inclusive='left')),('islands_100000km2_plus',d.land_type.eq('island')&d.area_km2.ge(100000)),('north_arctic_circle',d.latitude.ge(66.562)),('high_elevation_2500m_plus',d.elevation.ge(2500)),('low_elevation_0_500m',d.elevation.between(0,500,inclusive='left'))]:raw.append(summarize(d[mask],label))
    pd.DataFrame(raw).to_csv(out/'photo_hue_geography.csv',index=False)
    cells=community_table(d);cells.to_csv(out/'species_equal_1degree_landmass_cells.csv.gz',index=False)
    complete=community_table(d.dropna(subset=['heat','rain','solar','elevation']))
    fits,support,matched=fit_contrast(complete)
    pd.DataFrame(fits).to_csv(out/'island_mainland_adjusted.csv',index=False);matched.to_csv(out/'matched_island_mainland_cells.csv',index=False)
    cdesc=[]
    for gate in [1,5,10]:
        for key,g in cells[cells.n_taxa>=gate].groupby('land_type'):cdesc.append({'group':key,'minimum_taxa':gate,'cells':len(g),'species_cell_white_mean':float(g.p_white.mean()),'photos':int(g.n_photos.sum()),'taxon_cells':int(g.n_taxa.sum())})
    pd.DataFrame(cdesc).to_csv(out/'cell_equal_composition.csv',index=False)
    elev=[]
    for left,right in [(0,500),(500,1500),(1500,2500),(2500,10000)]:elev.append(summarize(d[d.elevation.between(left,right,inclusive='left')],f'elevation_{left}_{right}m'))
    pd.DataFrame(elev).to_csv(out/'altitude_photo_hues.csv',index=False)
    d[['photo_id','land_type','landmass','area_km2']].to_csv(out/'photo_landmass_cache.csv.gz',index=False)
    decision={'schema':'mois_island_mainland_macro_photo_geography_exploratory_v1','source_photo_ids':313535,'gshhg_source_sha256':SHA_GSHHG,'core_photo_sha256':source['files_sha256'][P],'mainland_sibling_ids':sorted(str(x) for x in continents),'raw':raw,'community':cdesc,'adjusted':fits,'matched_support':support,'inferences_not_identified':['actual pollinator scarcity','island-specific pigment selection','within-species changes','universal alpine response','flowering-individual abundance','UV or biochemical pigment cost','complete biogeographic realm matching'],'design':'Exploratory, frozen core. Species-equal 1degree landmass cells. Adjusted model uses photo-level complete cases first, and shared 10deg latitude x 30deg longitude strata; nominal 10deg spatial-cluster t intervals. No phylogeny or observer adjustment yet. Ocean/boundary/unresolved locations not defaulted to mainland. High-elevation threshold is not alpine-biome membership.','output_sha256':{p.name:digest(p) for p in out.glob('*.csv*')}}
    (out/'result.json').write_text(json.dumps(decision,indent=2,ensure_ascii=False,allow_nan=False))
    print('MOIS_ISLAND_MAINLAND_RESULT',json.dumps(decision,ensure_ascii=False),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--core',type=Path,required=True);ap.add_argument('--coast',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();run(a.core,a.coast,a.output)
