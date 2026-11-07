"""Exact 2020 subset of MOIS original global GSHHG point/cell assignments.

Source procedure is copied from the independently verified 2026-11-06
global all-land inventory. Island area and latitude are NEVER filters.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import make_valid

def require(ok, message):
    if not bool(ok):raise ValueError(message)

def area_class(area):
    require(np.isfinite(area) and area>0,'invalid landmass area')
    return 'mainland' if area>7000000 else 'island'

def valid_coordinates(d):
    return d.latitude.between(-90,90)&d.longitude.between(-180,180)

def prepare_land(raw):
    require(raw.crs is not None,'geography CRS missing')
    x=raw[['sibling_id','geometry']].to_crs(4326).copy()
    x['land_id']=pd.to_numeric(x.sibling_id,errors='raise').astype('int64').astype(str)
    bad=~x.is_valid
    if bad.any():x.loc[bad,'geometry']=x.loc[bad].geometry.map(make_valid)
    dup=x.land_id.duplicated(keep=False)
    one=x.loc[~dup,['land_id','geometry']]
    if dup.any():
        multi=x.loc[dup,['land_id','geometry']].dissolve(by='land_id',as_index=False)
        one=gpd.GeoDataFrame(pd.concat([one,multi],ignore_index=True),crs=4326)
    one=one.sort_values('land_id').reset_index(drop=True)
    one['area_km2']=one.to_crs(6933).area/1e6
    one['land_class']=one.area_km2.map(area_class)
    one['area_below_5km2']=one.area_km2.lt(5)
    require(one.land_id.is_unique,'duplicate prepared landmass ID')
    return one[['land_id','area_km2','land_class','area_below_5km2','geometry']]

def assign_points(frame,land,lakes):
    require(frame.unit_id.notna().all() and frame.unit_id.is_unique,'duplicate source unit ID')
    a=frame[['unit_id','longitude','latitude']].reset_index(drop=True)
    good=valid_coordinates(a)
    result=pd.DataFrame({'unit_id':a.unit_id,'land_id':None,'area_km2':np.nan,
         'land_class':'coordinate_missing_or_invalid'})
    points=gpd.GeoDataFrame(a.loc[good],geometry=gpd.points_from_xy(a.loc[good,'longitude'],a.loc[good,'latitude']),crs=4326)
    joined=gpd.sjoin(points,land[['land_id','area_km2','land_class','geometry']],how='left',predicate='within')
    require(not joined.index.duplicated().any(),'overlapping L1 land matches')
    result.loc[joined.index,['land_id','area_km2','land_class']]=joined[['land_id','area_km2','land_class']]
    result.loc[good&result.land_id.isna(),'land_class']='not_inside_L1'
    no=joined.loc[joined.land_id.isna(),points.columns]
    if len(no):
        boundary=gpd.sjoin(no,land[['land_id','geometry']],predicate='intersects',how='inner')
        result.loc[boundary.index.unique(),'land_class']='boundary_unassigned'
    if lakes is not None and len(lakes):
        water=gpd.sjoin(points,lakes[['geometry']],how='inner',predicate='intersects')
        ix=water.index.unique()
        result.loc[ix,['land_id','area_km2']]=[None,np.nan]
        result.loc[ix,'land_class']='inland_water_or_inland_island_unresolved'
    require(result.unit_id.is_unique and len(result)==len(frame),'geography join changed IDs')
    return result

def aggregate_targets(d,grid=.05):
    require(d.photo_id.is_unique,'duplicate photo ID in acquisition source')
    x=d.loc[valid_coordinates(d)&pd.to_numeric(d.year,errors='coerce').notna()].copy()
    x['year']=x.year.astype(int)
    x['gx']=np.minimum(np.floor((x.longitude+180)/grid).astype('int64'),int(360/grid)-1)
    x['gy']=np.minimum(np.floor((x.latitude+90)/grid).astype('int64'),int(180/grid)-1)
    x['cell_id']=x.gx.astype(str)+'_'+x.gy.astype(str)
    x['land_key']=x.land_id.fillna('UNRESOLVED:'+x.land_class)
    keys=['year','cell_id','land_key','land_class']
    anchors=x.sort_values('photo_id',kind='stable').drop_duplicates(keys)
    ans=anchors[keys+['photo_id','longitude','latitude']].rename(columns={
        'photo_id':'photo_id_anchor','longitude':'anchor_longitude','latitude':'anchor_latitude'})
    counts=x.groupby(keys,dropna=False,observed=True).agg(
        n_attempted=('photo_id','size'),n_strict9=('strict9','sum'),
        source_subarms=('source_subarm',lambda z:'|'.join(sorted(set(z.astype(str)))))).reset_index()
    ans=ans.merge(counts,on=keys,validate='one_to_one')
    ans['land_id']=ans.land_key.where(~ans.land_key.str.startswith('UNRESOLVED:'),None)
    ans['target_id']=ans.year.astype(str)+'|'+ans.cell_id+'|'+ans.land_key
    require(ans.target_id.is_unique,'duplicate request key')
    require(ans.n_attempted.sum()==len(x),'source attempts lost from environment queue')
    return ans.sort_values('target_id').reset_index(drop=True)
