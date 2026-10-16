"""Synthetic-only 12-raster geospatial and wet/dry tests, no original pixels."""
import csv
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    'mois_monthly_normals',ROOT/'scripts/analysis/mois_monthly_rainfall_flower_colour_20261016.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def raster_fixture(tmp):
    monthly=[]
    transform=from_origin(-180,90,60,45)
    for month in range(1,13):
        arr=np.full((4,6),float(10 if month<=6 else 100),dtype=np.float32)
        arr[:,1]=55
        if month==6:
            arr[3,0]=-9999
        path=tmp/f'wc2.1_10m_prec_{month:02}.tif'
        with rasterio.open(path,'w',driver='GTiff',height=4,width=6,count=1,
                           dtype='float32',crs='EPSG:4326',transform=transform,nodata=-9999) as ds:
            ds.write(arr,1)
        monthly.append(path)
    manifest=tmp/'monthly.json'
    manifest.write_text(json.dumps({
        'schema':m.MANIFEST_SCHEMA,'source':m.NORMALS,
        'units':'mm/month','resolution':'10m','months':{
            f'{i:02}':{'path':path.name,'sha256':m.sha256_file(path)}
            for i,path in enumerate(monthly,1)}}))
    return m.MonthlyNormals(manifest,m.sha256_file(manifest),test_small=True),manifest,monthly


def season_fixture(tmp):
    photo=tmp/'source_season.csv'
    ids=[]
    data=[
        ('1','Alpha','40','35','04','DRY','white',3),
        ('1','Alpha','40','35','10','WET','yellow_orange',3),
        ('2','Beta','-40','35','04','DRY','blue_purple',3),
        ('2','Beta','-40','35','10','WET','white',3),
        ('3','Gamma','10','-90','04','ASEASONAL','white',1)
    ]
    for tax,genus,lat,lon,month,label,hue,n in data:
        for j in range(n):
            ids.append({
                'measurement_id':f'MEASURE-{len(ids)}',
                'inat_taxon_id':tax,'genus':genus,'observation_id':f'uuid-{len(ids)}',
                'observer_id':f'observer_{j}',
                'latitude':lat,'longitude':lon,
                'bio12_mm':'660',
                'month':month,
                'domain':('NORTH_EXTRATROPICAL' if float(lat)>23.5 else
                          'SOUTH_EXTRATROPICAL' if float(lat)<-23.5 else 'TROPICAL'),
                'local_season':'SPRING','annotation_status':'NOT_VERIFIED',
                'region_stratum':'G_N1','precip_stratum':'P_MEDIUM',
                'sample_observation_weight':'1','sample_photo_weight':'1',
                'morph':hue,'measurement_status':'classified_four_state_morph'
            })
    with photo.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=m.FIELDS)
        w.writeheader()
        w.writerows(ids)
    receipt=tmp/'source_receipt.json'
    receipt.write_text('{"source":"synthetic"}')
    private=tmp/'private_season_support.json'
    private.write_text(json.dumps({
        'schema':m.SOURCE_SEASON_SCHEMA,
        'private_derived_photo_rows_sha256':m.sha256_file(photo),
        'source_sample_receipt_sha256':m.sha256_file(receipt),
        'source_join_sha256':'x'*64}))
    public=tmp/'public_season.json'
    public.write_text(json.dumps({
        'schema':m.SOURCE_PUBLIC_SCHEMA,
        'source_sealed_roi_join_sha256':'x'*64,
        'synthetic_test_only':True,
        'old_genus_13_environment_estimates_updated':False,
        'all_13_photo_site_environment_predictors_available_here':False,
        'selected_source_observations':len(ids),
        'source_classifiable_photographs':len(ids)}))
    return photo,private,public,receipt,len(ids)


def run(tmp,monthly,**kwargs):
    src,priv,pub,receipt,_=season_fixture(tmp)
    return m.audit(
        src,priv,pub,receipt,monthly,
        tmp/'monthly_photos.csv',tmp/'private_monthly.json',tmp/'public_monthly.json',
        min_classifiable=3,min_observers=2,
        min_public=kwargs.get('min_public',1),synthetic=True)


def test_wet_dry_is_local_site_not_arbitrary_northern_calendar():
    profile=[10.]*6+[100.]*6
    assert m.rainfall_label(profile,4)['label']=='DRY'
    assert m.rainfall_label(profile,10)['label']=='WET'
    assert m.rainfall_label(profile,None)['label']=='UNKNOWN_MONTH'
    assert m.rainfall_label([55.]*12,10)['label']=='ASEASONAL'
    assert m.rainfall_label([0.]*12,10)['label']=='ASEASONAL'
    assert m.rainfall_label(None,10)['label']=='CLIMATE_UNAVAILABLE'
    assert m.rainfall_label(profile,10)['annual_mm']==660
    with pytest.raises(ValueError,match='twelve'):
        m.rainfall_label([100.]*11,3)


def test_raster_hash_and_same_species_same_climate_cell_match(tmp_path):
    normals,_,_=raster_fixture(tmp_path)
    doc=run(tmp_path,normals)
    assert doc['n_source_selected_original_observations']==13
    assert doc['n_sites_with_comparable_annual_bio12']==13
    assert doc['local_label_totals']=={'ASEASONAL':1,'DRY':6,'WET':6}
    assert doc['n_classifiable_source_photos']==13
    pairs={(z['comparison'],z['domain']):z for z in doc['matched_same_species_and_site_wet_dry']}
    assert len(pairs)==4
    assert pairs[('same_species','NORTH_EXTRATROPICAL')][
        'equal_unit_mean_white_difference_wet_minus_dry']==-1
    assert pairs[('same_species_same_worldclim_cell','SOUTH_EXTRATROPICAL')][
        'equal_unit_mean_white_difference_wet_minus_dry']==1
    assert doc['global_13_environment_regression_not_run'] is True
    assert doc['normals_are_not_observed_weather_at_photo_date'] is True
    public=(tmp_path/'public_monthly.json').read_text()
    assert 'uuid-' not in public and 'observation_id' not in public
    assert 'inat_taxon_id' not in public
    assert doc['synthetic_test_only'] is True
    with pytest.raises(ValueError,match='Refuse overwrite'):
        run(tmp_path,normals)


def test_annual_bio12_raster_version_mismatch_aborts_before_public_results(tmp_path):
    normal,_,_=raster_fixture(tmp_path)
    src,priv,pub,receipt,_=season_fixture(tmp_path)
    src.write_text(src.read_text().replace(',660,',',950,'))
    data=json.loads(priv.read_text())
    data['private_derived_photo_rows_sha256']=m.sha256_file(src)
    priv.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='do not reconstruct original site BIO12'):
        m.audit(src,priv,pub,receipt,normal,
                tmp_path/'a.csv',tmp_path/'b.json',tmp_path/'c.json',
                synthetic=True)


def test_changed_monthly_raster_digest_is_rejected(tmp_path):
    _,manifest,paths=raster_fixture(tmp_path)
    paths[0].write_bytes(paths[0].read_bytes()+b'tamper')
    with pytest.raises(ValueError,match='GeoTIFF SHA256'):
        m.MonthlyNormals(manifest,m.sha256_file(manifest),test_small=True)


def test_suppressed_small_cells_and_no_fake_seasonal_result(tmp_path):
    norms,_,_=raster_fixture(tmp_path)
    r=run(tmp_path,norms,min_public=20)
    assert all('suppressed_below' in z for z in r['domain_monthly_rainfall_conditioned_photo_white'])
    assert all('privacy_suppressed_matching_units_below' in z for z in r['matched_same_species_and_site_wet_dry'])
    assert 'n_matched_units_private_before_suppression' not in r
