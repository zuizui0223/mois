"""Same-year coarse climate covariates from NOAA PSL NCEP/DOE Reanalysis II.

Nearest native Gaussian-grid location, not a 0.25-degree climate estimate.
Monthly mean temperature is weighted by days; precipitation rate integrated.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd
import xarray as xr

BASE='https://psl.noaa.gov/thredds/fileServer/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid/'

def annual(values,dates,kind):
    dates=pd.DatetimeIndex(dates)
    if len(dates)!=12 or set(dates.month)!=set(range(1,13)) or len(set(dates.year))!=1:raise ValueError('Annual climate requires exactly twelve months of one year')
    days=dates.days_in_month.to_numpy()[:,None]
    valid=np.isfinite(values).all(axis=0)
    if kind=='temperature':out=(values*days).sum(axis=0)/days.sum()-273.15
    else:
        valid&=(values>=0).all(axis=0)
        out=(values*days*86400).sum(axis=0)
    return np.where(valid,out,np.nan)

def run(queue,root,out):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    q=pd.read_csv(queue);assert not q.duplicated(['cell025','year']).any()
    xy=q.cell025.str.split('_',expand=True).astype(int)
    lon=((xy[0]+.5)*.25-180).to_numpy();lat=((xy[1]+.5)*.25-90).to_numpy()
    receipt={'product':'NCEP/DOE Reanalysis II provided by NOAA PSL, Boulder, Colorado, USA','documentation':'https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html','sampling':'nearest native Gaussian-grid location; no spatial downscaling','sources':{}}
    for file,var,kind,column in [('air.2m.mon.mean.nc','air','temperature','temperature_c'),('prate.sfc.mon.mean.nc','prate','precipitation','precipitation_mm')]:
        path=root/file;ds=xr.open_dataset(path);a=ds[var].squeeze(drop=True)
        units=a.attrs.get('units','')
        if kind=='temperature':assert units in ['K','degK']
        else:assert units in ['kg/m^2/s','kg m-2 s-1','Kg/m^2/s']
        yl=np.argmin(np.abs(ds.lat.values[:,None]-lat[None,:]),axis=0)
        xl=np.argmin(np.abs((ds.lon.values[:,None]-lon[None,:]+180)%360-180),axis=0)
        values=np.full(len(q),np.nan)
        for year in sorted(q.year.unique()):
            rows=np.flatnonzero(q.year.to_numpy()==year)
            subset=a.sel(time=str(year)).transpose('time','lat','lon')
            grid=subset.to_numpy();v=grid[:,yl[rows],xl[rows]]
            values[rows]=annual(v,subset.time.values,kind)
        q[column]=values
        q[kind+'_native_latitude']=ds.lat.values[yl];q[kind+'_native_longitude']=ds.lon.values[xl]
        receipt['sources'][file]={'url':BASE+file,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'units':units,'n_latitudes':len(ds.lat),'n_longitudes':len(ds.lon),'last_month':str(ds.time.values[-1])}
        ds.close()
    q.to_csv(out/'sameyear_coarse_climate.csv',index=False)
    receipt['rows']=len(q);receipt['complete_rows']=int(q[['temperature_c','precipitation_mm']].notna().all(axis=1).sum())
    (out/'climate_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--queue',required=True);p.add_argument('--root',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.queue,a.root,a.out)
