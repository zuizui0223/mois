#!/usr/bin/env python3
"""DEM-only terrain metrics for later mois validation; never infer from point elevation.

Metric-projected, north-up raster DEM must be supplied independently. Original
FCP source has photo-site point elevation only, not a DEM neighborhood.
This tool aggregates numeric terrain summaries without exporting coordinates.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np

def sha(path):
 h=hashlib.sha256()
 with open(path,"rb") as f:
  for b in iter(lambda:f.read(1<<20),b""):h.update(b)
 return h.hexdigest()

def terrain3x3(z,dx_m,dy_m):
 """Horn 3x3 gradient; geographic degrees are not valid metric cell sizes.

 Aspect is the DOWNSLOPE azimuth clockwise from north, undefined on flat cells.
 TRI follows sqrt(sum of squared center–eight-neighbor elevation differences).
 TPI is center minus mean of eight surrounding cells. Both depend on DEM grain.
 """
 q=np.asarray(z,float)
 if q.shape!=(3,3) or not np.isfinite(q).all():raise ValueError("Finite three-by-three DEM elevation neighbourhood required")
 if not (dx_m>0 and dy_m>0):raise ValueError("Projected, positive metric pixel scales required")
 dx=((q[0,2]+2*q[1,2]+q[2,2])-(q[0,0]+2*q[1,0]+q[2,0]))/(8*dx_m)
 dy=((q[2,0]+2*q[2,1]+q[2,2])-(q[0,0]+2*q[0,1]+q[0,2]))/(8*dy_m)
 slope=math.degrees(math.atan(math.hypot(dx,dy)))
 around=np.delete(q.reshape(-1),4)
 tpi=float(q[1,1]-around.mean())
 tri=float(np.sqrt(np.square(around-q[1,1]).sum()))
 aspect=None if math.hypot(dx,dy)<1e-12 else math.degrees(math.atan2(-dx,dy))%360
 return {"slope_deg":slope,"aspect_downhill_cw_from_north_deg":aspect,"TPI_m":tpi,"TRI_Riley_m":tri}

def self_test():
 flat=terrain3x3(np.ones((3,3))*100,5,5)
 assert flat["slope_deg"]==0 and flat["TPI_m"]==0 and flat["TRI_Riley_m"]==0 and flat["aspect_downhill_cw_from_north_deg"] is None
 ramp=np.array([[0,1,2],[0,1,2],[0,1,2]],float)
 a=terrain3x3(ramp,1,1)
 assert abs(a["slope_deg"]-45)<1e-12 and abs(a["aspect_downhill_cw_from_north_deg"]-270)<1e-12 and abs(a["TPI_m"])<1e-12
 summit=np.zeros((3,3));summit[1,1]=10
 h=terrain3x3(summit,1,1)
 assert h["TPI_m"]==10 and abs(h["TRI_Riley_m"]-math.sqrt(800))<1e-12 and h["slope_deg"]==0
 valley=-summit;v=terrain3x3(valley,1,1)
 assert v["TPI_m"]==-10 and v["TRI_Riley_m"]==h["TRI_Riley_m"]
 try:terrain3x3(np.zeros((3,3)),0,1)
 except ValueError:pass
 else:raise AssertionError("Nonmetric DEM pixel spacing was erroneously accepted")
 print("PASS synthetic projected DEM: flat, aspect, gradient, summit, valley, TRI and TPI")

def inspect(dem,sites):
 import rasterio
 from rasterio.windows import Window
 from pyproj import Transformer
 with rasterio.open(dem) as ds:
  if not ds.crs or not ds.crs.is_projected:
   raise ValueError("Reproject DEM to a distance-preserving LOCAL projected metre CRS; unprojected EPSG:4326 is not accepted")
  unit=(ds.crs.linear_units or "").lower()
  if unit not in ("metre","meter","metres","meters","m"):raise ValueError("Projected DEM CRS horizontal units are not metres")
  t=ds.transform
  if abs(t.b)>1e-10 or abs(t.d)>1e-10:raise ValueError("Rectilinear projected raster without rotated/sheared pixels required")
  dx,dy=abs(t.a),abs(t.e)
  if min(dx,dy)<1 or max(dx,dy)>10000:raise ValueError("DEM pixel grain outside supported 1 m–10 km range")
  to_local=Transformer.from_crs("EPSG:4326",ds.crs,always_xy=True)
  out=[];requested=0;reasons={"outside_or_edge":0,"nodata_or_nonfinite":0,"missing_coordinates":0}
  with open(sites,newline="",encoding="utf8") as f:
   reader=csv.DictReader(f)
   if not {"latitude","longitude"}.issubset(reader.fieldnames or []):raise ValueError("Private site CSV latitude and longitude required")
   for row in reader:
    requested+=1
    try:lat=float(row["latitude"]);lon=float(row["longitude"])
    except (KeyError,ValueError,TypeError):
     reasons["missing_coordinates"]+=1;continue
    if not -90<=lat<=90 or not -180<=lon<=180:
     reasons["missing_coordinates"]+=1;continue
    x,y=to_local.transform(lon,lat)
    if not np.isfinite(x) or not np.isfinite(y):
     reasons["missing_coordinates"]+=1;continue
    rr,cc=ds.index(x,y)
    if rr<1 or cc<1 or rr>=ds.height-1 or cc>=ds.width-1:
     reasons["outside_or_edge"]+=1;continue
    window=ds.read(1,window=Window(cc-1,rr-1,3,3),masked=True)
    if window.shape!=(3,3) or np.any(np.ma.getmaskarray(window)) or not np.isfinite(np.asarray(window)).all():
     reasons["nodata_or_nonfinite"]+=1;continue
    out.append(terrain3x3(np.asarray(window),dx,dy))
  if not out:raise ValueError("No valid 3x3 sites in independent DEM: cannot claim terrain was measured")
  def summary(key):
   v=np.array([row[key] for row in out])
   return {"median":float(np.median(v)),"p05":float(np.quantile(v,.05)),"p95":float(np.quantile(v,.95))}
  return {"status":"DEM_TERRAIN_MEASURED_FOR_PRIVATE_POINTS_NOT_FCP_GENUS_ASSOCIATION",
   "dem_sha256":sha(dem),"input_site_file_sha256":sha(sites),
   "dem_crs":str(ds.crs),"dem_grain_m":[dx,dy],
   "n_requested_sites":requested,"n_usable_site_neighbourhoods":len(out),"skipped_by_reason":reasons,
   "slope_deg":summary("slope_deg"),"TPI_m":summary("TPI_m"),"TRI_Riley_m":summary("TRI_Riley_m"),
   "n_defined_aspects":sum(row["aspect_downhill_cw_from_north_deg"] is not None for row in out),
   "boundary":"Original photo/site identities and precise positions not exported. This describes supplied terrain only; tests of colour×terrain with genus/species requires an explicitly linked and licensed source sample plus multiscale sensitivity."}

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--self-test",action="store_true")
 p.add_argument("--dem",type=Path)
 p.add_argument("--sites-csv",type=Path)
 p.add_argument("--output",type=Path)
 a=p.parse_args()
 if a.self_test:self_test();return
 if not (a.dem and a.sites_csv and a.output):
  p.error("Real terrain requires all three: --dem metric-projected GeoTIFF, --sites-csv private WGS84 points, --output aggregate JSON")
 r=inspect(a.dem,a.sites_csv)
 a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
 print("WROTE source-safe DEM-only terrain aggregate; original FCP colour/terrain relationship remains UNTESTED")

if __name__=="__main__":main()
