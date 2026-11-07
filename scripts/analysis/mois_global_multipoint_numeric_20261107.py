"""Global 2020 VIIRS/WorldPop bulk numeric sampling. No images, ROI, or imputation."""
from __future__ import annotations
import argparse,calendar,csv,datetime as dt,hashlib,json,math,re,statistics
from pathlib import Path
import requests

VIIRS="https://di-lawimagery1.img.arcgis.com/arcgis/rest/services/NighttimeLightsMDNB/ImageServer/getSamples"
POP="https://worldpop.arcgis.com/arcgis/rest/services/WorldPop_Population_Density_1km/ImageServer/getSamples"
MONTHS=(1,4,7,10)

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_BULK_GIS: "+msg)

def digest(data):return hashlib.sha256(data).hexdigest()

def daterange(y,m):
    a=dt.datetime(y,m,1,tzinfo=dt.timezone.utc)
    b=dt.datetime(y,m,calendar.monthrange(y,m)[1],23,59,59,tzinfo=dt.timezone.utc)
    return f'{int(a.timestamp()*1000)},{int(b.timestamp()*1000)}'

def api_params(product,month,rows,year=2020):
    require(product in ("viirs","population") and 1<=len(rows)<=128,"invalid multi-point source request")
    require(year==2020,"WorldPop Global1 has not been verified for years after 2020")
    require(product!="viirs" or month in MONTHS,"bad VIIRS month")
    points=[]
    for r in rows:
        x,y=float(r["query_longitude"]),float(r["query_latitude"])
        require(-180<=x<=180 and -90<=y<=90,"coordinate outside globe")
        points.append([x,y])
    q={"f":"json","geometryType":"esriGeometryMultipoint",
       "geometry":json.dumps({"points":points,"spatialReference":{"wkid":4326}},separators=(',',':')),
       "returnFirstValueOnly":"true","outFields":"*","interpolation":"RSP_NearestNeighbor"}
    q["time"]=daterange(year,month) if product=="viirs" else (
       daterange(year,1).split(",")[0]+","+daterange(year,12).split(",")[1])
    return q

def extract_samples(product,month,rows,response,year=2020):
    require(isinstance(response,dict) and not response.get("error"),"source-level provider error")
    samples=response.get("samples")
    require(isinstance(samples,list) and len(samples)<=len(rows),"invalid sample array")
    indexed={}
    for item in samples:
        i=item.get("locationId")
        require(isinstance(i,int) and 0<=i<len(rows) and i not in indexed,
                "duplicate/invalid returned point identity")
        loc=item.get("location") or {}
        expected=rows[i]
        require(abs(float(loc["x"])-float(expected["query_longitude"]))<=1e-7 and
                abs(float(loc["y"])-float(expected["query_latitude"]))<=1e-7,
                "provider location does not match target")
        source_value=item.get("value")
        token=re.split(r'[,\s]+',str(source_value).strip()) if source_value is not None else [""]
        if token in ([""],["NoData"]):
            indexed[i]={"status":"NO_DATA","value":None,"cloud":None}
            continue
        require(len(token)==(2 if product=="viirs" else 1),"unexpected raw VIIRS/cloud band dimension")
        nums=[float(x) for x in token]
        require(all(math.isfinite(x) for x in nums),"nonfinite band")
        if product=="viirs":
            meta=item.get("attributes") or {}
            require(meta.get("year")==year and meta.get("month")==month,
                    "raster has wrong source year/month")
            cloud=nums[1]
            require(cloud>=0 and abs(cloud-round(cloud))<1e-6,"bad cloud-free count")
            indexed[i]={"status":"OK" if cloud>=2 else "LOW_CLOUD_COUNT",
                        "value":nums[0] if cloud>=2 else None,"cloud":int(round(cloud))}
        else:
            indexed[i]={"status":"OK" if nums[0]>=0 else "NO_DATA",
                        "value":nums[0] if nums[0]>=0 else None,"cloud":None}
    return [indexed.get(i,{"status":"NO_DATA","value":None,"cloud":None})
            for i in range(len(rows))]

def validate_targets(path):
    rows=list(csv.DictReader(Path(path).open(newline="")))
    require(1<=len(rows)<=128 and len({r["target_id"] for r in rows})==len(rows),
            "unexpected target cohort or duplicate targets")
    require(all(r["target_id"].startswith("2020|") and
                r["land_class"] in ("mainland","island") for r in rows),
            "target has wrong source year or land category")
    require(len({(float(r["query_longitude"]),float(r["query_latitude"])) for r in rows})==len(rows),
            "two targets refer to the same GIS pixel request coordinates")
    return rows

def run(targets,out,offline=False):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rows=validate_targets(targets)
    merged={r["target_id"]:dict(r) for r in rows}
    receipt={"schema":"MOIS_2020_GIS_multipoint_strict_noarea_v1",
             "source_manifest_sha256":digest(Path(targets).read_bytes()),
             "target_rows":len(rows),"request_count_new":0,
             "response_sha256":{},"offline_source_replay":offline,
             "no_original_photo_downloads_or_ROI_inference":True,
             "no_area_or_latitude_cutoff":True,
             "population2020_API_response_year_not_echoed":True,
             "VIIRS_is_four_months_not_official_annual":True}
    for product,month in [("viirs",m) for m in MONTHS]+[("population",0)]:
        params=api_params(product,month,rows)
        filename=f'{product}_2020_{month:02d}.response'
        file=out/filename
        if offline:
            require(file.is_file(),"offline replay lacks original provider bytes")
            payload=file.read_bytes()
        else:
            url=VIIRS if product=="viirs" else POP
            reply=requests.get(url,params=params,timeout=65,
                  headers={"User-Agent":"MOIS global numerical bulk source QA/1.0"})
            payload=reply.content
            file.write_bytes(payload)
            require(reply.status_code==200,"provider HTTP "+str(reply.status_code))
            receipt["request_count_new"]+=1
        receipt["response_sha256"][filename]=digest(payload)
        try:body=json.loads(payload)
        except ValueError as exc:raise ValueError("invalid provider JSON") from exc
        values=extract_samples(product,month,rows,body)
        for row,val in zip(rows,values):
            x=merged[row["target_id"]]
            col="viirs_"+str(month) if product=="viirs" else "population_global1_0"
            x[col+"_status"]=val["status"]
            x[col+"_value"]=val["value"]
            if product=="viirs":x[col+"_cloudfree"]=val["cloud"]
    output=[]
    for row in rows:
        x=merged[row["target_id"]]
        good=[x[f"viirs_{m}_value"] for m in MONTHS if x[f"viirs_{m}_status"]=="OK"]
        pop=x["population_global1_0_value"] if x["population_global1_0_status"]=="OK" else None
        x["viirs_valid_months"]=len(good)
        x["viirs_median"]=statistics.median(good) if len(good)>=3 else None
        x["population_density"]=pop
        x["complete"]=len(good)>=3 and pop is not None
        output.append(x)
    from io import StringIO
    f=StringIO(newline='')
    writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for x in output for k in x)),
                          lineterminator='\n')
    writer.writeheader();writer.writerows(output)
    raw=f.getvalue().encode()
    filename="replay_values.csv" if offline else "values.csv"
    (out/filename).write_bytes(raw)
    receipt.update({"targets_complete":sum(x["complete"] for x in output),
        "targets_quality_incomplete":sum(not x["complete"] for x in output),
        "complete_island":sum(x["complete"] and x["land_class"]=="island" for x in output),
        "complete_mainland":sum(x["complete"] and x["land_class"]=="mainland" for x in output),
        "valid_numeric_values":sum(x[c+"_status"]=="OK" for x in output for c in
             [f"viirs_{m}" for m in MONTHS]+["population_global1_0"]),
        "output_values_sha256":digest(raw)})
    (out/("replay_receipt.json" if offline else "receipt.json")).write_text(
         json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_MULTIPOINT_2020",json.dumps({
         k:receipt[k] for k in ("target_rows","targets_complete","targets_quality_incomplete",
         "complete_island","complete_mainland","valid_numeric_values","request_count_new",
         "offline_source_replay")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--targets",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--offline",action="store_true")
    args=p.parse_args()
    run(args.targets,args.out,args.offline)
