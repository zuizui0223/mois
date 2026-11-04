#!/usr/bin/env python3
"""Source-locked ROI9 prediction of independently described green/brown species colour.

Out-of-genus, out-of-family and nonshared-species region holdouts. Outputs only
aggregates, never third-party rights-restricted trait rows or photo pixels.
Not a botanically verified flower contrast measure or MNAR sensitivity.
"""
import argparse,io,json,zipfile,hashlib
from pathlib import Path
import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss

FILES={
 "regional":("taxon_region_palette.csv.gz","17b76825f9fa8958a850e2308bb3181d092335ae8a15a64ad72dcebaf5f64963","5f9dee3b0ab41e42d0b1b1dd06650eda2b954b89cb118ddc8bcdd1a813da1479"),
 "taxonomy":("mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz","c10d8fc7bf39b0cb404d5eaa09e3e520596b8d1a9284e643dd007d5b485a93cf","dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2"),
 "accepted":("accepted_family_species_region.csv.gz","3eb62f10607c13893e410dc1b7abdf9f90d4f1e730b08edecb245669823d1c63","a138ab56c0233b769f41c56c7eee7df82f09365d052983c1460f74f8cbf881f0"),
 "direct":("direct_species_trait_ledger.csv.gz","5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd","7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577")
}
COMP=tuple("strict9_"+s for s in ("white","yellow","orange","red","pink","magenta","purple","blue","bronze"))
LABELS={"green_brown_inconspicuous":1,"white":0,"yellow_orange":0,"red_pink":0,"blue_purple":0}
ORDER={s:i for i,s in enumerate(("arctic","north_subpolar","north_temperate","south_temperate","island_cached","mainland_cached"))}

def ensure(ok,message):
    if not bool(ok):raise ValueError("MOIS_VALIDITY: "+message)

def sha(b):return hashlib.sha256(b).hexdigest()

def read_zip(path,key):
    filename,expected,inside=FILES[key]
    path=Path(path)
    ensure(sha(path.read_bytes())==expected,"ZIP source hash changed: "+key)
    with zipfile.ZipFile(path) as z:data=z.read(filename)
    ensure(sha(data)==inside,"Internal source hash changed: "+key)
    return pd.read_csv(io.BytesIO(data),compression="gzip",dtype=str,low_memory=False)

def build_tables(r,t,a,d):
    ensure(len(r)==30634 and not r.duplicated(["sector","inat_taxon_id"]).any(),"ROI9 regional source changed")
    ensure(set(r.sector)==set(ORDER),"ROI9 region list changed")
    ensure(t.inat_taxon_id.is_unique,"Taxon crosswalk identity drift")
    tax=t.loc[t.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
              ["inat_taxon_id","wcvp_accepted_plant_name_id"]]
    acc=a[["accepted_id","accepted_name","genus","family"]].drop_duplicates()
    ensure(acc.accepted_id.is_unique and acc.accepted_name.is_unique and len(acc)==22712,"WCVP accepted species drift")
    j=r.merge(tax,on="inat_taxon_id",validate="many_to_one").merge(
        acc,left_on="wcvp_accepted_plant_name_id",right_on="accepted_id",validate="many_to_one")
    ensure(len(j)==30158 and j.accepted_id.nunique()==22712,"Taxonomic join support drift")
    j["n"]=pd.to_numeric(j.n).astype(int)
    for k in COMP:j[k]=pd.to_numeric(j[k])
    v=j[list(COMP)].to_numpy(float)
    ensure((j.n>0).all() and np.isfinite(v).all() and ((v>=0)&(v<=1)).all(),"Invalid original ROI9 pixels")
    ensure(np.allclose(v.sum(axis=1),1,atol=1e-5),"Nine component total changed")
    for k in COMP:j[k]*=j.n
    agg=j.groupby(["sector","accepted_id","accepted_name","genus","family"],as_index=False,observed=True)[["n",*COMP]].sum()
    for k in COMP:agg[k]/=agg.n
    ensure(len(agg)==30054 and not agg.duplicated(["sector","accepted_id"]).any(),"Synonym collapse drift")
    original=a[["sector","accepted_id","n","strict9_white","strict9_bronze"]].copy()
    for k in ("n","strict9_white","strict9_bronze"):original[k]=pd.to_numeric(original[k])
    chk=agg.merge(original,on=["sector","accepted_id"],validate="one_to_one",suffixes=("","_prior"))
    ensure(len(chk)==len(agg) and np.array_equal(chk.n,chk.n_prior)
           and np.allclose(chk.strict9_white,chk.strict9_white_prior,atol=1e-6)
           and np.allclose(chk.strict9_bronze,chk.strict9_bronze_prior,atol=1e-6),
           "Independent accepted species photo aggregate drift")
    agg["rank"]=agg.sector.map(ORDER)
    reps=agg.sort_values(["accepted_id","n","rank"],ascending=[True,False,True],kind="mergesort").drop_duplicates("accepted_id")
    ensure(len(reps)==22712 and reps.accepted_id.is_unique,"Repeated accepted species")
    ensure(not d.duplicated(["accepted_species","trait_name"]).any(),"Nonunique external direct records")
    known=d.loc[(d.trait_name=="flower_primary_color")&d.quality.isin(["high","medium"])
        &d.resolution_status.eq("resolved")&~d.classification.eq("cultivar_contamination")
        &d.normalized_value.isin(LABELS),["accepted_species","normalized_value"]].copy()
    known["label"]=known.normalized_value.map(LABELS).astype(int)
    f=reps.merge(known,left_on="accepted_name",right_on="accepted_species",validate="one_to_one")
    ensure(len(f)==6019 and int(f.label.sum())==538,"Independent trait label support drift")
    return f,agg,known

def predict(trainx,trainy,testx):
    ensure(set(trainy)=={0,1},"Train split has one phenotype")
    model=make_pipeline(StandardScaler(),LogisticRegression(max_iter=3000,C=1.,solver="lbfgs"))
    model.fit(trainx,trainy)
    return model.predict_proba(testx)[:,1]

def measures(y,p):
    ensure(set(y)=={0,1},"Test group has one phenotype")
    return {"n":len(y),"positive":int(y.sum()),"prevalence":float(y.mean()),
        "auc":float(roc_auc_score(y,p)),"ap":float(average_precision_score(y,p)),
        "brier":float(brier_score_loss(y,p)),
        "baseline_brier":float(brier_score_loss(y,np.full(len(y),y.mean())))}

def groupfold(f,rank,feature):
    x=f[list(feature)].to_numpy(float);y=f.label.to_numpy(int);groups=f[rank].to_numpy()
    pred=np.full(len(f),np.nan);subs=[];assigned=np.zeros(len(f),dtype=int)
    for i,(tr,te) in enumerate(StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=20261105).split(x,y,groups)):
        ensure(not(set(groups[tr])&set(groups[te])),"Taxonomy leakage")
        pred[te]=predict(x[tr],y[tr],x[te]);assigned[te]+=1
        subs.append({"rank":rank,"features":("full9" if len(feature)==9 else "bronze"),"fold":i,
                     "test_groups":len(set(groups[te])),**measures(y[te],pred[te])})
    ensure((assigned==1).all() and np.isfinite(pred).all(),"OOF incomplete")
    return {"rank":rank,"features":("full9" if len(feature)==9 else "bronze"),
            "unique_groups":len(set(groups)),**measures(y,pred)},subs

def transfer(agg,known,source,target,feature):
    a=agg.loc[agg.sector.eq(source)].merge(known,left_on="accepted_name",right_on="accepted_species",validate="one_to_one")
    b=agg.loc[agg.sector.eq(target)].merge(known,left_on="accepted_name",right_on="accepted_species",validate="one_to_one")
    b=b.loc[~b.accepted_id.isin(set(a.accepted_id))]
    ensure(not(set(a.accepted_id)&set(b.accepted_id)),"Photo species leaked between regions")
    scores=predict(a[list(feature)].to_numpy(float),a.label.to_numpy(int),b[list(feature)].to_numpy(float))
    return {"source":source,"target":target,"features":("full9" if len(feature)==9 else "bronze"),
      "train_n":len(a),"train_positive":int(a.label.sum()),"test_unseen_genera":int((~b.genus.isin(set(a.genus))).sum()),
      "shared_species":0,**measures(b.label.to_numpy(int),scores)}

def run(inputs,folder):
    folder=Path(folder);ensure(not folder.exists(),"Do not overwrite analysis")
    tables={k:read_zip(inputs[k],k) for k in FILES}
    f,agg,known=build_tables(tables["regional"],tables["taxonomy"],tables["accepted"],tables["direct"])
    summaries=[];folds=[]
    for rank in ("genus","family"):
        for cols in ((COMP[-1],),COMP):
            s,p=groupfold(f,rank,cols);summaries.append(s);folds.extend(p)
    regions=[transfer(agg,known,a,b,cols) for a,b in
        (("north_temperate","south_temperate"),("south_temperate","north_temperate"))
        for cols in ((COMP[-1],),COMP)]
    folder.mkdir(parents=True)
    for filename,records in (("rank_holdout.csv",summaries),
                             ("fold_holdout.csv",folds),
                             ("north_south_species_transfer.csv",regions)):
        pd.DataFrame(records).to_csv(folder/filename,index=False,float_format="%.12g")
    rec={"schema":"mois_20261105_original_roi9_species_direct_label_group_holdout_v1",
         "original_source_zip_SHA256":{k:v[1] for k,v in FILES.items()},
         "input_inner_sha256":{k:v[2] for k,v in FILES.items()},
         "matched_distinct_accepted_species":int(agg.accepted_id.nunique()),
         "eligible_species_with_direct_single_colour":len(f),
         "green_brown_documented_species":int(f.label.sum()),
         "other_single_colour_documented_species":int((f.label==0).sum()),
         "genera":int(f.genus.nunique()),"families":int(f.family.nunique()),
         "one_measurement_sector_per_species":"largest strict9 n, fixed sector order for ties",
         "independent_source_trait_not_photo_visual_truth":True,
         "not_background_contrast_or_calibrated_green":True,
         "all_predictions_outside_training_group_or_species":True,
         "no_MNAR_scenarios_or_human_labels_or_image_downloads":True,
         "restricted_species_trait_rows_published":0,
         "library_versions":{"sklearn":sklearn.__version__,"pandas":pd.__version__,"numpy":np.__version__},
         "output_sha256":{p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir())}}
    (folder/"source_receipt.json").write_text(json.dumps(rec,sort_keys=True,indent=2)+"\n")
    print("MOIS_SOURCE_VERIFIED_INDEPENDENT_COLOUR_LABEL_HOLDOUT",json.dumps({
        "species":len(f),"green_brown":int(f.label.sum()),"rank_holdout":summaries,"regional_transfer":regions},sort_keys=True))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in FILES:p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args();run({name:getattr(a,name) for name in FILES},a.output)
