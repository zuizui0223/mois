#!/usr/bin/env python3
"""Reconstruct MOIS ROI9 counts from initial and later signed independent-photo ledgers.

Non-destructive correction of after_strict9 only; preserve source colour and status.
"""
from __future__ import annotations
import argparse, hashlib, io, json
from pathlib import Path
from zipfile import ZipFile
import pandas as pd

SOURCES = {
 "baseline": ("6cc7689938dcf13baf6498b954286924aa97f2b63b3418c41181b18e346c0d2b",
              "mois_species_first_42111_nine_palette_catalogue.csv.gz",
              "5f51b9115719ff3c94a46851f9a3972742f7c52192ceb31b3da110248028fad8"),
 "current": ("17da6386cc2970f8f55f6ee5373e0565027cac7d6c9f2f4d2fc84c532f91df6d",
             "mois_42111_taxa_postsnap_new1909_remaining18_original_photo9.csv.gz",
             "588f74124ff2982af5ada2fe8c0258b0e28f9d91010e0451b0b3023e7c67fbab"),
 "colxr": ("7b8d366d803c9f21d328495bc2517526e1ade6b448bb0c66fa731062e6ae4425",
           "current_exact_accepted_COLXR_numeric9_traits.csv.gz",
           "bc31db8b6b3123e3deb900ddcc0faa709a8ec3a7a04387fc8f0a7610897c4843")
}

def require(condition, message):
    if not bool(condition):
        raise ValueError("MOIS_STRICT9_REPAIR: "+message)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def read_signed(path, kind):
    zsha, name, csvsha = SOURCES[kind]
    raw = Path(path).read_bytes()
    require(sha(raw)==zsha, kind+" source ZIP digest changed")
    with ZipFile(io.BytesIO(raw)) as z:
        csv = z.read(name)
        receipt = json.loads(z.read("mois_species_first_receipt.json" if kind=="baseline" else "receipt.json"))
    require(sha(csv)==csvsha, kind+" internal CSV digest changed")
    if kind=="current":
        require(receipt["output_SHA256"][name]==csvsha, "current calculation receipt mismatch")
    if kind=="colxr":
        require(receipt["source_photo_atlas_SHA256"]==SOURCES["current"][2], "COLXR source atlas mismatch")
    result=pd.read_csv(io.BytesIO(csv),compression="gzip",
                       dtype={"inat_taxon_id":str},low_memory=False)
    require(result.inat_taxon_id.is_unique,kind+" duplicate source taxon")
    return result,receipt

def reconcile(baseline, current, *, new_attempts, new_numeric):
    require(len(baseline)==len(current) and
            baseline.inat_taxon_id.is_unique and current.inat_taxon_id.is_unique and
            set(baseline.inat_taxon_id)==set(current.inat_taxon_id),
            "taxon frame changed")
    prior=baseline.set_index("inat_taxon_id").loc[current.inat_taxon_id].reset_index()
    new_obs=(current.independent_roi9_observations.astype("int64")-
             prior.independent_roi9_observations.astype("int64"))
    attempts=(current.after_photo_attempts.astype("int64")-
              prior.after_photo_attempts.astype("int64"))
    require(new_obs.ge(0).all() and attempts.ge(0).all(),"nonmonotone observation history")
    require((new_obs<=attempts).all(),"numeric observations exceed photo attempts")
    require(int(new_obs.sum())==new_numeric,"signed new original ROI success count mismatch")
    require(int(attempts.sum())==new_attempts,"signed new photo-attempt count mismatch")
    corrected=prior.after_strict9.astype("int64")+new_obs
    as_signed=current.after_strict9.astype("int64")
    require((corrected>=as_signed).all() and
            (corrected>=current.independent_roi9_observations).all(),
            "post-reconciliation counter violates original evidence")
    require(corrected.gt(0).equals(current.hue_status.ne("no_strict9")),"strict9 status mismatch")
    result=current.copy()
    result.insert(result.columns.get_loc("after_strict9")+1,
                  "after_strict9_original_signed",as_signed.to_numpy())
    result["after_strict9"]=corrected.to_numpy()
    result["after_strict9_repaired"]=corrected.ne(as_signed).to_numpy()
    for column in ("hue_status","dominant_hue9","observed_hue9_set","after_photo_attempts",
                   "independent_roi9_observations",
                   "observed_mean_white","observed_mean_yellow","observed_mean_orange",
                   "observed_mean_red","observed_mean_pink","observed_mean_magenta",
                   "observed_mean_purple","observed_mean_blue","observed_mean_bronze"):
        require(result[column].equals(current[column]),"photo hue changed: "+column)
    return result

def run(baseline_zip, current_zip, colxr_zip, output):
    out=Path(output)
    require(not out.exists(),"refuse overwrite of corrected output directory")
    before,b=read_signed(baseline_zip,"baseline")
    current,c=read_signed(current_zip,"current")
    accepted,g=read_signed(colxr_zip,"colxr")
    require(len(before)==len(current)==42111,"42,111 source taxa changed")
    require(b["source_photo_attempts"]==1181146 and
            b["strict_nine_palette_original_photos"]==451643 and
            c["cumulative_original_photo_attempts"]==1231337 and
            c["cumulative_original_nine_ROI_photos"]==467565,
            "all-original source totals changed")
    fixed=reconcile(before,current,new_attempts=1231337-1181146,
                    new_numeric=467565-451643)
    dif=fixed.after_strict9-fixed.after_strict9_original_signed
    require(int(dif.gt(0).sum())==3812 and int(dif.sum())==7409 and
            int(fixed.after_strict9.sum())==457415,
            "3,812-taxon discrepancy was not exactly reconciled")
    require(467565-int(fixed.after_strict9.sum())==10150,
            "initial outside-taxon ROI9 totals not preserved")
    require(len(accepted)==37103 and accepted.inat_taxon_id.isin(fixed.inat_taxon_id).all(),
            "COLXR accepted taxa changed")
    matching=fixed.set_index("inat_taxon_id").loc[accepted.inat_taxon_id].reset_index()
    require((accepted.after_strict9.to_numpy()==
             matching.after_strict9_original_signed.to_numpy()).all() and
            accepted.photo_evidence_status.fillna("").equals(matching.hue_status.fillna("")),
            "COLXR trait bridge no longer matches signed source")
    accepted=accepted.copy()
    accepted.insert(accepted.columns.get_loc("after_strict9")+1,
                    "after_strict9_original_signed",accepted.after_strict9)
    accepted["after_strict9"]=matching.after_strict9.to_numpy()
    accepted["after_strict9_repaired"]=matching.after_strict9_repaired.to_numpy()
    out.mkdir(parents=True)
    reports={
       "mois_42111_taxa_corrected_strict9_counts.csv.gz":fixed,
       "mois_37103_COLXR_corrected_strict9_counts.csv.gz":accepted,
       "mois_3812_repaired_taxa.csv.gz":fixed.loc[fixed.after_strict9_repaired,
         ["inat_taxon_id","species","after_strict9_original_signed","after_strict9",
          "independent_roi9_observations","after_photo_attempts","hue_status","dominant_hue9"]]
    }
    output_sha={}
    for name,table in reports.items():
        p=out/name
        table.to_csv(p,index=False,compression={"method":"gzip","mtime":0})
        output_sha[name]=sha(p.read_bytes())
    receipt={
       "schema":"mois_signed_source_taxon_strict9_count_reconciliation_v1",
       "input_SHA256":{kind:{"zip":v[0],"csv":v[2]} for kind,v in SOURCES.items()},
       "attempted_original_photos":1231337,"numeric_original_ROI9":467565,
       "source_taxa":42111,"source_strict9_count_before":450006,
       "source_strict9_count_after":457415,"restored_photos":7409,
       "taxa_repaired":3812,"out_of_frame_numerical_ROI9_unchanged":10150,
       "COLXR_exact_accepted_species":37103,
       "COLXR_species_with_repaired_counts":int(accepted.after_strict9_repaired.sum()),
       "hue_vectors_and_status_unchanged":True,
       "botanical_petal_truth_proven":False,
       "global_GBIF_occurrences_downloaded":0,
       "output_SHA256":output_sha
    }
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_STRICT9_COUNT_REPAIR_SUCCESS",json.dumps(receipt,sort_keys=True))
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for v in ("baseline","current","colxr","output"):p.add_argument("--"+v,required=True)
    a=p.parse_args()
    run(a.baseline,a.current,a.colxr,a.output)
