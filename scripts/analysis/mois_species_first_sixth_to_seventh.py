#!/usr/bin/env python3
"""Verify all sixteen sixth-wave original RGB9 workers; append 1406 genuine photo
nine-fraction taxa; select next distinct photos breadth first. Signed data only.
No prior photo bytes are fetched, inferred, or silently labelled.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_species_first_palette9 import PALETTE,COLUMNS,class_labels,sha256
from scripts.analysis.mois_species_first_incremental_first256 import require,parse_boolean
from scripts.analysis.mois_species_first_wave3_incremental import add_to_catalogue
from scripts.analysis.mois_species_first_real1280_cumulative import choose_third_wave

SCHEMA="mois_species_first_sixth4096_authentic_append_seventh4096_v1"
ATLAS_SHA="650fa24c27d972b7c67fac99734a873b489a2a51be1a0d34c59005d68dbe357c"
PRIOR_SOURCE=1185498
PRIOR_STRICT=453219
PRIOR_SUPPORTED=31618
PRIOR_LABELLED=29168
HUES={"white","yellow_orange","red_pink","blue_purple"}
PIX=["flower_fraction_"+p for p in PALETTE]
AUDIT_SCHEMA="mois_species_first_sixth4096_real_original9_audit_v1"

def read_sixth(source,workers,audited):
    source,workers,audited=map(Path,(source,workers,audited))
    manifest=json.loads((source/"manifest.json").read_text())
    require(manifest.get("schema")=="mois_full_source_observation_sampling_cohort_v1" and
            manifest.get("source_arm")=="MOIS_SPECIES_FIRST_SIXTH_NEW4096" and
            manifest.get("colour_measured") is False and
            manifest.get("work_units")==16 and
            manifest.get("worker_ids")==list(range(16)) and
            manifest.get("new_photo_attempts_planned")==4096 and
            manifest.get("source_previous_attempts")==PRIOR_SOURCE,
            "Sixth frozen source manifest is not the distinct original-photo 4096 cohort")
    report=json.loads((audited/"receipt.json").read_text())
    require(report.get("schema")==AUDIT_SCHEMA and
            report.get("actual_independent_original_photo_attempts")==4096 and
            report.get("actual_genuine_strict9_numeric_success")==1406 and
            report.get("actual_clear_nine_palette_hues")==1087 and
            report.get("actual_technical_unknowns")==2690 and
            report.get("photo_observation_source_taxon_ID_collisions")==0,
            "Full sixth-wave audit success/unknown receipt mismatch")
    raw=[]
    for i in range(16):
        shard=source/"shards"/f"shard-{i:03d}.csv"
        folder=workers/f"mois-sixth4096-genuine-original-worker-{i//4}"/f"shard_{i:03d}"
        terminal=folder/"photos.csv"
        receipt=folder/"receipt.json"
        require(shard.is_file() and terminal.is_file() and receipt.is_file(),
                "Missing signed sixth original worker shard")
        require(sha256(shard)==manifest["photo_shard_sha256"][str(i)] and
                sha256(receipt)==report["sixteen_original_photo_worker_receipt_sha256"][str(i)],
                "Sixth original source or worker receipt SHA changed")
        sig=json.loads(receipt.read_text())
        require(sig.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1" and
                sig.get("shard")==i and sig.get("attempted")==256 and
                sig.get("source_sha256")==sha256(shard) and
                sig.get("photos_sha256")==sha256(terminal) and
                sig.get("photo_colours_measured") is True,
                "Sixth worker original photo pixel terminal signature failed")
        src=pd.read_csv(shard,dtype=str,keep_default_na=False)
        original=pd.read_csv(terminal,dtype=str,keep_default_na=False)
        require(len(src)==len(original)==256 and
                src.photo_id.is_unique and original.photo_id.is_unique and
                set(src.photo_id)==set(original.photo_id),
                "Sixth worker photo IDs do not equal signed source shard")
        for col in ("observation_id","inat_taxon_id"):
            require(src.set_index("photo_id").sort_index()[col].equals(
                    original.set_index("photo_id").sort_index()[col]),
                    "Sixth worker original photo taxon/observation identity drift")
        raw.append(original)
    data=pd.concat(raw,ignore_index=True)
    require(len(data)==4096 and data.photo_id.is_unique and
            data.observation_id.is_unique and data.inat_taxon_id.is_unique,
            "Sixth worker returned duplicate original photo/observation/taxon")
    actual=pd.read_csv(audited/"mois_sixth4096_actual_original_nine_roi.csv.gz",
                       dtype=str,keep_default_na=False)
    require(len(actual)==4096 and actual.photo_id.is_unique and
            set(actual.photo_id)==set(data.photo_id),
            "Sixth audit derivative differs from signed photo worker union")
    orig=data.set_index("photo_id").sort_index()
    der=actual.set_index("photo_id").sort_index()
    for col in ("inat_taxon_id","observation_id","morph","measurement_status","image_sha256"):
        if col in orig.columns and col in der.columns:
            require(orig[col].equals(der[col]),"Sixth audited original identity/provenance modified: "+col)
    num=orig[PIX].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    strict=(orig.morph.isin(HUES).to_numpy(bool) &
            np.isfinite(num).all(axis=1) & (num>=0).all(axis=1) &
            (num<=1.00001).all(axis=1) &
            (np.abs(np.nan_to_num(num).sum(axis=1)-1)<=0.001))
    flags=parse_boolean(der.strict9_original_photo_numeric).to_numpy(bool)
    require(np.array_equal(flags,strict) and int(strict.sum())==1406,
            "Sixth numeric nine-fraction strict flags changed after original ROI-v4 workers")
    labelled=np.full(len(orig),"",dtype=object)
    labelled[strict]=class_labels(num[strict])
    require((labelled==der.dominant_hue9_label.to_numpy()).all() and
            int((labelled!="").sum())==1087,
            "Sixth clear hue label not reproducible from genuine original pixel fractions")
    require(der.loc[~flags,PIX].eq("").all().all(),
            "Technical unknown must retain blank nine-fraction phenotype")
    observed=der.loc[flags,PIX].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    require(np.allclose(num[flags],observed,atol=1e-6,rtol=0) and
            np.allclose(num[flags].sum(axis=1),1,atol=.001),
            "Numerical genuine nine-palette fractions changed in audit derivative")
    frame=pd.DataFrame({"photo_id":orig.index.to_numpy(),
                        "observation_id":orig.observation_id.to_numpy(),
                        "inat_taxon_id":orig.inat_taxon_id.to_numpy(),
                        "strict9":strict,"dominant_label9":labelled})
    for i,col in enumerate(COLUMNS):
        frame[col]=np.where(strict,num[:,i],np.nan)
    return frame,{"audit_receipt_sha256":sha256(audited/"receipt.json"),
                  "audit_derived_photo_sha256":sha256(audited/"mois_sixth4096_actual_original_nine_roi.csv.gz"),
                  "sixteen_original_photo_worker_receipt_sha256":report["sixteen_original_photo_worker_receipt_sha256"],
                  "original_photo_source_manifest_sha256":sha256(source/"manifest.json")}

def append_and_select(atlas,sixth,previous,options,quota=4096):
    require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique and
            int(atlas.after_strict9.gt(0).sum())==PRIOR_SUPPORTED and
            int(atlas.dominant_hue9.fillna("").ne("").sum())==PRIOR_LABELLED,
            "Verified fifth 42111 nominal taxon atlas totals drift")
    require(len(previous)==PRIOR_SOURCE and previous.photo_id.is_unique,
            "Signed previous original photo attempt registry length/IDs changed")
    require(len(sixth)==4096 and sixth.photo_id.is_unique and
            sixth.observation_id.is_unique and sixth.inat_taxon_id.is_unique,
            "Sixth original photo/observation/taxon uniqueness lost")
    prev_id=set(previous.photo_id.astype(str))
    prev_obs=set(previous.observation_id.fillna("").astype(str))-{""}
    require(not (prev_id & set(sixth.photo_id.astype(str))) and
            not (prev_obs & set(sixth.observation_id.astype(str))),
            "Sixth new original photo/observation reused previously acquired original")
    updated=add_to_catalogue(atlas,sixth,expected_photos=4096)
    require(int(updated.after_strict9.gt(0).sum())==PRIOR_SUPPORTED+1406 and
            int(updated.dominant_hue9.fillna("").ne("").sum())==PRIOR_LABELLED+1087,
            "Sixth original ROI outcome did not add 1406 numerical/1087 labelled source taxa")
    used=pd.concat([previous[["photo_id","observation_id","inat_taxon_id"]],
                    sixth[["photo_id","observation_id","inat_taxon_id"]]],ignore_index=True)
    require(len(used)==PRIOR_SOURCE+4096 and used.photo_id.is_unique,
            "Full original photo-ID inventory changed after sixth")
    queue=choose_third_wave(updated,used,options,quota=quota)
    require(len(queue)==quota and queue.photo_id.is_unique and queue.observation_id.is_unique and
            queue.inat_taxon_id.is_unique and
            queue.acquisition_priority.isin(["first_unmeasured","retry_unclear"]).all(),
            "Seventh distinct source photo queue cannot be authenticated")
    q=queue.merge(updated[["inat_taxon_id","dominant_hue9"]],on="inat_taxon_id",validate="one_to_one",
                  suffixes=("","_atlas"))
    require(q.dominant_hue9_atlas.fillna("").eq("").all(),
            "Seventh queue remeasures a species with an already assigned photographed hue")
    queue["measurement_id"]="MOIS-SPECIES-FIRST-W7-"+queue.photo_id.astype(str)
    queue["image_filename"]=queue.measurement_id+".jpg"
    queue["shard"]=(np.arange(len(queue))//256).astype(str)
    return updated,queue

def run(atlas,source,workers,audited,previous,options,out,*,quota=4096):
    out=Path(out)
    require(not out.exists(),"Do not overwrite completed signed sixth original photo append")
    require(sha256(atlas)==ATLAS_SHA,"Fifth 42111 source atlas SHA changed")
    df=pd.read_csv(atlas,compression="infer",dtype={"inat_taxon_id":str})
    new,proof=read_sixth(source,workers,audited)
    old=pd.concat([pd.read_parquet(p,columns=["photo_id","observation_id","inat_taxon_id"])
                   for p in previous.values()],ignore_index=True)
    choices=pd.read_csv(options,compression="infer",dtype=str)
    updated,queue=append_and_select(df,new,old,choices,quota=quota)
    out.mkdir(parents=True)
    paths={"atlas":out/"mois_species_first_42111_after_sixth.csv.gz",
           "delta":out/"mois_sixth4096_real_original9_photo_delta.parquet",
           "next":out/"mois_seventh4096_prioritized_original_photo_queue.csv"}
    updated.to_csv(paths["atlas"],index=False,compression="gzip")
    new.to_parquet(paths["delta"],index=False)
    queue.to_csv(paths["next"],index=False)
    receipt={
        "schema":SCHEMA,"source_taxa_total":42111,
        "signed_previous_original_photo_attempts":PRIOR_SOURCE,
        "signed_previous_original_nineROI":PRIOR_STRICT,
        "sixth_new_original_photo_attempts":4096,
        "sixth_new_strict_numeric_nineROI":int(new.strict9.sum()),
        "sixth_new_clear_dominant_palette":int(new.dominant_label9.ne("").sum()),
        "sixth_new_technical_unknown":int((~new.strict9).sum()),
        "cumulative_original_photo_attempts":PRIOR_SOURCE+4096,
        "cumulative_original_numerical_nineROI":PRIOR_STRICT+int(new.strict9.sum()),
        "source_taxa_with_numerical_nineROI":int(updated.after_strict9.gt(0).sum()),
        "source_taxa_with_clear_dominant_nine_palette":int(updated.dominant_hue9.fillna("").ne("").sum()),
        "source_taxa_without_nineROI":int(updated.after_strict9.eq(0).sum()),
        "seventh_source_distinct_original_photos":len(queue),
        "seventh_source_first_unmeasured":int(queue.acquisition_priority.eq("first_unmeasured").sum()),
        "seventh_source_retry_unclear_taxa":int(queue.acquisition_priority.eq("retry_unclear").sum()),
        "previous_photo_or_observation_reused":False,
        "sixth_original_photo_provenance":proof,
        "output_sha256":{key:sha256(path) for key,path in paths.items()},
        "GBIF_species_occurrences_downloaded":0,
        "species_wide_colour_fixation_proven":False}
    (out/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_SIXTH4096_APPEND_SEVENTH4096",json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("atlas","source","workers","audited","baseline","first","second",
              "third","fourth","fifth","global-options","output"):
        p.add_argument("--"+k,required=True)
    p.add_argument("--quota",type=int,default=4096)
    a=p.parse_args()
    run(a.atlas,a.source,a.workers,a.audited,
        {k:getattr(a,k) for k in ("baseline","first","second","third","fourth","fifth")},
        a.global_options,a.output,quota=a.quota)

if __name__=="__main__":
    main()
