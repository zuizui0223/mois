#!/usr/bin/env python3
"""Post-hoc ecological-context diagnostic of out-of-genus soil vs water prediction.

Replays SHA256-verified original FCP photo/environment ZIPs through the existing
source loader, independently refits all seven heldout genus models with
training-fold standardisation, and reports ONLY original-source aggregate
coefficients / environmental-median strata. No taxon, photo, site, coordinate,
genus identifiers or copyrighted original photographs leave source memory.

The eight high-minus-low contrasts are only this post-hoc figure's selected
family. Neither they nor within-half 95% intervals correct the full search.
"""
from __future__ import annotations
import argparse, importlib.util, json
from pathlib import Path
import numpy as np

MODELS = {
 "geography": ["abs_lat","lon_sin","lon_cos"],
 "water2": ["abs_lat","lon_sin","lon_cos","wc_bio12","wc_vapr_annual_mean"],
 "soil5": ["abs_lat","lon_sin","lon_cos","soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"],
 "other24": ["abs_lat","lon_sin","lon_cos",*[f"wc_bio{i}" for i in range(1,20) if i!=12],
             "wc_srad_annual_mean","wc_srad_monthly_cv","wc_wind_annual_mean","wc_wind_monthly_cv",
             "wc_vapr_monthly_cv","wc_elevation_m"],
 "climate26": ["abs_lat","lon_sin","lon_cos",*[f"wc_bio{i}" for i in range(1,20)],
               "wc_srad_annual_mean","wc_srad_monthly_cv","wc_wind_annual_mean","wc_wind_monthly_cv",
               "wc_vapr_annual_mean","wc_vapr_monthly_cv","wc_elevation_m"],
}
MODELS["climate26_soil5"] = MODELS["climate26"] + MODELS["soil5"][3:]
MODELS["climate26_soil10"] = MODELS["climate26"] + MODELS["soil5"][3:] + [
 "soil_cec_0_30cm_source_raw","soil_sand_0_30cm_source_raw",
 "soil_silt_0_30cm_source_raw","soil_bdod_0_30cm_source_raw",
 "soil_cfvo_0_30cm_source_raw",
]
CONTEXT = {"annual_rainfall":"wc_bio12","peak_temperature":"wc_bio5",
           "elevation":"wc_elevation_m","latitude_abs":"abs_lat"}
TARGET = ("water_gain_vs_geography","soil5_increment_vs_climate26")
N_BOOT = 1999
SEED_HALVES = 202610151755
SEED_PAIRED = 202610151801
SHA_PHOTO = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
SHA_ENV = "6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"

def original_loader():
 path = Path(__file__).with_name("mois_full_biogeography_climate_soil_audit_20261015.py")
 spec = importlib.util.spec_from_file_location("mois_fcp_owner_input",path)
 module = importlib.util.module_from_spec(spec)
 spec.loader.exec_module(module)
 return module

def out_of_genus_losses(species,source_module):
 genera,gi=np.unique(species.genus.to_numpy(),return_inverse=True)
 sizes=np.bincount(gi).astype(float)
 weight=1.0/sizes[gi]
 fold=species.genus.map(source_module.fold).to_numpy(int)
 y=species.white.to_numpy(float)
 losses=np.full((len(species),len(MODELS)),np.nan)
 for f in range(5):
  train=fold!=f
  test=~train
  ww=weight[train]/weight[train].sum()
  mean_y=float(ww@y[train])
  for j,features in enumerate(MODELS.values()):
   X=species[features].to_numpy(float)
   mu=(X[train]*ww[:,None]).sum(axis=0)
   sd=np.maximum(np.sqrt(((X[train]-mu)**2*ww[:,None]).sum(axis=0)),1e-9)
   Z=(X[train]-mu)/sd;T=(X[test]-mu)/sd
   b=np.linalg.solve((Z.T*ww)@Z+.01*np.eye(len(features)),
                     (Z.T*ww)@(y[train]-mean_y))
   pred=np.clip(mean_y+T@b,0,1)
   losses[test,j]=(y[test]-pred)**2
 if not np.isfinite(losses).all():raise AssertionError("heldout predictions incomplete")
 return gi,genera,losses

def calc(photo_zip,env_zip):
 module=original_loader()
 assert (module.PHOTO_SHA,module.ENV_SHA)==(SHA_PHOTO,SHA_ENV)
 _,_,_,s,_=module.original(photo_zip,env_zip)
 if (len(s),s.genus.nunique(),int(s.n_photos.sum()))!=(16801,1665,26163):
  raise AssertionError("Original qualified soil36 taxonomy/photos drift")
 gi,genera,losses=out_of_genus_losses(s,module)
 total_sizes=np.bincount(gi).astype(float)
 names=list(MODELS)
 cols={name:i for i,name in enumerate(names)}
 genus_avg=np.stack([np.bincount(gi,weights=losses[:,j])/total_sizes
                     for j in range(len(names))],axis=1)
 mse={name:float(genus_avg[:,j].mean()) for j,name in enumerate(names)}
 if not (abs(mse["climate26"]-0.22813288536129273)<1e-10 and
         abs(mse["climate26_soil5"]-0.22853734681634522)<1e-10 and
         abs(mse["water2"]-0.22745535180555312)<1e-10):
  raise AssertionError("Independent original-source blocked CV does not reproduce earlier full26/soil/water ledgers")

 gain={
  "water_gain_vs_geography":losses[:,cols["geography"]]-losses[:,cols["water2"]],
  "soil5_increment_vs_climate26":losses[:,cols["climate26"]]-losses[:,cols["climate26_soil5"]],
 }
 source_context={};compare_boot=[];compare_idx=[]
 rng_common=np.random.default_rng(SEED_PAIRED)
 draw_genus=rng_common.multinomial(len(genera),np.ones(len(genera))/len(genera),size=N_BOOT)
 for i,(label,feature) in enumerate(CONTEXT.items()):
  levels=s[feature].to_numpy(float)
  cut=float(np.median(levels))
  halves={}
  saved={}
  for which,mask in (("lower",levels<cut),("higher",levels>=cut)):
   count=np.bincount(gi[mask],minlength=len(genera))
   used=count>0
   v={}
   for name,values in gain.items():
    by_genus=np.divide(np.bincount(gi[mask],weights=values[mask],minlength=len(genera)),
                       count,out=np.zeros(len(genera)),where=used)
    point=float(by_genus.sum()/used.sum())
    seed=SEED_HALVES+i*3+(which=="higher")
    rng=np.random.default_rng(seed)
    draws=rng.integers(0,used.sum(),size=(N_BOOT,used.sum()),dtype=np.int32)
    boots=by_genus[used][draws].mean(axis=1)
    saved[(which,name)]=(by_genus,used.astype(float))
    v[name]={"gain":point,"ordinary95":np.quantile(boots,[.025,.975]).tolist()}
   halves[which]={"n_species":int(mask.sum()),"n_genera_present":int(used.sum()),"metrics":v}
  differences={}
  for name in TARGET:
   lo,avail_lo=saved[("lower",name)]
   hi,avail_hi=saved[("higher",name)]
   lo_point=float(lo.sum()/avail_lo.sum());hi_point=float(hi.sum()/avail_hi.sum())
   lo_bs=(draw_genus@lo)/np.maximum(1,draw_genus@avail_lo)
   hi_bs=(draw_genus@hi)/np.maximum(1,draw_genus@avail_hi)
   d=hi_point-lo_point
   db=hi_bs-lo_bs
   differences[name]={"higher_minus_lower":float(d),
        "paired_genus_ordinary95":np.quantile(db,[.025,.975]).tolist()}
   compare_boot.append(db-d);compare_idx.append((label,name,differences[name]))
  source_context[label]={"source_median":cut,"halves":halves,"high_minus_low":differences}
 # Compare exactly eight post-hoc selected high-vs-low contrasts at once.
 E=np.column_stack(compare_boot)
 sd=np.std(E,axis=0,ddof=1)
 critical=float(np.quantile(np.max(np.abs(E/sd),axis=1),.95))
 for j,(label,name,entry) in enumerate(compare_idx):
  entry["eight_selected_contrast_maxT95"]=[float(entry["higher_minus_lower"]-critical*sd[j]),
                                            float(entry["higher_minus_lower"]+critical*sd[j])]
 main_soil_genus=genus_avg[:,cols["climate26"]]-genus_avg[:,cols["climate26_soil5"]]
 main_boot=np.random.default_rng(SEED_HALVES).integers(0,len(genera),size=(N_BOOT,len(genera)),dtype=np.int32)
 return {
  "schema":"mois_source_soil_moisture_context_eight_contrast_v1",
  "source_sha256":{"photo":SHA_PHOTO,"environment":SHA_ENV},
  "original_support":{"species":16801,"genera":1665,"photos_in_eligible_taxa":26163},
  "models":names,
  "model_heldout_equal_genus_mse":mse,
  "soil5_addition_to_full26":{"point_gain":float(main_soil_genus.mean()),
                             "ordinary95":np.quantile(main_soil_genus[main_boot].mean(axis=1),[.025,.975]).tolist()},
  "context":source_context,
  "uncertainty":{"bootstrap_draws":N_BOOT,"half_seed_base":SEED_HALVES,"paired_seed":SEED_PAIRED,
                 "eight_selected_contrast_maxT_critical":critical,
                 "maxT_only_eight_posthoc_contrasts":True,
                 "no_study_wide_multiple_testing_or_nested_model_refitting":True},
  "interpretation":"Water-related climate adds predictive information relative to source geography alone. Adding soil5 to all26 does not gain predictive accuracy globally and shows no positive point increment in any evaluated environmental median half. Eight climate-context differences of water2 or soil5 gain have intervals spanning zero after a joint eight-selected-contrast maxT adjustment; no stable evidence that soil increment or water2 gain is modified by these studied climate/altitude/latitude context splits.",
  "limitations":[
   "Source photo labels and SoilGrids eligibility nonrandom, and no actual terrain geometry measured.",
   "These four splits across same species, with genera present in both halves, are retrospective diagnostics of a frozen predictor family; not independent biome replications.",
   "No genetic morph frequencies, lifetime fitness, exact local soil water or causal pigment mechanisms.",
   "No source photo ID, coordinates, taxon names or original image/soil rows published."
  ]
 }

def check_numbers(observed,frozen,path="root"):
 if isinstance(observed,dict):
  assert isinstance(frozen,dict) and set(observed)==set(frozen),path
  for k,x in observed.items():check_numbers(x,frozen[k],path+"."+k)
 elif isinstance(observed,list):
  assert len(observed)==len(frozen),path
  for i,x in enumerate(observed):check_numbers(x,frozen[i],path+"["+str(i)+"]")
 elif isinstance(observed,(float,int)) and not isinstance(observed,bool):
  assert abs(observed-frozen)<2e-6,(path,observed,frozen)
 else:assert observed==frozen,(path,observed,frozen)

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--photo-zip",required=True,type=Path)
 p.add_argument("--env-zip",required=True,type=Path)
 p.add_argument("--output",required=True,type=Path)
 opt=p.add_mutually_exclusive_group(required=True)
 opt.add_argument("--write",action="store_true")
 opt.add_argument("--check",action="store_true")
 a=p.parse_args();result=calc(a.photo_zip,a.env_zip)
 if a.write:
  a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
  print("WROTE independent source climate-soil context result without original locations")
 else:
  prior=json.loads(a.output.read_text(encoding="utf8"))
  check_numbers(result,prior)
  print("PASS source SHA256, fixed genus-heldout CV and eight climate-context inference checks")

if __name__=="__main__":
 main()
