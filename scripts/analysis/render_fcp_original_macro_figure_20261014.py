#!/usr/bin/env python3
"""Render a descriptive Figure 6 for the ONE FCP New Phytologist manuscript.

The companion JSON is an exact, human-reviewed display ledger copied from
verified FCP PR #166, #191, #193 and #194 source-run results. This renderer
never refits or promotes selection hypotheses.
"""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EXPECTED_SCHEMA = "fcp_np_unified_figure6_source_20261014_v1"
PANELS = ("a", "b", "c", "d")

def make_figure(src, outdir):
    data=json.loads(Path(src).read_text(encoding="utf-8"))
    assert data["schema"]==EXPECTED_SCHEMA
    assert set(data["panels"])==set(PANELS)
    fig, axes=plt.subplots(2,2,figsize=(12.8,8.0),layout="constrained")
    fig.set_facecolor("white")
    ink="#244766"       # Colour-blind-safe blue for observations
    pale="#8c969d"      # Gray for zero and framing
    for key,ax in zip(PANELS, axes.ravel()):
        p=data["panels"][key]
        items=p["items"]
        pts=np.array([x["point"] for x in items],float)
        lows=np.array([x["lo"] for x in items],float)
        highs=np.array([x["hi"] for x in items],float)
        assert np.all(np.isfinite(pts)) and np.all(lows <= pts) and np.all(highs >= pts)
        y=np.arange(len(items))[::-1]
        xerr=np.vstack([pts-lows,highs-pts])
        ax.errorbar(pts,y,xerr=xerr,fmt="o",markersize=6.5,elinewidth=2.0,
                    capsize=4.5,capthick=1.4,color=ink,ecolor=ink,zorder=3)
        ax.axvline(0,color=pale,linestyle="--",lw=1.0,zorder=1)
        ax.set_yticks(y,labels=[x["label"] for x in items])
        ax.tick_params(axis="y",labelsize=10,length=0)
        ax.set_xlim(*p["xlim"])
        ax.set_ylim(-0.6,len(items)-0.4)
        ax.set_title(p["title"],loc="left",fontsize=12.5,fontweight="bold",pad=10)
        ax.set_xlabel(p["xlabel"],fontsize=10.5)
        ax.grid(axis="x",alpha=.17)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for i,(yy,xx) in enumerate(zip(y,pts)):
            # Sparse scientific panel; leave numerical intervals to caption and data ledger.
            if key=="a": ax.annotate(f"{xx:.4f}",(xx,yy),xytext=(4,9),textcoords="offset points",fontsize=8.5)
    fig.suptitle("Geographic scale and phenotype source bound global flower-colour–environment associations",
                 fontsize=15.0,fontweight="semibold")
    fig.supxlabel("Points show original estimates; whiskers denote the source-defined 95% uncertainty intervals. Different panels have different estimands.",
                  fontsize=9.0,color="#444444")
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/"figure6_unified_new_phytologist.svg",bbox_inches="tight")
    fig.savefig(out/"figure6_unified_new_phytologist.png",dpi=240,bbox_inches="tight")
    plt.close(fig)
    return [out/"figure6_unified_new_phytologist.svg",out/"figure6_unified_new_phytologist.png"]

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--source-json",type=Path,default=Path(__file__).resolve().parents[2]/"results/np_unified_20261014/figure6_source_data.json");p.add_argument("--outdir",type=Path,default=Path(__file__).parent)
    a=p.parse_args();files=make_figure(a.source_json,a.outdir)
    for item in files: print(item.name,item.stat().st_size)
