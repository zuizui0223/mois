"""Scientific figures from reviewed model coefficient tables, preserving estimands."""
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def forest(ax,d,label,title,xlabel):
    positions=np.arange(len(d))[::-1]
    ax.axvline(0,color='#777777',lw=.9,ls='--')
    for y,(_,r) in zip(positions,d.iterrows()):
        ax.plot([r.ci95_low,r.ci95_high],[y,y],color='#246387',lw=1.8)
        ax.plot(r.estimate,y,'o',color='#246387',ms=5)
    ax.set_yticks(positions);ax.set_yticklabels([label.get(x,x) for x in d['label']])
    ax.set_xlabel(xlabel);ax.set_title(title,loc='left',fontweight='bold',pad=12)
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
    ax.grid(axis='x',alpha=.15);ax.set_ylim(-.6,len(d)-.4)

def run(broad,nested,roi,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    b=pd.read_csv(broad);b=b.loc[b.term.eq('log_alan')].copy();b['label']=b.model
    names={'primary_region_year':'Primary: region + year','north':'Northern hemisphere','south':'Southern hemisphere','within_cell_year':'Within-cell changes','within_genus_region_year':'Within-genus comparison'}
    fig,ax=plt.subplots(figsize=(9,4.5));forest(ax,b,names,'Documented white-flowered species composition','Percentage-point difference per unit log1p(ALAN)')
    fig.text(.03,.02,'Points: conditional association estimates; lines: 95% cluster-robust confidence intervals.\nSpecies composition, not abundance or causal change. Spatial clusters; genus model also clusters genera.',fontsize=9)
    fig.subplots_adjust(left=.29,bottom=.22,top=.88,right=.96);fig.savefig(out/'broad_macro_coefficients.png',dpi=180);plt.close(fig)
    r=pd.read_csv(roi);r=r.loc[r.term.eq('log_alan')].copy();r['label']=r.model
    names={'original_ROI9_taxon_balanced':'ROI white area: equal taxa','ROI9_photo_weighted_sensitivity':'ROI white area: equal photos','technical_ROI9_success':'Technical ROI success rate'}
    fig,ax=plt.subplots(figsize=(9,4));forest(ax,r,names,'Original ROI9 image supplement: distinct outcomes','Percentage-point difference per unit log1p(ALAN)')
    fig.text(.03,.02,'Original measurement version; failed ROIs remain unknown. Technical success is not flower colour.\nRegion/year adjustments and spatial clustering; canonical coordinate precision is unavailable.',fontsize=9)
    fig.subplots_adjust(left=.33,bottom=.25,top=.86,right=.96);fig.savefig(out/'original_roi_sameyear_coefficients.png',dpi=180);plt.close(fig)
    n=pd.read_csv(nested);n=n.loc[n.term.isin(['white','night_centered','white_night'])].copy();n['label']=n.term
    names={'white':'White colour at night ratio 0.5','night_centered':'Night ratio among nonwhite','white_night':'White × night-ratio interaction'}
    fig,ax=plt.subplots(figsize=(9,4));forest(ax,n,names,'Experimental-species sampled light associations','Difference in log1p(ALAN)')
    fig.text(.03,.02,'Species-balanced; spatial/species clusters. Night ratio = night/(night+day) reproductive outcome.\nLiterature ratios treated as given; no local pollinator identification or full phylogenetic correction.',fontsize=9)
    fig.subplots_adjust(left=.36,bottom=.25,top=.86,right=.96);fig.savefig(out/'experimental_species_coefficients.png',dpi=180);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--broad',required=True);p.add_argument('--nested',required=True);p.add_argument('--roi',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.broad,a.nested,a.roi,a.out)
