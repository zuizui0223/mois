"""Plot all prelisted sector omissions and matched phylogenetic interactions."""
from pathlib import Path
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.analysis.mois_macro_figures import forest

def run(regional,phylogenetic,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    r=pd.read_csv(regional);r=r.loc[r.term.eq('log_alan')&(r.analysis.eq('southern_baseline')|r.analysis.str.startswith('omit_sector_'))].copy();r['label']=r.analysis
    labels={'southern_baseline':'All southern observations','omit_sector_east_of_60':'Omit longitude ≥60°E','omit_sector_minus30_to_60':'Omit longitude −30° to <60°','omit_sector_west_of_minus30':'Omit longitude <−30°'}
    fig,ax=plt.subplots(figsize=(10,4.5));forest(ax,r,labels,'Southern association depends on geographic support','Percentage-point difference per unit log1p(ALAN)')
    fig.text(.03,.025,'All three sector omissions shown. Lines: individual 95% spatial-cluster intervals.\nAttenuation after omission does not prove zero association outside that sector; these are post-primary checks.',fontsize=9)
    fig.subplots_adjust(left=.34,right=.97,bottom=.23,top=.87);fig.savefig(out/'southern_geographic_sensitivity.png',dpi=180);plt.close(fig)
    p=pd.read_csv(phylogenetic);p=p.loc[p.term.eq('white:ratio_centered')].copy();p['label']=p.analysis
    labels={}
    for definition,short in [('night_share','Original ratio (65 species)'),('raw_control_matched','Control-matched raw (46)'),('net_night_relative_outcome','Control-adjusted net (46)')]:
        for method,suffix in [('independent_species','independent'),('phylogenetic_plus_species','phylogeny')]:labels[definition+'__'+method]=short+'; '+suffix
    fig,ax=plt.subplots(figsize=(11,5.5));forest(ax,p,labels,'White-colour × diel-index interaction remains uncertain','Interaction coefficient for adjusted log1p(ALAN)')
    fig.text(.03,.025,'Same species and full first-stage sampling covariance within each method pair. Lines: 95% intervals.\nNet index = (night − exclusion)/(night + day), retaining negative contrasts; one published tree, ratios treated as given.',fontsize=9)
    fig.subplots_adjust(left=.43,right=.97,bottom=.20,top=.89);fig.savefig(out/'control_phylogeny_interaction_sensitivity.png',dpi=180);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--regional',required=True);p.add_argument('--phylogenetic',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.regional,a.phylogenetic,a.out)
