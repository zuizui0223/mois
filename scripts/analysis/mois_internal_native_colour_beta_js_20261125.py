"""MOIS internal native-colour spatial composition: no external trait joins.
Source artifacts: MOIS 14710329666 and 14400177816. Outputs aggregates only.
"""
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np,pandas as pd
HUES='white yellow orange red pink magenta purple blue bronze'.split()
SOURCES={'native999':('c031cf8fe3f80d6cf868ee20477fdb4c550b836435b79ae16d75e771c564f951','native_photo9_georegion_cells.csv.gz','e07c5906c10bc6944189053ee11ef35da6b76f0f93365b74ac0d13f1d71b327e'), 'geo':('cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd','global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz','449521d176ef8b1e029b5236ed464224d20b7538de1fb8d1f1552127d40c34ec')}
def sha(b):return hashlib.sha256(b).hexdigest()
def load(p,k):
    outer,member,inner=SOURCES[k];raw=Path(p).read_bytes();assert sha(raw)==outer
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        r=json.loads(z.read('receipt.json')); data=z.read(member)
        h=r.get('output_SHA256',r.get('output_sha256'))[member]
        assert sha(data)==h and h==inner
    return pd.read_csv(io.BytesIO(data),compression='gzip'),h
def entropy(p):
    return -(p*np.log(np.maximum(p,1e-300))).sum(axis=1)
def split(a,b):
    p=a/a.sum(axis=1)[:,None]; q=b/b.sum(axis=1)[:,None]; m=(p+q)/2
    total=entropy(m)-(entropy(p)+entropy(q))/2
    bp=np.column_stack([p[:,0],p[:,1:].sum(axis=1)]); bq=np.column_stack([q[:,0],q[:,1:].sum(axis=1)])
    binary=entropy((bp+bq)/2)-(entropy(bp)+entropy(bq))/2
    w=p[:,1:].sum(axis=1); v=q[:,1:].sum(axis=1); u=(w+v)/2
    cp=p[:,1:]/np.maximum(w[:,None],1e-300); cq=q[:,1:]/np.maximum(v[:,None],1e-300)
    cm=m[:,1:]/np.maximum(u[:,None],1e-300)
    conditional=u*entropy(cm)-(w*entropy(cp)+v*entropy(cq))/2
    assert np.max(np.abs(total-binary-conditional),initial=0)<1e-12
    assert np.min(np.c_[total,binary,conditional],initial=0)>-1e-12
    return np.maximum(0,np.c_[total,binary,conditional])
def pairs(d,cohort,threshold,drop_bronze):
    cols=[f'{cohort}_observed_{c}_species' for c in HUES]
    a=d[cols].to_numpy(int); assert (a>=0).all() and np.array_equal(a.sum(1),d[f'{cohort}_richness'])
    if drop_bronze:a=a[:,:-1]
    index={(x,y):i for i,(x,y) in enumerate(d[['grid_x','grid_y']].itertuples(index=False,name=None)) if a[i].sum()>=threshold}
    edges=[(i,index[(x+dx,y+dy)]) for (x,y),i in index.items() for dx,dy in [(1,0),(0,1)] if (x+dx,y+dy) in index]
    i,j=np.array(edges).T; assert len(set(edges))==len(edges)
    v=split(a[i],a[j]);x=d.grid_x.to_numpy();y=d.grid_y.to_numpy(); elev=d.sampled_land_median_elevation_m.to_numpy(float)
    bx=np.floor_divide(x,5);by=np.floor_divide(y,5)
    assert np.isfinite(v).all()
    return pd.DataFrame({'x1':x[i],'y1':y[i],'x2':x[j],'y2':y[j],'block_x':np.floor((x[i]+x[j])/10).astype(int),'block_y':np.floor((y[i]+y[j])/10).astype(int),'within_block':(bx[i]==bx[j])&(by[i]==by[j]),'north':(y[i]>=0)&(y[j]>=0),'south':(y[i]<0)&(y[j]<0),'elev_delta':abs(elev[i]-elev[j]),'elev_mean':(elev[i]+elev[j])/2,'total':v[:,0],'binary':v[:,1],'chromatic':v[:,2]})
def summary(f,seed=20261125):
    if not len(f):return {}
    a=f[['total','binary','chromatic']].to_numpy(); m=a.mean(0)
    blocks=f.groupby(['block_x','block_y'])[['total','binary','chromatic']].agg(['sum','count']); s=blocks.xs('sum',axis=1,level=1).to_numpy(); count=blocks[('total','count')].to_numpy(); k=len(s)
    equal=s/count[:,None]; be=equal.mean(0)
    rng=np.random.default_rng(seed); ratios=[]
    for _ in range(20):
        ii=rng.integers(k,size=(1000,k)); sums=s[ii].sum(1); ratios.extend(np.divide(sums[:,2],sums[:,0],out=np.full(len(sums),np.nan),where=sums[:,0]>0))
    return {'edges':len(f),'blocks':k,'mean_JS':m[0],'mean_binary_JS':m[1],'mean_chromatic_JS':m[2],'chromatic_share':m[2]/m[0],'chromatic_share_block_equal':be[2]/be[0],'chromatic_share_spatial_resampling_95':np.nanpercentile(ratios,[2.5,97.5]).tolist()}
def diagnostics(d):
    f=pairs(d,'strict_native',100,False);f2=pairs(d,'strict_native',100,True)
    keys=['x1','y1','x2','y2']; z=f.merge(f2[keys+['total','binary','chromatic']],on=keys,validate='one_to_one',suffixes=('_all','_no_bronze')); assert len(z)==len(f2)
    a=d[[f'strict_native_observed_{c}_species' for c in HUES]].to_numpy(int); idx={(x,y):i for i,(x,y) in enumerate(d[['grid_x','grid_y']].itertuples(index=False,name=None))}
    i=np.array([idx[(x,y)] for x,y in f[['x1','y1']].itertuples(index=False,name=None)]);j=np.array([idx[(x,y)] for x,y in f[['x2','y2']].itertuples(index=False,name=None)])
    p=a[i]/a[i].sum(1)[:,None];q=a[j]/a[j].sum(1)[:,None];m=(p+q)/2
    term=.5*(p*np.log(np.maximum(p,1e-300)/np.maximum(m,1e-300))+q*np.log(np.maximum(q,1e-300)/np.maximum(m,1e-300)))
    assert np.max(abs(term.sum(1)-f.total))<1e-12
    binary_one=[]
    for k in range(9):
        b=np.column_stack([a[i,k],a[i].sum(1)-a[i,k]]);c=np.column_stack([a[j,k],a[j].sum(1)-a[j,k]])
        binary_one.append(split(b,c)[:,0].mean())
    return {'n_same_edges':len(z),'same_edge_full9_chromatic_share':z.chromatic_all.sum()/z.total_all.sum(),'same_edge_exclude_bronze_chromatic_share':z.chromatic_no_bronze.sum()/z.total_no_bronze.sum(),'mean_category_js_contributions':dict(zip(HUES,term.mean(0).tolist())),'binary_one_vs_rest_js':dict(zip(HUES,binary_one)),'max_independent_kl_identity_error':float(np.max(abs(term.sum(1)-f.total))),'note':'One-vs-rest projections overlap; not additive.'}
def main(args):
    d,h1=load(args.native999,'native999');geo,h2=load(args.geo,'geo'); assert len(d)==len(geo)==18577 and not d.duplicated(['grid_x','grid_y']).any()
    d=d.merge(geo[['grid_x','grid_y','sampled_land_median_elevation_m','strict_native_n_photo9_species','moderate_native_n_photo9_species']],on=['grid_x','grid_y'],validate='one_to_one')
    for c,n in [('strict_native',2271852),('moderate_native',2604870)]:
        assert d[f'{c}_richness'].sum()==n and np.array_equal(d[f'{c}_richness'],d[f'{c}_n_photo9_species'])
    out=[]
    for c,t,drop in [('strict_native',30,False),('strict_native',100,False),('moderate_native',100,False),('strict_native',100,True)]:
        f=pairs(d,c,t,drop)
        for name,mask in [('all',np.ones(len(f),bool)),('north',f.north),('south',f.south),('within_block',f.within_block),('steep_elevation_delta_ge500',f.elev_delta>=500),('flat_elevation_delta_le100',f.elev_delta<=100)]:
            s=summary(f.loc[mask]);out.append(dict(cohort=c,threshold=t,exclude_bronze=drop,group=name,**s))
    result={'sources':{k:{'zip':v[0],'inner':h} for (k,v),h in zip(SOURCES.items(),[h1,h2])},'interpretation':'Adjacent cell composition, species counts not abundance, descriptive spatial resampling only. No Island DB. No permutation of species identity or causal inference.','rows':out,'diagnostics':diagnostics(d)}
    target=Path(args.output)
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native999',required=True);p.add_argument('--geo',required=True);p.add_argument('--output',required=True);main(p.parse_args())
