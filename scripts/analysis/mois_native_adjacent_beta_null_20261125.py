#!/usr/bin/env python3
"""Synchronous native species-label nulls for adjacent MOIS colour composition.
Read only a reconstructed source-validated native species-cell NPZ.
Required keys: cell, species, hue, genus, coords, zone. Order of the species
catalogue is the original GBIF coverage CSV order; grid order is signed
native_photo9_georegion_cells.csv.gz order. Run with exact source replay.
No image, external trait, new occurrence download, or per-cell independent shuffle.
"""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

ZONES=['Americas_N','Americas_S','Africa_Eurasia_W_N','Africa_Eurasia_W_S','Asia_Oceania_E_N','Asia_Oceania_E_S']
MEASURES=['total','binary','conditional_nonwhite']

def entropy(p):
    return -(p*np.log(np.maximum(p,1e-300))).sum(1)

def js_parts(a,b):
    p=a/a.sum(1)[:,None]; q=b/b.sum(1)[:,None]
    full=entropy((p+q)/2)-(entropy(p)+entropy(q))/2
    bp=np.column_stack([p[:,0],p[:,1:].sum(1)])
    bq=np.column_stack([q[:,0],q[:,1:].sum(1)])
    binary=entropy((bp+bq)/2)-(entropy(bp)+entropy(bq))/2
    conditional=full-binary
    assert np.min(np.c_[full,binary,conditional])>=-1e-12
    return np.column_stack([full,binary,conditional])

def groups_for(code,h,eligible,whitefixed):
    avail=np.flatnonzero(eligible & ((h>0) if whitefixed else True))
    order=avail[np.argsort(code[avail],kind='stable')]
    bounds=np.r_[0,1+np.flatnonzero(np.diff(code[order])),len(order)]
    groups=[order[a:b] for a,b in zip(bounds[:-1],bounds[1:])
            if b-a>1 and len(set(h[order[a:b]]))>1]
    return groups

def static_counter(index,species,h,mobile,nrows):
    mobmap=np.full(len(h),-1,np.int32); mobmap[mobile]=np.arange(len(mobile))
    m=mobmap[species]; keep=m>=0
    raw=np.bincount(index+h[species],minlength=nrows*9).reshape(nrows,9)
    i=index[keep]; j=m[keep]
    static=raw-np.bincount(i+h[mobile][j],minlength=nrows*9).reshape(nrows,9)
    assert static.min()>=0
    def draw(shuffled):
        c=static+np.bincount(i+shuffled[j],minlength=nrows*9).reshape(nrows,9)
        assert np.array_equal(c.sum(1),raw.sum(1))
        return c
    return raw,draw,int(keep.sum())

def summary(a,o):
    n=len(a); m=a.mean()
    pl=(1+(a<=o+1e-12).sum())/(n+1); pu=(1+(a>=o-1e-12).sum())/(n+1)
    return {'observed':float(o),'null_mean':float(m),'delta':float(o-m),
            'relative_delta':float((o-m)/m) if m else None,
            'null_interval_95':np.quantile(a,[.025,.975]).tolist(),
            'p_two':float(min(1,2*min(pl,pu))),'p_lower':float(pl),'p_upper':float(pu)}

def source(path,reference_csv):
    f=np.load(path,allow_pickle=False)
    cell=f['cell']; spp=f['species']; h=f['hue']; genus=f['genus']; coords=f['coords']; zn=f['zone']
    n=len(h); c=len(coords)
    assert n==37104 and c==18577 and len(cell)==len(spp)==2271852
    assert ((cell>=0)&(cell<c)&(spp>=0)&(spp<n)).all()
    assert len(np.unique(cell.astype(np.int64)*n+spp))==len(cell)
    assert ((h>=0)&(h<9)).all()
    ref=pd.read_csv(reference_csv)
    hues='white yellow orange red pink magenta purple blue bronze'.split()
    expected=ref[[f'strict_native_observed_{k}_species' for k in hues]].to_numpy()
    assert np.array_equal(coords,ref[['grid_x','grid_y']].to_numpy())
    assert np.array_equal(zn,ref.coordinate_sector.to_numpy(str))
    raw=np.bincount(cell.astype(np.int64)*9+h[spp],minlength=c*9).reshape(c,9)
    assert np.array_equal(raw,expected)
    return cell,spp,h,genus,coords,zn,raw

def run(path,reference_csv,out,seed=20261125,B=999):
    cell,spp,h,genus,coords,zn,raw=source(path,reference_csv)
    N=len(h); C=len(coords); zone=pd.Index(ZONES).get_indexer(zn); assert (zone>=0).all()
    gc=pd.factorize(genus,sort=True)[0]; occ=np.bincount(spp,minlength=N); eligible=occ>0
    breadth=np.digitize(occ,[2,6,21,101])
    cross=np.bincount(spp.astype(np.int64)*6+zone[cell],minlength=N*6).reshape(N,6)
    aff=np.where(cross.max(1)/np.maximum(occ,1)>=.65,cross.argmax(1),6)
    base=gc*5+breadth; strict=base*8+aff
    rich=raw.sum(1); index={tuple(x):i for i,x in enumerate(coords) if rich[i]>=100}
    edges=np.array([(i,index[(x+dx,y+dy)]) for (x,y),i in index.items()
                    for dx,dy in [(1,0),(0,1)] if (x+dx,y+dy) in index])
    assert len(edges)==6962
    I,J=edges.T
    masks={'global':np.ones(len(I),bool),
           'north':(coords[I,1]>=0)&(coords[J,1]>=0),
           'south':(coords[I,1]<0)&(coords[J,1]<0),
           'within_500km':(np.floor_divide(coords[I],5)==np.floor_divide(coords[J],5)).all(1)}
    def edge_stats(c):
        a=js_parts(c[I],c[J]);return np.array([a[m].mean(0) for m in masks.values()])
    obs=edge_stats(raw); assert abs(obs[0,0]-.004095388025525331)<1e-12
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result={'n_native_species':int(eligible.sum()),'n_native_pairs':len(cell),'n_cells':C,
            'n_edges':{k:int(v.sum()) for k,v in masks.items()},'seed':seed,'draws':B,'models':{}}
    # Same full species catalogue and strata as the local computation.
    for name,code,wf in [('genus',gc,False),('genus_native_width',base,False),
                         ('genus_native_width_zone',strict,False),
                         ('genus_native_width_zone_whitefixed',strict,True)]:
        groups=groups_for(code,h,eligible,wf);mobile=np.concatenate(groups)
        gid=np.repeat(np.arange(len(groups)),[len(g) for g in groups]);h0=h[mobile]
        rng=np.random.default_rng(seed)
        c0,counter,nmob=static_counter(cell.astype(np.int64)*9,spp,h,mobile,C)
        assert np.array_equal(raw,c0)
        draws=np.empty((B,len(masks),3)); alpha=np.empty(B)
        for b in range(B):
            sh=h0[np.lexsort((rng.random(len(mobile)),gid))]
            c=counter(sh)
            if wf:assert np.array_equal(c[:,0],raw[:,0])
            draws[b]=edge_stats(c)
            alpha[b]=entropy(c[rich>=30]/rich[rich>=30,None]).mean()
        stats={k:{v:summary(draws[:,i,j],obs[i,j]) for j,v in enumerate(MEASURES)}
               for i,k in enumerate(masks)}
        result['models'][name]={'groups':stats,'mobile_species':len(mobile),
                                'variable_strata':len(groups),'mobile_pairs':nmob,
                                'alpha_H9_ge30_null_mean':float(alpha.mean())}
        np.savez_compressed(out/(name+'_draws.npz'),draw=draws,alpha=alpha)
    # Secondary diagnostic: remove every shared accepted species from both endpoints.
    order=np.lexsort((spp,cell)); sp=spp[order]; ends=np.r_[0,rich.cumsum()]
    xa=[];xb=[];nn=[];sharing=[]
    for i,j in edges:
        a=sp[ends[i]:ends[i+1]]; b=sp[ends[j]:ends[j+1]]
        same=np.intersect1d(a,b,assume_unique=True)
        aa=np.setdiff1d(a,b,assume_unique=True); bb=np.setdiff1d(b,a,assume_unique=True)
        assert len(aa)+len(same)==len(a) and len(bb)+len(same)==len(b)
        sharing.append(len(same)/(len(same)+len(aa)+len(bb)))
        if min(len(aa),len(bb))>=30:
            xa.append(aa);xb.append(bb);nn.append([len(aa),len(bb),len(same)])
    nn=np.array(nn);E=len(nn);assert E==6024
    a=np.concatenate(xa);b=np.concatenate(xb)
    ia=np.repeat(np.arange(E,dtype=np.int64)*9,nn[:,0])
    ib=np.repeat(np.arange(E,dtype=np.int64)*9,nn[:,1])
    result['exclusive']={'n_edges':E,'mean_jaccard_similarity_all_edges':float(np.mean(sharing)),
                         'mean_exclusive_left_right_shared':nn.mean(0).tolist(),'models':{}}
    for wf in [False,True]:
        groups=groups_for(strict,h,eligible,wf);mobile=np.concatenate(groups)
        gid=np.repeat(np.arange(len(groups)),[len(g) for g in groups]);h0=h[mobile]
        ca,da,_=static_counter(ia,a,h,mobile,E);cb,db,_=static_counter(ib,b,h,mobile,E)
        observed=js_parts(ca,cb).mean(0); rng=np.random.default_rng(seed);draws=np.empty((B,3))
        for k in range(B):
            sh=h0[np.lexsort((rng.random(len(mobile)),gid))]
            A=da(sh);B_=db(sh)
            if wf:assert np.array_equal(A[:,0],ca[:,0]) and np.array_equal(B_[:,0],cb[:,0])
            draws[k]=js_parts(A,B_).mean(0)
        key='whitefixed' if wf else 'full9'
        result['exclusive']['models'][key]={k:summary(draws[:,j],observed[j]) for j,k in enumerate(MEASURES)}
        np.savez_compressed(out/f'exclusive_{key}_draws.npz',draw=draws)
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--native-npz',required=True)
    p.add_argument('--reference-csv',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args();run(a.native_npz,a.reference_csv,a.out)
