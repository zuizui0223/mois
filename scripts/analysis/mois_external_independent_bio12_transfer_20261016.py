#!/usr/bin/env python3
"""Run a new, source-independent Mois white–precipitation test.

NOT a completed empirical replication until independent GBIF specimen records,
a licensed/reproducible WorldClim BIO12 raster, and original Dryad colours
have been supplied. This replayer does not request remote data and never
exports individual taxon names or georeferenced specimen points.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

TRAIT_SHA='a253308785e4cbd0e361b3ca04cdfdae29c523eefa843375cebf8030ee0874af'
DOI='10.5061/dryad.r4xgxd2sc'
EXCLUDES={'INTRODUCED','INVASIVE','MANAGED','CULTIVATED','FOSSIL_SPECIMEN'}


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
    return h.hexdigest()


def binomial(x):
    """Keep only canonical Genus species, never fuzzy matching or synonyms."""
    terms=x.strip().replace('_',' ').split()
    if len(terms)<2 or not re.fullmatch(r'[A-Za-z][A-Za-z-]+',terms[0]) or not re.fullmatch(r'[A-Za-z][A-Za-z-]+',terms[1]):
        return ''
    return (terms[0]+' '+terms[1]).casefold()


def load_traits(path,validate=True):
    if validate and sha256(path)!=TRAIT_SHA:raise ValueError('Dryad species/colour checksum differs from frozen analysis')
    traits={}
    with open(path,encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f)
        if list(reader.fieldnames or [])!=['clade','species','flower_color','fruit_color']:
            raise ValueError('Unexpected source traits columns')
        for row in reader:
            taxon=binomial(row['species'])
            if not taxon or taxon in traits:raise ValueError('Invalid or duplicate species name after binomial normalization')
            traits[taxon]=(row['clade'],int(row['flower_color'].strip().casefold()=='white'))
    return traits


def points_from_gbif(path, traits, *, min_year=2010, max_uncertainty=5000, preserve_only=True):
    """Read source-identified GBIF occurrence CSV/TSV; never claim wild/native from blank metadata."""
    found=[];c=Counter();seen=set();cseen=set()
    with open(path,encoding='utf-8-sig',newline='') as f:
        head=f.read(10000);f.seek(0)
        delimiter='\t' if '\t' in head.splitlines()[0] else ','
        rdr=csv.DictReader(f,delimiter=delimiter)
        cols=set(rdr.fieldnames or [])
        required={'gbifID','scientificName','decimalLatitude','decimalLongitude','basisOfRecord','occurrenceStatus','year','datasetKey'}
        if not required.issubset(cols):raise ValueError(f'GBIF download lacks columns {sorted(required-cols)}')
        for row in rdr:
            c['rows_raw']+=1
            source_id=row['gbifID'].strip()
            if not source_id or source_id in seen:c['duplicate_or_missing_gbifid']+=1;continue
            seen.add(source_id)
            species=binomial(row['scientificName'])
            if species not in traits:c['taxon_not_matched']+=1;continue
            if row['occurrenceStatus'].upper() not in ('PRESENT',''):c['not_present']+=1;continue
            if preserve_only and row['basisOfRecord'].upper()!='PRESERVED_SPECIMEN':
                c['not_preserved_specimen']+=1;continue
            est=row.get('establishmentMeans','').strip().upper()
            if any(x in est for x in EXCLUDES):c['managed_or_introduced']+=1;continue
            issues=(row.get('issue') or '').upper()
            if 'COORDINATE' in issues or 'COUNTRY_COORDINATE_MISMATCH' in issues:c['coordinate_issues']+=1;continue
            try:
                lat=float(row['decimalLatitude']);lon=float(row['decimalLongitude'])
                yr=int(row['year'])
                unc=float(row.get('coordinateUncertaintyInMeters') or 'nan')
            except (ValueError,TypeError):c['malformed']+=1;continue
            if not (-90<=lat<=90 and -180<=lon<=180) or (abs(lat)<1e-6 and abs(lon)<1e-6):
                c['impossible_coords']+=1;continue
            if yr<min_year or yr>2026:c['out_of_year_range']+=1;continue
            if np.isfinite(unc) and unc>max_uncertainty:c['uncertain_location']+=1;continue
            if not row['datasetKey']:c['missing_dataset_key']+=1;continue
            # 0.25 degree lon/lat climate exposure cells, *not equal area*.
            key=(species,int(math.floor((lat+90)*4)),int(math.floor((lon+180)*4)))
            if key in cseen:c['duplicate_species_grid']+=1;continue
            cseen.add(key)
            found.append((species,lat,lon))
            c['accepted_unique_species_grid']+=1
    return found, dict(c)


def sample_bio12(points,raster_path, multiplier=1.0):
    # Import only at execution time: CI's pure numerical smoke test needs no GIS stack.
    import rasterio
    from rasterio.warp import transform
    result=defaultdict(list);bad=Counter()
    with rasterio.open(raster_path) as ds:
        if ds.count!=1:raise ValueError('Expected single-band annual precipitation GeoTIFF')
        if ds.crs is None:raise ValueError('Raster must declare a projected geographic CRS')
        nodata=ds.nodata
        lon=[p[2] for p in points];lat=[p[1] for p in points]
        xx,yy=transform('EPSG:4326',ds.crs,lon,lat)
        for (species,lat,lon),pix in zip(points,ds.sample(zip(xx,yy),indexes=1,masked=True)):
            if np.ma.is_masked(pix):bad['nodata']+=1;continue
            v=float(pix)*multiplier
            if (not np.isfinite(v)) or (nodata is not None and math.isclose(v,nodata*multiplier)) or v<0 or v>25000:
                bad['outside_plausible_mm']+=1;continue
            result[species].append((lat,lon,v))
    return result,dict(bad)


def fit_white_bio12(traits, site_climate, min_cells=5, min_clade_species=3, min_total_species=200, min_total_clades=15, min_variable_clades=10, bootstraps=1999,seed=202610162151):
    taxa=[]
    for species,(clade,white) in traits.items():
        records=site_climate.get(species,[])
        if len(records)<min_cells:continue
        median=float(np.median([z[2] for z in records]))
        latitude=float(np.mean([abs(z[0]) for z in records]))
        lsin=float(np.mean([math.sin(math.radians(z[1])) for z in records]))
        lcos=float(np.mean([math.cos(math.radians(z[1])) for z in records]))
        taxa.append((clade,white,median,latitude,lsin,lcos))
    eligible_counts=Counter(row[0] for row in taxa)
    selected=[row for row in taxa if eligible_counts[row[0]]>=min_clade_species]
    clade_white=defaultdict(set)
    for row in selected:clade_white[row[0]].add(row[1])
    variable=sum(len(vals)==2 for vals in clade_white.values())
    if len(selected)<min_total_species or len(clade_white)<min_total_clades or variable<min_variable_clades:
        return {'status':'insufficient_independent_species_and_clades',
                'n_taxa_candidate':len(taxa),'n_taxa_eligible':len(selected),
                'n_clades_eligible':len(clade_white),'n_clades_with_white_variation':variable,
                'n_eligible_white':sum(row[1] for row in selected),
                'predeclared_gate':{'min_species':min_total_species,'min_clades':min_total_clades,'min_variable_clades':min_variable_clades}}
    climate_sd=float(np.std([row[2] for row in selected]))
    if climate_sd<=0:raise ValueError('No source climatic gradient')
    groups=defaultdict(list)
    for clade,white,rain,alat,lsin,lcos in selected:
        groups[clade].append([rain/climate_sd,alat/10,lsin,lcos,float(white)])
    mom=[]
    for rows in groups.values():
        rows=np.asarray(rows,dtype=float)
        x=rows[:,:4]-rows[:,:4].mean(axis=0)
        y=rows[:,4]-rows[:,4].mean()
        w=1/len(rows)
        mom.append((w*x.T@x, w*x.T@y))
    xx=sum(t[0] for t in mom);xy=sum(t[1] for t in mom)
    if np.linalg.matrix_rank(xx)<4:raise ValueError('Rank-deficient source coordinates')
    beta=float(np.linalg.solve(xx,xy)[0])*100
    rng=np.random.default_rng(seed);bs=[]
    for _ in range(bootstraps):
        ix=rng.integers(0,len(mom),size=len(mom))
        A=sum((mom[i][0] for i in ix));b=sum((mom[i][1] for i in ix))
        if np.linalg.matrix_rank(A)==4:
            bs.append(float(np.linalg.solve(A,b)[0])*100)
    if len(bs)<bootstraps*.95:raise ValueError('Clade bootstrap too many non-estimable designs')
    lo,hi=map(float,np.quantile(bs,[.025,.975]))
    return {'status':'source_independent_climate_analysis_completed',
            'n_taxa_candidate':len(taxa),'n_taxa_eligible':len(selected),
            'n_clades_eligible':len(groups),'n_clades_with_white_variation':variable,
            'n_eligible_white':sum(row[1] for row in selected),
            'rainfall_median_species_source_sd_mm':climate_sd,
            'white_bio12_pp_per_source_sd':beta,'ordinary_clade_bootstrap95_pp':[lo,hi],
            'n_clade_resamples_estimable':len(bs),'bootstraps':bootstraps,
            'warning':'Clade-mean centered equal-clade WLS on recorded preserved specimen environments. Not causal trait selection; not equal to Mois photo-site BIO12 scaling, original clade composition, or fixed familywise maxT tests.'}


def self_test():
    assert binomial('Annona_reticulata')=='annona reticulata'
    assert binomial('Annona reticulata L.')=='annona reticulata'
    assert binomial('Invalid')==''
    traits={f'genus{g} species{i}':(f'clade{g}',int(i%2==0)) for g in range(6) for i in range(4)}
    climate={}
    for idx,k in enumerate(traits):
        g,i=divmod(idx,4)
        climate[k]=[(15+g*3+i*2+g*i*.12+z*.01,
                     -120+g*18+i*i+g*i*.23+z*.01,
                     500+g*35+37*math.sin(2*i+g)+z) for z in range(5)]
    out=fit_white_bio12(traits,climate,min_total_species=20,min_total_clades=5,min_variable_clades=5,bootstraps=99)
    assert out['status']=='source_independent_climate_analysis_completed'
    assert out['n_taxa_eligible']==24 and out['n_clades_eligible']==6
    assert len(out['ordinary_clade_bootstrap95_pp'])==2
    print('PASS synthetic independent occurrence BIO12 estimator; NO external GBIF/raster climate values tested')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--traits',type=Path)
    p.add_argument('--gbif-occurrences',type=Path)
    p.add_argument('--gbif-download-doi')
    p.add_argument('--worldclim-bio12',type=Path)
    p.add_argument('--worldclim-citation')
    p.add_argument('--raster-factor',type=float,default=1)
    p.add_argument('--output',type=Path)
    p.add_argument('--self-test',action='store_true')
    args=p.parse_args()
    if args.self_test:self_test();return
    if not all([args.traits,args.gbif_occurrences,args.worldclim_bio12,args.gbif_download_doi,args.worldclim_citation]):
        p.error('Need verified traits, independently sourced GBIF occurrence export and DOI, and credited WorldClim raster')
    if not re.fullmatch(r'10\.15468/dl\.[\w-]+',args.gbif_download_doi.strip(),re.I):
        p.error('A citable GBIF occurrence-download DOI is mandatory')
    traits=load_traits(args.traits)
    points,counts=points_from_gbif(args.gbif_occurrences,traits)
    bio12,bad=sample_bio12(points,args.worldclim_bio12,args.raster_factor)
    out=fit_white_bio12(traits,bio12)
    out.update({'source_trait_sha256':sha256(args.traits),'traits_doi':DOI,
                'gbif_occurrence_download_doi':args.gbif_download_doi,
                'gbif_download_sha256':sha256(args.gbif_occurrences),
                'worldclim_raster_sha256':sha256(args.worldclim_bio12),
                'worldclim_citation':args.worldclim_citation,
                'preserved_specimen_only':True,
                'geography_grid':'0.25 degree latitude-longitude source spatial thinning (not equal area)',
                'filter_counts':counts,'raster_invalid_counts':bad})
    output=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output)
    print(output)


if __name__=='__main__':main()