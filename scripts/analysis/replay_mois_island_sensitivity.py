"""Replay exploratory island-size and influential-landmass sensitivity.
Usage: python scripts/analysis/replay_mois_island_sensitivity.py EXTRACTED_RESULTS
No image download, photo inference, or coastline reclassification is required.
"""
import argparse,json
from pathlib import Path
import pandas as pd
import statsmodels.formula.api as smf

def run(root):
    c=pd.read_csv(root/'matched_island_mainland_cells.csv')
    reference=pd.read_csv(root/'island_mainland_adjusted.csv')
    f='p_white ~ island + z_heat + z_rain + z_solar + z_elevation + C(region)'
    ranked=c[c.island.eq(1)].groupby('landmass').size().sort_values(ascending=False)
    cases=[('replay_primary',c),('minimum_10_taxa',c[c.n_taxa>=10]),
           ('islands_at_least_5km2',c[c.island.eq(0)|c.area_km2.ge(5)]),
           ('exclude_islands_100000km2_plus',c[c.island.eq(0)|c.area_km2.lt(100000)])]
    cases += [(f'drop_island_{lm}',c[c.landmass.ne(lm)]) for lm in ranked.head(5).index]
    results=[]
    for label,q in cases:
        both=q.groupby('region').island.nunique()
        q=q[q.region.isin(both[both==2].index)]
        m=smf.ols(f,q).fit(cov_type='cluster',cov_kwds={'groups':q.block},use_t=True)
        ci=m.conf_int().loc['island']
        results.append({'analysis':label,'effect_pp':float(100*m.params['island']),
                        'low_pp':float(100*ci.iloc[0]),'high_pp':float(100*ci.iloc[1]),
                        'nominal_p':float(m.pvalues['island']),
                        'cells':len(q),'island_cells':int(q.island.sum()),
                        'clusters':int(q.block.nunique())})
    ref=reference[reference.model.eq('geography_climate_altitude')].iloc[0]
    if abs(results[0]['effect_pp']-ref.island_minus_mainland_pp)>=1e-8:
        raise ValueError('Independent primary replay differs from source coefficient')
    pd.DataFrame(results).to_csv(root/'replayed_sensitivity.csv',index=False)
    return results
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path)
    a=p.parse_args();print(json.dumps(run(a.directory),indent=2))
