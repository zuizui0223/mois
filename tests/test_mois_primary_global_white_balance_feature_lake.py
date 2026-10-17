"""Original-photo weighted and taxon-weighted near-half are not equivalent."""
import pandas as pd
from scripts.analysis import mois_primary_global_white_balance_feature_lake as m


def mock(ids,hues,taxa):
    x=pd.DataFrame({"photo_id":ids,"morph":hues,"inat_taxon_id":taxa,
                    "species":["Campanula a"]*len(ids)})
    x["classified"]=x.morph.isin(m.HUES)
    x["white"]=x.morph.eq("white")
    return x


def test_white_vs_three_colours_unknown_never_assumed_chromatic():
    a=mock(["1","2","3","4","5"],[
        "white","yellow_orange","red_pink","blue_purple","mixed_uncertain"],
        ["1","2","3","4","5"])
    r=m.fractions(a)
    assert (r["photo_attempts"],r["classifiable"],r["white"],
            r["chromatic"],r["unknown"])==(5,4,1,3,1)


def test_species_equal_photo_nearhalf_different_from_photo_weighted():
    a=mock(["1","2","3","4","5","6"],[
        "white","white","white","white","blue_purple","mixed_uncertain"],
        ["1","1","1","1","2","3"])
    r=m.balances(a,draws=29)
    v=r["near_half_rather_than_exact_50_50_fraction"]
    assert v["photos"]==.8
    assert v["species_equal"]==.5
    assert v["genus_equal"]==.5
    assert r["balance_is_not_maintained_50_50_equilibrium"]


def test_known_measured_population_and_solar_conversion_guard():
    assert (m.N,m.CLASS,m.DATED,m.PARTIAL)==(309439,149662,149653,248312)
    assert m.SHA_PHOTO!=m.SHA_CLIMATE!=m.SHA_STATIC
    assert len(m.HUES)==4


def test_decile_denominator_includes_unknown_photo_attempts():
    a=mock([str(i) for i in range(40)],
           ["white","mixed_uncertain","yellow_orange","red_pink"]*10,
           [str(i) for i in range(40)])
    for k in (*m.ENV[:-1],"annual_solar_Wm2"):
        a[k]=range(40)
    d=m.deciles(a)
    q=d[d.variable.eq("wc2_bio12_10m")]
    assert q.photo_attempts.sum()==40
    assert q.classifiable.sum()==30
    assert q.white.sum()==10


def test_genus_depth_restrictions_can_change_near_half_fraction():
    a=mock(["1","2","3","4","5","6"],[
       "white","white","red_pink","red_pink","white","mixed_uncertain"],
       ["1","1","2","2","3","4"])
    a["species"]=["Campanula a","Campanula a",
                  "Campanula b","Campanula b","Cirsium c","Cirsium d"]
    row=m.balance_sampling_sensitivity(a)
    initial=next(x for x in row if x["min_classified_original_photos_per_taxon"]==1
                 and x["min_named_taxa_per_genus"]==1)
    restricted=next(x for x in row if x["min_classified_original_photos_per_taxon"]==1
                 and x["min_named_taxa_per_genus"]==2)
    assert initial["eligible_genera"]==2
    assert restricted["eligible_genera"]==1
    assert restricted["genus_equal_white_fraction"]==0.5


def test_equal_geographic_cell_balance_not_identical_to_classified_count():
    a=mock([str(k) for k in range(90)],
           ["white"]*30+["red_pink"]*60,[str(k) for k in range(90)])
    a["latitude"]=[20.]*30+[35.]*60
    a["longitude"]=40.
    r=m.geographic_cell_sensitivity(a)
    assert r["eligible_cells"]==2
    assert r["photos_in_eligible_cells"]==90
    assert r["equal_cell_mean_white_fraction"]==.5
    assert m.fractions(a)["white_fraction"]==1/3
