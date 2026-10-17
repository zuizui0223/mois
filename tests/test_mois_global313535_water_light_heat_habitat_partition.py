"""MOIS original-photo hue denominators and matched environmental analysis guards."""
import numpy as np
import pandas as pd
from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as m

def source(seed=11,genera=74,species=3,photos=3):
    rng=np.random.default_rng(seed);rows=[];number=0
    for i in range(genera):
        genus="G"+chr(65+i//26)+chr(65+i%26)
        for j in range(species):
            rain=220+30*i+rng.normal(0,320)
            solar=150+rng.normal(0,45)
            heat=22+rng.normal(0,9)
            dry=max(0,rain*.04+rng.normal(0,20))
            seas=55+rng.normal(0,22)
            for k in range(photos):
                number+=1
                colour=("white" if rng.random()<.25+.00009*rain+
                        .0006*solar else "red_pink")
                if rng.random()<.08:colour="mixed_uncertain"
                rows.append({
                    "photo_id":str(number),
                    "inat_taxon_id":str(i*species+j),
                    "species":genus+" s"+str(j),
                    "morph":colour,
                    "classified":colour in m.HUES,
                    "white":colour=="white",
                    "old_mois_or_new_fcp":"historical" if k==0 else "fcp",
                    "latitude":(-44 if i%2 else 44)+rng.normal(0,2),
                    "longitude":121+i*.06+rng.normal(0,2),
                    m.RAIN:rain,m.DRY:dry,m.SEASON:seas,
                    m.HEAT:heat,m.SOLAR:solar,
                    "wc2_elev_2p5m":200+rng.normal(0,90),
                    "terrain_slope_2p5m_degree":abs(rng.normal(8,3)),
                    "terrain_TRI_2p5m_m":abs(rng.normal(50,14)),
                    "soil_pH":5+rng.normal(0,.9) if i<65 else np.nan,
                    "soil_SOC":15+rng.normal(0,3) if i<65 else np.nan,
                    "soil_N":8+rng.normal(0,2) if i<65 else np.nan,
                    "soil_clay":20+rng.normal(0,5) if i<65 else np.nan,
                    "soil_available_water_proxy":90+rng.normal(0,15)
                        if i<65 else np.nan,
                })
    d=pd.DataFrame(rows)
    d["genus"]=d.species.map(m.genus)
    return d

def test_photo_hues_not_equivalent_to_genetically_maintained_half():
    assert (m.N_PHOTOS,m.N_HUES,m.N_DATED,m.WHITE_N)==(
        313535,151239,151230,71427)
    assert len(m.HUES)==4
    assert "mixed_uncertain" not in m.HUES
    assert m.SOLAR=="annual_solar_Wm2"

def test_classification_denominator_includes_unknown_original_photos():
    d=source()
    r=m.rates(d)
    assert r["photo_attempts"]==len(d)
    assert r["photo_white_ROIs"]+r["photo_chromatic_ROIs"]==r[
        "classified_original_photo_ROIs"]
    assert r["unclassified_original_photo_ROIs"]==len(d)-int(d.classified.sum())
    assert r["hue_count_white"]==r["photo_white_ROIs"]

def test_species_and_genus_weights_differ_from_photographer_photo_volume():
    d=pd.DataFrame({
        "photo_id":["1","2","3","4","5"],
        "inat_taxon_id":["1","1","1","1","2"],
        "genus":["Campanula"]*5,
        "morph":["white"]*4+["red_pink"],
        "white":[True]*4+[False],
        "classified":[True]*5,
    })
    r=m.rates(d)
    assert r["photo_white_fraction"]==.8
    assert r["equal_named_taxon_white_fraction"]==.5
    assert r["equal_genus_white_fraction"]==.5

def test_tertile_cutoffs_formed_before_hue_filter_and_two_dimensional_cells():
    d=source()
    t,cuts=m.mark_tertiles(d)
    assert set(cuts)=={m.RAIN,m.HEAT,m.SOLAR}
    v=m.two_way_tables(t)
    for k,g in v.groupby("conditioning_axis"):
        assert len(g)==9
        assert int(g.photo_attempts.sum())==len(d)
        assert int(g.classified_original_photo_ROIs.sum())==int(d.classified.sum())

def test_latitude_conditional_contrast_remains_observational():
    d=source()
    t,_=m.mark_tertiles(d)
    tbl,r=m.latitude_strata_wet_dry(t,min_photo=2,min_taxa=2,draws=99)
    assert len(tbl)==r["matched_strata"]
    if not tbl.empty:
        assert r["not_causal_moisture_response"]
        assert tbl.within_latitude_solar_only_not_full_geographic_matching.all()

def test_terrain_soil_models_match_same_species_source_per_nested_block():
    d=source()
    terms,blocks,r=m.fit_habitat_blocks(d)
    assert len(terms)>30
    assert set(blocks.block)=={"landscape_terrain","partially_covered_soil"}
    assert r["all_annual_climate_source"]["qualifying_genera"]>=50
    assert r["landscape_terrain"]["matched_original_photo_support_within_block"]
    assert r["partially_covered_soil"]["qualifying_genera"]>=50
    assert r["partially_covered_soil"]["matched_original_photo_support_within_block"]
    assert any(terms.variable.isin(("rain_x_solar","rain_x_heat")))
