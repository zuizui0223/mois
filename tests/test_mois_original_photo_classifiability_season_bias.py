"""Synthetic photo classification/missingness tests; no simulated ecological claims."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scripts.analysis.mois_seasonal_267k_climate_normals import labelled_calendar
from scripts.analysis.mois_original_photo_classifiability_season_bias import (
    symmetric_white_detection, globe_by_zone, site_year_selection,
    geocell_bootstrap,
)


def original_fake():
    rows=[]
    counter=1
    for lat in (40.2,-40.2):
        for month in (1,7):
            for i in range(8):
                # six independently recorded/colour-classified photo observations
                # and two source ROI unknowns per local meteorological season.
                cc = "white" if i%2==0 else "red_pink"
                rows.append({
                    "photo_id":str(counter),"observation_id":str(counter),
                    "inat_taxon_id":str(300+i%2),
                    "latitude":lat,"longitude":138.05,
                    "observed_on":f"2025-{month:02d}-17",
                    "morph":cc if i<6 else "mixed_uncertain"})
                counter+=1
    return labelled_calendar(pd.DataFrame(rows))


def test_exact_algebra_of_detected_white_yield_across_attempts():
    x=symmetric_white_detection(0.7,0.4,0.5,0.25)
    assert x["summer_minus_winter_white_label_yield_per_all_attempts_pp"]==pytest.approx(25)
    assert x["source_white_hue_share_conditional_on_classified_component_pp"]==pytest.approx(13.75)
    assert x["source_four_hue_classifiability_component_pp"]==pytest.approx(11.25)
    assert sum([x["source_white_hue_share_conditional_on_classified_component_pp"],
                x["source_four_hue_classifiability_component_pp"]]) ==pytest.approx(25)
    with pytest.raises(ValueError,match="probabilities"):
        symmetric_white_detection(1.2,0.4,0.5,0.25)


def test_classification_failure_is_not_counted_as_nonwhite():
    d=original_fake()
    all_rows=globe_by_zone(d)
    assert len(all_rows)==2
    for x in all_rows:
        assert x["summer_actual_original_photo_attempts"]==8
        assert x["summer_four_hue_classified_original_photos"]==6
        assert x["summer_white_labeled_original_photos"]==3
        assert x["summer_classifiability_fraction"]==0.75
        assert x["summer_white_fraction_among_classified"]==0.5


def test_matching_classified_source_cell_preserves_all_missing_attempts():
    d=original_fake()
    # Change two summer-known observations to hue unknown in Northern latitudes.
    # The Southern local season is reversed. Prefer explicit Northern check.
    selected=d[(d.zone=="north_extratropical") & d.season.eq("summer")]
    extra_unknown=selected[selected.photo_id.isin(selected.photo_id.tail(2))]
    # Artificially add two new original source photo attempts with technical ROI
    # unknown. Their known original observation IDs remain independent.
    add=[]
    for _,r in extra_unknown.iterrows():
        row=r.to_dict()
        row["photo_id"]="added_"+str(r["photo_id"])
        row["observation_id"]="added_"+str(r["observation_id"])
        row["morph"]="mixed_uncertain"
        add.append(row)
    new=pd.concat([d,pd.DataFrame(add)],ignore_index=True)
    new["classified"]=new.morph.isin(("white","yellow_orange","red_pink","blue_purple"))
    new["white"]=new.morph.eq("white")
    site=site_year_selection(new,0.25,min_obs=5,min_taxa=2)
    north=next(x for x in site if x["hemisphere"]=="north_extratropical")
    assert north["n_selected_geocell_season_years"]==1
    assert north["classified_photo_records_in_paired_seasons"]==12
    assert north["classified_and_missing_photo_attempts_in_paired_seasons"]==18
    assert north["classifiability_season_source_diagnostic"]["estimate_pp"]==pytest.approx(-15,abs=1e-4)
    assert north["classifiability_season_source_diagnostic"]["cluster_bootstrap95_pp"] is None


def test_geocell_bootstrap_explicitly_holds_scarce_site_data():
    one=pd.DataFrame([{"latitude_cell":12,"longitude_cell":8,"rate":-10.0}])
    x=geocell_bootstrap(one,"rate",draws=10)
    assert x["cluster_bootstrap95_pp"] is None
    grid=pd.DataFrame([{"latitude_cell":i,"longitude_cell":8,"rate":-3.0}
                       for i in range(12)])
    r=geocell_bootstrap(grid,"rate",draws=50)
    assert r["estimate_pp"]==-3.0
    assert r["cluster_bootstrap95_pp"]==[-3.,-3.]


def test_no_real_season_data_is_not_missing_photo_hue_count():
    d=original_fake()
    d=d[d.zone.eq("north_extratropical")].copy()
    d["morph"]="mixed_uncertain"
    d["classified"]=False
    d["white"]=False
    out=site_year_selection(d,.25)
    assert all(x["n_selected_geocell_season_years"]==0 for x in out)
