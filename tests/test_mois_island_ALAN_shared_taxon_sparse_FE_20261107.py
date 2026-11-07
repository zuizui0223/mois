"""No-real-source/no-network sparse WLS species fixed-effect identity checks."""
import unittest
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scripts.analysis import mois_island_ALAN_shared_taxon_sparse_FE_20261107 as m

class SharedSourceTaxonFixedEffects(unittest.TestCase):
    def test_pure_same_taxon_across_land_classes_and_one_cell_equal_weight(self):
        sp="wcvp_accepted_plant_name_id"
        data=pd.DataFrame([
         {"target_id":"2019|a|1",sp:"A","year":2019,"is_island":1,"region5":"X",
          "photo_white_fraction":.3},
         {"target_id":"2019|a|1",sp:"B","year":2019,"is_island":1,"region5":"X",
          "photo_white_fraction":.5},
         {"target_id":"2019|b|2",sp:"A","year":2019,"is_island":0,"region5":"X",
          "photo_white_fraction":.1},
         {"target_id":"2019|c|1",sp:"C","year":2019,"is_island":1,"region5":"Y",
          "photo_white_fraction":.8},
        ])
        selected,c=m.true_crossland_accepted_taxa(data,
              pd.DataFrame({"target_id":data.target_id.unique()}))
        self.assertEqual(set(selected[sp]),{"A"})
        self.assertEqual(len(selected),2)
        self.assertEqual(c["source_accepted_taxa_both_land_classes"],1)
        self.assertTrue(np.allclose(
            selected.groupby("target_id").photo_cell_equal_species_weight.sum().to_numpy(),1))

    def test_sparse_weighted_species_FE_matches_dense_statsmodels_WLS(self):
        rng=np.random.default_rng(4519)
        species_ids=[f"WCVP{i}" for i in range(24)]
        records=[]
        for cell in range(480):
            year=[2019,2021,2023,2024,2025][cell%5]
            region=f"G{cell%16}"
            island=int((cell//5)%2)
            spp=rng.choice(species_ids,2,replace=False)
            light=float(rng.uniform(.0,2.1))
            pop=float(rng.uniform(.0,3))
            bio5=float(rng.normal(20,4))
            rain=float(rng.uniform(3,7))
            elev=float(rng.uniform(0,1500))
            effort=float(rng.uniform(1.2,3.3))
            for sid in spp:
                white=(45 + .65*light - 1.7*light*island + 2*island
                    - .8*pop + .11*bio5 + .16*rain + .003*elev
                    + 1.1*effort
                    + float(sid[4:])*.03 + float(rng.normal(0,3)))
                records.append(dict(
                    target_id=f"{year}|{cell}|{'i' if island else 'm'}",
                    year=year,region5=region,land_id="L"+str(cell%7),
                    wcvp_accepted_plant_name_id=sid,
                    is_island=island,photo_white_fraction=white/100.,
                    log_light=light,log_pop=pop,bio5=bio5,log_rain=rain,
                    elevation=elev,log_attempts=effort,
                    photo_cell_equal_species_weight=.5))
        d=pd.DataFrame(records)
        result=m.fit_sparse_taxon_FE(d)
        self.assertTrue(result["supported"])
        self.assertEqual(result["source"]["original_photo_cell_years"],480)
        self.assertEqual(result["source"]["observed_accepted_taxa"],24)
        beta=result["interaction_geographic_spatial_cluster"]["beta"]
        d["white_pp"]=100*d.photo_white_fraction
        full=smf.wls(
            "white_pp ~ log_light*is_island + log_pop"
            " + bio5 + log_rain + elevation + log_attempts"
            " + C(wcvp_accepted_plant_name_id) + C(year) + C(region5)",
            data=d,weights=d.photo_cell_equal_species_weight,
        ).fit()
        self.assertAlmostEqual(beta,float(full.params["log_light:is_island"]),places=6)
        self.assertTrue(result["interaction_geographic_spatial_cluster"]["estimable"])
        self.assertEqual(result["source"]["source_accepted_taxon_clusters"],24)
        without=m.fit_sparse_taxon_FE(d,effect_cols=("region5","year"))
        no_taxon=smf.wls(
            "white_pp ~ log_light*is_island + log_pop"
            " + bio5 + log_rain + elevation + log_attempts"
            " + C(year) + C(region5)",
            data=d,weights=d.photo_cell_equal_species_weight,
        ).fit()
        self.assertAlmostEqual(
            without["interaction_geographic_spatial_cluster"]["beta"],
            float(no_taxon.params["log_light:is_island"]),places=6)
        self.assertNotIn("species",without["source"]["fixed_effect_cardinalities"])
        self.assertEqual(without["source"]["original_taxon_photo_cell_years"],len(d))

    def test_interaction_ci_does_not_claim_precision_from_invalid_variance(self):
        x=m.interaction_stats(-2.5,-.1,12,500)
        self.assertFalse(x["estimable"])
        self.assertEqual(x["reason"],"nonpositive_or_nonfinite_cluster_sandwich_variance")

if __name__=="__main__":
    unittest.main()
