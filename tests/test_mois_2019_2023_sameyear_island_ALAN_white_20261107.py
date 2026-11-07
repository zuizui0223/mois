"""Offline contract checks for original 2019–2023 island–ALAN source refit."""
import unittest
import pandas as pd
from scripts.analysis import mois_2019_2023_sameyear_island_ALAN_white_20261107 as m

class ActualYearSourceContracts(unittest.TestCase):
    def test_only_real_available_original_years_allowed(self):
        self.assertEqual(m.YEAR_SET,(2019,2021,2022,2023))
        self.assertTrue(all(m.YEAR[y]["known"]>10000 for y in m.YEAR_SET))
        self.assertEqual(len(m.PHOTO_SHA),64)
        self.assertEqual(len(m.WCVP_SHA),64)
    def test_interaction_rejects_tiny_island_group_even_many_mainland_cells(self):
        rows=[]
        for k in range(80):
            rows.append(dict(region5=str(k),land_id="i"+str(k) if k<7 else "main",
                 land_class="island" if k<7 else "mainland",white_pp=40+k/5))
        d=pd.DataFrame(rows)
        result=m.fit_if_supported(d,"white_pp ~ 1","region5",interaction=True)
        self.assertFalse(result["estimable"])
        self.assertEqual(result["reason"],"too_few_island_cells_or_islands")
        self.assertEqual(result["island_cells"],7)
    def test_within_island_inference_requires_independent_islands(self):
        d=pd.DataFrame([dict(land_id=str(k),white_pp=40) for k in range(9) for _ in range(9)])
        result=m.fit_if_supported(d,"white_pp~1","land_id",within=True)
        self.assertFalse(result["estimable"])
        self.assertEqual(result["reason"],"too_few_independent_within_islands")
    def test_geography_ambiguity_exclusion_is_source_constrained(self):
        # A raw cell ID by itself is NOT a safe island/mainland identity.
        land=pd.DataFrame([dict(cell_id="1_2",land_id="4",land_class="island"),
             dict(cell_id="1_2",land_id="5",land_class="mainland"),
             dict(cell_id="3_4",land_id="8",land_class="island")])
        mask=land.cell_id.duplicated(keep=False)
        self.assertEqual(land.loc[~mask].cell_id.tolist(),["3_4"])
    def test_no_2022_october_substitution_is_part_of_model_contract(self):
        self.assertIn(2022,m.YEAR_SET)
        import inspect
        code=inspect.getsource(m.run)
        self.assertIn("PROVIDER_RASTER_400",code)
        self.assertIn("four_month_comparability_to_other_years",code)
if __name__=="__main__":
    unittest.main()
