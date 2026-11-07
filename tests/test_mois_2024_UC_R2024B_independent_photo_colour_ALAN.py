"""Independent R2024B UC source is never pasted into R2025A CN missing pixels."""
import unittest
import numpy as np,pandas as pd
from scripts.analysis import mois_2024_UC_R2024B_independent_photo_colour_ALAN as m

class UCPhotoLightIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        n=66456
        cls.cn=pd.DataFrame({
            "target_id":["2024|"+str(i) for i in range(n)],
            "viirs_median":np.ones(n),
            "viirs_valid_months":np.array([4]*64646+[2]*(n-64646)),
            "land_class":["island"]*8908+["mainland"]*(n-8908),
            "n_attempted":np.ones(n,dtype=int),"n_strict9":np.ones(n,dtype=int),
            "population_source_status":["OK"]*44667+["RASTER_OUTSIDE_OR_NODATA"]*(n-44667),
            "complete":np.array([True]*43595+[False]*(n-43595))})
        cls.uc=cls.cn.copy()
        cls.uc["population_UC_R2024B_status"]=["OK"]*50000+["NODATA"]*(n-50000)
        cls.uc["population_UC_R2024B_density_per_km2"]=[0.0]*50000+[np.nan]*(n-50000)
        cls.uc["UC_source_complete"]=cls.uc.viirs_valid_months.ge(3)&cls.uc.population_UC_R2024B_status.eq("OK")
        cls.receipt={"source_releases_NOT_mixed_in_one_covariate":True,
           "original_2024_constrained_2024_values_SHA256":m.INPUT_SHA["CN_2024_environment"]}
    def test_full_original_photo_source_passes_without_mixing(self):
        self.assertTrue(m.validate_sources(self.cn,self.uc,self.receipt))
    def test_duplicate_target_id_rejected(self):
        x=self.uc.copy()
        x.loc[0,"target_id"]=x.target_id.iloc[1]
        with self.assertRaisesRegex(ValueError,"duplicated"):
            m.validate_sources(self.cn,x,self.receipt)
    def test_changing_2024_viirs_original_cell_rejected(self):
        x=self.uc.copy()
        x.loc[1,"viirs_median"]=1000
        with self.assertRaisesRegex(ValueError,"same-year source changed"):
            m.validate_sources(self.cn,x,self.receipt)
    def test_wrong_release_provenance_rejected(self):
        bad=self.receipt.copy();bad["source_releases_NOT_mixed_in_one_covariate"]=False
        with self.assertRaisesRegex(ValueError,"release cannot be mixed"):
            m.validate_sources(self.cn,self.uc,bad)
    def test_uc_missing_source_never_imputed_to_zero(self):
        x=self.uc.copy()
        x.loc[x.population_UC_R2024B_status.ne("OK"),"population_UC_R2024B_density_per_km2"]=0
        with self.assertRaisesRegex(ValueError,"NoData has a numerical value"):
            m.validate_sources(self.cn,x,self.receipt)
    def test_uc_pop_must_be_nonnegative(self):
        x=self.uc.copy();x.loc[0,"population_UC_R2024B_density_per_km2"]=-1
        with self.assertRaisesRegex(ValueError,"Negative"):
            m.validate_sources(self.cn,x,self.receipt)

if __name__=="__main__":unittest.main()
