"""SHA-frozen source photo geography recovery without touching strict9 image pixels."""
import csv,hashlib,io,json,tarfile,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.mois_restore_bulk128_photo_geography_from_original_source import (
    run,restore,source_photo_roster,sha
)

class SourcePhotoGeoRestorationTests(unittest.TestCase):
    def fake(self,folder):
        root=Path(folder)
        src=[]
        for i in range(4):
            src.append({
                "photo_id":str(90100+i),
                "observation_id":f"photo-obs-{i}",
                "inat_taxon_id":str(700+i),
                "observed_on":"2025-12-15",
                "latitude":"-12.4",
                "longitude":"125.8",
                "source_arm":("mois_bulk128_tier10_photo_reference_v1" if i==0 else
                              "mois_bulk128_tier10_taxon_region_season_gap_v1")
            })
        workerbytes={}
        for unit in (0,1):
            stream=io.StringIO()
            writer=csv.DictWriter(stream,fieldnames=list(src[0]),
                                  lineterminator="\n")
            writer.writeheader();writer.writerows(src[unit*2:unit*2+2])
            workerbytes[f"shards/shard-{unit:03d}.csv"]=stream.getvalue().encode()
        manifest={
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "sampling_schema":"mois_original93m_tier10_bulk128_probability_reference_plus_gap_v1",
            "worker_ids":[0,1],
            "photo_shard_sha256":{str(i):hashlib.sha256(
                workerbytes[f"shards/shard-{i:03d}.csv"]).hexdigest() for i in (0,1)}
        }
        man=root/"source.json";man.write_text(json.dumps(manifest))
        tar=root/"source.tar.gz"
        with tarfile.open(tar,"w:gz") as archive:
            for name,data in workerbytes.items():
                info=tarfile.TarInfo(name);info.size=len(data)
                archive.addfile(info,io.BytesIO(data))
        world=[]
        for i in range(7):
            bulk=i<4
            world.append({
                "photo_id":str(90100+i) if bulk else str(10000+i),
                "observation_id":f"photo-obs-{i}",
                "inat_taxon_id":str(700+i),
                "source_design":"bulk128_131072" if bulk else
                                "historical_K1_FCP_source_heterogeneous",
                "source_subarm":("reference" if i==0 else "targeted") if bulk else "historical",
                "observed_on":"" if bulk else "2023-05-10",
                "latitude":np.nan if bulk else 40.0,
                "longitude":np.nan if bulk else 140.0,
                "strict9":i%2==0,
                "strict9_white":1.0 if i%2==0 else np.nan,
            })
        photo=root/"world.parquet";pd.DataFrame(world).to_parquet(photo,index=False)
        return photo,man,tar

    def test_original_obs_geo_and_date_restored_only_on_bulk_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            p,m,t=self.fake(directory)
            a=run(p,m,t,Path(directory)/"result",verify_sha=False,
                  expected_world=7,expected_bulk=4,expected_units=2)
            self.assertEqual(a["new_bulk_original_photo_IDs_restored"],4)
            self.assertEqual(a["prior_bulk_source_geopoints_known"],0)
            self.assertEqual(a["source_original_bulk_photo_geopoints_available"],4)
            self.assertEqual(a["bulk_original_photo_observation_days_nonempty"],4)
            self.assertEqual(a["real_original_strict9_numerical_images"],4)
            d=pd.read_parquet(Path(directory)/"result"/
                "mois_1181146_true_original_strict9_geography_restored_v2.parquet")
            self.assertTrue(d.loc[:3,"latitude"].eq(-12.4).all())
            self.assertTrue(d.loc[:3,"longitude"].eq(125.8).all())
            self.assertTrue(d.loc[:3,"observed_on"].eq("2025-12-15").all())
            self.assertTrue(d.loc[4:,"latitude"].eq(40.0).all())
            self.assertTrue(d.loc[4:,"observed_on"].eq("2023-05-10").all())
            self.assertTrue(d.strict9.equals(pd.read_parquet(p).strict9))
            self.assertTrue(d.strict9_white.equals(pd.read_parquet(p).strict9_white))

    def test_source_photo_taxon_conflict_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            p,m,t=self.fake(folder)
            df=pd.read_parquet(p)
            src,_=source_photo_roster(m,t,unit_count=2)
            df.loc[0,"inat_taxon_id"]="WRONG_PLANT_TAXON"
            with self.assertRaisesRegex(ValueError,"inconsistent source observation or plant taxon"):
                restore(df,src,actual_n=4)
            df=pd.read_parquet(p)
            df.loc[0,"source_subarm"]="targeted"
            with self.assertRaisesRegex(ValueError,"probability-reference"):
                restore(df,src,actual_n=4)

if __name__=="__main__":
    unittest.main()
