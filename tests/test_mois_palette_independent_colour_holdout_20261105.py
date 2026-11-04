"""Evidence and leakage checks for independent species-level colour validation."""
import tempfile,unittest
from pathlib import Path
import numpy as np,pandas as pd
from scripts.analysis import mois_palette_independent_colour_holdout_20261105 as m

class PaletteHoldoutTests(unittest.TestCase):
    def test_no_white_or_multistate_assigned_to_inconspicuous(self):
        self.assertEqual(m.LABELS["white"],0)
        self.assertEqual(m.LABELS["green_brown_inconspicuous"],1)
        for n in ("green_brown_inconspicuous|white","mixed_or_variable","unresolved","other_described"):
            self.assertNotIn(n,m.LABELS)
    def test_bronze_is_only_one_of_nine_components(self):
        self.assertEqual(len(m.COMP),9)
        self.assertEqual(m.COMP[-1],"strict9_bronze")
    def test_bad_archive_refused_before_reading_rows(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"fake.zip";f.write_bytes(b"not authenticated")
            with self.assertRaisesRegex(ValueError,"ZIP source hash"):
                m.read_zip(f,"direct")
    def test_metrics_reject_single_label(self):
        with self.assertRaisesRegex(ValueError,"one phenotype"):
            m.measures(np.array([0,0]),np.array([0.1,0.3]))
    def test_perfect_auc(self):
        x=m.measures(np.array([1,0,1,0]),np.array([.9,.1,.8,.2]))
        self.assertEqual(x["n"],4)
        self.assertAlmostEqual(x["auc"],1.)
    def test_predictor_does_not_write_species_to_public_metrics(self):
        rng=np.random.default_rng(20261105)
        n=300
        f=pd.DataFrame({"genus":["G"+str(i//3) for i in range(n)],
          "family":["F"+str(i//15) for i in range(n)],
          "label":np.repeat([0,0,1],100),
          **{q:rng.random(n) for q in m.COMP}})
        s,folds=m.groupfold(f,"genus",(m.COMP[-1],))
        self.assertEqual(len(folds),5)
        self.assertEqual(sum(q["n"] for q in folds),n)
        self.assertNotIn("accepted_id",str(s))
    def test_region_transfer_requires_species_nonoverlap(self):
        all=[]
        for i in range(260):
            sides=("north_temperate","south_temperate") if i<10 else (
              ("north_temperate",) if i<135 else ("south_temperate",))
            for j in sides:
                all.append({"sector":j,"accepted_id":str(i),"accepted_name":f"Genus species{i}",
                 "genus":"Genus"+str(i//4),"strict9_bronze":i/260})
        frame=pd.DataFrame(all)
        y=pd.DataFrame({"accepted_species":[f"Genus species{i}" for i in range(260)],
                        "label":[i%2 for i in range(260)]})
        res=m.transfer(frame,y,"north_temperate","south_temperate",(m.COMP[-1],))
        self.assertEqual(res["shared_species"],0)
        self.assertEqual(res["train_n"],135)
        self.assertEqual(res["n"],125)
    def test_pins_exact_four_independent_artifacts(self):
        self.assertEqual(len(m.FILES),4)
        self.assertTrue(all(len(a[1])==64 and len(a[2])==64 for a in m.FILES.values()))

if __name__=="__main__":unittest.main()
