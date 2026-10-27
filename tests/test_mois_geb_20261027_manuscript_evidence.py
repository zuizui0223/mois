"""Source evidence / manuscript publication-nonclaim regression.

All numerical assertions use only immutable, source-derived, pre-fitted
empirical results; no new S3 photos, image inference, model fit or global white
prevalence estimation. Do not allow mistaking photo opportunities for images.
"""
import json,math,re,unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"manuscript"
ART=DIR/"MOIS_GEB_RESEARCH_ARTICLE_20261027_MAIN_TEXT.md"
EVIDENCE=DIR/"MOIS_GEB_20261027_EVIDENCE_NUMBERS.json"
FIG=DIR/"MOIS_GEB_20261027_FIGURES_TABLES_AND_SOURCE_BOUNDARIES.md"
GATE=DIR/"MOIS_GEB_20261027_SUBMISSION_GATES_AND_ONE_PAPER_SCOPE.md"

class SourceOriginalFloralEcologyPaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paper=ART.read_text(encoding="utf-8")
        cls.ledger=json.loads(EVIDENCE.read_text(encoding="utf-8"))
        cls.fig=FIG.read_text(encoding="utf-8")
        cls.gates=GATE.read_text(encoding="utf-8")
    def test_journal_geb_abstract_and_manuscript_narrative(self):
        abstract=self.paper.split("## Structured abstract\n",1)[1].split("\n---",1)[0]
        w=len(re.sub(r"\*\*","",abstract).split())
        self.assertLessEqual(w,300,("GEB structured abstract >300 words",w))
        for part in ("Aim","Location","Time period","Major taxa studied",
                     "Methods","Results","Main conclusions"):
            self.assertIn("**"+part+":**",abstract)
        for heading in ("## 1 | Introduction","## 2 | Materials and methods",
                        "## 3 | Results","## 4 | Discussion","## 5 | Conclusions"):
            self.assertIn(heading,self.paper)
        self.assertTrue(all("Figure "+str(i) in self.fig for i in range(1,6)))
        self.assertIn("not submitted",self.paper)
        self.assertIn("different eligibility criterion",self.fig.lower())

    def test_photo_hue_counts_match_real_source_instead_of_opportunity_ids(self):
        d=self.ledger
        self.assertEqual(d["study_frames"][0]["upstream_photo_opportunities"],93209351)
        c=d["visible_four_hue_350399"]
        self.assertEqual(c["image_attempts"],350399)
        self.assertEqual(sum(c["hue_counts"].values()),164297)
        self.assertEqual(c["hue_unknown"]+c["hue_known"],350399)
        self.assertAlmostEqual(100*c["hue_counts"]["white"]/c["hue_known"],
                               c["weighting_white_fraction_percent"]["photo_equal"],places=3)
        frames={x["id"]:x for x in d["study_frames"]}
        self.assertEqual(frames["climate309439"]["measured_original_photo_ids"],309439)
        self.assertEqual(frames["geo313535"]["classifiable_four_hue"],151239)
        self.assertEqual(frames["composition350399"]["classifiable_four_hue"],164297)

    def test_strict9_union_is_deduplicated_not_fitted_climate_1181146(self):
        d=self.ledger
        self.assertEqual(1043930+2048+4096+131072,1181146)
        self.assertEqual(400710+682+1467+48784,451643)
        q=d["strict9_quality1181146"]
        self.assertEqual(q["strict9_valid"]+q["strict9_unknown"],1181146)
        self.assertEqual(q["source_taxa_with_valid9"]+
                         q["source_taxa_with_zero9"],42111)
        self.assertEqual(sum(q["remaining_unmeasured_source_taxa"].values()),12063)
        self.assertEqual(q["original_coordinate_unknown"],195)
        self.assertEqual(q["original_local_season_unknown"],197)
        self.assertFalse(d["inference_gates"]["new_annual_climate_regression_on_1181146_pictures_performed"])
        self.assertIn("not updated environmental regression coefficients",self.paper)
        self.assertIn("451,643",self.paper)

    def test_climate_interval_meaning_no_falsely_identified_local_effect(self):
        d=self.ledger["annual_moisture_genus_equal_309439"]
        self.assertEqual(d["resamples"],999)
        self.assertEqual(d["within_model_simultaneous_family"],24)
        estimates=d["comparisons"]
        self.assertEqual(len(estimates),4)
        for r in estimates:
            self.assertLessEqual(r["ci_low"],r["effect"])
            self.assertGreaterEqual(r["ci_high"],r["effect"])
        self.assertTrue(all(r["ci_low"]>0 for r in estimates[:2]))
        self.assertTrue(all(r["ci_low"]<0<r["ci_high"] for r in estimates[2:]))
        self.assertIn("geographically conditioned",self.paper)

    def test_environment_matching_and_exact_photo_seasons_keep_limits(self):
        d=self.ledger
        s=d["seasonal_350399"]["summer_winter_by_hemisphere"]
        for hemi in ("north","south"):
            row=s[hemi]
            self.assertAlmostEqual(sum(row["four_term_components_pp"].values()),
                                   row["sum_white_difference_pp"],places=4)
        self.assertEqual(d["seasonal_350399"]
          ["taxon_cell_season_year_strict_at_least_two_distinct_obs_each_season"],
          {"north":2,"south":2})
        geography=d["elevation_and_island"]
        self.assertLess(geography["elevation"]["nominal_95"][1],0)
        lo,hi=geography["island"]["region_cluster_nominal_95"]
        self.assertLess(lo,0)
        self.assertGreater(hi,0)
        self.assertIn("ordinary whole-genus",self.paper)
        self.assertIn("at least",self.paper)

    def test_missingness_and_unverified_labels_never_upgraded(self):
        d=self.ledger["strict9_quality1181146"]
        p=d["extra_sampled_all_prior_roi9_failed"]
        self.assertEqual(p["YOLO_500px_flower_candidate"]+
                         p["negative"]+p["error"],p["n"])
        self.assertTrue(p["anatomical_flower_petal_truth_unknown"])
        self.assertTrue(self.ledger["inference_gates"]["technical_unknown_imputed_as_nonwhite"] is False)
        self.assertTrue(self.ledger["inference_gates"]["anthocyanin_cost_or_fitness_tested"] is False)
        self.assertIn("not",self.gates.lower())
        self.assertIn("not",self.fig.lower())
        self.assertIn("128-photo human",self.paper)

if __name__=="__main__":
    unittest.main()
