"""End-to-end test of the onset-age-vs-score analysis (onset_age_vs_score.py).

tests/reference_onset.json holds, for both score columns, the (alias,
onset_age, score) points that survive the original 33_onset_age_vs_score.py
filtering - computed with the ORIGINAL MySQL-based code - and the resulting
Spearman/Pearson correlations. The test scores a temp copy of the committed db
with score_n_store.py (both dosage modes), collects the points with the ported
code, and checks that points and correlations match the reference. It also
runs the script end to end and checks that it writes a plot file.

Run with:  python3 -m unittest discover tests    (or pytest, if installed)
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scipy.stats import spearmanr, pearsonr

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from abca4_score.sqlite_utils import connect  # noqa: E402
from onset_age_vs_score import collect_onset_age_vs_score  # noqa: E402

DATA_DB = REPO_ROOT / "data" / "abca4_pub70_test.db"
REFERENCE = REPO_ROOT / "tests" / "reference_onset.json"
PUBLICATION_ID = "70"
TOLERANCE = 1.e-9


class TestOnsetAgeVsScore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(REFERENCE) as f:
            cls.reference = json.load(f)

        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.db_path = str(Path(cls.tmpdir.name) / "abca4_pub70_test.db")
        shutil.copy(DATA_DB, cls.db_path)

        script = str(REPO_ROOT / "score_n_store.py")
        for extra_args in ([], ["--assume_dosage_compensation"]):
            run = subprocess.run([sys.executable, script, "--db", cls.db_path, "--quiet"] + extra_args,
                                 cwd=REPO_ROOT, capture_output=True, text=True)
            if run.returncode != 0:
                raise RuntimeError(f"score_n_store.py {' '.join(extra_args)} failed:\n"
                                   f"{run.stdout}\n{run.stderr}")

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def collect_points(self, score_column):
        db, cursor = connect(self.db_path)
        scores, onset_ages, _, aliases = \
            collect_onset_age_vs_score(cursor, PUBLICATION_ID, score_column)
        cursor.close()
        db.close()
        points = sorted(zip(aliases, onset_ages, scores), key=lambda p: (p[0] or "", p[1], p[2]))
        return points, scores, onset_ages

    def assert_matches_reference(self, score_column):
        expected = self.reference[score_column]
        points, scores, onset_ages = self.collect_points(score_column)

        self.assertEqual(len(points), len(expected["points"]),
                         f"{score_column}: {len(points)} points, reference has {len(expected['points'])}")
        for (alias, onset_age, score), (ref_alias, ref_onset_age, ref_score) in \
                zip(points, expected["points"]):
            self.assertEqual(alias, ref_alias)
            self.assertEqual(onset_age, ref_onset_age, f"case {alias}: onset age differs")
            self.assertAlmostEqual(score, ref_score, delta=TOLERANCE, msg=f"case {alias}: score differs")

        spearman = spearmanr(scores, onset_ages)
        pearson = pearsonr(scores, onset_ages)
        self.assertAlmostEqual(spearman.statistic, expected["spearman"], delta=TOLERANCE)
        self.assertAlmostEqual(spearman.pvalue, expected["spearman_pvalue"], delta=TOLERANCE)
        self.assertAlmostEqual(pearson.statistic, expected["pearson"], delta=TOLERANCE)
        self.assertAlmostEqual(pearson.pvalue, expected["pearson_pvalue"], delta=TOLERANCE)

    def test_points_and_correlations_match_original_pipeline(self):
        self.assert_matches_reference("score")

    def test_points_and_correlations_match_original_pipeline_w_dosage_compensation(self):
        self.assert_matches_reference("score_w_dosage_compensation")

    def test_reference_is_not_vacuous(self):
        self.assertGreater(len(self.reference["score"]["points"]), 50)

    def test_script_end_to_end_writes_plot(self):
        plot_path = Path(self.tmpdir.name) / "onset_vs_score.png"
        script = str(REPO_ROOT / "onset_age_vs_score.py")
        env = dict(os.environ, MPLBACKEND="Agg")
        run = subprocess.run([sys.executable, script, "--db", self.db_path,
                              "-p", PUBLICATION_ID, "--out", str(plot_path)],
                             cwd=REPO_ROOT, capture_output=True, text=True, env=env)
        self.assertEqual(run.returncode, 0, f"onset_age_vs_score.py failed:\n{run.stdout}\n{run.stderr}")
        self.assertIn("Spearman correlation", run.stdout)
        self.assertIn("Pearson  correlation", run.stdout)
        self.assertTrue(plot_path.exists() and plot_path.stat().st_size > 0)


if __name__ == "__main__":
    unittest.main()
