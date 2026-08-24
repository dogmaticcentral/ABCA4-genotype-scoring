"""End-to-end test of the SQLite scoring pipeline.

The committed data/abca4_pub70_test.db holds the (unscored) publication-70
cohort extracted from the original MySQL database. tests/reference_scores.json
holds, for each genotype, the score computed by the ORIGINAL MySQL-based code,
in both dosage-compensation modes, together with the per-allele scoring
provenance. The test copies the db to a temp file, runs score_n_store.py on it
(both modes), and checks that the port reproduces the reference exactly.

Run with:  python3 -m unittest discover tests    (or pytest, if installed)
"""

import json
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DB = REPO_ROOT / "data" / "abca4_pub70_test.db"
REFERENCE = REPO_ROOT / "tests" / "reference_scores.json"
TOLERANCE = 1.e-9


class TestScoreNStore(unittest.TestCase):

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

        db = sqlite3.connect(cls.db_path)
        cls.scores = {str(row[0]): (row[1], row[2]) for row in
                      db.execute("select id, score, score_w_dosage_compensation from genotypes")}
        cls.scoring_source = {str(row[0]): json.loads(row[1]) for row in
                              db.execute("select allele_id, source from allele_scoring_source")}
        db.close()

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def assert_score_matches(self, genotype_id, computed, expected, label):
        if expected is None:
            self.assertIsNone(computed, f"genotype {genotype_id}: {label} should be NULL, got {computed}")
        else:
            self.assertIsNotNone(computed, f"genotype {genotype_id}: {label} is NULL, expected {expected}")
            self.assertAlmostEqual(computed, expected, delta=TOLERANCE,
                                   msg=f"genotype {genotype_id}: {label} {computed} != reference {expected}")

    def test_all_reference_genotypes_present(self):
        self.assertEqual(set(self.scores.keys()), set(self.reference.keys()))

    def test_scores_match_original_pipeline(self):
        for genotype_id, entry in self.reference.items():
            computed_score, _ = self.scores[genotype_id]
            self.assert_score_matches(genotype_id, computed_score, entry["score"], "score")

    def test_scores_with_dosage_compensation_match_original_pipeline(self):
        for genotype_id, entry in self.reference.items():
            _, computed_score = self.scores[genotype_id]
            self.assert_score_matches(genotype_id, computed_score,
                                      entry["score_w_dosage_compensation"], "score_w_dosage_compensation")

    def test_some_genotypes_actually_scored(self):
        # guard against the vacuous pass where nothing gets scored on either side
        scored = [g for g, (s, _) in self.scores.items() if s is not None]
        self.assertGreater(len(scored), 50)

    def test_scoring_provenance_matches_original_pipeline(self):
        for genotype_id, entry in self.reference.items():
            for allele_id, expected_source in entry["scoring_source"].items():
                self.assertIn(allele_id, self.scoring_source,
                              f"allele {allele_id} missing from allele_scoring_source")
                self.assertEqual(self.scoring_source[allele_id], expected_source,
                                 f"allele {allele_id}: scoring source differs from the original pipeline")


if __name__ == "__main__":
    unittest.main()
