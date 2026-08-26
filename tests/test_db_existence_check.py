"""Checks that the entry-point scripts refuse to run against a nonexistent --db path.

sqlite3.connect() silently creates an empty database file for a path that
doesn't exist, which would otherwise mask a typo'd --db argument as a
(vacuous) successful run. score_n_store.py and onset_age_vs_score.py both
guard against that before connecting; this test exercises that guard as a
subprocess, the same way test_scoring.py exercises the happy path.

Run with:  python3 -m unittest discover tests    (or pytest, if installed)
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = [str(REPO_ROOT / "score_n_store.py"), str(REPO_ROOT / "onset_age_vs_score.py")]


class TestDbExistenceCheck(unittest.TestCase):

    def run_script(self, script, db_path):
        return subprocess.run([sys.executable, script, "--db", db_path],
                              cwd=REPO_ROOT, capture_output=True, text=True)

    def test_nonexistent_db_exits_nonzero(self):
        for script in SCRIPTS:
            with self.subTest(script=script), tempfile.TemporaryDirectory() as tmpdir:
                missing_path = str(Path(tmpdir) / "does_not_exist.db")
                run = self.run_script(script, missing_path)
                self.assertNotEqual(run.returncode, 0,
                                   f"{script} should not exit successfully on a missing db")

    def test_nonexistent_db_reports_informative_message(self):
        for script in SCRIPTS:
            with self.subTest(script=script), tempfile.TemporaryDirectory() as tmpdir:
                missing_path = str(Path(tmpdir) / "does_not_exist.db")
                run = self.run_script(script, missing_path)
                self.assertIn(missing_path, run.stderr,
                             f"{script}: error message should name the missing db path")
                self.assertNotIn("Traceback", run.stderr,
                                f"{script}: missing db should be handled gracefully, not raise a traceback")

    def test_nonexistent_db_is_not_created(self):
        for script in SCRIPTS:
            with self.subTest(script=script), tempfile.TemporaryDirectory() as tmpdir:
                missing_path = str(Path(tmpdir) / "does_not_exist.db")
                self.run_script(script, missing_path)
                self.assertFalse(Path(missing_path).exists(),
                                f"{script} must not silently create the db file")


if __name__ == "__main__":
    unittest.main()
