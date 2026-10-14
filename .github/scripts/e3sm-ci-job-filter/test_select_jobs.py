"""Unit tests for select_jobs.py. Run with: python3 -m unittest discover .github/scripts/e3sm-ci-job-filter"""

import unittest

from select_jobs import parse_names, path_matches, select_jobs

JOBS = ["gnu", "gnugpu", "intel"]
PR = {"pull_request": {"number": 1, "draft": False}}
DRAFT = {"pull_request": {"number": 1, "draft": True}}

class TestParse(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_names("a, b\n c"), ["a", "b", "c"])
        with self.assertRaises(ValueError):
            parse_names("a a")

class TestPaths(unittest.TestCase):
    def test_globs(self):
        self.assertTrue(path_matches("a/b/c.txt", ["a/**"]))
        self.assertTrue(path_matches("a/c.txt", ["a/*"]))
        self.assertFalse(path_matches("a/b/c.txt", ["a/*"]))
        self.assertTrue(path_matches("x/y/z.F90", ["**/*.F90"]))
        self.assertTrue(path_matches("z.F90", ["**/*.F90"]))
        self.assertTrue(path_matches("externals/ekat", ["externals/ekat"]))
        self.assertFalse(path_matches("externals/ekat2", ["externals/ekat"]))
        self.assertFalse(path_matches("a.b", ["a?c"]))

    def test_negation_last_match_wins(self):
        pats = ["components/eamxx/**", "!components/eamxx/docs/**"]
        self.assertTrue(path_matches("components/eamxx/src/a.cpp", pats))
        self.assertFalse(path_matches("components/eamxx/docs/a.md", pats))
        # A later positive pattern re-includes
        self.assertTrue(path_matches("components/eamxx/docs/a.md", pats + ["**/a.md"]))

class TestSelect(unittest.TestCase):
    def sel(self, name, event, **kw):
        return select_jobs(JOBS, name, event, **kw)[0]

    def test_events(self):
        self.assertEqual(self.sel("schedule", {}), JOBS)
        self.assertEqual(self.sel("push", {}), JOBS)
        self.assertEqual(self.sel("pull_request", PR), JOBS)

    def test_draft(self):
        self.assertEqual(self.sel("pull_request", DRAFT), [])
        self.assertEqual(self.sel("pull_request", DRAFT, skip_drafts=False), JOBS)
        # drafts do not affect other events
        self.assertEqual(self.sel("schedule", DRAFT), JOBS)

    def test_dispatch(self):
        kw = dict(dispatch_input="job_to_run")
        self.assertEqual(self.sel("workflow_dispatch", {"inputs": {"job_to_run": "all"}}, **kw), JOBS)
        self.assertEqual(self.sel("workflow_dispatch", {"inputs": {}}, **kw), JOBS)
        self.assertEqual(self.sel("workflow_dispatch", {}, **kw), JOBS)
        # exact match: 'gnu' must not select 'gnugpu'
        self.assertEqual(self.sel("workflow_dispatch", {"inputs": {"job_to_run": "gnu"}}, **kw), ["gnu"])
        self.assertEqual(self.sel("workflow_dispatch", {"inputs": {"job_to_run": "intel"}}, **kw), ["intel"])
        with self.assertRaises(ValueError):
            self.sel("workflow_dispatch", {"inputs": {"job_to_run": "nope"}}, **kw)
        # No dispatch-input configured: run all
        self.assertEqual(self.sel("workflow_dispatch", {"inputs": {"job_to_run": "gnu"}}), JOBS)

    def test_paths(self):
        pats = ["components/eamxx/**", "!components/eamxx/docs/**"]
        def run(files, **kw):
            return self.sel("pull_request", PR, patterns=pats, get_files=lambda: files, **kw)
        self.assertEqual(run(["components/eamxx/src/a.cpp"]), JOBS)
        self.assertEqual(run(["README.md", "components/eamxx/docs/a.md"]), [])
        self.assertEqual(run(["README.md", "components/eamxx/src/a.cpp"]), JOBS)
        # Self paths and the caller workflow are always relevant
        self.assertEqual(run([".github/workflows/e3sm-ci-job-filter.yml"]), JOBS)
        self.assertEqual(run([".github/workflows/me.yml"]), [])
        self.assertEqual(run([".github/workflows/me.yml"], caller_path=".github/workflows/me.yml"), JOBS)
        # Unknown file list (too many files): fail open
        self.assertEqual(self.sel("pull_request", PR, patterns=pats, get_files=lambda: None), JOBS)
        # Paths do not apply to non-PR events
        self.assertEqual(self.sel("schedule", {}, patterns=pats), JOBS)
        # Draft wins over paths
        self.assertEqual(self.sel("pull_request", DRAFT, patterns=pats,
                                  get_files=lambda: ["components/eamxx/a"]), [])

if __name__ == "__main__":
    unittest.main()
