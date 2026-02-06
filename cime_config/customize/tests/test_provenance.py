#!/usr/bin/env python3

import os
import sys
import tempfile
import unittest
from unittest import mock

from CIME import utils

import provenance


# pylint: disable=protected-access
class TestProvenance(unittest.TestCase):
    def test_parse_dot_git_path_error(self):
        with self.assertRaises(utils.CIMEError):
            provenance._parse_dot_git_path("/src/CIME")

    def test_parse_dot_git_path(self):
        value = provenance._parse_dot_git_path("/src/CIME/.git/worktrees/test")

        assert value == "/src/CIME/.git"

    def test_read_gitdir(self):
        with tempfile.TemporaryDirectory() as tempdir:
            dot_git_path = os.path.join(tempdir, ".git")

            with open(dot_git_path, "w") as fd:
                fd.write("gitdir:     /src/CIME/.git/worktrees/test")

            value = provenance._read_gitdir(dot_git_path)

            assert value == "/src/CIME/.git/worktrees/test"

            with open(dot_git_path, "w") as fd:
                fd.write("gitdir:/src/CIME/.git/worktrees/test")

            value = provenance._read_gitdir(dot_git_path)

            assert value == "/src/CIME/.git/worktrees/test"

    def test_read_gitdir_not_file(self):
        with tempfile.TemporaryDirectory() as tempdir:
            dot_git_path = os.path.join(tempdir, ".git")

            os.makedirs(dot_git_path)

            with self.assertRaises(utils.CIMEError):
                provenance._read_gitdir(dot_git_path)

    def test_read_gitdir_bad_contents(self):
        with tempfile.TemporaryDirectory() as tempdir:
            dot_git_path = os.path.join(tempdir, ".git")

            with open(dot_git_path, "w") as fd:
                fd.write("some value: /src/CIME/.git/worktrees/test")

            with self.assertRaises(utils.CIMEError):
                provenance._read_gitdir(dot_git_path)

    def test_find_git_root(self):
        with tempfile.TemporaryDirectory() as tempdir:
            os.makedirs(os.path.join(tempdir, ".git"))

            value = provenance._find_git_root(tempdir)

            assert value == f"{tempdir}/.git"

    def test_find_git_root_does_not_exist(self):
        with tempfile.TemporaryDirectory() as tempdir:
            with self.assertRaises(utils.CIMEError):
                provenance._find_git_root(tempdir)

    def test_find_git_root_submodule(self):
        with tempfile.TemporaryDirectory() as tempdir:
            cime_path = os.path.join(tempdir, "cime")

            os.makedirs(cime_path)

            cime_git_dot_file = os.path.join(cime_path, ".git")

            with open(cime_git_dot_file, "w") as fd:
                fd.write(f"gitdir: {tempdir}/.git/modules/cime")

            temp_dot_git_path = os.path.join(tempdir, ".git", "modules", "cime")

            os.makedirs(temp_dot_git_path)

            temp_config = os.path.join(temp_dot_git_path, "config")

            with open(temp_config, "w") as fd:
                fd.write("")

            value = provenance._find_git_root(cime_path)

            assert value == f"{tempdir}/.git/modules/cime"

            # relative path
            with open(cime_git_dot_file, "w") as fd:
                fd.write(f"gitdir: ../.git/modules/cime")

            value = provenance._find_git_root(cime_path)

            assert value == f"{tempdir}/.git/modules/cime"

    def test_find_git_root_worktree(self):
        with tempfile.TemporaryDirectory() as tempdir:
            git_dot_file = os.path.join(tempdir, ".git")

            with open(git_dot_file, "w") as fd:
                fd.write("gitdir: /src/CIME/.git/worktrees/test")

            value = provenance._find_git_root(tempdir)

            assert value == "/src/CIME/.git"

    def test_find_git_root_worktree_bad_contents(self):
        with tempfile.TemporaryDirectory() as tempdir:
            with open(os.path.join(tempdir, ".git"), "w") as fd:
                fd.write("some value: /src/CIME/.git/worktrees/test")

            with self.assertRaises(utils.CIMEError):
                provenance._find_git_root(tempdir)

    @mock.patch("CIME.utils.run_cmd")
    def test_run_git_cmd_recursively(self, run_cmd):
        run_cmd.return_value = (0, "data", None)

        with mock.patch("provenance.open", mock.mock_open()) as m:
            provenance._run_git_cmd_recursively(
                "status", "/srcroot", "/output.txt"
            )  # pylint: disable=protected-access

        m.assert_called_with("/output.txt", "w")

        write = m.return_value.__enter__.return_value.write

        write.assert_any_call("data\n\n")
        write.assert_any_call("data\n")

        run_cmd.assert_any_call("git status", from_dir="/srcroot")
        run_cmd.assert_any_call(
            'git submodule foreach --recursive "git status; echo"', from_dir="/srcroot"
        )

    @mock.patch("CIME.utils.run_cmd")
    def test_run_git_cmd_recursively_error(self, run_cmd):
        run_cmd.return_value = (1, "data", "error")

        with mock.patch("provenance.open", mock.mock_open()) as m:
            provenance._run_git_cmd_recursively(
                "status", "/srcroot", "/output.txt"
            )  # pylint: disable=protected-access

        m.assert_called_with("/output.txt", "w")

        write = m.return_value.__enter__.return_value.write

        write.assert_any_call("error\n\n")
        write.assert_any_call("error\n")

        run_cmd.assert_any_call("git status", from_dir="/srcroot")
        run_cmd.assert_any_call(
            'git submodule foreach --recursive "git status; echo"', from_dir="/srcroot"
        )

    @mock.patch("CIME.utils.safe_copy")
    @mock.patch("CIME.utils.run_cmd")
    def test_record_git_provenance(self, run_cmd, safe_copy):
        run_cmd.return_value = (0, "data", None)

        with mock.patch("provenance.open", mock.mock_open()) as m:
            with tempfile.TemporaryDirectory() as tempdir:
                os.makedirs(os.path.join(tempdir, ".git"))

                provenance._record_git_provenance(
                    tempdir, "/output", "5"
                )  # pylint: disable=protected-access

        m.assert_any_call("/output/GIT_STATUS.5", "w")
        m.assert_any_call("/output/GIT_DIFF.5", "w")
        m.assert_any_call("/output/GIT_LOG.5", "w")
        m.assert_any_call("/output/GIT_REMOTE.5", "w")

        write = m.return_value.__enter__.return_value.write

        write.assert_any_call("data\n\n")
        write.assert_any_call("data\n")

        run_cmd.assert_any_call("git status", from_dir=tempdir)
        run_cmd.assert_any_call(
            'git submodule foreach --recursive "git status; echo"', from_dir=tempdir
        )
        run_cmd.assert_any_call("git diff", from_dir=tempdir)
        run_cmd.assert_any_call(
            'git submodule foreach --recursive "git diff; echo"', from_dir=tempdir
        )
        run_cmd.assert_any_call(
            "git log --first-parent --pretty=oneline -n 5", from_dir=tempdir
        )
        run_cmd.assert_any_call(
            'git submodule foreach --recursive "git log --first-parent'
            ' --pretty=oneline -n 5; echo"',
            from_dir=tempdir,
        )
        run_cmd.assert_any_call("git remote -v", from_dir=tempdir)
        run_cmd.assert_any_call(
            'git submodule foreach --recursive "git remote -v; echo"', from_dir=tempdir
        )

        safe_copy.assert_any_call(
            f"{tempdir}/.git/config", "/output/GIT_CONFIG.5", preserve_meta=False
        )

    @mock.patch("provenance.logger")
    def test_record_timing_invalid_timing_dir(self, mock_logger):
        """Test that logger.warning properly handles RuntimeError without TypeError"""
        mock_case = mock.Mock()
        mock_case.get_value.side_effect = lambda key, **kwargs: {
            "PROJECT": "test_project",
            "RUNDIR": "/run",
            "EXEROOT": "/exe",
            "CASEROOT": "/case",
            "SRCROOT": "/src",
            "CASE": "test_case",
            "SAVE_TIMING_DIR": "/invalid/path"
        }.get(key, None)
        mock_case.is_save_timing_dir_project.return_value = True
        
        # This should not raise TypeError anymore
        provenance._record_timing(mock_case, "test_lid")
        
        # Verify that logger.warning was called with the error message
        mock_logger.warning.assert_called_once()
        # Get the actual warning message
        call_args = mock_logger.warning.call_args[0][0]
        # Should contain the error message about invalid SAVE_TIMING_DIR
        assert "SAVE_TIMING_DIR" in call_args
        assert "not valid" in call_args


if __name__ == "__main__":
    unittest.main()
