#!/usr/bin/env python

import unittest
import shutil
import tempfile
import os
from CIME.user_mod_support import apply_user_mods
from CIME.utils import CIMEError
import six

# ========================================================================
# Define some parameters
# ========================================================================

_SOURCEMODS = os.path.join("SourceMods", "src.drv")

class TestUserModSupport(unittest.TestCase):

    # ========================================================================
    # Test helper functions
    # ========================================================================

    def setUp(self):
        self._caseroot = tempfile.mkdtemp()
        self._caseroot_sourcemods = os.path.join(self._caseroot, _SOURCEMODS)
        os.makedirs(self._caseroot_sourcemods)
        self._user_mods_parent_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._caseroot, ignore_errors=True)
        shutil.rmtree(self._user_mods_parent_dir, ignore_errors=True)

    def createUserMod(self, name, include_dirs=None):
        """Create a user_mods directory with the given name.

        This directory is created within self._user_mods_parent_dir

        For name='foo', it will contain:

        - A user_nl_cpl file with contents:
          foo

        - A shell_commands file with contents:
          echo foo >> /PATH/TO/CASEROOT/shell_commands_result

        - A file in _SOURCEMODS named myfile.F90 with contents:
          foo

        If include_dirs is given, it should be a list of strings, giving names
        of other user_mods directories to include. e.g., if include_dirs is
        ['foo1', 'foo2'], then this will create a file 'include_user_mods' that
        contains paths to the 'foo1' and 'foo2' user_mods directories, one per
        line.
        """

        mod_dir = os.path.join(self._user_mods_parent_dir, name)
        os.makedirs(mod_dir)
        mod_dir_sourcemods = os.path.join(mod_dir, _SOURCEMODS)
        os.makedirs(mod_dir_sourcemods)

        with open(os.path.join(mod_dir, "user_nl_cpl"), "w") as user_nl_cpl:
            user_nl_cpl.write(name + "\n")
        with open(os.path.join(mod_dir, "shell_commands"), "w") as shell_commands:
            command = "echo {} >> {}/shell_commands_result\n".format(name, self._caseroot)
            shell_commands.write(command)
        with open(os.path.join(mod_dir_sourcemods, "myfile.F90"), "w") as f90_file:
            f90_file.write(name + "\n")

        if include_dirs:
            with open(os.path.join(mod_dir, "include_user_mods"), "w") as include_user_mods:
                for one_include in include_dirs:
                    include_user_mods.write(os.path.join(self._user_mods_parent_dir, one_include) + "\n")

    def assertResults(self, expected_user_nl_cpl,
                      expected_shell_commands_result,
                      expected_sourcemod,
                      msg = ""):
        """Asserts that the contents of the files in self._caseroot match expectations

        If msg is provided, it is printed for some failing assertions
        """

        path_to_user_nl_cpl = os.path.join(self._caseroot, "user_nl_cpl")
        self.assertTrue(os.path.isfile(path_to_user_nl_cpl),
                        msg = msg + ": user_nl_cpl does not exist")
        with open(path_to_user_nl_cpl, "r") as user_nl_cpl:
            contents = user_nl_cpl.read()
            self.assertEqual(expected_user_nl_cpl, contents)

        path_to_shell_commands_result = os.path.join(self._caseroot, "shell_commands_result")
        self.assertTrue(os.path.isfile(path_to_shell_commands_result),
                        msg = msg + ": shell_commands_result does not exist")
        with open(path_to_shell_commands_result, "r") as shell_commands_result:
            contents = shell_commands_result.read()
            self.assertEqual(expected_shell_commands_result, contents)

        path_to_sourcemod = os.path.join(self._caseroot_sourcemods, "myfile.F90")
        self.assertTrue(os.path.isfile(path_to_sourcemod),
                        msg = msg + ": sourcemod file does not exist")
        with open(path_to_sourcemod, "r") as sourcemod:
            contents = sourcemod.read()
            self.assertEqual(expected_sourcemod, contents)

    # ========================================================================
    # Begin actual tests
    # ========================================================================

    def test_basic(self):
        self.createUserMod("foo")
        apply_user_mods(self._caseroot,
                        os.path.join(self._user_mods_parent_dir, "foo"))
        self.assertResults(expected_user_nl_cpl = "foo\n",
                           expected_shell_commands_result = "foo\n",
                           expected_sourcemod = "foo\n",
                           msg = "test_basic")

    def test_keepexe(self):
        self.createUserMod("foo")
        with six.assertRaisesRegex(self, CIMEError, "cannot have any source mods"):
            apply_user_mods(self._caseroot,
                            os.path.join(self._user_mods_parent_dir, "foo"), keepexe=True)

    def test_two_applications(self):
        """If apply_user_mods is called twice, the second should appear after the first so that it takes precedence."""

        self.createUserMod("foo1")
        self.createUserMod("foo2")
        apply_user_mods(self._caseroot,
                        os.path.join(self._user_mods_parent_dir, "foo1"))
        apply_user_mods(self._caseroot,
                        os.path.join(self._user_mods_parent_dir, "foo2"))
        self.assertResults(expected_user_nl_cpl = "foo1\nfoo2\n",
                           expected_shell_commands_result = "foo1\nfoo2\n",
                           expected_sourcemod = "foo2\n",
                           msg = "test_two_applications")

    def test_include(self):
        """If there is an included mod, the main one should appear after the included one so that it takes precedence."""

        self.createUserMod("base")
        self.createUserMod("derived", include_dirs=["base"])

        apply_user_mods(self._caseroot,
                        os.path.join(self._user_mods_parent_dir, "derived"))

        self.assertResults(expected_user_nl_cpl = "base\nderived\n",
                           expected_shell_commands_result = "base\nderived\n",
                           expected_sourcemod = "derived\n",
                           msg = "test_include")

    def test_duplicate_includes(self):
        """Test multiple includes, where both include the same base mod.

        The base mod should only be included once.
        """

        self.createUserMod("base")
        self.createUserMod("derived1", include_dirs=["base"])
        self.createUserMod("derived2", include_dirs=["base"])
        self.createUserMod("derived_combo",
                           include_dirs = ["derived1", "derived2"])

        apply_user_mods(self._caseroot,
                        os.path.join(self._user_mods_parent_dir, "derived_combo"))

        # NOTE(wjs, 2017-04-15) The ordering of derived1 vs. derived2 is not
        # critical here: If this aspect of the behavior changes, the
        # expected_contents can be changed to match the new behavior in this
        # respect.
        expected_contents = """base
derived2
derived1
derived_combo
"""
        self.assertResults(expected_user_nl_cpl = expected_contents,
                           expected_shell_commands_result = expected_contents,
                           expected_sourcemod = "derived_combo\n",
                           msg = "test_duplicate_includes")
