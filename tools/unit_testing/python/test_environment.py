#!/usr/bin/env python
"""Unit tests for the environment module.

Public classes:
TestEnvSys - Test the EnvSystemInterface abstract base class.
TestNoMod - NoModuleInterface tests.
TestMod - ModuleInterface tests.
TestSoft - SoftEnvInterface tests.
TestExpandEnv - expand_env tests.

To test the actual module system, you must specify these environment
variables when running the test:
MODULE_SYSTEM - Type of module system available locally.
TEST_MODULE - Name of a module which is *not* loaded, to use for testing
              purposes.
MODULE_FILE - File that must be executed to load the interface to python
              (if any).
"""

from os import environ
import unittest

from environment import *

__all__ = ("TestEnvSys", "TestNoMod", "TestMod", "TestSoft",
           "TestExpandEnv")

# Test data gleaned from the environment.
module_system = "none"
test_module = "foo"
module_file = None
if "MODULE_SYSTEM" in environ:
    module_system = environ["MODULE_SYSTEM"]
    assert "TEST_MODULE" in environ, \
        "MODULE_SYSTEM is set, but TEST_MODULE is not defined."
    test_module = environ["TEST_MODULE"]
    if module_system == "module":
        assert "MODULE_FILE" in environ, \
            "MODULE_SYSTEM is module, but MODULE_FILE is not defined."
        module_file = environ["MODULE_FILE"]


class TestEnvSys(unittest.TestCase):

    """Tests for the EnvSystemInterface class.

    This tests that each method is NotImplemented. Other test classes in
    this module inherit from this one in order to detect when a method has
    been implemented but no test is defined for the new method.

    However, figuring out which class causes the error is a bit of a
    guessing game. Addressing that issue is difficult in Python 2.6 because
    of how assertRaises works.
    """

    test_class = EnvSystemInterface

    def setUp(self):
        """Set up by creating an instance of the class under test."""
        self.test_obj = self.test_class()

    def test_is_loaded(self):
        """is_loaded raises NotImplementedError."""
        self.assertRaises(NotImplementedError,
                          self.test_obj.is_loaded, test_module)

    def test_purge(self):
        """purge raises NotImplementedError."""
        self.assertRaises(NotImplementedError,
                          self.test_obj.purge)

    def test_load(self):
        """load raises NotImplementedError."""
        self.assertRaises(NotImplementedError,
                          self.test_obj.load, test_module)

    def test_purge_str(self):
        """purge_str raises NotImplementedError."""
        self.assertRaises(NotImplementedError,
                          self.test_obj.purge_str)

    def test_load_str(self):
        """load_str raises NotImplementedError."""
        self.assertRaises(NotImplementedError,
                          self.test_obj.load_str, test_module)

    def test_unload_str(self):
        """unload_str raises NotImplementedError."""
        self.assertRaises(NotImplementedError,
                          self.test_obj.unload_str, test_module)


class TestNoMod(TestEnvSys):

    """Tests for the NoModuleInterface class."""

    test_class = NoModuleInterface

    def test_purge_str(self):
        """User can call NoModuleInterface.purge_str and get a ":"."""
        self.assertEqual(self.test_obj.purge_str(), ":")


class TestMod(TestEnvSys):

    """Tests for the ModuleInterface class.

    It's not easy to verify that this behaves correctly, so in many cases
    we do some kind of check for self consistency only. If module_system
    is not "module", then most of the tests look for an exception.
    """

    test_class = ModuleInterface

    def setUp(self):
        """Set up module system if present, then call parent setUp."""
        if module_system == "module":
            ModuleInterface.python_init(module_file)
        super(TestMod, self).setUp()

    def tearDown(self):
        """Remove the test module from the current environment."""
        if module_system == "module":
            self.test_obj.unload(test_module)

    def test_is_loaded(self):
        """ModuleInterface.is_loaded returns false for unloaded module."""
        if module_system == "module":
            self.assertFalse(self.test_obj.is_loaded(test_module))
        else:
            self.assertRaises(AssertionError,
                              self.test_obj.is_loaded, test_module)

    def test_purge(self):
        """ModuleInterface.purge removes the test module."""
        if module_system == "module":
            self.test_obj.load(test_module)
            self.test_obj.purge()
            self.assertFalse(self.test_obj.is_loaded(test_module))
        else:
            self.assertRaises(AssertionError,
                              self.test_obj.purge)

    def test_load(self):
        """ModuleInterface.load loads a module."""
        if module_system == "module":
            self.test_obj.load(test_module)
            self.assertTrue(self.test_obj.is_loaded(test_module))
        else:
            self.assertRaises(AssertionError,
                              self.test_obj.load, test_module)

    def test_unload(self):
        """ModuleInterface.unload unloads a module."""
        if module_system == "module":
            self.test_obj.load(test_module)
            self.test_obj.unload(test_module)
            self.assertFalse(self.test_obj.is_loaded(test_module))
        else:
            self.assertRaises(AssertionError,
                              self.test_obj.unload, test_module)

    def test_purge_str(self):
        """User gets non-null string from ModuleInterface.purge_str."""
        self.assertNotEqual(self.test_obj.purge_str().strip(),
                            "")

    def test_load_str(self):
        """User gets non-null string from ModuleInterface.load_str."""
        self.assertNotEqual(self.test_obj.load_str(test_module).strip(),
                            "")

    def test_unload_str(self):
        """User gets non-null string from ModuleInterface.unload_str."""
        self.assertNotEqual(self.test_obj.unload_str(test_module).strip(),
                            "")


class TestSoft(TestEnvSys):

    """Tests for the SoftEnvInterface class.

    It's not easy to verify that this behaves correctly without a specific
    test machine in mind. As with the module tests, we really just verify
    that no exceptions are raised and the return isn't the null string.
    """

    test_class = SoftEnvInterface

    def test_purge_str(self):
        """User gets non-null string from SoftEnvInterface.purge_str."""
        self.assertNotEqual(self.test_obj.purge_str().strip(), "")

    def test_load_str(self):
        """User gets non-null string from SoftEnvInterface.load_str."""
        self.assertNotEqual(self.test_obj.load_str(test_module).strip(),
                            "")

    def test_unload_str(self):
        """User gets non-null string from SoftEnvInterface.unload_str."""
        self.assertNotEqual(self.test_obj.unload_str(test_module).strip(),
                            "")


class TestExpandEnv(unittest.TestCase):

    """Tests for the expand_env function."""

    def test_no_variable(self):
        """With no variables, expand_env is a no-op."""
        self.assertEqual(expand_env("foo", {"UNUSED": "not used"}), "foo")

    def test_variable_missing(self):
        """With variables not present, expand_env is a no-op."""
        self.assertEqual(expand_env("${NOT_HERE}", {}), "${NOT_HERE}")

    def test_brace_expansion(self):
        """Test that an expansion works with curly braces."""
        self.assertEqual(expand_env("${FOO}", {"FOO": "bar"}), "bar")

    def test_bare_expansion(self):
        """Test that an expansion works with no braces."""
        self.assertEqual(expand_env("$FOO", {"FOO": "bar"}), "bar")

    def test_recursive_expansion(self):
        """Test that expansion is done recursively."""
        self.assertEqual(expand_env("${FOO}",
                                    {"FOO": "${FOO2}", "FOO2": "bar"}),
                         "bar")

    def test_brace_closing(self):
        """Test that braces must be closed for expansion to occur."""
        self.assertEqual(expand_env("${FOO", {"FOO": "bar"}), "${FOO")


if __name__ == "__main__":
    unittest.main()
