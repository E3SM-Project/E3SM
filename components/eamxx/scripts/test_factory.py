# This module contains classes that describe a test type for test-all-scream
# Each test type (here represented by a TestProperty object) can have different
# flags, build type, profiling, cmake options, etc.
# The function "test_factory" can be used to get a list of test types from
# their string representation.

from collections import OrderedDict
from utils import expect

###############################################################################
class TestProperty(object):
###############################################################################

    """
    Parent class of predefined test types for SCREAM standalone. test-all-scream
    offers a number of customization points, but you may need to just use
    cmake if you need maximal customization. You can run test-all-scream --dry-run
    to get the corresponding cmake command which can then be used as a starting
    point for making your own cmake command.
    """

    def __init__(self, longname, description, cmake_args,
                 uses_baselines=True, on_by_default=True, default_test_len=None):
        # What the user uses to select tests via test-all-scream CLI.
        # Should also match the class name when converted to caps
        self.shortname      = type(self).__name__.lower()

        # A longer name used to name baseline and test directories for a test.
        # Also used in output/error messages to refer to the test
        self.longname       = longname

        # A longer decription of the test
        self.description    = description

        # Cmake config args for this test. Check that quoting is done with
        # single quotes.
        self.cmake_args     = cmake_args
        for name, arg in self.cmake_args:
            expect('"' not in arg,
                   f"In test definition for {longname}, found cmake args with double quotes {name}='{arg}'"
                   "Please use single quotes if quotes are needed.")

        # Does the test do baseline testing
        self.uses_baselines = uses_baselines

        # Should this test be run if the user did not specify tests at all?
        self.on_by_default  = on_by_default

        # Should this test have a default test size
        self.default_test_len = default_test_len

        #
        # Properties not set by constructor (Set by the main TestAllScream object)
        #

        # Resources used by this test.
        self.compile_res_count = None
        self.testing_res_count = None

        # Does this test need baselines
        self.baselines_missing = False
        self.baselines_expired = False

        #
        # Common
        #

        if not self.uses_baselines:
            self.cmake_args += [("SCREAM_ENABLE_BASELINE_TESTS", "False")]

    def disable_baselines(self):
        if self.uses_baselines:
            self.uses_baselines = False
            self.cmake_args += [("SCREAM_ENABLE_BASELINE_TESTS", "False")]

    # Tests will generally be referred to via their longname
    def __str__(self):
        return self.longname

###############################################################################
class DBG(TestProperty):
###############################################################################

    def __init__(self, _):
        TestProperty.__init__(
            self,
            "full_debug",
            "debug",
            [("CMAKE_BUILD_TYPE", "Debug"), ("EKAT_DEFAULT_BFB", "True"),
             ("Kokkos_ENABLE_DEBUG_BOUNDS_CHECK", "True")]
        )

###############################################################################
class SP(TestProperty):
###############################################################################

    def __init__(self, _):
        TestProperty.__init__(
            self,
            "full_sp_debug",
            "debug single precision",
            [("CMAKE_BUILD_TYPE", "Debug"), ("EKAT_DEFAULT_BFB", "True"),
             ("SCREAM_DOUBLE_PRECISION", "False")],
        )

###############################################################################
class FPE(TestProperty):
###############################################################################

    def __init__(self, tas):
        TestProperty.__init__(
            self,
            "debug_nopack_fpe",
            "debug pksize=1 floating point exceptions on",
            [("CMAKE_BUILD_TYPE", "Debug"), ("EKAT_DEFAULT_BFB", "True"),
             ("SCREAM_PACK_SIZE", "1"), ("SCREAM_FPE","True")],
            uses_baselines=False,
            on_by_default=(tas is not None and not tas._machine.uses_gpu())
        )

###############################################################################
class OPT(TestProperty):
###############################################################################

    def __init__(self, _):
        TestProperty.__init__(
            self,
            "release",
            "release",
            [("CMAKE_BUILD_TYPE", "Release")],
        )

###############################################################################
class COV(TestProperty):
###############################################################################

    def __init__(self, _):
        TestProperty.__init__(
            self,
            "coverage",
            "debug coverage",
            [("CMAKE_BUILD_TYPE", "Debug"), ("EKAT_ENABLE_COVERAGE", "True")],
            uses_baselines=False,
            on_by_default=False,
            default_test_len="short"
        )

###############################################################################
class VALG(TestProperty):
###############################################################################

    def __init__(self, tas):
        TestProperty.__init__(
            self,
            "valgrind",
            "Release build where tests run through valgrind",
            [("CMAKE_BUILD_TYPE", "RelWithDebInfo"),
             ("EKAT_ENABLE_VALGRIND", "True"),
             ("SCREAM_TEST_MAX_THREADS", "2")],
            uses_baselines=False,
            on_by_default=False,
            default_test_len="short"
        )
        if tas is not None:
            # If a stored suppression file exists for this machine, use it
            persistent_supp_file = tas.get_root_dir() / "scripts" / "jenkins" / "valgrind" / f"{tas.get_machine().name}.supp"
            if persistent_supp_file.exists():
                self.cmake_args.append( ("EKAT_VALGRIND_SUPPRESSION_FILE", str(persistent_supp_file)) )

###############################################################################
class CSM(TestProperty):
###############################################################################

    def __init__(self, _):
        TestProperty.__init__(
            self,
            "compute_sanitizer_memcheck",
            "debug with compute sanitizer memcheck",
            [("CMAKE_BUILD_TYPE", "Debug"),
             ("EKAT_ENABLE_COMPUTE_SANITIZER", "True"),
             ("EKAT_COMPUTE_SANITIZER_OPTIONS", "'--tool=memcheck'")],
            uses_baselines=False,
            on_by_default=False,
            default_test_len="short"
        )

###############################################################################
class CSR(TestProperty):
###############################################################################

    def __init__(self, _):
        TestProperty.__init__(
            self,
            "compute_sanitizer_racecheck",
            "debug with compute sanitizer racecheck",
            [("CMAKE_BUILD_TYPE", "Debug"),
             ("EKAT_ENABLE_COMPUTE_SANITIZER", "True"),
             ("EKAT_COMPUTE_SANITIZER_OPTIONS", "'--tool=racecheck --racecheck-detect-level=error'")],
            uses_baselines=False,
            on_by_default=False,
            default_test_len="short"
        )

###############################################################################
class CSI(TestProperty):
###############################################################################

    def __init__(self, _):
        TestProperty.__init__(
            self,
            "compute_sanitizer_initcheck",
            "debug with compute sanitizer initcheck",
            [("CMAKE_BUILD_TYPE", "Debug"),
             ("EKAT_ENABLE_COMPUTE_SANITIZER", "True"),
             ("EKAT_COMPUTE_SANITIZER_OPTIONS", "'--tool=initcheck'")],
            uses_baselines=False,
            on_by_default=False,
            default_test_len="short"
        )

###############################################################################
class CSS(TestProperty):
###############################################################################

    def __init__(self, _):
        TestProperty.__init__(
            self,
            "compute_sanitizer_synccheck",
            "debug with compute sanitizer synccheck",
            [("CMAKE_BUILD_TYPE", "Debug"),
             ("EKAT_ENABLE_COMPUTE_SANITIZER", "True"),
             ("EKAT_COMPUTE_SANITIZER_OPTIONS", "'--tool=synccheck'")],
            uses_baselines=False,
            on_by_default=False,
            default_test_len="short"
        )

###############################################################################
def create_tests(user_req_tests, tas):
###############################################################################
    testclasses = TestProperty.__subclasses__()
    if not user_req_tests:
        result = [testclass(tas) for testclass in testclasses
                  if testclass(tas).on_by_default]
    else:
        valid_names = [testclass(tas).shortname for testclass in testclasses]
        for user_req_test in user_req_tests:
            expect(user_req_test in valid_names, f"'{user_req_test}' is not a known test")

        result = [testclass(tas) for testclass in testclasses if testclass(tas).shortname in user_req_tests]

    return result

###########################################################################
def get_test_name_dict():
###########################################################################
    """
    Returns a dict mapping short test names to full names
    """
    testclasses = TestProperty.__subclasses__()
    return OrderedDict([(testc(None).shortname, testc(None).description) for testc in testclasses])
