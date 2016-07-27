"""
Base class for CIME system tests that involve doing two runs and comparing their
output.

If your two runs can use the same executable, then you should inherit directly
from this class. If your two runs need separate executables, then you should
inherit from SystemTestsCompareTwoDiffferentBuilds.

Classes that inherit from this are REQUIRED to implement the following methods:

(1) _setup_first_phase
    This method will be called to set up the run for the first phase of the
    two-phase test

(2) _setup_second_phase
    This method will be called to set up the run for the second phase of the
    two-phase test

In addition, the __init__ method in subclasses MAY set various attributes AFTER
calling SystemTestsCompareTwo.__init__, via calls to the following. However,
this is only for the sake of enhancing communication via logging, file names,
etc. - these don't affect the functionality of the test:

(1) self.set_test_suffix(test_suffix)

(2) self.set_description_first_phase(description)

(3) self.set_description_second_phase(description)

"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class SystemTestsCompareTwo(SystemTestsCommon):

    def __init__(self, case, expected=None):
        """
        Initialize a SystemTestsCompareTwo object.

        Individual test cases that inherit from SystemTestsCompareTwo should
        call this __init__ method first. Afterwards, they may then want to call
        methods to customize the running of the test:

        (1) self.set_test_suffix(test_suffix): sets the test suffix appended to
            netcdf files output by the second run of the test

        (2) self.set_description_first_phase(description): sets the description
        printed to the log file when starting the first run

        (3) self.set_description_second_phase(description): sets the description
        printed to the log file when starting the second run
        """
        SystemTestsCommon.__init__(self, case, expected)

        # Set some default values, which can be overridden by subclasses
        self.set_test_suffix('test')
        self.set_description_first_phase('')
        self.set_description_second_phase('')

        # Initialize results
        # TODO(wjs, 2016-07-27) Currently these results of the individual pieces
        # aren't used anywhere, but I'm storing them because I think it would be
        # helpful to use them in the test reporting
        self._status_run1 = "NOT RUN"
        self._status_run2 = "NOT RUN"
        self._status_compare = "NOT RUN"


    # ========================================================================
    # Methods that can be called by subclasses in their init methods to
    # customize various aspects of the test.
    #
    # These should be called AFTER calling SystemTestsCompareTwo.__init__
    # ========================================================================

    def set_test_suffix(self, test_suffix):
        """
        Sets the test suffix appended to netcdf files that are output by the
        second run of the test.

        Typically, this will be called in the __init__ method of subclasses to
        set the test suffix to something reasonable.

        Arguments:

        test_suffix: string
        """

        self._test_suffix = test_suffix

    def set_description_first_phase(self, description):
        """
        Sets the description written to the log file when starting the first
        run.

        Typically, this will be called in the __init__ method of subclasses to
        set the description to something reasonable.

        Arguments:

        description: string
        """

        self._description_first_phase = description

    def set_description_second_phase(self, description):
        """
        Sets the description written to the log file when starting the second
        run.

        Typically, this will be called in the __init__ method of subclasses to
        set the description to something reasonable.

        Arguments:

        description: string
        """

        self._description_second_phase = description

    # ========================================================================
    # Methods that must be implemented by specific tests that inherit from this
    # base class
    # ========================================================================

    def _setup_first_phase(self):
        """
        Sets up the run for the first phase of the two-phase test.

        All tests inheriting from this base class MUST implement this method.
        """
        raise NotImplementedError

    def _setup_second_phase(self):
        """
        Sets up the run for the second phase of the two-phase test.

        All tests inheriting from this base class MUST implement this method.
        """
        raise NotImplementedError


    # ========================================================================
    # Main public methods
    # ========================================================================

    def run(self):
        """
        Runs both phases of the two-phase test and compares their results
        """

        self._setup_first_phase()
        self._case.flush()
        logger.info('Doing first run: ' + self._description_first_phase)
        success = self._run()
        if success:
            self._status_run1 = "PASS"
        else:
            self._status_run1 = "FAIL"
            return False

        self._setup_second_phase()
        self._case.flush()
        logger.info('Doing second run: ' + self._description_second_phase)
        success = self._run(self._test_suffix)
        if success:
            self._status_run2 = "PASS"
        else:
            self._status_run2 = "FAIL"
            return False

        success = self._component_compare_test("base", self._test_suffix)
        if success:
            self._status_compare = "PASS"
        else:
            self._status_compare = "FAIL"
            return False

        return success



