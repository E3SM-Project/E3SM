"""
Base class for CIME system tests that involve doing two runs and comparing their
output.

In the __init__ method for your test, you MUST call
    SystemTestsCompareTwo.__init__
See the documentation of that method for details.

Classes that inherit from this are REQUIRED to implement the following methods:

(1) _run_common_setup
    This method will be called before both the first and second phase of the
    two-phase test, before either _run_one_setup or
    _run_two_setup. This should contain settings needed for both phases of
    the run, such as setting CONTINUE_RUN to False. In principle, these settings
    could just be made once, but for robustness and communication purposes, this
    is executed before both run phases.

(2) _run_one_setup
    This method will be called to set up the run for the first phase of the
    two-phase test

(3) _run_two_setup
    This method will be called to set up the run for the second phase of the
    two-phase test

Classes that inherit from this MAY implement the following methods, if they have
any work to be done in these phases:

(1) _pre_build
    This method will be called immediately before the build. This can contain
    work like copying user_nl files to some saved location
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class SystemTestsCompareTwo(SystemTestsCommon):

    def __init__(self,
                 case,
                 run_two_suffix = 'test',
                 run_one_description = '',
                 run_two_description = ''):
        """
        Initialize a SystemTestsCompareTwo object. Individual test cases that
        inherit from SystemTestsCompareTwo MUST call this __init__ method.

        Args:
            case: case object passsed to __init__ method of individual test
            run_two_suffix (str, optional): Suffix appended to files output by
                the second run. Defaults to 'test'. This can be anything other
                than 'base' (which is the suffix used for the first run).
            run_one_description (str, optional): Description printed to log file
                when starting the first run. Defaults to ''.
            run_two_description (str, optional): Description printed to log file
                when starting the second run. Defaults to ''.
        """
        SystemTestsCommon.__init__(self, case)

        # NOTE(wjs, 2016-08-03) It is currently CRITICAL for run_one_suffix to
        # be 'base', because this is assumed for baseline comparison and
        # generation. Once that assumption is relaxed, then run_one_suffix can
        # be set in the call to the constructor just like run_two_suffix
        # currently is.
        self._run_one_suffix = 'base'
        self._run_two_suffix = run_two_suffix.rstrip()
        expect(self._run_two_suffix != self._run_one_suffix,
               "ERROR: Must have different suffixes for run one and run two")

        self._run_one_description = run_one_description
        self._run_two_description = run_two_description

        # Initialize results
        # TODO(wjs, 2016-07-27) Currently these results of the individual pieces
        # aren't used anywhere, but I'm storing them because I think it would be
        # helpful to use them in the test reporting
        self._status_run1 = "NOT RUN"
        self._status_run2 = "NOT RUN"
        self._status_compare = "NOT RUN"


    # ========================================================================
    # Methods that MUST be implemented by specific tests that inherit from this
    # base class
    # ========================================================================

    def _run_common_setup(self):
        """
        This method will be called before both the first and second phase of the
        two-phase test, before either _run_one_setup or
        _run_two_setup. This should contain settings needed for both phases of
        the run, such as setting CONTINUE_RUN to False. In principle, these settings
        could just be made once, but for robustness and communication purposes, this
        is executed before both run phases.
        """
        raise NotImplementedError


    def _run_one_setup(self):
        """
        Sets up the run for the first phase of the two-phase test.

        All tests inheriting from this base class MUST implement this method.
        """
        raise NotImplementedError

    def _run_two_setup(self):
        """
        Sets up the run for the second phase of the two-phase test.

        All tests inheriting from this base class MUST implement this method.
        """
        raise NotImplementedError

    # ========================================================================
    # Methods that MAY be implemented by specific tests that inherit from this
    # base class, if they have any work that needs to be done during these
    # phases
    # ========================================================================

    def _pre_build(self):
        """
        This method will be called immediately before the build. This can
        contain work like copying user_nl files to some saved location.

        Note that this may be called multiple times: once for the sharedlib
        build, once for the model build, and potentially more times if the case
        is rebuilt. Thus, anything that is done in this method must be written
        so that it is robust to being called multiple times.
        """
        pass

    # ========================================================================
    # Main public methods
    # ========================================================================

    def build(self, sharedlib_only=False, model_only=False):
        self._pre_build()
        SystemTestsCommon.build(self, sharedlib_only=sharedlib_only, model_only=model_only)

    def run(self):
        """
        Runs both phases of the two-phase test and compares their results
        """

        # First run
        self._run_common_setup()
        self._run_one_setup()
        logger.info('Doing first run: ' + self._run_one_description)
        success = self._run(self._run_one_suffix)
        if success:
            self._status_run1 = "PASS"
        else:
            self._status_run1 = "FAIL"
            return False

        # Second run
        self._run_common_setup()
        self._run_two_setup()
        logger.info('Doing second run: ' + self._run_two_description)
        success = self._run(self._run_two_suffix)
        if success:
            self._status_run2 = "PASS"
        else:
            self._status_run2 = "FAIL"
            return False

        # Compare results
        success = self._component_compare_test(self._run_one_suffix, self._run_two_suffix)
        if success:
            self._status_compare = "PASS"
        else:
            self._status_compare = "FAIL"
            return False

        return success

    def get_run_one_status(self):
        """
        Returns a string specifying the status of run 1
        """
        return self._status_run1

    def get_run_two_status(self):
        """
        Returns a string specifying the status of run 2
        """
        return self._status_run2

    def get_compare_status(self):
        """
        Returns a string specifying the status of the comparison between run 1
        and run 2
        """
        return self._status_compare

