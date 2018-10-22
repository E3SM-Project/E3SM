"""
Base class for CIME system tests that involve doing two runs and comparing their
output.

In the __init__ method for your test, you MUST call
    SystemTestsCompareTwo.__init__
See the documentation of that method for details.

Classes that inherit from this are REQUIRED to implement the following methods:

(1) _case_one_setup
    This method will be called to set up case 1, the "base" case

(2) _case_two_setup
    This method will be called to set up case 2, the "test" case

In addition, they MAY require the following methods:

(1) _common_setup
    This method will be called to set up both cases. It should contain any setup
    that's needed in both cases. This is called before _case_one_setup or
    _case_two_setup.

(2) _case_one_custom_prerun_action(self):
    Use this to do arbitrary actions immediately before running case one

(3) _case_two_custom_prerun_action(self):
    Use this to do arbitrary actions immediately before running case two

(4) _case_one_custom_postrun_action(self):
    Use this to do arbitrary actions immediately after running case one

(5) _case_two_custom_postrun_action(self):
    Use this to do arbitrary actions immediately after running case two
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case import Case
from CIME.utils import get_model
from CIME.test_status import *

import shutil, os, glob

logger = logging.getLogger(__name__)

class SystemTestsCompareTwo(SystemTestsCommon):

    def __init__(self,
                 case,
                 separate_builds = False,
                 run_two_suffix = 'test',
                 run_one_description = '',
                 run_two_description = '',
                 multisubmit = False):
        """
        Initialize a SystemTestsCompareTwo object. Individual test cases that
        inherit from SystemTestsCompareTwo MUST call this __init__ method.

        Args:
            case: case object passsed to __init__ method of individual
                test. This is the main case associated with the test.
            separate_builds (bool): Whether separate builds are needed for the
                two cases. If False, case2 uses the case1 executable.
            run_two_suffix (str, optional): Suffix appended to the case name for
                the second run. Defaults to 'test'. This can be anything other
                than 'base'.
            run_one_description (str, optional): Description printed to log file
                when starting the first run. Defaults to ''.
            run_two_description (str, optional): Description printed to log file
                when starting the second run. Defaults to ''.
            multisubmit (bool): Do first and second runs as different submissions.
                Designed for tests with RESUBMIT=1
        """
        SystemTestsCommon.__init__(self, case)

        self._separate_builds = separate_builds

        # run_one_suffix is just used as the suffix for the netcdf files
        # produced by the first case; we may eventually remove this, but for now
        # it is needed by the various component_*.sh scripts. run_two_suffix is
        # also used as the suffix for netcdf files, but more importantly is used
        # to create the case name for the clone case.
        #
        # NOTE(wjs, 2016-08-03) It is currently CRITICAL for run_one_suffix to
        # be 'base', because this is assumed for baseline comparison and
        # generation. Once that assumption is relaxed, then run_one_suffix can
        # be set in the call to the constructor just like run_two_suffix
        # currently is. Or, if these tools are rewritten to work without any
        # suffix, then run_one_suffix can be removed entirely.
        self._run_one_suffix = 'base'
        self._run_two_suffix = run_two_suffix.rstrip()
        expect(self._run_two_suffix != self._run_one_suffix,
               "ERROR: Must have different suffixes for run one and run two")

        self._run_one_description = run_one_description
        self._run_two_description = run_two_description

        # Save case for first run so we can return to it if we switch self._case
        # to point to self._case2
        self._case1 = self._case
        self._caseroot1 = self._get_caseroot()

        self._caseroot2 = self._get_caseroot2()
        # Initialize self._case2; it will get set to its true value in
        # _setup_cases_if_not_yet_done
        self._case2 = None

        self._setup_cases_if_not_yet_done()

        self._multisubmit = multisubmit and self._case1.get_value("BATCH_SYSTEM") != "none"

    # ========================================================================
    # Methods that MUST be implemented by specific tests that inherit from this
    # base class
    # ========================================================================

    def _case_one_setup(self):
        """
        This method will be called to set up case 1, the "base" case.

        This should be written to refer to self._case: this object will point to
        case1 at the point that this is called.
        """
        raise NotImplementedError

    def _case_two_setup(self):
        """
        This method will be called to set up case 2, the "test" case

        This should be written to refer to self._case: this object will point to
        case2 at the point that this is called.
        """
        raise NotImplementedError

    # ========================================================================
    # Methods that MAY be implemented by specific tests that inherit from this
    # base class, if they have any work to do in these methods
    # ========================================================================

    def _common_setup(self):
        """
        This method will be called to set up both cases. It should contain any setup
        that's needed in both cases. This is called before _case_one_setup or
        _case_two_setup.

        This should be written to refer to self._case: It will be called once with
        self._case pointing to case1, and once with self._case pointing to case2.
        """
        pass

    def _case_one_custom_prerun_action(self):
        """
        Use to do arbitrary actions immediately before running case one
        """
        pass

    def _case_two_custom_prerun_action(self):
        """
        Use to do arbitrary actions immediately before running case two
        """
        pass

    def _case_one_custom_postrun_action(self):
        """
        Use to do arbitrary actions immediately after running case one
        """
        pass

    def _case_two_custom_postrun_action(self):
        """
        Use to do arbitrary actions immediately after running case two
        """
        pass

    # ========================================================================
    # Main public methods
    # ========================================================================

    def build_phase(self, sharedlib_only=False, model_only=False):
        if self._separate_builds:
            self._activate_case1()
            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
            self._activate_case2()
            # Although we're doing separate builds, it still makes sense
            # to share the sharedlibroot area with case1 so we can reuse
            # pieces of the build from there.
            if get_model() != "e3sm":
                # We need to turn off this change for E3SM because it breaks
                # the MPAS build system
                self._case2.set_value("SHAREDLIBROOT",
                                      self._case1.get_value("SHAREDLIBROOT"))

            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
        else:
            self._activate_case1()
            self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
            # pio_typename may be changed during the build if the default is not a
            # valid value for this build, update case2 to reflect this change
            for comp in self._case1.get_values("COMP_CLASSES"):
                comp_pio_typename = "{}_PIO_TYPENAME".format(comp)
                self._case2.set_value(comp_pio_typename, self._case1.get_value(comp_pio_typename))

            # The following is needed when _case_two_setup has a case_setup call
            # despite sharing the build (e.g., to change NTHRDS)
            self._case2.set_value("BUILD_COMPLETE",True)
            self._case2.flush()

    def run_phase(self, success_change=False):  # pylint: disable=arguments-differ
        """
        Runs both phases of the two-phase test and compares their results
        If success_change is True, success requires some files to be different
        """
        first_phase = self._case1.get_value("RESUBMIT") == 1 # Only relevant for multi-submit tests
        run_type = self._case1.get_value("RUN_TYPE")

        # First run
        if not self._multisubmit or first_phase:
            logger.info('Doing first run: ' + self._run_one_description)

            # Add a PENDing compare phase so that we'll notice if the second part of compare two
            # doesn't run.
            with self._test_status:
                self._test_status.set_status("{}_{}_{}".format(COMPARE_PHASE, self._run_one_suffix, self._run_two_suffix), TEST_PEND_STATUS)

            self._activate_case1()
            self._case_one_custom_prerun_action()
            self.run_indv(suffix = self._run_one_suffix)
            self._case_one_custom_postrun_action()

        # Second run
        logger.info("_multisubmit {} first phase {}".format(self._multisubmit, first_phase))
        if not self._multisubmit or not first_phase:
            # Subtle issue: case1 is already in a writeable state since it tends to be opened
            # with a with statement in all the API entrances in CIME. case2 was created via clone,
            # not a with statement, so it's not in a writeable state, so we need to use a with
            # statement here to put it in a writeable state.
            with self._case2:
                logger.info('Doing second run: ' + self._run_two_description)
                self._activate_case2()
                # This assures that case two namelists are populated
                self._skip_pnl = False
                # we need to make sure run2 is properly staged.
                if run_type != "startup":
                    self._case2.check_case()

                self._case_two_custom_prerun_action()
                self.run_indv(suffix = self._run_two_suffix)
                self._case_two_custom_postrun_action()
            # Compare results
            # Case1 is the "main" case, and we need to do the comparisons from there
            self._activate_case1()
            self._link_to_case2_output()
            self._component_compare_test(self._run_one_suffix, self._run_two_suffix, success_change=success_change)

    def copy_case1_restarts_to_case2(self):
        """
        Makes a copy (or symlink) of restart files and related files
        (necessary history files, rpointer files) from case1 to case2.

        This is not done automatically, but can be called by individual
        tests where case2 does a continue_run using case1's restart
        files.
        """
        rundir2 = self._case2.get_value("RUNDIR")
        self._case1.archive_last_restarts(archive_restdir = rundir2,
                                          rundir=self._case1.get_value("RUNDIR"),
                                          link_to_restart_files = True)

    # ========================================================================
    # Private methods
    # ========================================================================

    def _get_caseroot2(self):
        """
        Determines and returns caseroot for case2

        Assumes that self._case1 is already set to point to the case1 object
        """
        casename2 = self._case1.get_value("CASE")
        caseroot1 = self._case1.get_value("CASEROOT")

        # Nest the case directory for case2 inside the case directory for case1
        caseroot2 = os.path.join(caseroot1, "case2", casename2)

        return caseroot2

    def _get_output_root2(self):
        """
        Determines and returns cime_output_root for case2

        Assumes that self._case1 is already set to point to the case1 object
        """
        # Since case2 has the same name as case1, its CIME_OUTPUT_ROOT
        # must also be different, so that anything put in
        # $CIME_OUTPUT_ROOT/$CASE/ is not accidentally shared between
        # case1 and case2. (Currently nothing is placed here, but this
        # helps prevent future problems.)
        output_root2 = os.path.join(self._case1.get_value("CIME_OUTPUT_ROOT"),
                                    self._case1.get_value("CASE"), "case2_output_root")
        return output_root2

    def _get_case2_exeroot(self):
        """
        Gets exeroot for case2.

        Returns None if we should use the default value of exeroot.
        """
        if self._separate_builds:
            # case2's EXEROOT needs to be somewhere that (1) is unique
            # to this case (considering that case1 and case2 have the
            # same case name), and (2) does not have too long of a path
            # name (because too-long paths can make some compilers
            # fail).
            case1_exeroot = self._case1.get_value("EXEROOT")
            case2_exeroot = os.path.join(case1_exeroot, "case2bld")
        else:
            # Use default exeroot
            case2_exeroot = None
        return case2_exeroot

    def _get_case2_rundir(self):
        """
        Gets rundir for case2.
        """
        # case2's RUNDIR needs to be somewhere that is unique to this
        # case (considering that case1 and case2 have the same case
        # name). Note that the location below is symmetrical to the
        # location of case2's EXEROOT set in _get_case2_exeroot.
        case1_rundir = self._case1.get_value("RUNDIR")
        case2_rundir = os.path.join(case1_rundir, "case2run")
        return case2_rundir

    def _setup_cases_if_not_yet_done(self):
        """
        Determines if case2 already exists on disk. If it does, this method
        creates the self._case2 object pointing to the case directory. If it
        doesn't exist, then this method creates case2 as a clone of case1, and
        sets the self._case2 object appropriately.

        This also does the setup for both case1 and case2.

        Assumes that the following variables are already set in self:
            _caseroot1
            _caseroot2
            _case1

        Sets self._case2
        """

        # Use the existence of the case2 directory to signal whether we have
        # done the necessary test setup for this test: When we initially create
        # the case2 directory, we set up both test cases; then, if we find that
        # the case2 directory already exists, we assume that the setup has
        # already been done. (In some cases it could be problematic to redo the
        # test setup when it's not needed - e.g., by appending things to user_nl
        # files multiple times. This is why we want to make sure to just do the
        # test setup once.)
        if os.path.exists(self._caseroot2):
            self._case2 = self._case_from_existing_caseroot(self._caseroot2)
        else:
            try:
                self._case2 = self._case1.create_clone(
                    self._caseroot2,
                    keepexe = not self._separate_builds,
                    cime_output_root = self._get_output_root2(),
                    exeroot = self._get_case2_exeroot(),
                    rundir = self._get_case2_rundir())
                self._write_info_to_case2_output_root()
                self._setup_cases()
            except:
                # If a problem occurred in setting up the test cases, it's
                # important to remove the case2 directory: If it's kept around,
                # that would signal that test setup was done successfully, and
                # thus doesn't need to be redone - which is not the case. Of
                # course, we'll likely be left in an inconsistent state in this
                # case, but if we didn't remove the case2 directory, the next
                # re-build of the test would think, "okay, setup is done, I can
                # move on to the build", which would be wrong.
                if os.path.isdir(self._caseroot2):
                    shutil.rmtree(self._caseroot2)
                self._activate_case1()
                logger.warning("WARNING: Test case setup failed. Case2 has been removed, "
                               "but the main case may be in an inconsistent state. "
                               "If you want to rerun this test, you should create "
                               "a new test rather than trying to rerun this one.")
                raise

    def _case_from_existing_caseroot(self, caseroot):
        """
        Returns a Case object from an existing caseroot directory

        Args:
            caseroot (str): path to existing caseroot
        """
        return Case(case_root=caseroot, read_only=False)

    def _activate_case1(self):
        """
        Make case 1 active for upcoming calls
        """
        os.chdir(self._caseroot1)
        self._set_active_case(self._case1)

    def _activate_case2(self):
        """
        Make case 2 active for upcoming calls
        """
        os.chdir(self._caseroot2)
        self._set_active_case(self._case2)

    def _write_info_to_case2_output_root(self):
        """
        Writes a file with some helpful information to case2's
        output_root.

        The motivation here is two-fold:

        (1) Currently, case2's output_root directory is empty. This
            could be confusing.

        (2) For users who don't know where to look, it could be hard to
            find case2's bld and run directories. It is somewhat easier
            to stumble upon case2's output_root, so we put a file there
            pointing them to the right place.
        """

        readme_path = os.path.join(self._get_output_root2(), "README")
        try:
            with open(readme_path, "w") as fd:
                fd.write("This directory is typically empty.\n\n")
                fd.write("case2's run dir is here: {}\n\n".format(
                    self._case2.get_value("RUNDIR")))
                fd.write("case2's bld dir is here: {}\n".format(
                    self._case2.get_value("EXEROOT")))
        except IOError:
            # It's not a big deal if we can't write the README file
            # (e.g., because the directory doesn't exist or isn't
            # writeable; note that the former may be the case in unit
            # tests). So just continue merrily on our way if there was a
            # problem.
            pass

    def _setup_cases(self):
        """
        Does all test-specific set up for the two test cases.
        """

        # Set up case 1
        self._activate_case1()
        self._common_setup()
        self._case_one_setup()
        # Flush the case so that, if errors occur later, then at least case 1 is
        # in a correct, post-setup state. This is important because the mere
        # existence of a case 2 directory signals that setup is done. So if the
        # build fails and the user rebuilds, setup won't be redone - so it's
        # important to ensure that the results of setup are flushed to disk.
        #
        # Note that case 1 will be in its post-setup state even if case 2 setup
        # fails. Putting the case1 flush after case 2 setup doesn't seem to help
        # with that (presumably some flush is called automatically), and anyway
        # wouldn't help with things like appending to user_nl files (which don't
        # rely on flush). So we just have to live with that possibility (but
        # note that we print a warning to the log file if that happens, in the
        # caller of this method).
        self._case.flush()
        # This assures that case one namelists are populated
        # and creates the case.test script
        self._case.case_setup(test_mode=False, reset=True)

        # Set up case 2
        self._activate_case2()
        self._common_setup()
        self._case_two_setup()
        # Flush the case so that, if errors occur later, then at least case2 is
        # in a correct, post-setup state
        self._case.flush()

        # Go back to case 1 to ensure that's where we are for any following code
        self._activate_case1()

    def _link_to_case2_output(self):
        """
        Looks for all files in rundir2 matching the pattern casename2*.nc.run2suffix

        For each file found, makes a link in rundir1 pointing to this file; the
        link is renamed so that the original occurrence of casename2 is replaced
        with casename1.

        For example:

        /glade/scratch/sacks/somecase/run/somecase.clm2.h0.nc.run2 ->
        /glade/scratch/sacks/somecase.run2/run/somecase.run2.clm2.h0.nc.run2

        If the destination link already exists and points to the correct
        location, it is maintained as is. However, an exception will be raised
        if the destination link is not exactly as it should be: we avoid
        overwriting some existing file or link.
        """

        casename1 = self._case1.get_value("CASE")
        casename2 = self._case2.get_value("CASE")
        rundir1 = self._case1.get_value("RUNDIR")
        rundir2 = self._case2.get_value("RUNDIR")
        run2suffix = self._run_two_suffix

        pattern = '{}*.nc.{}'.format(casename2, run2suffix)
        case2_files = glob.glob(os.path.join(rundir2, pattern))
        for one_file in case2_files:
            file_basename = os.path.basename(one_file)
            modified_basename = file_basename.replace(casename2, casename1, 1)
            one_link = os.path.join(rundir1, modified_basename)
            if (os.path.islink(one_link) and
                os.readlink(one_link) == one_file):
                # Link is already set up correctly: do nothing
                # (os.symlink raises an exception if you try to replace an
                # existing file)
                pass
            else:
                os.symlink(one_file, one_link)
