"""
Implementation of the CIME LII test.  This class inherits from SystemTestsCommon

This is a CLM specific test:
Verifies that namelist variable 'use_init_interp' works correctly
(1) do a run with use_init_interp false (suffix base)
(2) do a run with use_init_interp true (suffix init_interp_on)
"""

import shutil, glob
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class LII(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):

        # Make copies of the namelist files for each part of the test. Enclose the
        # copies in conditionals so that we only do this namelist setup the first time
        # the build script is invoked - otherwise, if the build is rerun, the namelist
        # files would build up repeated instances of the setting of force_init_intep.
        #
        # Note the use of shell wildcards to make sure we apply these mods to
        # multi-instance versions

        if not os.path.exists("user_nl_nointerp"):
            os.makedirs("user_nl_nointerp")
            for filename in glob.glob(r'user_nl_clm*'):
                shutil.copy(filename, os.path.join("user_nl_nointerp",filename))
                with open(os.path.join("user_nl_nointerp",filename), "a") as newfile:
                    newfile.write("use_init_interp = .false.")

        if not os.path.exists("user_nl_interp"):
            os.makedirs("user_nl_interp")
            for filename in glob.glob(r'user_nl_clm*'):
                shutil.copy(filename, os.path.join("user_nl_interp",filename))
                with open(os.path.join("user_nl_interp",filename), "a") as newfile:
                    newfile.write("use_init_interp = .true.")

        self.clean_build()
        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

    def run_phase(self):
        '''
        Do a run with init_interp false, a run with init_interp true and
        compare
        '''
        caseroot = self._case.get_value("CASEROOT")

        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")
        self._case.flush()
        for user_nl_dir in ("nointerp", "interp"):
            for filename in glob.glob(r'user_nl_%s/*'%user_nl_dir):
                shutil.copy(filename,
                            os.path.join(caseroot,os.path.basename(filename)))

            stop_n = self._case.get_value("STOP_N")
            stop_option = self._case.get_value("STOP_OPTION")
            logger.info("doing a %d %s initial test with init_interp set to %s, no restarts written"
                        % (stop_n, stop_option, user_nl_dir == "interp"))

            self.run_indv(suffix=user_nl_dir)

        self._component_compare_test("nointerp", "interp")
