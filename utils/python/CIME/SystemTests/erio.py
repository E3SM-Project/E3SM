"""
ERIO tests restart with different PIO methods

This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
import shutil

logger = logging.getLogger(__name__)

class ERIO(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        SystemTestsCommon.__init__(self, case, expected=["TEST"])

        self._pio_types = self._case.get_env("run").get_valid_values("PIO_TYPENAME")

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        Build one executable per PIO method
        """
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")

        for pio_type in self._pio_types:
            if pio_type == "default":
                continue
            else:
                self._case.set_value("PIO_TYPENAME", pio_type)
                self._case.flush()

                self.clean_build()
                self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)
                if sharedlib_only:
                    continue

                logger.info("doing a build with pio type %s" % pio_type)

                shutil.move("%s/%s.exe" % (exeroot, cime_model),
                            "%s/%s.exe.%s" % (exeroot, cime_model, pio_type))

    def _full_run(self, pio_type):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")

        base_exe_path = "%s/%s.exe" % (exeroot, cime_model)
        pio_exe_path  = "%s/%s.exe.%s" % (exeroot, cime_model, pio_type)

        if os.path.exists(base_exe_path):
            os.remove(base_exe_path)

        shutil.copy(pio_exe_path, base_exe_path)

        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        expect(stop_n > 0, "Bad STOP_N: %d" % stop_n)

        # Move to config_tests.xml once that's ready
        rest_n = stop_n/2 + 1
        self._case.set_value("REST_N", rest_n)
        self._case.set_value("REST_OPTION", stop_option)
        self._case.set_value("HIST_N", stop_n)
        self._case.set_value("HIST_OPTION", stop_option)
        self._case.set_value("CONTINUE_RUN", False)
        self._case.flush()

        expect(stop_n > 2, "ERROR: stop_n value %d too short"%stop_n)
        logger.info("doing an %s %s initial test with restart file at %s %s with pio type %s"
                    %(str(stop_n), stop_option, str(rest_n), stop_option, pio_type))
        self.run_indv(suffix=pio_type)

    def _restart_run(self, pio_type, other_pio_type):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = self._case.get_value("MODEL")

        base_exe_path = "%s/%s.exe" % (exeroot, cime_model)
        pio_exe_path  = "%s/%s.exe.%s" % (exeroot, cime_model, pio_type)

        if os.path.exists(base_exe_path):
            os.remove(base_exe_path)

        shutil.copy(pio_exe_path, base_exe_path)

        stop_n      = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")

        rest_n = stop_n/2 + 1
        stop_new = stop_n - rest_n
        expect(stop_new > 0, "ERROR: stop_n value %d too short %d %d"%(stop_new,stop_n,rest_n))

        self._case.set_value("STOP_N", stop_new)
        self._case.set_value("CONTINUE_RUN", True)
        self._case.set_value("REST_OPTION","never")
        self._case.flush()
        logger.info("doing an %s %s restart test with %s against %s"
                    %(str(stop_n), stop_option, pio_type, other_pio_type))

        suffix = "%s.%s" % (other_pio_type, pio_type)
        self.run_indv(suffix=suffix)

        # Compare restart file
        self._component_compare_test(other_pio_type, suffix)

    def run_phase(self):

        for idx, pio_type1 in enumerate(self._pio_types):
            if pio_type1 != "default":
                self._full_run(pio_type1)
                for pio_type2 in self._pio_types[idx+1:]:
                    if pio_type2 != "default":
                        self._restart_run(pio_type2, pio_type1)
