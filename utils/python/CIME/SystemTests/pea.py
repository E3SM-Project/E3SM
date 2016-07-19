"""
Implementation of the CIME PEA test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case_setup import case_setup
import CIME.utils
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)

class PEA(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, case)

    def build(self, sharedlib_only=False, model_only=False):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        # first set all component mpi tasks to 1
        for comp in ['ATM','CPL','OCN','WAV','GLC','ICE','ROF','LND']:
            self._case.set_value("NTASKS_%s"%comp, 1)

        build1 = os.path.join("LockedFiles","env_build.PEA1.xml")
        if ( os.path.isfile(build1) ):
            shutil.copy(build1,"env_build.xml")

        mpilib = self._case.get_value("MPILIB")
        for mpilib in [mpilib, "mpi-serial"]:
            logging.warn("Starting bld for %s"%mpilib)
            self._case.set_value("MPILIB",mpilib)
            self._case.flush()
            case_setup(self._case, reset=True)
            self.clean_build()
            SystemTestsCommon.build(self, sharedlib_only=sharedlib_only, model_only=model_only)
            if (not sharedlib_only):
                shutil.move("%s/%s.exe"%(exeroot,cime_model),
                            "%s/%s.exe.PEA_%s"%(exeroot,cime_model,mpilib))
            shutil.copy("env_build.xml",os.path.join("LockedFiles",
                                                     "env_build_PEA_%s.xml"%mpilib))

    def _pea_first_phase(self):

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        exefile  = "%s/%s.exe"%(exeroot,cime_model)
        exefile1 = "%s/%s.exe.PEA1"%(exeroot,cime_model)
        if (os.path.isfile(exefile)):
            os.remove(exefile)
        shutil.copy(exefile1, exefile)

        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")
        self._case.flush()

        stop_n = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        logger.info("doing an %d %s initial test with 1pe and mpi, no restarts written"
                    % (stop_n, stop_option))

        return SystemTestsCommon._run(self)

    def _pea_second_phase(self):

        expect(os.path.isfile("env_mach_pes.xml.2"),
               "ERROR: env_mach_pes.xml.2 does not exist, run case.build" )

        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        exefile  = "%s/%s.exe"%(exeroot,cime_model)
        exefile2 = "%s/%s.exe.PEA2"%(exeroot,cime_model)
        if (os.path.isfile(exefile)):
            os.remove(exefile)
        shutil.copy(exefile2, exefile)

        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")
        self._case.flush()

        stop_n = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        logger.info("doing an %d %s initial test with 1pe and serial mpi, no restarts written"
                    % (stop_n, stop_option))

        success = SystemTestsCommon._run(self, "mpiserial")

        # Compare restart file
        if success:
            return self._component_compare_test("base", "mpiserial")
        else:
            return False

    def run(self):
        success = self._pea_first_phase()

        if success:
            return self._pea_second_phase()
        else:
            return False

    def report(self):
        SystemTestsCommon.report(self)
