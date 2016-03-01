"""
Implementation of the CIME PEA test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
import CIME.utils
from system_tests_common import SystemTestsCommon

class PEA(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def build(self, sharedlib_only=False, model_only=False):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        for comp in ['ATM','CPL','OCN','WAV','GLC','ICE','ROF','LND']:
            self._case.set_value("NTASKS_%s"%comp,"1")

        build1 = os.path.join("LockedFiles","env_build.PEA1.xml")
        if ( os.path.isfile(build1) ):
            shutil.copy(build1,"env_build.xml")

        mpilib = self._case.get_value("MPILIB")
        for mpilib in [mpilib, "mpi-serial"]:
            logging.warn("Starting bld for %s"%mpilib)
            self._case.set_value("MPILIB",mpilib)
            self._case.flush()
            run_cmd("case.setup -clean ")
            run_cmd("case.setup")
            run_cmd('case.clean_build')
            SystemTestsCommon.build(self, sharedlib_only=sharedlib_only, model_only=model_only)
            if (not sharedlib_only):
                shutil.move("%s/%s.exe"%(exeroot,cime_model),
                            "%s/%s.exe.PEA_%s"%(exeroot,cime_model,mpilib))
            shutil.copy("env_build.xml",os.path.join("LockedFiles",
                                                     "env_build_PEA_%s.xml"%mpilib))

    def run(self):
        SystemTestsCommon.run(self)

    def report(self):
        SystemTestsCommon.report(self)
