"""
Implementation of the CIME ERP test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
import CIME.utils
from system_tests_common import SystemTestsCommon

class ERP(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def build(self, sharedlib_only=False, model_only=False):
        """
        Build two cases.  Case one uses defaults, case2 uses half the number of threads
        and tasks.   This test will fail for components (pop) that do not reproduce exactly
        with different numbers of mpi tasks.
        """
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        machpes1 = os.path.join("LockedFiles","env_mach_pes.ERP1.xml")
        if ( os.path.isfile(machpes1) ):
            shutil.copy(machpes1,"env_mach_pes.xml")
        else:
            logging.warn("Copying env_mach_pes.xml to %s"%(machpes1))
            shutil.copy("env_mach_pes.xml", machpes1)

        for bld in range(1,3):
            logging.warn("Starting bld %s"%bld)
            if(bld == 2):
                # halve the number of tasks and threads
                for comp in ['ATM','CPL','OCN','WAV','GLC','ICE','ROF','LND']:
                    ntasks      = int(self._case.get_value("NTASKS_%s"%comp))
                    nthreads  = int(self._case.get_value("NTHRDS_%s"%comp))
                    rootpe      = int(self._case.get_value("ROOTPE_%s"%comp))
                    if ( nthreads > 1 ):
                        self._case.set_value("BUILD_THREADED","TRUE")
                    if ( ntasks > 1 ):
                        self._case.set_value("NTASKS_%s"%comp, "%s"%int(ntasks/2))
                        self._case.set_value("ROOTPE_%s"%comp, "%s"%int(rootpe/2))
                    if ( nthreads > 1 ):
                        self._case.set_value("NTHRDS_%s"%comp, "%s"%int(nthreads/2))
            self._case.set_value("SMP_BUILD","0")
            self._case.flush()
            run_cmd("case.setup -clean -testmode")
            run_cmd("case.setup")
            run_cmd('case.clean_build')
            SystemTestsCommon.build(self, sharedlib_only=sharedlib_only, model_only=model_only)
            if (not sharedlib_only):
                shutil.move("%s/%s.exe"%(exeroot,cime_model),
                            "%s/%s.exe.ERP%s"%(exeroot,cime_model,bld))
            shutil.copy("env_build.xml",os.path.join("LockedFiles","env_build_ERP%s.xml"%bld))

        #
        # Because mira/cetus interprets its run script differently than
        # other systems we need to copy the original env_mach_pes.xml
        # back
        #
        shutil.copy(machpes1,"env_mach_pes.xml")
        shutil.copy("env_mach_pes.xml",
                    os.path.join("LockedFiles","env_mach_pes.xml"))

    def run(self):
        SystemTestsCommon.run(self)

    def report(self):
        SystemTestsCommon.report(self)
