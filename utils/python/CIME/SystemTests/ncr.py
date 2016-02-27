"""
Implementation of the CIME NCR test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
import CIME.utils
from system_tests_common import SystemTestsCommon


class NCR(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def build(self):
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        machpes1 = os.path.join("LockedFiles","env_mach_pes.NCR1.xml")
        if ( os.path.isfile(machpes1) ):
            shutil.copy(machpes1,"env_mach_pes.xml")

        for bld in range(1,3):
            """
            Build two exectuables for this test, the first is a default build
            the second halves the number of tasks and runs two instances
            for each component
            """
            logging.warn("Starting bld %s"%bld)
            machpes = os.path.join("LockedFiles","env_mach_pes.NCR%s.xml"%bld)
            ntasks_sum = 0
            for comp in ['ATM','OCN','WAV','GLC','ICE','ROF','LND']:
                self._case.set_value("NINST_%s"%comp,str(bld))
                ntasks      = int(self._case.get_value("NTASKS_%s"%comp))
                if(bld == 1):
                    self._case.set_value("ROOTPE_%s"%comp,"0")
                    if ( ntasks > 1 ):
                        self._case.set_value("NTASKS_%s"%comp, "%s"%int(ntasks/2))
                else:
                    self._case.set_value("ROOTPE_%s"%comp, "%s"% ntasks_sum)
                    ntasks_sum += ntasks*2
                    self._case.set_value("NTASKS_%s"%comp, "%s"%int(ntasks*2))
            self._case.flush()


            run_cmd("case.setup -clean -testmode")
            run_cmd("case.setup")
            run_cmd('case.clean_build')
            SystemTestsCommon.build(self)
            shutil.move("%s/%s.exe"%(exeroot,cime_model),
                        "%s/%s.exe.NCR%s"%(exeroot,cime_model,bld))
            shutil.copy("env_build.xml",os.path.join("LockedFiles","env_build.NCR%s.xml"%bld))
            shutil.copy("env_mach_pes.xml", machpes)
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
