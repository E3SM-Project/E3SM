"""
Implementation of the CIME PEM test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
import CIME.utils
from system_tests_common import SystemTestsCommon


class PEM(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def build(self, sharedlib_only=False, model_only=False):
        """
        build two cases, the first is default the second has halve the tasks per component
        """
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        machpes1 = os.path.join("LockedFiles","env_mach_pes.PEM1.xml")
        if ( os.path.isfile(machpes1) ):
            shutil.copy(machpes1,"env_mach_pes.xml")

        for bld in range(1,3):
            logging.warn("Starting bld %s"%bld)
            self._case.flush()
            run_cmd("case.setup -clean -testmode")
            run_cmd("case.setup")
            run_cmd('case.clean_build')
            SystemTestsCommon.build(self, sharedlib_only=sharedlib_only, model_only=model_only)
            machpes = os.path.join("LockedFiles","env_mach_pes.PEM%s.xml"%bld)

            if (not sharedlib_only):
                shutil.move("%s/%s.exe"%(exeroot,cime_model),
                            "%s/%s.exe.PEM%s"%(exeroot,cime_model,bld))

            shutil.copy("env_build.xml",os.path.join("LockedFiles","env_build_PEM%s.xml"%bld))
            shutil.copy("env_mach_pes.xml", machpes)

            if bld == 1:
                ntasks = int(self._case.get_value("NTASKS_%s"%comp))
                rootpe = int(self._case.get_value("ROOTPE_%s"%comp))
                if ntasks > 1:
                    self._case.set_value("NTASKS_%s"%comp, "%s"%int(ntasks/2))
                    self._case.set_value("ROOTPE_%s"%comp, "%s"%int(rootpe/2))


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
