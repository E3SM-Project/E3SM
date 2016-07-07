"""
Implementation of the CIME PEM test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.case_setup import case_setup
import CIME.utils
from system_tests_common import SystemTestsCommon

class PEM(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize a test object
        """
        SystemTestsCommon.__init__(self, case)

    def build(self, sharedlib_only=False, model_only=False):
        """
        build two cases, the first is default the second has halve the tasks per component
        """
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()

        machpes1 = os.path.join("LockedFiles","env_mach_pes.PEM1.xml")
        if ( os.path.isfile(machpes1) ):
            shutil.copy(machpes1,"env_mach_pes.xml")
        else:
            logging.warn("Copying env_mach_pes.xml to %s"%(machpes1))
            shutil.copy("env_mach_pes.xml", machpes1)

        for bld in range(1,3):
            logging.warn("Starting bld %s"%bld)
            self._case.flush()
            case_setup(self._case, reset=True)
            self.clean_build()
            SystemTestsCommon.build(self, sharedlib_only=sharedlib_only, model_only=model_only)
            machpes = os.path.join("LockedFiles","env_mach_pes.PEM%s.xml"%bld)

            if (not sharedlib_only):
                shutil.move("%s/%s.exe"%(exeroot,cime_model),
                            "%s/%s.exe.PEM%s"%(exeroot,cime_model,bld))

            shutil.copy("env_build.xml",os.path.join("LockedFiles","env_build_PEM%s.xml"%bld))
            shutil.copy("env_mach_pes.xml", machpes)

            if bld == 1:
                for comp in ["CPL", "ATM", "LND", "ICE", "OCN", "ROF", "GLC", "WAV"]:
                    ntasks = self._case.get_value("NTASKS_%s"%comp)
                    rootpe = self._case.get_value("ROOTPE_%s"%comp)
                    if ntasks > 1:
                        self._case.set_value("NTASKS_%s"%comp, ntasks/2)
                        self._case.set_value("ROOTPE_%s"%comp, rootpe/2)

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
