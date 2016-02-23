"""
Interface to the env_test.xml file.  This class inherits from EnvBase
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
import CIME.utils
from system_tests_common import SystemTestsCommon


class CME(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def build(self):
        self._case.set_value('USE_ESMF_LIB','TRUE')
        exeroot = self._case.get_value("EXEROOT")
        cime_model = CIME.utils.get_model()
        for CPL in ['MCT','ESMF']:
            self._case.set_value('COMP_INTERFACE',CPL)
            self._case.flush()
            run_cmd('case.clean_build')
            SystemTestsCommon.build(self)
            shutil.move("%s/%s.exe"%(exeroot,cime_model),
                        "%s/%s.exe.%s"%(exeroot,cime_model,CPL))
            shutil.copy("env_build.xml",os.path.join("LockedFiles","env_build_%s.xml"%CPL))

    def run(self):
        SystemTestsCommon.run(self)

    def report(self):
        SystemTestsCommon.report(self)
