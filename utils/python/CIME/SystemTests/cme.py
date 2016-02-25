"""
Interface to the CME system test.  This class inherits from SystemTestsCommon
"""
import shutil
from CIME.XML.standard_module_setup import *
from CIME.case import Case
import CIME.utils
from system_tests_common import SystemTestsCommon


class CME(SystemTestsCommon):
    def __init__(self, caseroot, case):
        """
        initialize an object interface to the CME test
        """
        SystemTestsCommon.__init__(self, caseroot, case)

    def build(self):
        """
        Build two exectuables for the CME test, one with ESMF interfaces
        the other with MCT interfaces and compare results - they should be exact
        """
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
