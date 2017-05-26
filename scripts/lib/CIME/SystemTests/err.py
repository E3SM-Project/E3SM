"""
CIME ERR test  This class inherits from ERS
ERR tests short term archiving and restart capabilities
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.ers import ERS

import shutil, glob

logger = logging.getLogger(__name__)

class ERR(ERS):

    def __init__(self, case):
        """
        initialize an object interface to the ERR system test
        """
        ERS.__init__(self, case)

    def run_phase(self):
        first_phase = self._case.get_value("RESUBMIT") == 1

        if first_phase:
            self._case.set_value("DOUT_S", True)
            self._case.flush()
            self._ers_first_phase()
        else:
            dout_s_root = self._case.get_value("DOUT_S_ROOT")
            rundir = self._case.get_value("RUNDIR")
            case = self._case.get_value("CASE")
            # First remove restart, history and rpointer files from the run directory
            for item in glob.iglob(os.path.join(rundir, "{}.*".format(case))):
                if not item.endswith("base"):
                    os.remove(item)
            for item in glob.iglob(os.path.join(rundir, "rpointer.*")):
                os.remove(item)
            # Then replace them from the restart directory
            for item in glob.iglob(os.path.join(dout_s_root,"rest","*","*")):
                shutil.copy(item, rundir)

            self._ers_second_phase()
