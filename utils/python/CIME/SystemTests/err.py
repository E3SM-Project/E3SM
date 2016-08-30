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
            logger.info("staging files from archive %s" % dout_s_root)
            for item in glob.glob(os.path.join(dout_s_root, "*", "hist", "*base")):
                shutil.copy(item, rundir)

            self._ers_second_phase()
