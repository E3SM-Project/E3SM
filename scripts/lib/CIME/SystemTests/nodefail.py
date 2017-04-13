"""
CIME restart upon failed node test.
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.ers import ERS
from CIME.utils import get_model

logger = logging.getLogger(__name__)

class NODEFAIL(ERS):

    def __init__(self, case):
        """
        initialize an object interface to the ERS system test
        """
        ERS.__init__(self, case)

        self._fail_sentinel = os.path.join(case.get_value("RUNDIR"), "FAIL_SENTINEL")
        self._fail_str      = case.get_value("NODE_FAIL_REGEX")

    def _restart_fake_phase(self):
        # Swap out model.exe for one that emits node failures
        rundir = self._case.get_value("RUNDIR")
        exeroot = self._case.get_value("EXEROOT")
        fake_exe = \
"""#!/bin/bash

fail_sentinel=%s
cpl_log=%s/cpl.log.$LID
model_log=%s/%s.log.$LID
touch $cpl_log
touch $fail_sentinel
declare -i num_fails=$(cat $fail_sentinel | wc -l)
declare -i times_to_fail=${NODEFAIL_NUM_FAILS:-3}

if ((num_fails < times_to_fail)); then
  echo FAKE FAIL >> $cpl_log
  echo FAIL >> $fail_sentinel
  echo '%s' >> $model_log
  sleep 1
  exit -1
else
  echo Insta pass
  echo SUCCESSFUL TERMINATION > $cpl_log
fi
""" % (self._fail_sentinel, rundir, rundir, get_model(), self._fail_str)

        fake_exe_file = os.path.join(exeroot, "fake.sh")
        with open(fake_exe_file, "w") as fd:
            fd.write(fake_exe)

        os.chmod(fake_exe_file, 0755)

        prev_run_exe = self._case.get_value("run_exe")
        env_mach_specific = self._case.get_env("mach_specific")
        env_mach_specific.set_value("run_exe", fake_exe_file)
        self._case.flush(flushall=True)

        self.run_indv(suffix=None)

        env_mach_specific = self._case.get_env("mach_specific")
        env_mach_specific.set_value("run_exe", prev_run_exe)
        self._case.flush(flushall=True)

    def run_phase(self):
        self._ers_first_phase()
        self._restart_fake_phase()
        self._ers_second_phase()
