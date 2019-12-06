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
        driver = self._case.get_value("COMP_INTERFACE")
        if driver == "nuopc":
            logname = "med"
        else:
            logname = "cpl"
        fake_exe = \
"""#!/bin/bash

fail_sentinel={0}
cpl_log={1}/{4}.log.$LID
model_log={1}/{2}.log.$LID
touch $cpl_log
touch $fail_sentinel
declare -i num_fails=$(cat $fail_sentinel | wc -l)
declare -i times_to_fail=${{NODEFAIL_NUM_FAILS:-3}}

if ((num_fails < times_to_fail)); then
  echo FAKE FAIL >> $cpl_log
  echo FAIL >> $fail_sentinel
  echo '{3}' >> $model_log
  sleep 1
  exit -1
else
  echo Insta pass
  echo SUCCESSFUL TERMINATION > $cpl_log
fi
""".format(self._fail_sentinel, rundir, get_model(), self._fail_str, logname)

        fake_exe_file = os.path.join(exeroot, "fake.sh")
        with open(fake_exe_file, "w") as fd:
            fd.write(fake_exe)

        os.chmod(fake_exe_file, 0o755)

        prev_run_exe = self._case.get_value("run_exe")
        env_mach_specific = self._case.get_env("mach_specific")
        env_mach_specific.set_value("run_exe", fake_exe_file)
        self._case.flush(flushall=True)

        # This flag is needed by mpt to run a script under mpiexec
        mpilib = self._case.get_value("MPILIB")
        if mpilib == "mpt":
            os.environ["MPI_SHEPHERD"] = "true"

        self.run_indv(suffix=None)

        if mpilib == "mpt":
            del os.environ["MPI_SHEPHERD"]

        env_mach_specific = self._case.get_env("mach_specific")
        env_mach_specific.set_value("run_exe", prev_run_exe)
        self._case.flush(flushall=True)

    def run_phase(self):
        self._ers_first_phase()
        self._restart_fake_phase()
        self._ers_second_phase()
