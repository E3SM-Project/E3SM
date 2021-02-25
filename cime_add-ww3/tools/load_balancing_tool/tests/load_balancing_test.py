#!/usr/bin/env python
"""
Test for pulp, glpk

IceLndAtmOcn:
* Test writing to json file, then running from it
* test running solve from timing dir files
* test for pulp
* test can save json file if pulp is not available
* test running submit using X
* test running solve from X
* test writing pes file
* test extended algorithm
"""
try:
    from Tools.standard_script_setup import *
except ImportError, e:
    print 'Error importing Tools.standard_script_setup'
    print 'May need to add cime/scripts to PYTHONPATH\n'
    raise ImportError(e)
try:
    import optimize_model
except ImportError, e:
    print 'Error importing optimize_model'
    print 'May need to add cime/tools/load_balancing_tool to PYTHONPATH\n'
    raise ImportError(e)



from CIME.utils import run_cmd_no_fail, get_full_test_name
from CIME.XML.machines import Machines
from CIME.XML import pes
import unittest, json, tempfile, sys, re, copy

SCRIPT_DIR  = CIME.utils.get_scripts_root()
MACHINE = Machines()
CODE_DIR = os.path.join(SCRIPT_DIR, "..", "tools", "load_balancing_tool")
TEST_DIR = os.path.join(SCRIPT_DIR, "..", "tools", "load_balancing_tool",
                        "tests")
X_OPTIONS = """
STOP_N=1
"""
JSON_DICT = {
  "description" : "Optimize using data available from original load balancing tool. The original tool solved the problem using a different model, so we do not expect exact replication: (Original solution: NTASKS_ATM: 1006 NTASKS_ICE:  889 NTASKS_LND:  117 NTASKS_OCN:   18 TOTAL_COST:  28.749 s/mday)",
  "layout" : "IceLndAtmOcn",
  "totaltasks" : 1024,
  "ATM" : {
      "ntasks" : [32,64,128,256,512],
      "blocksize" : 8,
      "nthrds" : [1],
      "cost" : [427.471, 223.332, 119.580, 66.182, 37.769]
  },
  "OCN" : {
      "ntasks" : [32,64,128,256,512],
      "blocksize" : 8,
      "nthrds" : [1],
      "cost" : [ 15.745, 7.782, 4.383, 3.181, 2.651]
  },
  "LND" : {
      "ntasks" : [32,64,128,256,512],
      "blocksize" : 8,
      "nthrds" : [1],
      "cost" : [  4.356, 2.191, 1.191, 0.705, 0.560]
  },
  "ICE" : {
      "ntasks" : [32,64,160,320,640],
      "blocksize" : 8,
      "nthrds" : [1],
      "cost" : [8.018, 4.921, 2.368, 1.557, 1.429]
  }
}

PES_XML = """
<config_pes>
  <grid name="any">
    <mach name="any">
      <pes compset="any" pesize="0">
        <comment>none</comment>
        <ntasks>
          <ntasks_atm>2</ntasks_atm>
          <ntasks_lnd>2</ntasks_lnd>
          <ntasks_rof>2</ntasks_rof>
          <ntasks_ice>2</ntasks_ice>
          <ntasks_ocn>2</ntasks_ocn>
          <ntasks_glc>2</ntasks_glc>
          <ntasks_wav>2</ntasks_wav>
          <ntasks_cpl>2</ntasks_cpl>
        </ntasks>
        <nthrds>
          <nthrds_atm>1</nthrds_atm>
          <nthrds_lnd>1</nthrds_lnd>
          <nthrds_rof>1</nthrds_rof>
          <nthrds_ice>1</nthrds_ice>
          <nthrds_ocn>1</nthrds_ocn>
          <nthrds_glc>1</nthrds_glc>
          <nthrds_wav>1</nthrds_wav>
          <nthrds_cpl>1</nthrds_cpl>
        </nthrds>
        <rootpe>
          <rootpe_atm>0</rootpe_atm>
          <rootpe_lnd>0</rootpe_lnd>
          <rootpe_rof>0</rootpe_rof>
          <rootpe_ice>0</rootpe_ice>
          <rootpe_ocn>0</rootpe_ocn>
          <rootpe_glc>0</rootpe_glc>
          <rootpe_wav>0</rootpe_wav>
          <rootpe_cpl>0</rootpe_cpl>
        </rootpe>
      </pes>
    </mach>
  </grid>

  <grid name="any">
    <mach name="any">
      <pes compset="any" pesize="1">
        <comment>none</comment>
        <ntasks>
          <ntasks_atm>4</ntasks_atm>
          <ntasks_lnd>4</ntasks_lnd>
          <ntasks_rof>4</ntasks_rof>
          <ntasks_ice>4</ntasks_ice>
          <ntasks_ocn>4</ntasks_ocn>
          <ntasks_glc>4</ntasks_glc>
          <ntasks_wav>4</ntasks_wav>
          <ntasks_cpl>4</ntasks_cpl>
        </ntasks>
        <nthrds>
          <nthrds_atm>1</nthrds_atm>
          <nthrds_lnd>1</nthrds_lnd>
          <nthrds_rof>1</nthrds_rof>
          <nthrds_ice>1</nthrds_ice>
          <nthrds_ocn>1</nthrds_ocn>
          <nthrds_glc>1</nthrds_glc>
          <nthrds_wav>1</nthrds_wav>
          <nthrds_cpl>1</nthrds_cpl>
        </nthrds>
        <rootpe>
          <rootpe_atm>0</rootpe_atm>
          <rootpe_lnd>0</rootpe_lnd>
          <rootpe_rof>0</rootpe_rof>
          <rootpe_ice>0</rootpe_ice>
          <rootpe_ocn>0</rootpe_ocn>
          <rootpe_glc>0</rootpe_glc>
          <rootpe_wav>0</rootpe_wav>
          <rootpe_cpl>0</rootpe_cpl>
        </rootpe>
      </pes>
    </mach>
  </grid>
</config_pes>
"""

###############################################################################
def _main_func(description):
###############################################################################

    unittest.main(verbosity=2, catchbreak=True)

###############################################################################

class LoadBalanceTests(unittest.TestCase):
    def _check_solution(self, output, var, val):
        """
        Utility function, checks output of milp solve to make sure solution
        value is expected
        """
        pattern = var + ' = (\d+)'
        m = re.search(pattern, output)
        if not m:
            self.fail("pattern '%s' not found in output" % (pattern))
        check = int(m.groups()[0])
        self.assertTrue(check == val, "%s = %d, expected %d" % (var, check, val))


    def test_pulp(self):
        try:
            import pulp
        except ImportError, e:
            self.fail("ERROR: pulp not found. Install or set PYTHONPATH")
        x = pulp.LpVariable('x')
        p = pulp.LpProblem('p', pulp.LpMinimize)
        p.solve()
        self.assertTrue(p.status == 1, "ERROR: simple pulp solve failed")


    def test_read_and_write_json(self):
        "Solve from json file, writing to new json file, solve from new file"
        with tempfile.NamedTemporaryFile('w+') as jsonfile1, tempfile.NamedTemporaryFile('w+') as jsonfile2:
            json.dump(JSON_DICT, jsonfile1)
            jsonfile1.flush()
            cmd = "./load_balancing_solve.py --json-input %s --json-output %s" % (jsonfile1.name, jsonfile2.name)
            output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)
            self._check_solution(output, "NTASKS_ATM", 992)
            cmd = "./load_balancing_solve.py --json-input %s" % jsonfile2.name
            output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)
            self._check_solution(output, "NTASKS_ATM", 992)


    def test_solve_from_timing_dir(self):
        cmd = "./load_balancing_solve.py --timing-dir %s --total-tasks 64 --blocksize 2 --layout IceLndAtmOcn" % os.path.join(TEST_DIR, "timing")
        output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)
        self._check_solution(output, "NTASKS_ATM", 62)

    def test_write_pes(self):
        with tempfile.NamedTemporaryFile('w+') as jsonfile1, tempfile.NamedTemporaryFile('w+') as pes_file:
            json.dump(JSON_DICT, jsonfile1)
            jsonfile1.flush()
            cmd = "./load_balancing_solve.py --json-input %s --pe-output %s" % (jsonfile1.name, pes_file.name)
            output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)

            self.assertTrue(os.access(pes_file.name, os.R_OK), "pesfile %s not written" % pes_file.name)
            pesobj = CIME.XML.pes.Pes(pes_file.name)

        pes_ntasks, pes_nthrds, pes_rootpe, _, _, _ = \
               pesobj.find_pes_layout('any', 'any', 'any', '')
        self.assertTrue(pes_ntasks['NTASKS_ATM']==992)
        

    def test_set_blocksize_atm(self):
        cmd = "./load_balancing_solve.py --timing-dir %s --total-tasks 64 --blocksize 2 --blocksize-atm 4 --layout IceLndAtmOcn" % os.path.join(TEST_DIR, "timing")
        output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)
        self._check_solution(output, "NTASKS_ATM", 60)
        self._check_solution(output, "NBLOCKS_ATM", 15)
        self._check_solution(output, "NTASKS_OCN", 4)
        self._check_solution(output, "NBLOCKS_OCN", 2)

    def test_graph_models(self):
        try:
            import matplotlib
        except ImportError, e:
            self.skipTest("matplotlib not found")

        with tempfile.NamedTemporaryFile('w+') as jsonfile:
            json.dump(JSON_DICT, jsonfile)
            jsonfile.flush()
            cmd = "./load_balancing_solve.py --json-input %s --graph-models" % (jsonfile.name)
            output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)
            self._check_solution(output, "NTASKS_ATM", 992)

    def test_xcase_submit(self):
        test_root = MACHINE.get_value("CIME_OUTPUT_ROOT")
        machine = MACHINE.get_machine_name()
        compiler = MACHINE.get_default_compiler()

        test_name = get_full_test_name("PFS_I0",grid="f19_g16", compset="X",
                                             machine=machine, compiler=compiler)
        expected_dir = os.path.join(test_root,
                                    "{}.test_lbt".format(test_name),
                                    "timing")
        if not os.path.isdir(expected_dir):
            with tempfile.NamedTemporaryFile('w+') as tfile, tempfile.NamedTemporaryFile('w+') as xfile:
                tfile.write(PES_XML)
                tfile.flush()
                xfile.write(X_OPTIONS)
                xfile.flush()
                cmd = "./load_balancing_submit.py --pesfile {} --res f19_g16 --compset X --test-id test_lbt  --extra-options-file {} --test-root {}".format(tfile.name, xfile.name, test_root)
                if MACHINE.has_batch_system():
                    sys.stdout.write("Jobs will be submitted to queue. Rerun "
                                     "load_balancing_test.py after jobs have "
                    "finished.")
                else:
                    cmd += " --force-purge"
                output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)

                self.assertTrue(output.find("Timing jobs submitted") >= 0,
                                "Expected 'Timing jobs submitted' in output")

        if os.path.isdir(expected_dir):

            cmd = "./load_balancing_solve.py --total-tasks 32 --blocksize 1 --test-id test_lbt --print-models --test-root {} --layout IceLndAtmOcn".format(test_root)
            output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)
            self.assertTrue(output.find("***ATM***") > 0,
                            "--print-models failed to print ATM data")
            self._check_solution(output, "NTASKS_ATM", 31)

    def test_use_atm_lnd(self):
        "Solve layout atm_lnd from json file"
        with tempfile.NamedTemporaryFile('w+') as jsonfile1:
            atmlnd_dict = copy.deepcopy(JSON_DICT)
            # Fake data for ROF, CPL
            atmlnd_dict['ROF'] = {"ntasks" : [32,64,128,256],
                                  "blocksize" : 8,
                                  "nthrds" : [1],
                                  "cost" : [8.0, 4.0, 2.0, 1.0]}
            atmlnd_dict['CPL'] = {"ntasks" : [32,64,128,256],
                                  "blocksize" : 8,
                                  "nthrds" : [1],
                                  "cost" : [8.0, 4.0, 2.0, 1.0]}
            json.dump(atmlnd_dict, jsonfile1)
            jsonfile1.flush()
            cmd = "./load_balancing_solve.py --json-input %s --print-models --layout tests.atm_lnd.AtmLnd" % (jsonfile1.name)
            output = run_cmd_no_fail(cmd, from_dir=CODE_DIR)
            self._check_solution(output, "Natm", 976)
            self._check_solution(output, "NBatm", 976/8)

if __name__ == '__main__':
    _main_func(__doc__)
