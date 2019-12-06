"""
CIME ERI test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from CIME.utils import safe_copy
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from stat import S_ISDIR, ST_CTIME, ST_MODE
import shutil, glob, os

logger = logging.getLogger(__name__)

def _get_rest_date(archive_root):
    restdir = os.path.join(archive_root,"rest")
    # get all entries in the directory w/ stats
    entries = (os.path.join(restdir, fn) for fn in os.listdir(restdir))
    entries = ((os.stat(path), path) for path in entries)
    entries = sorted((stat[ST_CTIME], path)
                     for stat, path in entries if S_ISDIR(stat[ST_MODE]))
    last_dir = os.path.basename(entries[-1][1])
    ref_sec = last_dir[-5:]
    ref_date = last_dir[:10]
    return ref_date, ref_sec


def _helper(dout_sr, refdate, refsec, rundir):
    rest_path = os.path.join(dout_sr, "rest", "{}-{}".format(refdate, refsec))

    for item in glob.glob("{}/*{}*".format(rest_path, refdate)):
        dst = os.path.join(rundir, os.path.basename(item))
        if os.path.exists(dst):
            os.remove(dst)
        os.symlink(item, dst)

    for item in glob.glob("{}/*rpointer*".format(rest_path)):
        safe_copy(item, rundir)

class ERI(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the ERI system test
        """
        SystemTestsCommon.__init__(self, case)
        self._testname = "ERI"

    def run_phase(self):
        caseroot = self._case.get_value("CASEROOT")
        clone1_path = "{}.ref1".format(caseroot)
        clone2_path = "{}.ref2".format(caseroot)
        #self._case.set_value("CHECK_TIMING", False)

        #
        # clone the main case to create ref1 and ref2 cases
        #
        for clone_path in [clone1_path, clone2_path]:
            if os.path.exists(clone_path):
                shutil.rmtree(clone_path)

        clone1, clone2 = [self._case.create_clone(clone_path, keepexe=True) for clone_path in [clone1_path, clone2_path]]
        orig_case = self._case
        orig_casevar = orig_case.get_value("CASE")
        #
        # determine run lengths needed below
        #
        stop_n = self._case.get_value("STOP_N")
        stop_option = self._case.get_value("STOP_OPTION")
        run_startdate = self._case.get_value("RUN_STARTDATE")
        start_tod = self._case.get_value("START_TOD")
        if start_tod == 0:
            start_tod="00000"

        stop_n1 = int(stop_n / 6)
        rest_n1 = stop_n1
        start_1 = run_startdate

        stop_n2 = stop_n - stop_n1
        rest_n2 = int(stop_n2 / 2 + 1)
        hist_n  = stop_n2

        start_1_year, start_1_month, start_1_day = [int(item) for item in start_1.split("-")]
        start_2_year = start_1_year + 2
        start_2 = "{:04d}-{:02d}-{:02d}".format(start_2_year, start_1_month, start_1_day)

        stop_n3 = stop_n2 - rest_n2
        rest_n3 = int(stop_n3 / 2 + 1)

        stop_n4 = stop_n3 - rest_n3

        expect(stop_n4 >= 1 and stop_n1 >= 1, "Run length too short")

        #
        # (1) Test run:
        # do an initial ref1 case run
        # cloned the case and running there
        # (NOTE: short term archiving is on)
        #

        os.chdir(clone1_path)
        self._set_active_case(clone1)

        logger.info("ref1 startup: doing a {} {} startup run from {} and {} seconds".format(stop_n1, stop_option, start_1, start_tod))
        logger.info("  writing restarts at {} {}".format(rest_n1, stop_option))
        logger.info("  short term archiving is on ")

        with clone1:
            clone1.set_value("CONTINUE_RUN", False)
            clone1.set_value("RUN_STARTDATE", start_1)
            clone1.set_value("STOP_N", stop_n1)
            clone1.set_value("REST_OPTION", stop_option)
            clone1.set_value("REST_N", rest_n1)
            clone1.set_value("HIST_OPTION", "never")

        dout_sr1 = clone1.get_value("DOUT_S_ROOT")

        # force cam namelist to write out initial file at end of run
        if os.path.exists("user_nl_cam"):
            if "inithist" not in open("user_nl_cam", "r").read():
                with open("user_nl_cam", "a") as fd:
                    fd.write("inithist = 'ENDOFRUN'\n")

        with clone1:
            clone1.case_setup(test_mode=True, reset=True)
            # if the initial case is hybrid this will put the reference data in the correct location
            clone1.check_all_input_data()

            self._skip_pnl = False
            self.run_indv(st_archive=True, suffix=None)

        #
        # (2) Test run:
        # do a hybrid ref2 case run
        # cloned the main case and running with ref1 restarts
        # (NOTE: short term archiving is on)
        #

        os.chdir(clone2_path)
        self._set_active_case(clone2)

        # Set startdate to start2, set ref date based on ref1 restart
        refdate_2, refsec_2 = _get_rest_date(dout_sr1)

        logger.info("ref2 hybrid: doing a {} {} startup hybrid run".format(stop_n2, stop_option))
        logger.info("  starting from {} and using ref1 {} and {} seconds".format(start_2, refdate_2, refsec_2))
        logger.info("  writing restarts at {} {}".format(rest_n2, stop_option))
        logger.info("  short term archiving is on ")

        # setup ref2 case
        with clone2:
            clone2.set_value("RUN_TYPE",      "hybrid")
            clone2.set_value("RUN_STARTDATE", start_2)
            clone2.set_value("RUN_REFCASE",   "{}.ref1".format(orig_casevar))
            clone2.set_value("RUN_REFDATE",   refdate_2)
            clone2.set_value("RUN_REFTOD",    refsec_2)
            clone2.set_value("GET_REFCASE",   False)
            clone2.set_value("CONTINUE_RUN",  False)
            clone2.set_value("STOP_N",        stop_n2)
            clone2.set_value("REST_OPTION",   stop_option)
            clone2.set_value("REST_N",        rest_n2)
            clone2.set_value("HIST_OPTION",   stop_option)
            clone2.set_value("HIST_N",        hist_n)

        rundir2 = clone2.get_value("RUNDIR")
        dout_sr2 = clone2.get_value("DOUT_S_ROOT")

        _helper(dout_sr1, refdate_2, refsec_2, rundir2)

        # run ref2 case (all component history files will go to short term archiving)
        with clone2:
            clone2.case_setup(test_mode=True, reset=True)

            self._skip_pnl = False
            self.run_indv(suffix="hybrid", st_archive=True)

        #
        # (3a) Test run:
        # do a branch run from ref2 restart (short term archiving is off)
        #

        os.chdir(caseroot)
        self._set_active_case(orig_case)
        refdate_3, refsec_3 = _get_rest_date(dout_sr2)


        logger.info("branch: doing a {} {} branch".format(stop_n3, stop_option))
        logger.info("  starting from ref2 {} and {} seconds restarts".format(refdate_3, refsec_3))
        logger.info("  writing restarts at {} {}".format(rest_n3, stop_option))
        logger.info("  short term archiving is off")

        self._case.set_value("RUN_TYPE"      , "branch")
        self._case.set_value("RUN_REFCASE"   , "{}.ref2".format(self._case.get_value("CASE")))
        self._case.set_value("RUN_REFDATE"   , refdate_3)
        self._case.set_value("RUN_REFTOD"    , refsec_3)
        self._case.set_value("GET_REFCASE"   , False)
        self._case.set_value("CONTINUE_RUN"  , False)
        self._case.set_value("STOP_N"        , stop_n3)
        self._case.set_value("REST_OPTION"   , stop_option)
        self._case.set_value("REST_N"        , rest_n3)
        self._case.set_value("HIST_OPTION"   , stop_option)
        self._case.set_value("HIST_N"        , stop_n2)
        self._case.set_value("DOUT_S"        , False)
        self._case.flush()

        rundir = self._case.get_value("RUNDIR")
        if not os.path.exists(rundir):
            os.makedirs(rundir)

        _helper(dout_sr2, refdate_3, refsec_3, rundir)

        # link the hybrid history files from ref2 to the run dir for comparison
        for item in glob.iglob("%s/*.hybrid"%rundir2):
            newfile = "{}".format(item.replace(".ref2", ""))
            newfile = os.path.basename(newfile)
            dst = os.path.join(rundir, newfile)
            if os.path.exists(dst):
                os.remove(dst)
            os.symlink(item, dst)

        self._skip_pnl = False
        # run branch case (short term archiving is off)
        self.run_indv()

        #
        # (3b) Test run:
        # do a restart continue from (3a) (short term archiving off)
        #

        logger.info("branch restart: doing a {} {} continue restart test".format(stop_n4, stop_option))

        self._case.set_value("CONTINUE_RUN",  True)
        self._case.set_value("STOP_N",        stop_n4)
        self._case.set_value("REST_OPTION",   "never")
        self._case.set_value("DOUT_S",        False)
        self._case.set_value("HIST_OPTION",   stop_option)
        self._case.set_value("HIST_N",        hist_n)
        self._case.flush()

        # do the restart run (short term archiving is off)
        self.run_indv(suffix="rest")

        self._component_compare_test("base", "hybrid")
        self._component_compare_test("base", "rest")
