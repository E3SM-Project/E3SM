from CIME.XML.standard_module_setup import *
from CIME.utils                     import expect, does_file_have_string, append_status
from CIME.XML.lt_archive            import LTArchive

import time

logger = logging.getLogger(__name__)

###############################################################################
def case_lt_archive(case):
###############################################################################
    caseroot = case.get_value("CASEROOT")

    # max number of threads needed by scripts
    os.environ["maxthrds"] = 1

    # document start
    append_status("lt_archive starting",caseroot=caseroot,sfile="CaseStatus")

    # determine status of run and short term archiving
    runComplete = does_file_have_string(os.path.join(caseroot, "CaseStatus"),
                                        "run SUCCESSFUL")
    staComplete = does_file_have_string(os.path.join(caseroot, "stArchiveStatus"),
                                        "st_archive_complete")

    # set up envrionment vars and call the lt_archive.sh script
    if runComplete and staComplete:
        os.environ["DOUT_S_ROOT"] = case.get_value("DOUT_S_ROOT")
        os.environ["DOUT_L_MSROOT"] = case.get_value("DOUT_L_MSROOT")
        os.environ["DOUT_L_HPSS_ACCNT"] = case.get_value("DOUT_L_HPSS_ACCNT")

        lid = time.strftime("%y%m%d-%H%M%S")
        lt_archive = LTArchive(case.get_value("MACH"))
        lt_archive_args = lt_archive.get_lt_archive_args()
        cmd = os.path.join(caseroot, "Tools/lt_archive.sh") \
            + lt_archive_args + "ltArchiveStatus." + lid + " 2>&1"
        run_cmd_no_fail(cmd, from_dir=caseroot)
    else:
        expect(False,
               "lt_archive: run or st_archive is not yet complete or was not successful."
               "Unable to perform long term archive...")

    # document completion
    append_status("lt_archive completed" ,caseroot=caseroot, sfile="CaseStatus")

    return True
