from CIME.XML.standard_module_setup import *
from CIME.case_submit               import submit
from CIME.case                      import Case
from CIME.utils                     import get_model, run_cmd, append_status

logger = logging.getLogger(__name__)

###############################################################################
def case_st_archive(case):
###############################################################################
    caseroot = case.get_value("CASEROOT")

    # do short-term archiving
    append_status("st_archiving starting",
                 caseroot=caseroot, sfile="CaseStatus")

    cmd = os.path.join(caseroot, "Tools/st_archive") + " >> stArchiveStatus 2>&1"
    run_cmd(cmd)

    append_status("st_archiving completed",
                 caseroot=caseroot, sfile="CaseStatus")

    # resubmit case if appropriate
    resubmit = case.get_value("RESUBMIT")
    if resubmit > 0:
        append_status("resubmitting from st_archive",
                      caseroot=caseroot, sfile="CaseStatus")
        submit(case, resubmit=True)

    return True
