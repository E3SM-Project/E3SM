from CIME.XML.standard_module_setup import *
from CIME.case_submit               import submit
from CIME.case                      import Case
from CIME.utils                     import get_model, run_cmd, append_status

logger = logging.getLogger(__name__)

###############################################################################
def case_st_archive(case):
###############################################################################
    caseroot = case.get_value("CASEROOT")
    logger.info("st_archive starting")
    # do short-term archiving
    append_status("st_archiving starting",
                 caseroot=caseroot, sfile="CaseStatus")

    cmd = os.path.join(caseroot, "Tools/st_archive") + " >> stArchiveStatus 2>&1"
    rc, out, err = run_cmd(cmd, ok_to_fail=True)
    if rc != 0:
        append_status("st_archive failed: %s \nerr = %s"%(out,err),sfile="CaseStatus")
        return False

    append_status("st_archiving completed",
                 caseroot=caseroot, sfile="CaseStatus")
    logger.info("st_archive completed")

    # resubmit case if appropriate
    resubmit = case.get_value("RESUBMIT")
    if resubmit > 0:
        append_status("resubmitting from st_archive",
                      caseroot=caseroot, sfile="CaseStatus")
        logger.info("resubmitting from st_archive, resubmit=%d"%resubmit)
        submit(case, resubmit=True)

    return True
