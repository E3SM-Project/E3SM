from CIME.XML.standard_module_setup import *

from CIME.utils                     import expect, run_cmd
from CIME.case import Case
from CIME.env_module        import EnvModule
from CIME.preview_namelists        import preview_namelists
from CIME.check_lockedfiles        import check_lockedfiles
from CIME.XML.batch                 import Batch
from CIME.batch_utils import BatchUtils

logger = logging.getLogger(__name__)

def submit(caseroot, job="case.run", resubmit=None):
    case = Case(caseroot)
#    if resubmit is None:
#        check_case(case, caseroot)
        #Load Modules
    env_module = EnvModule(case.get_value("MACH"), case.get_value("COMPILER"),
                           case.get_value("CIMEROOT"),caseroot, case.get_value("MPILIB"),
                           case.get_value("DEBUG"))
    env_module.load_env_for_case()
    batchobj = BatchUtils(job, case)
    case.set_value("RUN_WITH_SUBMIT",True)
    case.flush()
    batchobj.submit_single_job(job)



def check_case(case,caseroot):
    check_lockedfiles(caseroot)
    preview_namelists(dryrun=False, casedir=caseroot)
    expect(case.get_value("BUILD_COMPLETE"), "Build complete is "
           "not True please rebuild the model by calling case.build")
    logger.info("Check case OK")

