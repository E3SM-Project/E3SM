#!/usr/bin/env python2

"""
case.submit - Submit a cesm workflow to the queueing system or run it if there is no queueing system.   A cesm workflow may include multiple jobs.
"""

from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.utils import expect
from CIME.env_module        import EnvModule
from CIME.preview_namelists        import preview_namelists
from CIME.check_lockedfiles        import check_lockedfiles
from CIME.XML.batch                 import Batch
from CIME.batch_utils import BatchUtils

logger = logging.getLogger(__name__)

def submit(caseroot, job=None, resubmit=None, no_batch=False):
    case = Case(caseroot)
    if job is None:
        if case.get_value("TEST"):
            job = "case.test"
        else:
            job = "case.run"
    if resubmit:
        resub = case.get_value("RESUBMIT")
        case.set_value("RESUBMIT",resub-1)
        if case.get_value("RESUBMIT_SETS_CONTINUE_RUN"):
            case.set_value("CONTINUE_RUN",True)
    else:
        check_case(case, caseroot)
        check_DA_settings(case)

    cimeroot = case.get_value("CIMEROOT")
    os.environ["CIMEROOT"] = cimeroot
    #Load Modules
    env_module = EnvModule(case.get_value("MACH"), case.get_value("COMPILER"),
                           case.get_value("CIMEROOT"),caseroot, case.get_value("MPILIB"),
                           case.get_value("DEBUG"))
    env_module.load_env_for_case()
    batchobj = BatchUtils(job, case)
    case.set_value("RUN_WITH_SUBMIT",True)
    case.flush()
    batchobj.submit_jobs(no_batch=no_batch)

def check_case(case,caseroot):
    check_lockedfiles(caseroot)
    preview_namelists(dryrun=False, casedir=caseroot)
    expect(case.get_value("BUILD_COMPLETE"), "Build complete is "
           "not True please rebuild the model by calling case.build")
    logger.info("Check case OK")

def check_DA_settings(case):
    if case.get_value("DATA_ASSIMILATION"):
        script = case.get_value("DATA_ASSIMILATION_SCRIPT")
        cycles = case.get_value("DATA_ASSIMILATION_CYCLES")
        logger.info("Data Assimilation enabled using script %s with %d cycles"%(script,cycles))

