#!/usr/bin/env python2

"""
case.submit - Submit a cesm workflow to the queueing system or run it
if there is no queueing system.  A cesm workflow may include multiple
jobs.
"""

from CIME.XML.standard_module_setup import *
from CIME.case import Case
from CIME.utils import expect
from CIME.preview_namelists        import preview_namelists
from CIME.check_lockedfiles        import check_lockedfiles
from CIME.XML.batch                 import Batch
from CIME.XML.env_mach_specific  import EnvMachSpecific

logger = logging.getLogger(__name__)

def submit(case, job=None, resubmit=False, no_batch=False, prereq_jobid=None):
    caseroot = case.get_value("CASEROOT")
    if job is None:
        if case.get_value("TEST"):
            job = "case.test"
        else:
            job = "case.run"

    logger.info("Submitting job '%s', resubmit=%s" % (job, resubmit))

    if resubmit:
        resub = case.get_value("RESUBMIT")
        case.set_value("RESUBMIT",resub-1)
        if case.get_value("RESUBMIT_SETS_CONTINUE_RUN"):
            case.set_value("CONTINUE_RUN", True)
    else:
        check_case(case, caseroot)
        check_DA_settings(case)

    cimeroot = case.get_value("CIMEROOT")
    os.environ["CIMEROOT"] = cimeroot

    # if case.submit is called with the no_batch flag then we assume that this
    # flag will stay in effect for the duration of the RESUBMITs
    env_batch = case._get_env("batch")
    if not resubmit:
        case.set_value("IS_FIRST_RUN", True)
        if no_batch:
            batch_system = "none"
        else:
            batch_system = env_batch.get_batch_system_type()
        case.set_value("BATCH_SYSTEM", batch_system)
    else:
        if env_batch.get_batch_system_type() == "none":
            no_batch = True

        # This is a resubmission, do not reinitialize test values
        case.set_value("IS_FIRST_RUN", False)

    #Load Modules
    env_module = case._get_env("mach_specific")


    env_module.load_env_for_case(compiler=case.get_value("COMPILER"),
                                 debug=case.get_value("DEBUG"),
                                 mpilib=case.get_value("MPILIB"))

    case.set_value("RUN_WITH_SUBMIT",True)
    case.flush()
    case.submit_jobs(no_batch=no_batch)

def check_case(case, caseroot):
    check_lockedfiles(caseroot)
    preview_namelists(case, dryrun=False)
    expect(case.get_value("BUILD_COMPLETE"), "Build complete is "
           "not True please rebuild the model by calling case.build")
    logger.info("Check case OK")

def check_DA_settings(case):
    if case.get_value("DATA_ASSIMILATION"):
        script = case.get_value("DATA_ASSIMILATION_SCRIPT")
        cycles = case.get_value("DATA_ASSIMILATION_CYCLES")
        logger.info("Data Assimilation enabled using script %s with %d cycles"%(script,cycles))

