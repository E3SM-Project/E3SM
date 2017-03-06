#!/usr/bin/env python

"""
case.submit - Submit a cesm workflow to the queueing system or run it
if there is no queueing system.  A cesm workflow may include multiple
jobs.
"""
import socket
from CIME.XML.standard_module_setup import *
from CIME.utils                     import expect, append_status
from CIME.preview_namelists         import create_namelists
from CIME.check_lockedfiles         import check_lockedfiles
from CIME.check_input_data          import check_all_input_data
from CIME.case_cmpgen_namelists     import case_cmpgen_namelists

logger = logging.getLogger(__name__)

def submit(case, job=None, resubmit=False, no_batch=False):
    caseroot = case.get_value("CASEROOT")

    if job is None:
        if case.get_value("TEST"):
            job = "case.test"
        else:
            job = "case.run"

    if resubmit:
        resub = case.get_value("RESUBMIT")
        logger.info("Submitting job '%s', resubmit=%d" % (job, resub))
        case.set_value("RESUBMIT",resub-1)
        if case.get_value("RESUBMIT_SETS_CONTINUE_RUN"):
            case.set_value("CONTINUE_RUN", True)
    else:
        if job in ("case.test","case.run"):
            check_case(case, caseroot)
            check_DA_settings(case)
            if case.get_value("MACH") == "mira":
                with open(".original_host","w") as fd:
                    fd.write( socket.gethostname())

    # if case.submit is called with the no_batch flag then we assume that this
    # flag will stay in effect for the duration of the RESUBMITs
    env_batch = case.get_env("batch")
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
    case.load_env()

    case.set_value("RUN_WITH_SUBMIT",True)
    case.flush()

    logger.warn("submit_jobs %s"%job)
    job_ids = case.submit_jobs(no_batch=no_batch, job=job)
    msg = "Submitted jobs %s"%job_ids
    append_status(msg, caseroot=caseroot, sfile="CaseStatus")

def check_case(case, caseroot):
    check_lockedfiles(caseroot)
    create_namelists(case) # Must be called before check_all_input_data
    logger.info("Checking that inputdata is available as part of case submission")
    check_all_input_data(case)
    # Now that we have baselines, do baseline operations
    if case.get_value("TEST"):
        case_cmpgen_namelists(case)

    expect(case.get_value("BUILD_COMPLETE"), "Build complete is "
           "not True please rebuild the model by calling case.build")
    logger.info("Check case OK")

def check_DA_settings(case):
    if case.get_value("DATA_ASSIMILATION"):
        script = case.get_value("DATA_ASSIMILATION_SCRIPT")
        cycles = case.get_value("DATA_ASSIMILATION_CYCLES")
        logger.info("Data Assimilation enabled using script %s with %d cycles"%(script,cycles))

