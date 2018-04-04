#!/usr/bin/env python

"""
case.submit - Submit a cesm workflow to the queueing system or run it
if there is no queueing system.  A cesm workflow may include multiple
jobs.
submit, check_case and check_da_settings are members of class Case in file case.py
"""
import socket
from CIME.XML.standard_module_setup import *
from CIME.utils                     import expect, run_and_log_case_status, verbatim_success_msg
from CIME.locked_files              import unlock_file, lock_file
from CIME.test_status               import *

logger = logging.getLogger(__name__)

def _submit(case, job=None, no_batch=False, prereq=None, resubmit=False,
            skip_pnl=False, mail_user=None, mail_type=None, batch_args=None):
    if job is None:
        job = case.get_primary_job()

    rundir = case.get_value("RUNDIR")
    continue_run = case.get_value("CONTINUE_RUN")
    expect(os.path.isdir(rundir) or not continue_run,
           " CONTINUE_RUN is true but RUNDIR {} does not exist".format(rundir))

    # if case.submit is called with the no_batch flag then we assume that this
    # flag will stay in effect for the duration of the RESUBMITs
    env_batch = case.get_env("batch")
    if resubmit:
        if env_batch.get_batch_system_type() == "none":
            no_batch = True

        # This is a resubmission, do not reinitialize test values
        if job == "case.test":
            case.set_value("IS_FIRST_RUN", False)

        resub = case.get_value("RESUBMIT")
        logger.info("Submitting job '{}', resubmit={:d}".format(job, resub))
        case.set_value("RESUBMIT", resub-1)
        if case.get_value("RESUBMIT_SETS_CONTINUE_RUN"):
            case.set_value("CONTINUE_RUN", True)

    else:
        if job == "case.test":
            case.set_value("IS_FIRST_RUN", True)

        if no_batch:
            batch_system = "none"
        else:
            batch_system = env_batch.get_batch_system_type()

        case.set_value("BATCH_SYSTEM", batch_system)

        env_batch_has_changed = False
        try:
            case.check_lockedfile(os.path.basename(env_batch.filename))
        except SystemExit:
            env_batch_has_changed = True

        if env_batch.get_batch_system_type() != "none" and env_batch_has_changed:
            # May need to regen batch files if user made batch setting changes (e.g. walltime, queue, etc)
            logger.warning(\
"""
env_batch.xml appears to have changed, regenerating batch scripts
manual edits to these file will be lost!
""")
            env_batch.make_all_batch_files(case)

        unlock_file(os.path.basename(env_batch.filename))
        lock_file(os.path.basename(env_batch.filename))

        if job == case.get_primary_job():
            case.check_case()
            case.check_DA_settings()
            if case.get_value("MACH") == "mira":
                with open(".original_host", "w") as fd:
                    fd.write( socket.gethostname())

    #Load Modules
    case.load_env()

    case.flush()

    logger.warning("submit_jobs {}".format(job))
    job_ids = case.submit_jobs(no_batch=no_batch, job=job, skip_pnl=skip_pnl,
                               prereq=prereq, mail_user=mail_user,
                               mail_type=mail_type, batch_args=batch_args)

    xml_jobids = []
    for jobname, jobid in job_ids.items():
        logger.info("Submitted job {} with id {}".format(jobname, jobid))
        if jobid:
            xml_jobids.append("{}:{}".format(jobname, jobid))

    xml_jobid_text = ", ".join(xml_jobids)
    if xml_jobid_text:
        case.set_value("JOB_IDS", xml_jobid_text)

    return xml_jobid_text

def submit(self, job=None, no_batch=False, prereq=None, resubmit=False,
           skip_pnl=False, mail_user=None, mail_type=None, batch_args=None):
    if self.get_value("TEST"):
        caseroot = self.get_value("CASEROOT")
        casebaseid = self.get_value("CASEBASEID")
        # This should take care of the race condition where the submitted job
        # begins immediately and tries to set RUN phase. We proactively assume
        # a passed SUBMIT phase. If this state is already PASS, don't set it again
        # because then we'll lose RUN phase info if it's there. This info is important
        # for system_tests_common to know if it needs to reinitialize the test or not.
        with TestStatus(test_dir=caseroot, test_name=casebaseid) as ts:
            phase_status = ts.get_status(SUBMIT_PHASE)
            if phase_status != TEST_PASS_STATUS:
                ts.set_status(SUBMIT_PHASE, TEST_PASS_STATUS)

    try:
        functor = lambda: _submit(self, job=job, no_batch=no_batch, prereq=prereq,
                                  resubmit=resubmit, skip_pnl=skip_pnl,
                                  mail_user=mail_user, mail_type=mail_type,
                                  batch_args=batch_args)
        run_and_log_case_status(functor, "case.submit", caseroot=self.get_value("CASEROOT"),
                                custom_success_msg_functor=verbatim_success_msg)
    except:
        # If something failed in the batch system, make sure to mark
        # the test as failed if we are running a test.
        if self.get_value("TEST"):
            with TestStatus(test_dir=caseroot, test_name=casebaseid) as ts:
                ts.set_status(SUBMIT_PHASE, TEST_FAIL_STATUS)

        raise

def check_case(self):
    self.check_lockedfiles()
    self.create_namelists() # Must be called before check_all_input_data
    logger.info("Checking that inputdata is available as part of case submission")
    self.check_all_input_data()

    expect(self.get_value("BUILD_COMPLETE"), "Build complete is "
           "not True please rebuild the model by calling case.build")
    logger.info("Check case OK")

def check_DA_settings(self):
    script = self.get_value("DATA_ASSIMILATION_SCRIPT")
    cycles = self.get_value("DATA_ASSIMILATION_CYCLES")
    if len(script) > 0 and os.path.isfile(script) and cycles > 0:
        logger.info("Data Assimilation enabled using script {} with {:d} cycles".format(script,cycles))
