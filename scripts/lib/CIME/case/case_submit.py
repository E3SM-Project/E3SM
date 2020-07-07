#!/usr/bin/env python

"""
case.submit - Submit a cesm workflow to the queueing system or run it
if there is no queueing system.  A cesm workflow may include multiple
jobs.
submit, check_case and check_da_settings are members of class Case in file case.py
"""
from six.moves                      import configparser
from CIME.XML.standard_module_setup import *
from CIME.utils                     import expect, run_and_log_case_status, verbatim_success_msg, CIMEError
from CIME.locked_files              import unlock_file, lock_file
from CIME.test_status               import *

import socket

logger = logging.getLogger(__name__)

def _build_prereq_str(case, prev_job_ids):
    delimiter = case.get_value("depend_separator")
    prereq_str = ""
    for job_id in prev_job_ids.values():
        prereq_str += str(job_id) + delimiter
    return prereq_str[:-1]

def _submit(case, job=None, no_batch=False, prereq=None, allow_fail=False, resubmit=False,
            resubmit_immediate=False, skip_pnl=False, mail_user=None, mail_type=None,
            batch_args=None, workflow=True):
    if job is None:
        job = case.get_first_job()

    # Check mediator
    hasMediator = True
    comp_classes = case.get_values("COMP_CLASSES")
    if 'CPL' not in comp_classes:
        hasMediator = False

    # Check if CONTINUE_RUN value makes sense
    if job != "case.test" and case.get_value("CONTINUE_RUN") and hasMediator:
        rundir = case.get_value("RUNDIR")
        expect(os.path.isdir(rundir),
               "CONTINUE_RUN is true but RUNDIR {} does not exist".format(rundir))
        # only checks for the first instance in a multidriver case
        if case.get_value("COMP_INTERFACE") == "nuopc":
            rpointer = "rpointer.cpl"
        elif case.get_value("MULTI_DRIVER"):
            rpointer = "rpointer.drv_0001"
        else:
            rpointer = "rpointer.drv"
        expect(os.path.exists(os.path.join(rundir,rpointer)),
               "CONTINUE_RUN is true but this case does not appear to have restart files staged in {} {}".format(rundir,rpointer))
        # Finally we open the rpointer file and check that it's correct
        casename = case.get_value("CASE")
        with open(os.path.join(rundir,rpointer), "r") as fd:
            ncfile = fd.readline().strip()
            expect(ncfile.startswith(casename) and
                   os.path.exists(os.path.join(rundir,ncfile)),
                   "File {ncfile} not present or does not match case {casename}".
                   format(ncfile=os.path.join(rundir,ncfile),casename=casename))

    # if case.submit is called with the no_batch flag then we assume that this
    # flag will stay in effect for the duration of the RESUBMITs
    env_batch = case.get_env("batch")
    external_workflow = case.get_value("EXTERNAL_WORKFLOW")
    if env_batch.get_batch_system_type() == "none" or resubmit and external_workflow:
        no_batch = True

    if no_batch:
        batch_system = "none"
    else:
        batch_system = env_batch.get_batch_system_type()
    unlock_file(os.path.basename(env_batch.filename))
    case.set_value("BATCH_SYSTEM", batch_system)

    env_batch_has_changed = False
    if not external_workflow:
        try:
            case.check_lockedfile(os.path.basename(env_batch.filename))
        except:
            env_batch_has_changed = True

    if batch_system != "none" and env_batch_has_changed and not external_workflow:
        # May need to regen batch files if user made batch setting changes (e.g. walltime, queue, etc)
        logger.warning(\
"""
env_batch.xml appears to have changed, regenerating batch scripts
manual edits to these file will be lost!
""")
        env_batch.make_all_batch_files(case)
    case.flush()
    lock_file(os.path.basename(env_batch.filename))

    if resubmit:
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
        except CIMEError:
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

        case.check_case()
        if job == case.get_primary_job():
            case.check_DA_settings()
            if case.get_value("MACH") == "mira":
                with open(".original_host", "w") as fd:
                    fd.write( socket.gethostname())

    #Load Modules
    case.load_env()

    case.flush()

    logger.warning("submit_jobs {}".format(job))
    job_ids = case.submit_jobs(no_batch=no_batch, job=job, prereq=prereq,
                               skip_pnl=skip_pnl, resubmit_immediate=resubmit_immediate,
                               allow_fail=allow_fail, mail_user=mail_user,
                               mail_type=mail_type, batch_args=batch_args, workflow=workflow)

    xml_jobids = []
    for jobname, jobid in job_ids.items():
        logger.info("Submitted job {} with id {}".format(jobname, jobid))
        if jobid:
            xml_jobids.append("{}:{}".format(jobname, jobid))

    xml_jobid_text = ", ".join(xml_jobids)
    if xml_jobid_text:
        case.set_value("JOB_IDS", xml_jobid_text)

    return xml_jobid_text

def submit(self, job=None, no_batch=False, prereq=None, allow_fail=False, resubmit=False,
           resubmit_immediate=False, skip_pnl=False, mail_user=None, mail_type=None,
           batch_args=None, workflow=True):
    if resubmit_immediate and self.get_value("MACH") in ['mira', 'cetus']:
        logger.warning("resubmit_immediate does not work on Mira/Cetus, submitting normally")
        resubmit_immediate = False

    caseroot = self.get_value("CASEROOT")
    if self.get_value("TEST"):
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

    # If this is a resubmit check the hidden file .submit_options for
    # any submit options used on the original submit and use them again
    submit_options = os.path.join(caseroot, ".submit_options")
    if resubmit and os.path.exists(submit_options):
        config = configparser.RawConfigParser()
        config.read(submit_options)
        if not skip_pnl and config.has_option('SubmitOptions','skip_pnl'):
            skip_pnl = config.getboolean('SubmitOptions', 'skip_pnl')
        if mail_user is None and config.has_option('SubmitOptions', 'mail_user'):
            mail_user = config.get('SubmitOptions', 'mail_user')
        if mail_type is None and config.has_option('SubmitOptions', 'mail_type'):
            mail_type = str(config.get('SubmitOptions', 'mail_type')).split(',')
        if batch_args is None and config.has_option('SubmitOptions', 'batch_args'):
            batch_args = config.get('SubmitOptions', 'batch_args')

    try:
        functor = lambda: _submit(self, job=job, no_batch=no_batch, prereq=prereq,
                                  allow_fail=allow_fail, resubmit=resubmit,
                                  resubmit_immediate=resubmit_immediate, skip_pnl=skip_pnl,
                                  mail_user=mail_user, mail_type=mail_type,
                                  batch_args=batch_args, workflow=workflow)
        run_and_log_case_status(functor, "case.submit", caseroot=caseroot,
                                custom_success_msg_functor=verbatim_success_msg)
    except BaseException: # Want to catch KeyboardInterrupt too
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

    if self.get_value('COMP_WAV') == 'ww':
        # the ww3 buildnml has dependencies on inputdata so we must run it again
        self.create_namelists(component='WAV')

    expect(self.get_value("BUILD_COMPLETE"), "Build complete is "
           "not True please rebuild the model by calling case.build")
    logger.info("Check case OK")

def check_DA_settings(self):
    script = self.get_value("DATA_ASSIMILATION_SCRIPT")
    cycles = self.get_value("DATA_ASSIMILATION_CYCLES")
    if len(script) > 0 and os.path.isfile(script) and cycles > 0:
        logger.info("Data Assimilation enabled using script {} with {:d} cycles".format(script,
                                                                                        cycles))
