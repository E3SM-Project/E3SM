from CIME.XML.standard_module_setup import *
from CIME.case_submit               import submit
from CIME.utils                     import gzip_existing_file, new_lid, run_and_log_case_status, get_timestamp
from CIME.check_lockedfiles         import check_lockedfiles
from CIME.get_timing                import get_timing
from CIME.provenance                import save_prerun_provenance, save_postrun_provenance
from CIME.preview_namelists         import create_namelists
from CIME.case_st_archive           import case_st_archive, restore_from_archive

import shutil, time, sys, os, glob

logger = logging.getLogger(__name__)

###############################################################################
def pre_run_check(case, lid, skip_pnl=False):
###############################################################################

    # Pre run initialization code..
    caseroot = case.get_value("CASEROOT")
    din_loc_root = case.get_value("DIN_LOC_ROOT")
    batchsubmit = case.get_value("BATCHSUBMIT")
    mpilib = case.get_value("MPILIB")
    rundir = case.get_value("RUNDIR")
    build_complete = case.get_value("BUILD_COMPLETE")

    if case.get_value("TESTCASE") == "PFS":
        env_mach_pes = os.path.join(caseroot,"env_mach_pes.xml")
        shutil.copy(env_mach_pes,"{}.{}".format(env_mach_pes, lid))

    # check for locked files.
    check_lockedfiles(case)
    logger.debug("check_lockedfiles OK")

    # check that build is done
    expect(build_complete,
           "BUILD_COMPLETE is not true\nPlease rebuild the model interactively")
    logger.debug("build complete is {} ".format(build_complete))

    # load the module environment...
    case.load_env()

    # set environment variables
    # This is a requirement for yellowstone only
    if mpilib == "mpi-serial" and "MP_MPILIB" in os.environ:
        del os.environ["MP_MPILIB"]
    else:
        os.environ["MPILIB"] = mpilib

    if batchsubmit is None or len(batchsubmit) == 0:
        os.environ["LBQUERY"] = "FALSE"
        os.environ["BATCHQUERY"] = "undefined"
    elif batchsubmit == 'UNSET':
        os.environ["LBQUERY"] = "FALSE"
        os.environ["BATCHQUERY"] = "undefined"
    else:
        os.environ["LBQUERY"] = "TRUE"

    # create the timing directories, optionally cleaning them if needed.
    if not os.path.isdir(rundir):
        os.mkdir(rundir)

    if os.path.isdir(os.path.join(rundir, "timing")):
        shutil.rmtree(os.path.join(rundir, "timing"))

    os.makedirs(os.path.join(rundir, "timing", "checkpoints"))

    # This needs to be done everytime the LID changes in order for log files to be set up correctly
    # The following also needs to be called in case a user changes a user_nl_xxx file OR an env_run.xml
    # variable while the job is in the queue
    if not skip_pnl:
        create_namelists(case)

    logger.info("-------------------------------------------------------------------------")
    logger.info(" - Prestage required restarts into {}".format(rundir))
    logger.info(" - Case input data directory (DIN_LOC_ROOT) is {} ".format(din_loc_root))
    logger.info(" - Checking for required input datasets in DIN_LOC_ROOT")
    logger.info("-------------------------------------------------------------------------")

###############################################################################
def _run_model_impl(case, lid, skip_pnl=False):
###############################################################################

    pre_run_check(case, lid, skip_pnl=skip_pnl)

    model = case.get_value("MODEL")

    # Set OMP_NUM_THREADS
    env_mach_pes = case.get_env("mach_pes")
    comp_classes = case.get_values("COMP_CLASSES")
    thread_count = env_mach_pes.get_max_thread_count(comp_classes)
    os.environ["OMP_NUM_THREADS"] = str(thread_count)

    # Run the model
    logger.info("{} MODEL EXECUTION BEGINS HERE".format(time.strftime("%Y-%m-%d %H:%M:%S")))

    cmd = case.get_mpirun_cmd(job="case.run")
    cmd = case.get_resolved_value(cmd)
    logger.info("run command is {} ".format(cmd))

    rundir = case.get_value("RUNDIR")
    loop = True

    while loop:
        loop = False
        stat = run_cmd(cmd, from_dir=rundir)[0]
        model_logfile = os.path.join(rundir, model + ".log." + lid)
        # Determine if failure was due to a failed node, if so, try to restart
        if stat != 0:
            node_fail_re = case.get_value("NODE_FAIL_REGEX")
            if node_fail_re:
                node_fail_regex = re.compile(node_fail_re)
                model_logfile = os.path.join(rundir, model + ".log." + lid)
                if os.path.exists(model_logfile):
                    num_fails = len(node_fail_regex.findall(open(model_logfile, 'r').read()))
                    if num_fails > 0 and case.spare_nodes >= num_fails:
                        # We failed due to node failure!
                        logger.warning("Detected model run failed due to node failure, restarting")

                        # Archive the last consistent set of restart files and restore them
                        case_st_archive(case, no_resubmit=True)
                        restore_from_archive(case)

                        orig_cont = case.get_value("CONTINUE_RUN")
                        if not orig_cont:
                            case.set_value("CONTINUE_RUN", True)
                            create_namelists(case)

                        lid = new_lid()
                        loop = True

                        case.spare_nodes -= num_fails

            if not loop:
                # We failed and we're not restarting
                expect(False, "RUN FAIL: Command '{}' failed\nSee log file for details: {}".format(cmd, model_logfile))

    logger.info("{} MODEL EXECUTION HAS FINISHED".format(time.strftime("%Y-%m-%d %H:%M:%S")))

    post_run_check(case, lid)

    return lid

###############################################################################
def run_model(case, lid, skip_pnl=False):
###############################################################################
    functor = lambda: _run_model_impl(case, lid, skip_pnl=skip_pnl)
    return run_and_log_case_status(functor, "case.run", caseroot=case.get_value("CASEROOT"))

###############################################################################
def post_run_check(case, lid):
###############################################################################

    rundir = case.get_value("RUNDIR")
    model = case.get_value("MODEL")
    cpl_ninst = 1
    if case.get_value("MULTI_DRIVER"):
        cpl_ninst = case.get_value("NINST_MAX")
    cpl_logs = []
    if cpl_ninst > 1:
        for inst in range(cpl_ninst):
            cpl_logs.append(os.path.join(rundir, "cpl_%04d.log." % (inst+1) + lid))
    else:
        cpl_logs = [os.path.join(rundir, "cpl" + ".log." + lid)]
    cpl_logfile = cpl_logs[0]

    # find the last model.log and cpl.log
    model_logfile = os.path.join(rundir, model + ".log." + lid)

    if not os.path.isfile(model_logfile):
        expect(False, "Model did not complete, no {} log file ".format(model_logfile))
    elif os.stat(model_logfile).st_size == 0:
        expect(False, "Run FAILED")
    else:
        count_ok = 0
        for cpl_logfile in cpl_logs:
            if not os.path.isfile(cpl_logfile):
                break
            with open(cpl_logfile, 'r') as fd:
                if 'SUCCESSFUL TERMINATION' in fd.read():
                    count_ok += 1
        if count_ok != cpl_ninst:
            expect(False, "Model did not complete - see {} \n " .format(cpl_logfile))

###############################################################################
def save_logs(case, lid):
###############################################################################
    logdir = case.get_value("LOGDIR")
    if logdir is not None and len(logdir) > 0:
        if not os.path.isdir(logdir):
            os.makedirs(logdir)

        caseroot = case.get_value("CASEROOT")
        rundir = case.get_value("RUNDIR")
        logfiles = glob.glob(os.path.join(rundir, "*.log.{}".format(lid)))
        for logfile in logfiles:
            if os.path.isfile(logfile):
                logfile_gz = gzip_existing_file(logfile)
                shutil.copy(logfile_gz,
                            os.path.join(caseroot, logdir, os.path.basename(logfile_gz)))

###############################################################################
def resubmit_check(case):
###############################################################################

    # check to see if we need to do resubmission from this particular job,
    # Note that Mira requires special logic

    dout_s = case.get_value("DOUT_S")
    logger.warn("dout_s {} ".format(dout_s))
    mach = case.get_value("MACH")
    logger.warn("mach {} ".format(mach))
    testcase = case.get_value("TESTCASE")
    resubmit_num = case.get_value("RESUBMIT")
    logger.warn("resubmit_num {}".format(resubmit_num))
    # If dout_s is True than short-term archiving handles the resubmit
    # If dout_s is True and machine is mira submit the st_archive script
    resubmit = False
    if not dout_s and resubmit_num > 0:
        resubmit = True
    elif dout_s and mach == 'mira':
        caseroot = case.get_value("CASEROOT")
        cimeroot = case.get_value("CIMEROOT")
        cmd = "ssh cooleylogin1 'cd {}; CIMEROOT={} ./case.submit {} --job case.st_archive'".format(caseroot, cimeroot, caseroot)
        run_cmd(cmd, verbose=True)

    if resubmit:
        if testcase is not None and testcase in ['ERR']:
            job = "case.test"
        else:
            job = "case.run"
        submit(case, job=job, resubmit=True)

###############################################################################
def do_external(script_name, caseroot, rundir, lid, prefix):
###############################################################################
    filename = "{}.external.log.{}".format(prefix, lid)
    outfile = os.path.join(rundir, filename)
    cmd = script_name + " 1> {} {} 2>&1".format(outfile, caseroot)
    logger.info("running {}".format(script_name))
    run_cmd_no_fail(cmd)

###############################################################################
def do_data_assimilation(da_script, caseroot, cycle, lid, rundir):
###############################################################################
    filename = "da.log.{}".format(lid)
    outfile = os.path.join(rundir, filename)
    cmd = da_script + " 1> {} {} {:d} 2>&1".format(outfile, caseroot, cycle)
    logger.info("running {}".format(da_script))
    run_cmd_no_fail(cmd)

###############################################################################
def case_run(case, skip_pnl=False):
###############################################################################
    # Set up the run, run the model, do the postrun steps
    run_with_submit = case.get_value("RUN_WITH_SUBMIT")
    expect(run_with_submit,
           "You are not calling the run script via the submit script. "
           "As a result, short-term archiving will not be called automatically."
           "Please submit your run using the submit script like so:"
           " ./case.submit Time: {}".format(get_timestamp()))

    # Forces user to use case.submit if they re-submit
    if case.get_value("TESTCASE") is None:
        case.set_value("RUN_WITH_SUBMIT", False)

    prerun_script = case.get_value("PRERUN_SCRIPT")
    postrun_script = case.get_value("POSTRUN_SCRIPT")

    data_assimilation = case.get_value("DATA_ASSIMILATION")
    data_assimilation_cycles = case.get_value("DATA_ASSIMILATION_CYCLES")
    data_assimilation_script = case.get_value("DATA_ASSIMILATION_SCRIPT")

    # set up the LID
    lid = new_lid()

    save_prerun_provenance(case)

    for cycle in range(data_assimilation_cycles):
        # After the first DA cycle, runs are restart runs
        if cycle > 0:
            case.set_value("CONTINUE_RUN", "TRUE")
            lid = new_lid()

        if prerun_script:
            case.flush()
            do_external(prerun_script, case.get_value("CASEROOT"), case.get_value("RUNDIR"),
                        lid, prefix="prerun")
            case.read_xml()

        lid = run_model(case, lid, skip_pnl)
        save_logs(case, lid)       # Copy log files back to caseroot
        if case.get_value("CHECK_TIMING") or case.get_value("SAVE_TIMING"):
            get_timing(case, lid)     # Run the getTiming script

        if data_assimilation:
            case.flush()
            do_data_assimilation(data_assimilation_script, case.get_value("CASEROOT"), cycle, lid,
                                 case.get_value("RUNDIR"))
            case.read_xml()

        if postrun_script:
            case.flush()
            do_external(postrun_script, case.get_value("CASEROOT"), case.get_value("RUNDIR"),
                        lid, prefix="postrun")
            case.read_xml()

        save_postrun_provenance(case)

    logger.warn("check for resubmit")
    resubmit_check(case)

    return True
