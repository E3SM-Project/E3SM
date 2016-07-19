
from CIME.XML.standard_module_setup import *
from CIME.case_submit               import submit
from CIME.XML.files                 import Files
from CIME.XML.component             import Component
from CIME.XML.machines              import Machines
from CIME.utils                     import append_status, touch, gzip_existing_file
from CIME.check_lockedfiles         import check_lockedfiles
from CIME.preview_namelists         import preview_namelists
from CIME.task_maker                import TaskMaker

import shutil, time, sys, os, getpass, tarfile, glob, signal

logger = logging.getLogger(__name__)

###############################################################################
def pre_run_check(case):
###############################################################################

    # Pre run initialization code..
    caseroot = case.get_value("CASEROOT")
    din_loc_root = case.get_value("DIN_LOC_ROOT")
    batchsubmit = case.get_value("BATCHSUBMIT")
    mpilib = case.get_value("MPILIB")
    rundir = case.get_value("RUNDIR")
    build_complete = case.get_value("BUILD_COMPLETE")

    # check for locked files.
    check_lockedfiles(case.get_value("CASEROOT"))
    logger.debug("check_lockedfiles OK")

    # check that build is done
    expect (build_complete,
            "BUILD_COMPLETE is not true\nPlease rebuild the model interactively")
    logger.debug("build complete is %s " %build_complete)

    # load the module environment...
    env_module = case.get_env("mach_specific")
    env_module.load_env_for_case(compiler=case.get_value("COMPILER"),
                                 debug=case.get_value("DEBUG"),
                                 mpilib=case.get_value("MPILIB"))

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

    if os.path.isdir(os.path.join(rundir,"timing")):
        shutil.rmtree(os.path.join(rundir,"timing"))

    os.makedirs(os.path.join(rundir,"timing","checkpoints"))

    # run preview namelists
    preview_namelists(case)

    # document process
    append_status("Run started ",caseroot=caseroot,
                 sfile="CaseStatus")

    logger.info( "-------------------------------------------------------------------------")
    logger.info( " - To prestage required restarts, untar a restart.tar file into %s" %(rundir))
    logger.info( " - Case input data directory (DIN_LOC_ROOT) is %s " %(din_loc_root))
    logger.info( " - Checking for required input datasets in DIN_LOC_ROOT")
    logger.info( "-------------------------------------------------------------------------")

###############################################################################
def run_model(case):
###############################################################################

    # Set OMP_NUM_THREADS
    tm = TaskMaker(case)
    num_threads = tm.thread_count
    os.environ["OMP_NUM_THREADS"] = str(num_threads)

    # Run the model
    logger.info("%s MODEL EXECUTION BEGINS HERE" %(time.strftime("%Y-%m-%d %H:%M:%S")))

    machine = Machines(machine=case.get_value("MACH"))
    cmd = machine.get_full_mpirun(tm, case, "case.run")
    cmd = case.get_resolved_value(cmd)

    logger.info("run command is %s " %cmd)
    rundir = case.get_value("RUNDIR")
    run_cmd_no_fail(cmd, from_dir=rundir)
    logger.info( "%s MODEL EXECUTION HAS FINISHED" %(time.strftime("%Y-%m-%d %H:%M:%S")))

###############################################################################
def post_run_check(case, lid):
###############################################################################

    caseroot = case.get_value("CASEROOT")
    rundir = case.get_value("RUNDIR")
    model = case.get_value("MODEL")

    # find the last model.log and cpl.log
    model_logfile = os.path.join(rundir,model + ".log." + lid)
    cpl_logfile   = os.path.join(rundir,"cpl" + ".log." + lid)

    if not os.path.isfile(model_logfile):
        msg = "Model did not complete, no %s log file "%model_logfile
        append_status(msg, caseroot=caseroot, sfile="CaseStatus")
        expect(False, msg)
    elif not os.path.isfile(cpl_logfile):
        msg = "Model did not complete, no cpl log file"
        append_status(msg, caseroot=caseroot, sfile="CaseStatus")
        expect(False, msg)
    elif os.stat(model_logfile).st_size == 0:
        msg = " Run FAILED "
        append_status(msg, caseroot=caseroot, sfile="CaseStatus")
        expect (False, msg)
    else:
        with open(cpl_logfile, 'r') as fd:
            if 'SUCCESSFUL TERMINATION' in fd.read():
                msg = "Run SUCCESSFUL"
                append_status(msg, caseroot=caseroot, sfile="CaseStatus" )
            else:
                msg = "Model did not complete - see %s \n " %(cpl_logfile)
                append_status(msg, caseroot=caseroot, sfile="CaseStatus")
                expect (False, msg)

###############################################################################
def _get_batch_job_id(case):
###############################################################################
    mach = case.get_value("MACH")
    if mach == 'titan':
        return os.environ("PBS_JOBID")
    elif mach in ['edison', 'corip1']:
        return os.environ("SLURM_JOB_ID")
    elif mach == 'mira':
        return os.environ("COBALT_JOBID")
    else:
        return None

###############################################################################
def save_timing_setup_acme(case, lid):
###############################################################################
    if not case.get_value("SAVE_TIMING") or case.get_value("MODEL") != "acme":
        return

    timing_dir = case.get_value("SAVE_TIMING_DIR")
    if timing_dir is None or timing_dir == 'UNSET':
        logger.warning("ACME requires SAVE_TIMING_DIR to be set in order to save timings. Skipping save timings")
        return
    logger.warn("timing dir is %s"%timing_dir)
    rundir = case.get_value("RUNDIR")
    caseroot = case.get_value("CASEROOT")
    cimeroot = case.get_value("CIMEROOT")
    base_case = case.get_value("CASEBASEID")
    full_timing_dir = os.path.join(timing_dir, "performance_archive", getpass.getuser(), base_case, lid)
    expect(not os.path.exists(full_timing_dir), "%s already exists" % full_timing_dir)

    os.makedirs(full_timing_dir)
    mach = case.get_value("MACH")
    compiler = case.get_value("COMPILER")

    # For some batch machines save queue info
    job_id = _get_batch_job_id(case)
    if mach == "mira":
        for cmd, filename in [("qstat -lf", "qstatf"), ("qstat -lf %s" % job_id, "qstatf_jobid")]:
            run_cmd_no_fail("%s > %s.%s" % (cmd, filename, lid), from_dir=full_timing_dir)
            gzip_existing_file(os.path.join(full_timing_dir, filename))
    elif mach == ["corip1", "edison"]:
        for cmd, filename in [("sqs -f", "sqsf"), ("sqs -w -a", "sqsw"), ("sqs -f %s" % job_id, "sqsf_jobid"), ("squeue", "squeuef")]:
            run_cmd_no_fail("%s > %s.%s" % (cmd, filename, lid), from_dir=full_timing_dir)
            gzip_existing_file(os.path.join(full_timing_dir, filename))
    elif mach == "titan":
        for cmd, filename in [("xtdb2proc -f xtdb2proc", "xtdb2procf"),
                              ("qstat -f > qstat", "qstatf"),
                              ("qstat -f %s > qstatf_jobid" % job_id, "qstatf_jobid"),
                              ("xtnodestat > xtnodestat", "xtnodestatf"),
                              ("showq > showqf", "showqf")]:
            run_cmd_no_fail(cmd + "." + lid, from_dir=full_timing_dir)
            gzip_existing_file(os.path.join(full_timing_dir, filename + "." + lid))

        mdiag_reduce = os.path.join(full_timing_dir, "mdiag_reduce." + lid)
        run_cmd_no_fail("./mdiag_reduce.csh > %s" % mdiag_reduce, from_dir=os.path.join(caseroot, "Tools"))
        gzip_existing_file(mdiag_reduce)

    # copy/tar SourceModes
    source_mods_dir = os.path.join(caseroot, "SourceMods")
    if os.path.isdir(source_mods_dir):
        with tarfile.open(os.path.join(full_timing_dir, "SourceMods.%s.tar.gz" % lid), "w:gz") as tfd:
            tfd.add(source_mods_dir)

    # Save various case configuration items
    case_docs = os.path.join(full_timing_dir, "CaseDocs")
    os.mkdir(case_docs)
    globs_to_copy = [
        "CaseDocs/*",
        "*.run",
        "*.xml",
        "user_nl_*",
        "*env_mach_specific*",
        "Macros",
        "README.case",
        "Depends.%s" % mach,
        "Depends.%s" % compiler,
        "Depends.%s.%s" % (mach, compiler),
        "software_environment.txt"
        ]
    for glob_to_copy in globs_to_copy:
        for item in glob.glob(os.path.join(caseroot, glob_to_copy)):
            shutil.copy(item, os.path.join(case_docs, os.path.basename(item) + "." + lid))

    if job_id is not None:
        sample_interval = case.get_value("SYSLOG_N")
        if sample_interval > 0:
            archive_checkpoints = os.path.join(full_timing_dir, "checkpoints")
            os.mkdir(archive_checkpoints)
            touch("%s/acme.log.%s" % (rundir, lid))
            syslog_jobid = run_cmd_no_fail("./mach_syslog %d %s %s %s %s/timing/checkpoints %s/checkpoints >& /dev/null & echo $!" %
                                           (sample_interval, job_id, lid, rundir, rundir, archive_checkpoints),
                                           from_dir=os.path.join(caseroot, "Tools"))
            with open(os.path.join(rundir, "syslog_jobid", ".%s" % job_id), "w") as fd:
                fd.write("%s\n" % syslog_jobid)

    # Save state of repo
    run_cmd_no_fail("git describe > %s" % os.path.join(full_timing_dir, "GIT_DESCRIBE"), from_dir=cimeroot)

###############################################################################
def save_timing_cesm(case, lid):
###############################################################################
    rundir = case.get_value("RUNDIR")
    timing_dir = case.get_value("SAVE_TIMING_DIR")
    timing_dir = os.path.join(timing_dir, case.get_value("CASE"))
    shutil.move(os.path.join(rundir,"timing"),
                os.path.join(timing_dir,"timing."+lid))

###############################################################################
def save_timing_acme(case, lid):
###############################################################################
    rundir = case.get_value("RUNDIR")
    timing_dir = case.get_value("SAVE_TIMING_DIR")
    caseroot = case.get_value("CASEROOT")
    mach = case.get_value("MACH")
    base_case = case.get_value("CASEBASEID")
    full_timing_dir = os.path.join(timing_dir, "performance_archive", getpass.getuser(), base_case, lid)

    # Kill mach_syslog
    job_id = _get_batch_job_id(case)
    if job_id is not None:
        syslog_jobid_path = os.path.join(rundir, "syslog_jobid", ".%s" % job_id)
        if os.path.exists(syslog_jobid_path):
            try:
                with open(syslog_jobid_path, "r") as fd:
                    syslog_jobid = int(fd.read().strip())
                os.kill(syslog_jobid, signal.SIGTERM)
            except (ValueError, OSError) as e:
                logger.warning("Failed to kill syslog: %s" % e)
            finally:
                os.remove(syslog_jobid_path)

    # copy/tar timings
    with tarfile.open(os.path.join(full_timing_dir, "timing.%s.tar.gz" % lid), "w:gz") as tfd:
        tfd.add(os.path.join(rundir, "timing"))

    #
    # save output files and logs
    #
    globs_to_copy = []
    if mach == "titan":
        globs_to_copy.append("%s*OU" % job_id)
    elif mach == "mira":
        globs_to_copy.append("%s*output" % job_id)
        globs_to_copy.append("%s*cobaltlog" % job_id)
    elif mach in ["edison", "corip1"]:
        globs_to_copy.append("%s" % case.get_value("CASE"))

    globs_to_copy.append("logs/acme.log.%s.gz" % lid)
    globs_to_copy.append("logs/cpl.log.%s.gz" % lid)
    globs_to_copy.append("timing/*.%s" % lid)
    globs_to_copy.append("CaseStatus")

    for glob_to_copy in globs_to_copy:
        for item in glob.glob(os.path.join(caseroot, glob_to_copy)):
            shutil.copy(item, full_timing_dir)

###############################################################################
def get_timings(case, lid):
###############################################################################
    check_timing = case.get_value("CHECK_TIMING")
    if check_timing:
        caseroot = case.get_value("CASEROOT")
        timingDir = os.path.join(caseroot, "timing")
        if not os.path.isdir(timingDir):
            os.makedirs(timingDir)

        logger.info("Running timing script %s " %(os.path.join(caseroot, "Tools", "getTiming")))
        run_cmd_no_fail("%s -lid %s " % (os.path.join(caseroot, "Tools", "getTiming"), lid))

        # save the timing files if desired. Some of the details here are
        # model dependent.
        model = case.get_value("MODEL")
        save_timing = case.get_value("SAVE_TIMING")
        if save_timing:
            if model == "acme":
                save_timing_acme(case, lid)
            else:
                save_timing_cesm(case, lid)

        # compress relevant timing files
        logger.info( "gzipping timing stats.." )
        timingfile = os.path.join(timingDir, model + "_timing_stats." + lid)
        gzip_existing_file(timingfile)
        logger.info("Done with timings")

###############################################################################
def save_logs(case, lid):
###############################################################################
    logdir = case.get_value("LOGDIR")
    if logdir is not None and len(logdir) > 0:
        if not os.path.isdir(logdir):
            os.makedirs(logdir)

        caseroot = case.get_value("CASEROOT")
        rundir = case.get_value("RUNDIR")

        # get components
        files = Files()
        config_file = files.get_value("CONFIG_DRV_FILE")
        component = Component(config_file)
        comps = [x.lower() for x in component.get_valid_model_components()]
        comps = [x.replace('drv', 'cpl') for x in comps]
        model = [case.get_value("MODEL")]
        comps = comps + model

        # for each component, compress log files and copy to logdir
        for comp in comps:
            logfile = os.path.join(rundir, comp + '.log.' + lid)
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
    logger.warn("dout_s %s "%(dout_s))
    mach = case.get_value("MACH")
    logger.warn("mach %s "%(mach))
    testcase = case.get_value("TESTCASE")
    resubmit_num = case.get_value("RESUBMIT")
    logger.warn("resubmit_num %s"%(resubmit_num))
    # If dout_s is True than short-term archiving handles the resubmit
    # If dout_s is True and machine is mira submit the st_archive script
    resubmit = False
    if not dout_s and resubmit_num > 0:
        resubmit = True
    elif dout_s and mach == 'mira':
        caseroot = case.get_value("CASEROOT")
        cimeroot = case.get_value("CIMEROOT")
        cmd = "ssh cooleylogin1 'cd %s; CIMEROOT=%s ./case.submit %s --job case.st_archive' "%(caseroot, cimeroot, caseroot)
        run_cmd(cmd, verbose=True)

    if resubmit:
        if testcase is not None and testcase in ['ERR']:
            job = "case.test"
        else:
            job = "case.run"
        submit(case, job=job, resubmit=True)

###############################################################################
def do_data_assimilation(da_script, lid):
###############################################################################
    cmd = da_script + "1> da.log.%s 2>&1" %(lid)
    logger.debug("running %s" %da_script)
    run_cmd_no_fail(cmd)
    # disposeLog(case, 'da', lid)  THIS IS UNDEFINED!

###############################################################################
def case_run(case):
###############################################################################
    # Set up the run, run the model, do the postrun steps
    run_with_submit = case.get_value("RUN_WITH_SUBMIT")
    expect (run_with_submit,
            "You are not calling the run script via the submit script. "
            "As a result, short-term archiving will not be called automatically."
            "Please submit your run using the submit script like so:"
            " ./case.submit")

    data_assimilation = case.get_value("DATA_ASSIMILATION")
    data_assimilation_cycles = case.get_value("DATA_ASSIMILATION_CYCLES")
    data_assimilation_script = case.get_value("DATA_ASSIMILATION_SCRIPT")

    # set up the LID
    lid = time.strftime("%y%m%d-%H%M%S")
    os.environ["LID"] = lid

    save_timing_setup_acme(case, lid)

    for _ in range(data_assimilation_cycles):
        pre_run_check(case)
        run_model(case)
        post_run_check(case, lid)
        save_logs(case, lid)       # Copy log files back to caseroot
        get_timings(case, lid)     # Run the getTiming script
        if data_assimilation:
            do_data_assimilation(data_assimilation_script, lid)

    logger.warn("check for resubmit")
    resubmit_check(case)

    return True
