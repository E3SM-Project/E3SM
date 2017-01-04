#!/usr/bin/env python

"""
Library for saving build/run provenance.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import touch, gzip_existing_file, SharedArea, copy_umask

import tarfile, getpass, signal, glob, shutil

logger = logging.getLogger(__name__)

def _get_batch_job_id_for_syslog(case):
    """
    mach_syslog only works on certain machines
    """
    mach = case.get_value("MACH")
    if mach == 'titan':
        return os.environ["PBS_JOBID"]
    elif mach in ['edison', 'corip1']:
        return os.environ["SLURM_JOB_ID"]
    elif mach == 'mira':
        return os.environ["COBALT_JOBID"]
    else:
        return None

def save_build_provenance_acme(case, lid=None):
    cimeroot = case.get_value("CIMEROOT")
    exeroot = case.get_value("EXEROOT")
    caseroot = case.get_value("CASEROOT")
    lid = os.environ["LID"] if lid is None else lid

    # Save git describe
    describe_prov = os.path.join(exeroot, "GIT_DESCRIBE.%s" % lid)
    if os.path.exists(describe_prov):
        os.remove(describe_prov)
    run_cmd_no_fail("git describe > %s" % describe_prov, from_dir=cimeroot)

    # Save HEAD
    headfile = os.path.join(cimeroot, ".git", "logs", "HEAD")
    headfile_prov = os.path.join(exeroot, "GIT_LOGS_HEAD.%s" % lid)
    if os.path.exists(headfile_prov):
        os.remove(headfile_prov)
    if os.path.exists(headfile):
        copy_umask(headfile, headfile_prov)

    # Save SourceMods
    sourcemods = os.path.join(caseroot, "SourceMods")
    sourcemods_prov = os.path.join(exeroot, "SourceMods.%s.tar.gz" % lid)
    if os.path.exists(sourcemods_prov):
        os.remove(sourcemods_prov)
    if os.path.isdir(sourcemods):
        with tarfile.open(sourcemods_prov, "w:gz") as tfd:
            tfd.add(sourcemods, arcname="SourceMods")

    # Save build env
    env_prov = os.path.join(exeroot, "build_environment.%s.txt" % lid)
    if os.path.exists(env_prov):
        os.remove(env_prov)
    copy_umask(os.path.join(caseroot, "software_environment.txt"), env_prov)

    # For all the just-created post-build provenance files, symlink a generic name
    # to them to indicate that these are the most recent or active.
    for item in ["GIT_DESCRIBE", "GIT_LOGS_HEAD", "SourceMods", "build_environment"]:
        globstr = "%s/%s.%s*" % (exeroot, item, lid)
        matches = glob.glob(globstr)
        expect(len(matches) < 2, "Multiple matches for glob %s should not have happened" % globstr)
        if matches:
            the_match = matches[0]
            generic_name = the_match.replace(".%s" % lid, "")
            if os.path.exists(generic_name):
                os.remove(generic_name)
            os.symlink(the_match, generic_name)

def save_build_provenance_cesm(case, lid=None): # pylint: disable=unused-argument
    pass

def save_build_provenance(case, lid=None):
    with SharedArea():
        model = case.get_value("MODEL")
        if model == "acme":
            save_build_provenance_acme(case, lid=lid)
        elif model == "cesm":
            save_build_provenance_cesm(case, lid=lid)

def save_prerun_provenance_acme(case, lid=None):
    if not case.get_value("SAVE_TIMING"):
        return

    lid = os.environ["LID"] if lid is None else lid

    timing_dir = case.get_value("SAVE_TIMING_DIR")
    if timing_dir is None or timing_dir == 'UNSET':
        logger.warning("ACME requires SAVE_TIMING_DIR to be set in order to save timings. Skipping save timings")
        return

    logger.info("timing dir is %s" % timing_dir)
    rundir = case.get_value("RUNDIR")
    blddir = case.get_value("EXEROOT")
    caseroot = case.get_value("CASEROOT")
    cimeroot = case.get_value("CIMEROOT")
    base_case = case.get_value("CASE")
    full_timing_dir = os.path.join(timing_dir, "performance_archive", getpass.getuser(), base_case, lid)
    expect(not os.path.exists(full_timing_dir), "%s already exists" % full_timing_dir)

    os.makedirs(full_timing_dir)
    expect(os.path.exists(full_timing_dir), "%s does not exists" % full_timing_dir)
    mach = case.get_value("MACH")
    compiler = case.get_value("COMPILER")

    # For some batch machines save queue info
    job_id = _get_batch_job_id_for_syslog(case)
    if mach == "mira":
        for cmd, filename in [("qstat -lf", "qstatf"), ("qstat -lf %s" % job_id, "qstatf_jobid")]:
            filename = "%s.%s" % (filename, lid)
            run_cmd_no_fail("%s > %s" % (cmd, filename), from_dir=full_timing_dir)
            gzip_existing_file(os.path.join(full_timing_dir, filename))
    elif mach in ["corip1", "edison"]:
        for cmd, filename in [("sqs -f", "sqsf"), ("sqs -w -a", "sqsw"), ("sqs -f %s" % job_id, "sqsf_jobid"), ("squeue", "squeuef")]:
            filename = "%s.%s" % (filename, lid)
            run_cmd_no_fail("%s > %s" % (cmd, filename), from_dir=full_timing_dir)
            gzip_existing_file(os.path.join(full_timing_dir, filename))
    elif mach == "titan":
        for cmd, filename in [("xtdb2proc -f", "xtdb2proc"),
                              ("qstat -f >", "qstatf"),
                              ("qstat -f %s >" % job_id, "qstatf_jobid"),
                              ("xtnodestat >", "xtnodestat"),
                              ("showq >", "showq")]:
            full_cmd = cmd + " " + filename
            run_cmd_no_fail(full_cmd + "." + lid, from_dir=full_timing_dir)
            gzip_existing_file(os.path.join(full_timing_dir, filename + "." + lid))

        mdiag_reduce = os.path.join(full_timing_dir, "mdiag_reduce." + lid)
        run_cmd_no_fail("./mdiag_reduce.csh > %s" % mdiag_reduce, from_dir=os.path.join(caseroot, "Tools"))
        gzip_existing_file(mdiag_reduce)

    # copy/tar SourceModes
    source_mods_dir = os.path.join(caseroot, "SourceMods")
    if os.path.isdir(source_mods_dir):
        with tarfile.open(os.path.join(full_timing_dir, "SourceMods.%s.tar.gz" % lid), "w:gz") as tfd:
            tfd.add(source_mods_dir, arcname="SourceMods")

    # Save various case configuration items
    case_docs = os.path.join(full_timing_dir, "CaseDocs.%s" % lid)
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
            copy_umask(item, os.path.join(case_docs, os.path.basename(item) + "." + lid))

    # Copy some items from build provenance
    blddir_globs_to_copy = [
        "GIT_LOGS_HEAD",
        "build_environment.txt"
        ]
    for blddir_glob_to_copy in blddir_globs_to_copy:
        for item in glob.glob(os.path.join(blddir, blddir_glob_to_copy)):
            copy_umask(item, os.path.join(full_timing_dir, os.path.basename(item) + "." + lid))

    # What this block does is mysterious to me (JGF)
    if job_id is not None:
        sample_interval = case.get_value("SYSLOG_N")
        if sample_interval > 0:
            archive_checkpoints = os.path.join(full_timing_dir, "checkpoints.%s" % lid)
            os.mkdir(archive_checkpoints)
            touch("%s/acme.log.%s" % (rundir, lid))
            syslog_jobid = run_cmd_no_fail("./mach_syslog %d %s %s %s %s/timing/checkpoints %s >& /dev/null & echo $!" %
                                           (sample_interval, job_id, lid, rundir, rundir, archive_checkpoints),
                                           from_dir=os.path.join(caseroot, "Tools"))
            with open(os.path.join(rundir, "syslog_jobid.%s" % job_id), "w") as fd:
                fd.write("%s\n" % syslog_jobid)

    # Save state of repo
    run_cmd_no_fail("git describe > %s" % os.path.join(full_timing_dir, "GIT_DESCRIBE.%s" % lid), from_dir=os.path.dirname(cimeroot))

def save_prerun_provenance_cesm(case, lid=None): # pylint: disable=unused-argument
    pass

def save_prerun_provenance(case, lid=None):
    with SharedArea():
        model = case.get_value("MODEL")
        if model == "acme":
            save_prerun_provenance_acme(case, lid=lid)
        elif model == "cesm":
            save_prerun_provenance_cesm(case, lid=lid)

def save_postrun_provenance_cesm(case, lid=None):
    save_timing = case.get_value("SAVE_TIMING")
    if save_timing:
        lid = os.environ["LID"] if lid is None else lid
        rundir = case.get_value("RUNDIR")
        timing_dir = case.get_value("SAVE_TIMING_DIR")
        timing_dir = os.path.join(timing_dir, case.get_value("CASE"))
        shutil.move(os.path.join(rundir,"timing"),
                    os.path.join(timing_dir,"timing."+lid))

def save_postrun_provenance_acme(case, lid):
    save_timing = case.get_value("SAVE_TIMING")
    if not save_timing:
        return

    lid = os.environ["LID"] if lid is None else lid

    rundir = case.get_value("RUNDIR")
    timing_dir = case.get_value("SAVE_TIMING_DIR")
    caseroot = case.get_value("CASEROOT")
    mach = case.get_value("MACH")
    base_case = case.get_value("CASE")
    full_timing_dir = os.path.join(timing_dir, "performance_archive", getpass.getuser(), base_case, lid)

    # Kill mach_syslog
    job_id = _get_batch_job_id_for_syslog(case)
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
    rundir_timing_dir = os.path.join(rundir, "timing." + lid)
    shutil.move(os.path.join(rundir, "timing"), rundir_timing_dir)
    with tarfile.open("%s.tar.gz" % rundir_timing_dir, "w:gz") as tfd:
        tfd.add(rundir_timing_dir, arcname=os.path.basename(rundir_timing_dir))

    shutil.rmtree(rundir_timing_dir)
    copy_umask("%s.tar.gz" % rundir_timing_dir, full_timing_dir)

    gzip_existing_file(os.path.join(caseroot, "timing", "acme_timing_stats.%s" % lid))

    # JGF: not sure why we do this
    timing_saved_file = "timing.%s.saved" % lid
    touch(os.path.join(caseroot, "timing", timing_saved_file))

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
    globs_to_copy.append("timing/*.%s*" % lid)
    globs_to_copy.append("CaseStatus")

    for glob_to_copy in globs_to_copy:
        for item in glob.glob(os.path.join(caseroot, glob_to_copy)):
            basename = os.path.basename(item)
            if basename != timing_saved_file:
                if lid not in basename and not basename.endswith(".gz"):
                    copy_umask(item, os.path.join(full_timing_dir, "%s.%s" % (basename, lid)))
                else:
                    copy_umask(item, full_timing_dir)

    # zip everything
    for root, _, files in os.walk(full_timing_dir):
        for filename in files:
            if not filename.endswith(".gz"):
                gzip_existing_file(os.path.join(root, filename))

def save_postrun_provenance(case, lid=None):
    with SharedArea():
        model = case.get_value("MODEL")
        if model == "acme":
            save_postrun_provenance_acme(case, lid=lid)
        elif model == "cesm":
            save_postrun_provenance_cesm(case, lid=lid)
