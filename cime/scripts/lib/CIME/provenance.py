#!/usr/bin/env python

"""
Library for saving build/run provenance.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import touch, gzip_existing_file, SharedArea, copy_umask, convert_to_babylonian_time, get_current_commit, indent_string, run_cmd, run_cmd_no_fail

import tarfile, getpass, signal, glob, shutil, sys

logger = logging.getLogger(__name__)

def _get_batch_job_id_for_syslog(case):
    """
    mach_syslog only works on certain machines
    """
    mach = case.get_value("MACH")
    try:
        if mach in ['anvil', 'titan']:
            return os.environ["PBS_JOBID"]
        elif mach in ['edison', 'cori-haswell', 'cori-knl']:
            return os.environ["SLURM_JOB_ID"]
        elif mach in ['mira', 'theta']:
            return os.environ["COBALT_JOBID"]
        elif mach in ['summit']:
            return os.environ["LSB_JOBID"]
    except:
        pass

    return None

def _save_build_provenance_e3sm(case, lid):
    cimeroot = case.get_value("CIMEROOT")
    exeroot = case.get_value("EXEROOT")
    caseroot = case.get_value("CASEROOT")

    # Save git describe
    describe_prov = os.path.join(exeroot, "GIT_DESCRIBE.{}".format(lid))
    desc = get_current_commit(tag=True, repo=cimeroot)
    with open(describe_prov, "w") as fd:
        fd.write(desc)

    # Save HEAD
    headfile = os.path.join(cimeroot, ".git", "logs", "HEAD")
    headfile_prov = os.path.join(exeroot, "GIT_LOGS_HEAD.{}".format(lid))
    if os.path.exists(headfile_prov):
        os.remove(headfile_prov)
    if os.path.exists(headfile):
        copy_umask(headfile, headfile_prov)

    # Save SourceMods
    sourcemods = os.path.join(caseroot, "SourceMods")
    sourcemods_prov = os.path.join(exeroot, "SourceMods.{}.tar.gz".format(lid))
    if os.path.exists(sourcemods_prov):
        os.remove(sourcemods_prov)
    if os.path.isdir(sourcemods):
        with tarfile.open(sourcemods_prov, "w:gz") as tfd:
            tfd.add(sourcemods, arcname="SourceMods")

    # Save build env
    env_prov = os.path.join(exeroot, "build_environment.{}.txt".format(lid))
    if os.path.exists(env_prov):
        os.remove(env_prov)
    env_module = case.get_env("mach_specific")
    env_module.save_all_env_info(env_prov)

    # For all the just-created post-build provenance files, symlink a generic name
    # to them to indicate that these are the most recent or active.
    for item in ["GIT_DESCRIBE", "GIT_LOGS_HEAD", "SourceMods", "build_environment"]:
        globstr = "{}/{}.{}*".format(exeroot, item, lid)
        matches = glob.glob(globstr)
        expect(len(matches) < 2, "Multiple matches for glob {} should not have happened".format(globstr))
        if matches:
            the_match = matches[0]
            generic_name = the_match.replace(".{}".format(lid), "")
            if os.path.exists(generic_name):
                os.remove(generic_name)
            os.symlink(the_match, generic_name)

def _save_build_provenance_cesm(case, lid): # pylint: disable=unused-argument
    version = case.get_value("MODEL_VERSION")
    # version has already been recorded
    srcroot = case.get_value("SRCROOT")
    manic = os.path.join("manage_externals","checkout_externals")
    manic_full_path = os.path.join(srcroot, manic)
    out = None
    if os.path.exists(manic_full_path):
        args = " --status --verbose --no-logging"
        stat, out, err = run_cmd(manic_full_path + args, from_dir=srcroot)
        errmsg = """Error gathering provenance information from manage_externals.

manage_externals error message:
{err}

manage_externals output:
{out}

To solve this, either:

(1) Find and fix the problem: From {srcroot}, try to get this command to work:
    {manic}{args}

(2) If you don't need provenance information, rebuild with --skip-provenance-check
""".format(out=indent_string(out, 4), err=indent_string(err, 4),
           srcroot=srcroot, manic=manic, args=args)
        expect(stat==0,errmsg)

    caseroot = case.get_value("CASEROOT")
    with open(os.path.join(caseroot, "CaseStatus"), "a") as fd:
        if version is not None and version != "unknown":
            fd.write("CESM version is {}\n".format(version))
        if out is not None:
            fd.write("{}\n".format(out))

def save_build_provenance(case, lid=None):
    with SharedArea():
        model = case.get_value("MODEL")
        lid = os.environ["LID"] if lid is None else lid

        if model == "e3sm":
            _save_build_provenance_e3sm(case, lid)
        elif model == "cesm":
            _save_build_provenance_cesm(case, lid)

def _save_prerun_timing_e3sm(case, lid):
    project = case.get_value("PROJECT", subgroup=case.get_primary_job())
    if not case.is_save_timing_dir_project(project):
        return

    timing_dir = case.get_value("SAVE_TIMING_DIR")
    if timing_dir is None or not os.path.isdir(timing_dir):
        logger.warning("SAVE_TIMING_DIR {} is not valid. E3SM requires a valid SAVE_TIMING_DIR to archive timing data.".format(timing_dir))
        return

    logger.info("Archiving timing data and associated provenance in {}.".format(timing_dir))
    rundir = case.get_value("RUNDIR")
    blddir = case.get_value("EXEROOT")
    caseroot = case.get_value("CASEROOT")
    cimeroot = case.get_value("CIMEROOT")
    base_case = case.get_value("CASE")
    full_timing_dir = os.path.join(timing_dir, "performance_archive", getpass.getuser(), base_case, lid)
    if os.path.exists(full_timing_dir):
        logger.warning("{} already exists. Skipping archive of timing data and associated provenance.".format(full_timing_dir))
        return

    try:
        os.makedirs(full_timing_dir)
    except OSError:
        logger.warning("{} cannot be created. Skipping archive of timing data and associated provenance.".format(full_timing_dir))
        return

    mach = case.get_value("MACH")
    compiler = case.get_value("COMPILER")

    # For some batch machines save queue info
    job_id = _get_batch_job_id_for_syslog(case)
    if job_id is not None:
        if mach == "mira":
            for cmd, filename in [("qstat -f", "qstatf"), ("qstat -lf %s" % job_id, "qstatf_jobid")]:
                filename = "%s.%s" % (filename, lid)
                run_cmd_no_fail(cmd, arg_stdout=filename, from_dir=full_timing_dir)
                gzip_existing_file(os.path.join(full_timing_dir, filename))
        elif mach == "theta":
            for cmd, filename in [("qstat -l --header JobID:JobName:User:Project:WallTime:QueuedTime:Score:RunTime:TimeRemaining:Nodes:State:Location:Mode:Command:Args:Procs:Queue:StartTime:attrs:Geometry", "qstatf"),
                                  ("qstat -lf %s" % job_id, "qstatf_jobid"),
                                  ("xtnodestat", "xtnodestat"),
                                  ("xtprocadmin", "xtprocadmin")]:
                filename = "%s.%s" % (filename, lid)
                run_cmd_no_fail(cmd, arg_stdout=filename, from_dir=full_timing_dir)
                gzip_existing_file(os.path.join(full_timing_dir, filename))
        elif mach in ["edison", "cori-haswell", "cori-knl"]:
            for cmd, filename in [("sinfo -a -l", "sinfol"), ("sqs -f %s" % job_id, "sqsf_jobid"),
                                  # ("sqs -f", "sqsf"),
                                  ("squeue -o '%.10i %.15P %.20j %.10u %.7a %.2t %.6D %.8C %.10M %.10l %.20S %.20V'", "squeuef"),
                                  ("squeue -t R -o '%.10i %R'", "squeues")]:
                filename = "%s.%s" % (filename, lid)
                run_cmd_no_fail(cmd, arg_stdout=filename, from_dir=full_timing_dir)
                gzip_existing_file(os.path.join(full_timing_dir, filename))
        elif mach == "titan":
            for cmd, filename in [("qstat -f %s >" % job_id, "qstatf_jobid"),
                                  ("xtnodestat >", "xtnodestat"),
                                  # ("qstat -f >", "qstatf"),
                                  # ("xtdb2proc -f", "xtdb2proc"),
                                  ("showq >", "showq")]:
                full_cmd = cmd + " " + filename
                run_cmd_no_fail(full_cmd + "." + lid, from_dir=full_timing_dir)
                gzip_existing_file(os.path.join(full_timing_dir, filename + "." + lid))

            # mdiag_reduce = os.path.join(full_timing_dir, "mdiag_reduce." + lid)
            # run_cmd_no_fail("./mdiag_reduce.csh", arg_stdout=mdiag_reduce, from_dir=os.path.join(caseroot, "Tools"))
            # gzip_existing_file(mdiag_reduce)
        elif mach == "anvil":
            for cmd, filename in [("qstat -f -1 acme >", "qstatf"),
                                  ("qstat -f %s >" % job_id, "qstatf_jobid"),
                                  ("qstat -r acme >", "qstatr")]:
                full_cmd = cmd + " " + filename
                run_cmd_no_fail(full_cmd + "." + lid, from_dir=full_timing_dir)
                gzip_existing_file(os.path.join(full_timing_dir, filename + "." + lid))
        elif mach == "summit":
            for cmd, filename in [("bjobs -u all >", "bjobsu_all"),
                                  ("bjobs -r -u all -o 'jobid slots exec_host' >", "bjobsru_allo"),
                                  ("bjobs -l -UF %s >" % job_id, "bjobslUF_jobid")]:
                full_cmd = cmd + " " + filename
                run_cmd_no_fail(full_cmd + "." + lid, from_dir=full_timing_dir)
                gzip_existing_file(os.path.join(full_timing_dir, filename + "." + lid))

    # copy/tar SourceModes
    source_mods_dir = os.path.join(caseroot, "SourceMods")
    if os.path.isdir(source_mods_dir):
        with tarfile.open(os.path.join(full_timing_dir, "SourceMods.{}.tar.gz".format(lid)), "w:gz") as tfd:
            tfd.add(source_mods_dir, arcname="SourceMods")

    # Save various case configuration items
    case_docs = os.path.join(full_timing_dir, "CaseDocs.{}".format(lid))
    os.mkdir(case_docs)
    globs_to_copy = [
        "CaseDocs/*",
        "*.run",
        ".*.run",
        "*.xml",
        "user_nl_*",
        "*env_mach_specific*",
        "Macros*",
        "README.case",
        "Depends.{}".format(mach),
        "Depends.{}".format(compiler),
        "Depends.{}.{}".format(mach, compiler),
        "software_environment.txt"
        ]
    for glob_to_copy in globs_to_copy:
        for item in glob.glob(os.path.join(caseroot, glob_to_copy)):
            copy_umask(item, os.path.join(case_docs, "{}.{}".format(os.path.basename(item).lstrip("."), lid)))

    # Copy some items from build provenance
    blddir_globs_to_copy = [
        "GIT_LOGS_HEAD",
        "build_environment.txt"
        ]
    for blddir_glob_to_copy in blddir_globs_to_copy:
        for item in glob.glob(os.path.join(blddir, blddir_glob_to_copy)):
            copy_umask(item, os.path.join(full_timing_dir, os.path.basename(item) + "." + lid))

    # Save state of repo
    from_repo = cimeroot if os.path.exists(os.path.join(cimeroot, ".git")) else os.path.dirname(cimeroot)
    desc = get_current_commit(tag=True, repo=from_repo)
    with open(os.path.join(full_timing_dir, "GIT_DESCRIBE.{}".format(lid)), "w") as fd:
        fd.write(desc)

    # What this block does is mysterious to me (JGF)
    if job_id is not None:

        # Kill mach_syslog from previous run if one exists
        syslog_jobid_path = os.path.join(rundir, "syslog_jobid.{}".format(job_id))
        if os.path.exists(syslog_jobid_path):
            try:
                with open(syslog_jobid_path, "r") as fd:
                    syslog_jobid = int(fd.read().strip())
                os.kill(syslog_jobid, signal.SIGTERM)
            except (ValueError, OSError) as e:
                logger.warning("Failed to kill syslog: {}".format(e))
            finally:
                os.remove(syslog_jobid_path)

        # If requested, spawn a mach_syslog process to monitor job progress
        sample_interval = case.get_value("SYSLOG_N")
        if sample_interval > 0:
            archive_checkpoints = os.path.join(full_timing_dir, "checkpoints.{}".format(lid))
            os.mkdir(archive_checkpoints)
            touch("{}/e3sm.log.{}".format(rundir, lid))
            syslog_jobid = run_cmd_no_fail("./mach_syslog {si} {jobid} {lid} {rundir} {rundir}/timing/checkpoints {ac} >& /dev/null & echo $!".format(si=sample_interval, jobid=job_id, lid=lid, rundir=rundir, ac=archive_checkpoints),
                                           from_dir=os.path.join(caseroot, "Tools"))
            with open(os.path.join(rundir, "syslog_jobid.{}".format(job_id)), "w") as fd:
                fd.write("{}\n".format(syslog_jobid))

def _save_prerun_provenance_e3sm(case, lid):
    if case.get_value("SAVE_TIMING"):
        _save_prerun_timing_e3sm(case, lid)

def _save_prerun_provenance_cesm(case, lid): # pylint: disable=unused-argument
    pass

def save_prerun_provenance(case, lid=None):
    with SharedArea():
        # Always save env
        lid = os.environ["LID"] if lid is None else lid
        env_module = case.get_env("mach_specific")
        logdir = os.path.join(case.get_value("CASEROOT"), "logs")
        if not os.path.isdir(logdir):
            os.makedirs(logdir)
        env_module.save_all_env_info(os.path.join(logdir, "run_environment.txt.{}".format(lid)))

        model = case.get_value("MODEL")
        if model == "e3sm":
            _save_prerun_provenance_e3sm(case, lid)
        elif model == "cesm":
            _save_prerun_provenance_cesm(case, lid)

def _save_postrun_provenance_cesm(case, lid):
    save_timing = case.get_value("SAVE_TIMING")
    if save_timing:
        rundir = case.get_value("RUNDIR")
        timing_dir = os.path.join("timing", case.get_value("CASE"))
        shutil.move(os.path.join(rundir,"timing"),
                    os.path.join(timing_dir,"timing."+lid))

def _save_postrun_timing_e3sm(case, lid):
    caseroot = case.get_value("CASEROOT")
    rundir = case.get_value("RUNDIR")

    # tar timings
    rundir_timing_dir = os.path.join(rundir, "timing." + lid)
    shutil.move(os.path.join(rundir, "timing"), rundir_timing_dir)
    with tarfile.open("%s.tar.gz" % rundir_timing_dir, "w:gz") as tfd:
        tfd.add(rundir_timing_dir, arcname=os.path.basename(rundir_timing_dir))

    shutil.rmtree(rundir_timing_dir)

    gzip_existing_file(os.path.join(caseroot, "timing", "e3sm_timing_stats.%s" % lid))

    # JGF: not sure why we do this
    timing_saved_file = "timing.%s.saved" % lid
    touch(os.path.join(caseroot, "timing", timing_saved_file))

    project = case.get_value("PROJECT", subgroup=case.get_primary_job())
    if not case.is_save_timing_dir_project(project):
        return

    timing_dir = case.get_value("SAVE_TIMING_DIR")
    if timing_dir is None or not os.path.isdir(timing_dir):
        return

    mach = case.get_value("MACH")
    base_case = case.get_value("CASE")
    full_timing_dir = os.path.join(timing_dir, "performance_archive", getpass.getuser(), base_case, lid)

    if not os.path.isdir(full_timing_dir):
        return

    # Kill mach_syslog
    job_id = _get_batch_job_id_for_syslog(case)
    if job_id is not None:
        syslog_jobid_path = os.path.join(rundir, "syslog_jobid.{}".format(job_id))
        if os.path.exists(syslog_jobid_path):
            try:
                with open(syslog_jobid_path, "r") as fd:
                    syslog_jobid = int(fd.read().strip())
                os.kill(syslog_jobid, signal.SIGTERM)
            except (ValueError, OSError) as e:
                logger.warning("Failed to kill syslog: {}".format(e))
            finally:
                os.remove(syslog_jobid_path)

    # copy timings
    copy_umask("%s.tar.gz" % rundir_timing_dir, full_timing_dir)

    #
    # save output files and logs
    #
    globs_to_copy = []
    if job_id is not None:
        if mach == "titan":
            globs_to_copy.append("%s*OU" % job_id)
        elif mach == "anvil":
            globs_to_copy.append("/home/%s/%s*OU" % (getpass.getuser(), job_id))
        elif mach in ["mira", "theta"]:
            globs_to_copy.append("%s*error" % job_id)
            globs_to_copy.append("%s*output" % job_id)
            globs_to_copy.append("%s*cobaltlog" % job_id)
        elif mach in ["edison", "cori-haswell", "cori-knl"]:
            globs_to_copy.append("%s*run*%s" % (case.get_value("CASE"), job_id))
        elif mach == "summit":
            globs_to_copy.append("e3sm.stderr.%s" % job_id)
            globs_to_copy.append("e3sm.stdout.%s" % job_id)

    globs_to_copy.append("logs/run_environment.txt.{}".format(lid))
    globs_to_copy.append(os.path.join(rundir, "e3sm.log.{}.gz".format(lid)))
    globs_to_copy.append(os.path.join(rundir, "cpl.log.{}.gz".format(lid)))
    globs_to_copy.append("timing/*.{}*".format(lid))
    globs_to_copy.append("CaseStatus")

    for glob_to_copy in globs_to_copy:
        for item in glob.glob(os.path.join(caseroot, glob_to_copy)):
            basename = os.path.basename(item)
            if basename != timing_saved_file:
                if lid not in basename and not basename.endswith(".gz"):
                    copy_umask(item, os.path.join(full_timing_dir, "{}.{}".format(basename, lid)))
                else:
                    copy_umask(item, full_timing_dir)

    # zip everything
    for root, _, files in os.walk(full_timing_dir):
        for filename in files:
            if not filename.endswith(".gz"):
                gzip_existing_file(os.path.join(root, filename))

def _save_postrun_provenance_e3sm(case, lid):
    if case.get_value("SAVE_TIMING"):
        _save_postrun_timing_e3sm(case, lid)

def save_postrun_provenance(case, lid=None):
    with SharedArea():
        model = case.get_value("MODEL")
        lid = os.environ["LID"] if lid is None else lid

        if model == "e3sm":
            _save_postrun_provenance_e3sm(case, lid)
        elif model == "cesm":
            _save_postrun_provenance_cesm(case, lid)

_WALLTIME_BASELINE_NAME = "walltimes"
_WALLTIME_FILE_NAME     = "walltimes"
_GLOBAL_MINUMUM_TIME    = 900
_GLOBAL_WIGGLE          = 1000
_WALLTIME_TOLERANCE     = ( (600, 2.0), (1800, 1.5), (9999999999, 1.25) )

def get_recommended_test_time_based_on_past(baseline_root, test, raw=False):
    if baseline_root is not None:
        try:
            the_path = os.path.join(baseline_root, _WALLTIME_BASELINE_NAME, test, _WALLTIME_FILE_NAME)
            if os.path.exists(the_path):
                last_line = int(open(the_path, "r").readlines()[-1])
                if raw:
                    best_walltime = last_line
                else:
                    best_walltime = None
                    for cutoff, tolerance in _WALLTIME_TOLERANCE:
                        if last_line <= cutoff:
                            best_walltime = int(float(last_line) * tolerance)
                            break

                    if best_walltime < _GLOBAL_MINUMUM_TIME:
                        best_walltime = _GLOBAL_MINUMUM_TIME

                    best_walltime += _GLOBAL_WIGGLE

                return convert_to_babylonian_time(best_walltime)
        except:
            # We NEVER want a failure here to kill the run
            logger.warning("Failed to read test time: {}".format(sys.exc_info()[1]))

    return None

def save_test_time(baseline_root, test, time_seconds):
    if baseline_root is not None:
        try:
            the_dir  = os.path.join(baseline_root, _WALLTIME_BASELINE_NAME, test)
            if not os.path.exists(the_dir):
                os.makedirs(the_dir)

            the_path = os.path.join(the_dir, _WALLTIME_FILE_NAME)
            with open(the_path, "a") as fd:
                fd.write("{}\n".format(int(time_seconds)))
        except:
            # We NEVER want a failure here to kill the run
            logger.warning("Failed to store test time: {}".format(sys.exc_info()[1]))
