import os
import tarfile
import glob
import re
import logging
import getpass
import signal
import shutil
from contextlib import contextmanager

from CIME import utils

logger = logging.getLogger(__name__)


def save_build_provenance(case, lid=None):
    with utils.SharedArea():
        lid = os.environ["LID"] if lid is None else lid

        srcroot = case.get_value("SRCROOT")
        exeroot = case.get_value("EXEROOT")
        caseroot = case.get_value("CASEROOT")

        _record_git_provenance(srcroot, exeroot, lid)

        _archive_source_mods(lid, exeroot, caseroot)

        _archive_build_environment(lid, exeroot, case)

        _archive_build_times(lid, exeroot)

        _symlink_current(lid, exeroot)


def _record_git_provenance(srcroot, exeroot, lid):
    """Records git provenance

    Records git status, diff and logs for main repo and all submodules.
    """
    # Save git describe
    describe_prov = os.path.join(exeroot, "GIT_DESCRIBE.{}".format(lid))
    desc = utils.get_current_commit(tag=True, repo=srcroot)
    with open(describe_prov, "w") as fd:
        fd.write(desc)

    gitroot = _find_git_root(srcroot)

    # Save HEAD
    headfile = os.path.join(gitroot, "logs", "HEAD")
    headfile_prov = os.path.join(exeroot, "GIT_LOGS_HEAD.{}".format(lid))
    if os.path.exists(headfile_prov):
        os.remove(headfile_prov)
    if os.path.exists(headfile):
        utils.safe_copy(headfile, headfile_prov, preserve_meta=False)

    # Save git submodule status
    submodule_prov = os.path.join(exeroot, "GIT_SUBMODULE_STATUS.{}".format(lid))
    subm_status = utils.get_current_submodule_status(recursive=True, repo=srcroot)
    with open(submodule_prov, "w") as fd:
        fd.write(subm_status)

    # Git Status
    status_prov = os.path.join(exeroot, "GIT_STATUS.{}".format(lid))
    _run_git_cmd_recursively("status", srcroot, status_prov)

    # Git Diff
    diff_prov = os.path.join(exeroot, "GIT_DIFF.{}".format(lid))
    _run_git_cmd_recursively("diff", srcroot, diff_prov)

    # Git Log
    log_prov = os.path.join(exeroot, "GIT_LOG.{}".format(lid))
    cmd = "log --first-parent --pretty=oneline -n 5"
    _run_git_cmd_recursively(cmd, srcroot, log_prov)

    # Git remote
    remote_prov = os.path.join(exeroot, "GIT_REMOTE.{}".format(lid))
    _run_git_cmd_recursively("remote -v", srcroot, remote_prov)

    # Git config
    config_src = os.path.join(gitroot, "config")
    if os.path.exists(config_src):
        config_prov = os.path.join(exeroot, "GIT_CONFIG.{}".format(lid))
        utils.safe_copy(config_src, config_prov, preserve_meta=False)


def _find_git_root(srcroot):
    """Finds the `.git` directory.

    NOTICE: It's assumed `srcroot` is an absolute path.

    There are three scenarios to locate the correct `.git` directory:

    NOTICE: In the 3rd case git `.git` directory is not actually `.git`.

    1. In a standard git repository it will be located at `{srcroot}/.git`.
    2. In a git worktree the `{srcroot}/.git` is a file containing a path,
       `{gitdir}`, which the `{gitroot}` can be parsed from.
    3. In a git submodule the `{srcroot}/.git` is a file containing a path,
       `{gitdir}`, where the `{gitroot}` is `{gitdir}`.

    To aid in finding the correct `{gitroot}` the file `{gitroot}/config`
    is checked, this file will always be located in the correct directory.

    Args:
        srcroot (str): Path to the source root.

    Returns:
        str: Absolute path to `.git` directory.
    """
    gitroot = f"{srcroot}/.git"

    utils.expect(
        os.path.exists(gitroot),
        f"{srcroot!r} is not a git repository, failed to collect provenance",
    )

    # Handle 1st scenario
    if os.path.isdir(gitroot):
        return gitroot

    # ensure we're in correct directory so abspath works correctly
    with _swap_cwd(srcroot):
        gitdir = os.path.abspath(_read_gitdir(gitroot))

    # Handle 3rd scenario, gitdir is the `.git` directory
    config_path = os.path.join(gitdir, "config")

    if os.path.exists(config_path):
        return gitdir

    # Handle 2nd scenario, gitdir should already be absolute path
    parsed_gitroot = _parse_dot_git_path(gitdir)

    return parsed_gitroot


@contextmanager
def _swap_cwd(new_cwd):
    old_cwd = os.getcwd()

    os.chdir(new_cwd)

    try:
        yield
    finally:
        os.chdir(old_cwd)


def _run_git_cmd_recursively(cmd, srcroot, output):
    """Runs a git command recursively

    Runs the git command in srcroot then runs it on each submodule.
    Then output from both commands is written to the output file.
    """
    rc1, output1, err1 = utils.run_cmd("git {}".format(cmd), from_dir=srcroot)

    rc2, output2, err2 = utils.run_cmd(
        'git submodule foreach --recursive "git {}; echo"'.format(cmd), from_dir=srcroot
    )

    with open(output, "w") as fd:
        fd.write((output1 if rc1 == 0 else err1) + "\n\n")
        fd.write((output2 if rc2 == 0 else err2) + "\n")


def _read_gitdir(gitroot):
    """Read `gitdir` from `.git` file.

    Reads `.git` file in a worktree or submodule and parse `gitdir`.

    Args:
        gitroot (str): Path ending with `.git` file.

    Returns:
        str: Path contained in `.git` file.
    """
    utils.expect(os.path.isfile(gitroot), f"Expected {gitroot!r} to be a file")

    with open(gitroot) as fd:
        line = fd.readline()

    gitdir_pattern = r"^gitdir:\s*(.*)$"

    m = re.match(gitdir_pattern, line)

    utils.expect(m is not None, f"Could parse gitdir path from file {gitroot!r}")

    return m.group(1)


def _parse_dot_git_path(gitdir):
    """Parse `.git` path.

    Take a path e.g. `/storage/cime/.git/worktrees/cime` and parse `.git`
    directory e.g. `/storage/cime/.git`.

    Args:
        gitdir (str): Path containing `.git` directory.

    Returns:
        str: Path ending with `.git`
    """
    dot_git_pattern = r"^(.*/\.git).*"

    m = re.match(dot_git_pattern, gitdir)

    utils.expect(m is not None, f"Could not parse .git from path {gitdir!r}")

    return m.group(1)


def _archive_source_mods(lid, exeroot, caseroot):
    # Save SourceMods
    sourcemods = os.path.join(caseroot, "SourceMods")
    sourcemods_prov = os.path.join(exeroot, "SourceMods.{}.tar.gz".format(lid))
    if os.path.exists(sourcemods_prov):
        os.remove(sourcemods_prov)
    if os.path.isdir(sourcemods):
        with tarfile.open(sourcemods_prov, "w:gz") as tfd:
            tfd.add(sourcemods, arcname="SourceMods")


def _archive_build_environment(lid, exeroot, case):
    # Save build env
    env_prov = os.path.join(exeroot, "build_environment.{}.txt".format(lid))
    if os.path.exists(env_prov):
        os.remove(env_prov)
    env_module = case.get_env("mach_specific")
    env_module.save_all_env_info(env_prov)


def _archive_build_times(lid, exeroot):
    # Save build times
    build_times = os.path.join(exeroot, "build_times.{}.txt".format(lid))
    if os.path.exists(build_times):
        os.remove(build_times)
    globstr = "{}/*bldlog*{}.gz".format(exeroot, lid)
    matches = glob.glob(globstr)
    if matches:
        _extract_times(matches, build_times)


def _extract_times(zipfiles, target_file):
    contents = "Target Build_time\n"
    total_build_time = 0.0
    for zipfile in zipfiles:
        stat, output, _ = utils.run_cmd("zgrep 'built in' {}".format(zipfile))
        if stat == 0:
            for line in output.splitlines():
                line = line.strip()
                if line:
                    items = line.split()
                    target, the_time = items[1], items[-2]
                    contents += "{} {}\n".format(target, the_time)

        stat, output, _ = utils.run_cmd("zgrep -E '^real [0-9.]+$' {}".format(zipfile))
        if stat == 0:
            for line in output.splitlines():
                line = line.strip()
                if line:
                    total_build_time += float(line.split()[-1])

    with open(target_file, "w") as fd:
        fd.write(contents)
        fd.write("Total_Elapsed_Time {}".format(str(total_build_time)))


def _symlink_current(lid, exeroot):
    # For all the just-created post-build provenance files, symlink a generic name
    # to them to indicate that these are the most recent or active.
    for item in [
        "GIT_DESCRIBE",
        "GIT_LOGS_HEAD",
        "GIT_SUBMODULE_STATUS",
        "GIT_STATUS",
        "GIT_DIFF",
        "GIT_LOG",
        "GIT_CONFIG",
        "GIT_REMOTE",
        "SourceMods",
        "build_environment",
        "build_times",
    ]:
        globstr = "{}/{}.{}*".format(exeroot, item, lid)
        matches = glob.glob(globstr)
        utils.expect(
            len(matches) < 2,
            "Multiple matches for glob {} should not have happened".format(globstr),
        )
        if matches:
            the_match = matches[0]
            generic_name = the_match.replace(".{}".format(lid), "")
            if os.path.exists(generic_name):
                os.remove(generic_name)
            os.symlink(the_match, generic_name)


def save_prerun_provenance(case, lid=None):
    with utils.SharedArea():
        # Always save env
        lid = os.environ["LID"] if lid is None else lid
        env_module = case.get_env("mach_specific")
        logdir = os.path.join(case.get_value("CASEROOT"), "logs")
        if not os.path.isdir(logdir):
            os.makedirs(logdir)
        env_module.save_all_env_info(
            os.path.join(logdir, "run_environment.txt.{}".format(lid))
        )

        run_dir = case.get_value("RUNDIR")

        base_preview_run = os.path.join(run_dir, "preview_run.log")
        preview_run = f"{base_preview_run}.{lid}"

        if os.path.exists(base_preview_run):
            os.remove(base_preview_run)

        with open(base_preview_run, "w") as fd:
            case.preview_run(lambda x: fd.write("{}\n".format(x)), None)

            # Create copy rather than symlink, the log is automatically gzipped
            utils.safe_copy(base_preview_run, preview_run)

        _cleanup_spio_stats(case)

        if case.get_value("SAVE_TIMING"):
            _record_timing(case, lid)


def _check_timing_dir(timing_dir, base_case, lid):
    if timing_dir is None or not os.path.isdir(timing_dir):
        raise RuntimeError(
            f"SAVE_TIMING_DIR {timing_dir} is not valid. E3SM requires a valid SAVE_TIMING_DIR to archive timing data."
        )

    logger.info(
        "Archiving timing data and associated provenance in {}.".format(timing_dir)
    )

    full_timing_dir = os.path.join(
        timing_dir, "performance_archive", getpass.getuser(), base_case, lid
    )

    if os.path.exists(full_timing_dir):
        raise RuntimeError(
            f"{full_timing_dir} already exists. Skipping archive of timing data and associated provenance."
        )

    try:
        os.makedirs(full_timing_dir)
    except OSError:
        raise RuntimeError(
            f"{full_timing_dir} cannot be created. Skipping archive of timing data and associated provenance."
        )

    return full_timing_dir


def _record_timing(case, lid):
    project = case.get_value("PROJECT", subgroup=case.get_primary_job())
    if not case.is_save_timing_dir_project(project):
        return

    rundir = case.get_value("RUNDIR")
    blddir = case.get_value("EXEROOT")
    caseroot = case.get_value("CASEROOT")
    srcroot = case.get_value("SRCROOT")
    base_case = case.get_value("CASE")
    timing_dir = case.get_value("SAVE_TIMING_DIR")

    try:
        full_timing_dir = _check_timing_dir(timing_dir, base_case, lid)
    except RuntimeError as e:
        logger.warning("{!s}", e)

        return

    mach = case.get_value("MACH")
    compiler = case.get_value("COMPILER")

    # For some batch machines save queue info
    job_id = _get_batch_job_id_for_syslog(case)

    if job_id is not None:
        _record_queue_info(mach, job_id, lid, full_timing_dir)

    # copy/tar SourceModes
    source_mods_dir = os.path.join(caseroot, "SourceMods")
    if os.path.isdir(source_mods_dir):
        with tarfile.open(
            os.path.join(full_timing_dir, "SourceMods.{}.tar.gz".format(lid)), "w:gz"
        ) as tfd:
            tfd.add(source_mods_dir, arcname="SourceMods")

    # Save various case configuration items
    case_docs = os.path.join(full_timing_dir, "CaseDocs.{}".format(lid))
    os.mkdir(case_docs)

    _copy_caseroot_files(mach, compiler, caseroot, case_docs, lid)
    _copy_blddir_files(blddir, full_timing_dir, lid)
    _copy_rundir_files(rundir, full_timing_dir, lid)

    # Save state of repo
    from_repo = (
        srcroot
        if os.path.exists(os.path.join(srcroot, ".git"))
        else os.path.dirname(srcroot)
    )

    desc = utils.get_current_commit(tag=True, repo=from_repo)
    with open(os.path.join(full_timing_dir, "GIT_DESCRIBE.{}".format(lid)), "w") as fd:
        fd.write(desc)

    # Collect syslog if enabled on machine. (Sarat)
    # Ref: https://github.com/ESMCI/cime/blob/655de638182ba9381a5d6607cdbade3b0a70a040/CIME/case/case.py#L1696
    # If machines_dir has a syslog.<machine_name> script, it is copied as caseroot/Tools/mach_syslog and run.
    # Otherwise, syslog.noop is used and no output is produced.
    if job_id is not None:
        _record_syslog(case, lid, job_id, caseroot, rundir, full_timing_dir)


def _record_queue_info(mach, job_id, lid, full_timing_dir):
    if mach in ["anvil", "chrysalis", "compy"]:
        _record_anl_queue(job_id, lid, full_timing_dir)
#TODO: Add Perlmutter
    elif mach in ["frontier", "crusher"]:
        _record_slurm_queue(job_id, lid, full_timing_dir)
    elif mach == "summit":
        _record_olcf_queue(job_id, lid, full_timing_dir)

# TODO: Switch to generic Slurm routine
def _record_nersc_queue(job_id, lid, full_timing_dir):
    for cmd, filename in [
        ("sinfo -a -l", "sinfol"),
        ("scontrol show jobid %s" % job_id, "sqsf_jobid"),
        # ("sqs -f", "sqsf"),
        (
            "squeue -o '%.10i %.15P %.20j %.10u %.7a %.2t %.6D %.8C %.10M %.10l %.20S %.20V'",
            "squeuef",
        ),
        ("squeue -t R -o '%.10i %R'", "squeues"),
    ]:
        filename = "%s.%s" % (filename, lid)
        utils.run_cmd_no_fail(cmd, arg_stdout=filename, from_dir=full_timing_dir)
        utils.gzip_existing_file(os.path.join(full_timing_dir, filename))

# Generic Slurm queue info (Frontier)
# TODO: Consolidate _record routines based on batch system if generalization is adequate
def _record_slurm_queue(job_id, lid, full_timing_dir):
    for cmd, filename in [
        ("sinfo -a -l", "sinfol"),
        ("scontrol show jobid %s" % job_id, "scontrol_jobid"),
        (
            "squeue -o '%.10i %.15P %.20j %.10u %.7a %.2t %.6D %.8C %.10M %.10l %.20S %.20V'",
            "squeuef",
        ),
        ("squeue -t R -o '%.10i %R'", "squeues"),
    ]:
        filename = "%s.%s" % (filename, lid)
        utils.run_cmd_no_fail(cmd, arg_stdout=filename, from_dir=full_timing_dir)
        utils.gzip_existing_file(os.path.join(full_timing_dir, filename))


# TODO: Switch to generic Slurm routine
def _record_anl_queue(job_id, lid, full_timing_dir):
    for cmd, filename in [
        ("sinfo -l", "sinfol"),
        ("squeue -o '%all' --job {}".format(job_id), "squeueall_jobid"),
        (
            "squeue -o '%.10i %.10P %.15u %.20a %.2t %.6D %.8C %.12M %.12l %.20S %.20V %j'",
            "squeuef",
        ),
        ("squeue -t R -o '%.10i %R'", "squeues"),
    ]:
        filename = "%s.%s" % (filename, lid)
        utils.run_cmd_no_fail(cmd, arg_stdout=filename, from_dir=full_timing_dir)
        utils.gzip_existing_file(os.path.join(full_timing_dir, filename))


def _record_olcf_queue(job_id, lid, full_timing_dir):
    for cmd, filename in [
        ("bjobs -u all >", "bjobsu_all"),
        ("bjobs -r -u all -o 'jobid slots exec_host' >", "bjobsru_allo"),
        ("bjobs -l -UF %s >" % job_id, "bjobslUF_jobid"),
    ]:
        full_cmd = cmd + " " + filename
        utils.run_cmd_no_fail(full_cmd + "." + lid, from_dir=full_timing_dir)
        utils.gzip_existing_file(os.path.join(full_timing_dir, filename + "." + lid))


def _cleanup_spio_stats(case):
    rundir = case.get_value("RUNDIR")
    for item in glob.glob(os.path.join(rundir, "io_perf_summary*")):
        os.remove(item)

    spio_stats_dir = os.path.join(rundir, "spio_stats")
    if os.path.exists(spio_stats_dir):
        shutil.rmtree(spio_stats_dir)

    try:
        os.makedirs(spio_stats_dir)
    except OSError:
        logger.warning(
            "{} could not be created. Scorpio I/O statistics will be stored in the run directory.".format(
                spio_stats_dir
            )
        )


def _copy_caseroot_files(mach, compiler, caseroot, case_docs, lid):
    globs_to_copy = [
        "CaseDocs/*",
        "run_script_provenance/*",
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
        "software_environment.txt",
    ]

    utils.copy_globs([os.path.join(caseroot, x) for x in globs_to_copy], case_docs, lid)


def _copy_blddir_files(blddir, full_timing_dir, lid):
    # Copy some items from build provenance
    blddir_globs_to_copy = [
        "GIT_LOGS_HEAD",
        "GIT_STATUS",
        "GIT_DIFF",
        "GIT_LOG",
        "GIT_REMOTE",
        "GIT_CONFIG",
        "GIT_SUBMODULE_STATUS",
        "build_environment.txt",
        "build_times.txt",
    ]

    utils.copy_globs(
        [os.path.join(blddir, x) for x in blddir_globs_to_copy], full_timing_dir, lid
    )


def _copy_rundir_files(rundir, full_timing_dir, lid):
    rundir_globs_to_copy = [
        "preview_run.log",
    ]

    utils.copy_globs(
        [os.path.join(rundir, x) for x in rundir_globs_to_copy], full_timing_dir, lid
    )


def _record_syslog(case, lid, job_id, caseroot, rundir, full_timing_dir):
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
        archive_checkpoints = os.path.join(
            full_timing_dir, "checkpoints.{}".format(lid)
        )
        os.mkdir(archive_checkpoints)
        utils.touch("{}/e3sm.log.{}".format(rundir, lid))
        syslog_jobid = utils.run_cmd_no_fail(
            "./mach_syslog {si} {jobid} {lid} {rundir} {rundir}/timing/checkpoints {ac} >& /dev/null & echo $!".format(
                si=sample_interval,
                jobid=job_id,
                lid=lid,
                rundir=rundir,
                ac=archive_checkpoints,
            ),
            from_dir=os.path.join(caseroot, "Tools"),
        )
        with open(os.path.join(rundir, "syslog_jobid.{}".format(job_id)), "w") as fd:
            fd.write("{}\n".format(syslog_jobid))


def save_postrun_provenance(case, lid=None):
    with utils.SharedArea():
        lid = os.environ["LID"] if lid is None else lid

        if case.get_value("SAVE_TIMING"):
            caseroot = case.get_value("CASEROOT")
            rundir = case.get_value("RUNDIR")

            rundir_timing_dir = _archive_timings(lid, rundir)

            _archive_atm_costs(lid, rundir)

            _archive_memory_profile(lid, rundir)

            _archive_spio_stats(lid, rundir)

            utils.gzip_existing_file(
                os.path.join(caseroot, "timing", "e3sm_timing_stats.%s" % lid)
            )

            # JGF: not sure why we do this
            timing_saved_file = "timing.%s.saved" % lid
            utils.touch(os.path.join(caseroot, "timing", timing_saved_file))

            project = case.get_value("PROJECT", subgroup=case.get_primary_job())
            if not case.is_save_timing_dir_project(project):
                return

            timing_dir = case.get_value("SAVE_TIMING_DIR")
            if timing_dir is None or not os.path.isdir(timing_dir):
                return

            mach = case.get_value("MACH")
            base_case = case.get_value("CASE")
            full_timing_dir = os.path.join(
                timing_dir, "performance_archive", getpass.getuser(), base_case, lid
            )

            if not os.path.isdir(full_timing_dir):
                return

            # Kill mach_syslog
            job_id = _get_batch_job_id_for_syslog(case)
            if job_id is not None:
                _kill_mach_syslog(job_id, rundir)

            # copy timings
            utils.safe_copy(
                "%s.tar.gz" % rundir_timing_dir, full_timing_dir, preserve_meta=False
            )

            _copy_performance_archive_files(
                case,
                lid,
                job_id,
                mach,
                rundir,
                caseroot,
                full_timing_dir,
                timing_saved_file,
            )

            # zip everything
            for root, _, files in os.walk(full_timing_dir):
                for filename in files:
                    if not filename.endswith(".gz"):
                        utils.gzip_existing_file(os.path.join(root, filename))


def _archive_timings(lid, rundir):
    # tar timings
    rundir_timing_dir = os.path.join(rundir, "timing." + lid)
    shutil.move(os.path.join(rundir, "timing"), rundir_timing_dir)
    with tarfile.open("%s.tar.gz" % rundir_timing_dir, "w:gz") as tfd:
        tfd.add(rundir_timing_dir, arcname=os.path.basename(rundir_timing_dir))

    shutil.rmtree(rundir_timing_dir)

    return rundir_timing_dir


def _archive_atm_costs(lid, rundir):
    atm_chunk_costs_src_path = os.path.join(rundir, "atm_chunk_costs.txt")
    if os.path.exists(atm_chunk_costs_src_path):
        atm_chunk_costs_dst_path = os.path.join(
            rundir, "atm_chunk_costs.{}".format(lid)
        )
        shutil.move(atm_chunk_costs_src_path, atm_chunk_costs_dst_path)
        utils.gzip_existing_file(atm_chunk_costs_dst_path)


def _archive_memory_profile(lid, rundir):
    # gzip memory profile log
    glob_to_copy = "memory.[0-4].*.log"
    for item in glob.glob(os.path.join(rundir, glob_to_copy)):
        mprof_dst_path = os.path.join(
            os.path.dirname(item), (os.path.basename(item) + ".{}").format(lid)
        )
        shutil.move(item, mprof_dst_path)
        utils.gzip_existing_file(mprof_dst_path)


def _archive_spio_stats(lid, rundir):
    # Copy Scorpio I/O performance stats in "spio_stats" to "spio_stats.[LID]" + tar + compress
    spio_stats_dir = os.path.join(rundir, "spio_stats")
    if not os.path.exists(spio_stats_dir):
        os.mkdir(spio_stats_dir)

    for item in glob.glob(os.path.join(rundir, "io_perf_summary*")):
        utils.safe_copy(item, spio_stats_dir)

    spio_stats_job_dir = os.path.join(rundir, "spio_stats." + lid)
    shutil.copytree(spio_stats_dir, spio_stats_job_dir)
    with tarfile.open("%s.tar.gz" % spio_stats_job_dir, "w:gz") as tfd:
        tfd.add(spio_stats_job_dir, arcname=os.path.basename(spio_stats_job_dir))

    shutil.rmtree(spio_stats_job_dir)


def _kill_mach_syslog(job_id, rundir):
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


def _copy_performance_archive_files(
    case, lid, job_id, mach, rundir, caseroot, full_timing_dir, timing_saved_file
):
    globs_to_copy = []
    if job_id is not None:
        if mach in ["anvil", "chrysalis", "compy"]:
            globs_to_copy.append("run*%s*%s" % (case.get_value("CASE"), job_id))
        elif mach == "summit":
            globs_to_copy.append("e3sm.stderr.%s" % job_id)
            globs_to_copy.append("e3sm.stdout.%s" % job_id)

    globs_to_copy.append("logs/run_environment.txt.{}".format(lid))
    globs_to_copy.append(os.path.join(rundir, "e3sm.log.{}.gz".format(lid)))
    globs_to_copy.append(os.path.join(rundir, "cpl.log.{}.gz".format(lid)))
    globs_to_copy.append(os.path.join(rundir, "atm_chunk_costs.{}.gz".format(lid)))
    globs_to_copy.append(os.path.join(rundir, "memory.[0-4].*.log.{}.gz".format(lid)))
    globs_to_copy.append("timing/*.{}*".format(lid))
    globs_to_copy.append("CaseStatus")
    globs_to_copy.append(os.path.join(rundir, "spio_stats.{}.tar.gz".format(lid)))
    globs_to_copy.append(os.path.join(caseroot, "replay.sh"))

    for glob_to_copy in globs_to_copy:
        for item in glob.glob(os.path.join(caseroot, glob_to_copy)):
            basename = os.path.basename(item)
            if basename != timing_saved_file:
                if lid not in basename and not basename.endswith(".gz"):
                    utils.safe_copy(
                        item,
                        os.path.join(full_timing_dir, "{}.{}".format(basename, lid)),
                        preserve_meta=False,
                    )
                else:
                    utils.safe_copy(item, full_timing_dir, preserve_meta=False)


def _get_batch_job_id_for_syslog(case):
    """
    mach_syslog only works on certain machines
    """
    mach = case.get_value("MACH")
    try:
        if mach in ["anvil", "chrysalis", "compy", "pm-cpu", "pm-gpu", "muller-cpu", "muller-gpu", "alvarez","frontier"]:
            # Note: Besides, SLURM_JOB_ID, equivalent SLURM_JOBID is also present on some systems (Frontier).
            return os.environ["SLURM_JOB_ID"]
        elif mach in ["summit"]:
            return os.environ["LSB_JOBID"]
    except KeyError:
        pass

    return None
