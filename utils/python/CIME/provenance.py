#!/usr/bin/env python

"""
Library for saving build/run provenance.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import touch, gzip_existing_file

import tarfile, getpass, signal, glob, shutil

logger = logging.getLogger(__name__)

def _get_batch_job_id(case):
    mach = case.get_value("MACH")
    if mach == 'titan':
        return os.environ("PBS_JOBID")
    elif mach in ['edison', 'corip1']:
        return os.environ("SLURM_JOB_ID")
    elif mach == 'mira':
        return os.environ("COBALT_JOBID")
    else:
        return None

def save_build_provenance_acme(case, lid=None): # pylint: disable=unused-argument
    pass

def save_build_provenance_cesm(case, lid=None): # pylint: disable=unused-argument
    pass

def save_build_provenance(case, lid=None):
    model = case.get_value("MODEL")
    if model == "acme":
        save_build_provenance_acme(case, lid=lid)
    elif model == "cesm":
        save_build_provenance_cesm(case, lid=lid)

def save_prerun_provenance_acme(case, lid=None):
    if not case.get_value("SAVE_TIMING") or case.get_value("MODEL") != "acme":
        return

    lid = os.environ["LID"] if lid is None else lid
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

def save_prerun_provenance_cesm(case, lid=None): # pylint: disable=unused-argument
    pass

def save_prerun_provenance(case, lid=None):
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
    base_case = case.get_value("CASEBASEID")
    full_timing_dir = os.path.join(timing_dir, "performance_archive", getpass.getuser(), base_case, lid)

    if not os.path.isdir(full_timing_dir):
        os.makedirs(full_timing_dir)

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

def save_postrun_provenance(case, lid=None):
    model = case.get_value("MODEL")
    if model == "acme":
        save_postrun_provenance_acme(case, lid=lid)
    elif model == "cesm":
        save_postrun_provenance_cesm(case, lid=lid)

