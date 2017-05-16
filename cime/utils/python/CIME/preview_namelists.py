"""
API for preview namelist
"""

from CIME.XML.standard_module_setup import *

import glob, shutil
logger = logging.getLogger(__name__)

def preview_namelists(case, dryrun=False):
    # refresh case xml files from object
    case.flush()

    # Get data from XML
    exeroot = case.get_value("EXEROOT")
    libroot = case.get_value("LIBROOT")
    incroot = case.get_value("INCROOT")
    rundir = case.get_value("RUNDIR")
    caseroot = case.get_value("CASEROOT")
    casebuild = case.get_value("CASEBUILD")
    testcase = case.get_value("TESTCASE")

    logger.debug("LID is: '%s'" % os.getenv("LID", ""))
    logger.debug("caseroot is: '%s'" % caseroot)

    dryrun = True if (testcase == "SBN") else dryrun

    models = ["atm", "lnd", "ice", "ocn", "glc", "wav", "rof", "cpl"]
    docdir = os.path.join(caseroot, "CaseDocs")

    if (dryrun):
        # Only create rundir
        try:
            os.makedirs(rundir)
        except OSError:
            logger.warning("Not able to create $RUNDIR, trying a subdirectory of $CASEROOT")
            rundir = os.path.join(caseroot, rundir)
            try:
                os.makedirs(rundir)
                logger.info("Success! Setting RUNDIR=%s" % rundir)
                case.set_value("RUNDIR", rundir)
            except OSError:
                expect(False, "Could not create rundir")

    else:

        # Load modules
        env_module = case.get_env("mach_specific")
        env_module.load_env_for_case(compiler=case.get_value("COMPILER"),
                                     debug=case.get_value("DEBUG"),
                                     mpilib=case.get_value("MPILIB"))

        # Make necessary directories
        dirs_to_make = [os.path.join(exeroot, model, "obj") for model in models]
        dirs_to_make.extend([exeroot, libroot, incroot, rundir, docdir])

        for dir_to_make in dirs_to_make:
            if (not os.path.isdir(dir_to_make)):
                try:
                    logger.debug("Making dir '%s'" % dir_to_make)
                    os.makedirs(dir_to_make)
                except OSError as e:
                    expect(False, "Could not make directory '%s', error: %s" % (dir_to_make, e))

    # Create namelists
    for model in models:
        model_str = "drv" if model == "cpl" else model
        config_file = case.get_value("CONFIG_%s_FILE" % model_str.upper())
        config_dir = os.path.dirname(config_file)
        cmd = os.path.join(config_dir, "buildnml")
        logger.info("Running %s:"%cmd)
        if (logger.level == logging.DEBUG):
            rc, out, err = run_cmd("PREVIEW_NML=1 %s %s" % (cmd, caseroot))
            expect(rc==0,"Command %s failed rc=%d\nout=%s\nerr=%s"%(cmd,rc,out,err))
        else:
            rc, out, err = run_cmd("%s %s" % (cmd, caseroot))
            expect(rc==0,"Command %s failed rc=%d\nout=%s\nerr=%s"%(cmd,rc,out,err))
        if out is not None:
            logger.info("     %s"%out)
        if err is not None:
            logger.info("     %s"%err)
    # refresh case xml object from file
    case.read_xml()
    # Save namelists to docdir
    if (not os.path.isdir(docdir)):
        os.makedirs(docdir)
        try:
            with open(os.path.join(docdir, "README"), "w") as fd:
                fd.write(" CESM Resolved Namelist Files\n   For documentation only DO NOT MODIFY\n")
        except (OSError, IOError) as e:
            expect(False, "Failed to write %s/README: %s" % (docdir, e))


    for cpglob in ["*_in_[0-9]*", "*modelio*", "*_in",
                   "*streams*txt*", "*stxt", "*maps.rc", "*cism.config*"]:
        for file_to_copy in glob.glob(os.path.join(rundir, cpglob)):
            logger.debug("Copy file from '%s' to '%s'" % (file_to_copy, docdir))
            shutil.copy2(file_to_copy, docdir)

    # Copy over chemistry mechanism docs if they exist
    if (os.path.isdir(os.path.join(casebuild, "camconf"))):
        for file_to_copy in glob.glob(os.path.join(casebuild, "camconf", "*chem_mech*")):
            shutil.copy2(file_to_copy, docdir)
