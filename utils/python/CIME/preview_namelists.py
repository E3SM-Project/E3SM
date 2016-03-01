"""
API for preview namelist
"""

from XML.standard_module_setup import *
from CIME.utils import expect, run_cmd
from CIME.case import Case

import glob, shutil

def preview_namelists(dryrun=False, case=None, casedir=None):
    if (case is None):
        if (casedir is None):
            case = Case()
        else:
            case = Case(case_root=casedir)

    # Get data from XML
    exeroot = case.get_value("EXEROOT")
    libroot = case.get_value("LIBROOT")
    incroot = case.get_value("INCROOT")
    rundir = case.get_value("RUNDIR")
    sharedlibroot = case.get_value("SHAREDLIBROOT")
    caseroot = case.get_value("CASEROOT")
    casebuild = case.get_value("CASEBUILD")
    testcase = case.get_value("TESTCASE")

    dryrun = True if (testcase == "SBN") else dryrun

    models = ["atm", "lnd", "ice", "ocn", "glc", "wav", "rof", "cpl"]
    docdir = os.path.join(caseroot, "CaseDocs")

    if (dryrun):
        # Only create rundir
        try:
            os.makedirs(rundir)
        except OSError:
            logging.warning("Not able to create $RUNDIR, trying a subdirectory of $CASEROOT")
            rundir = os.path.join(caseroot, rundir)
            try:
                os.makedirs(rundir)
                logging.info("Success! Setting RUNDIR=%s" % rundir)
                case.set_value("RUNDIR", rundir)
            except OSError:
                expect(False, "Could not create rundir")

    else:
        # Make necessary directories
        dirs_to_make = [os.path.join(exeroot, model, "obj") for model in models]
        dirs_to_make.extend([exeroot, libroot, incroot, rundir, sharedlibroot, docdir])

        for dir_to_make in dirs_to_make:
            if (not os.path.isdir(dir_to_make)):
                try:
                    os.makedirs(dir_to_make)
                except OSError as e:
                    expect(False, "Could not make directory '%s', error: %s" % (dir_to_make, e))

    # Create namelists
    for model in models:
        model_str = "drv" if model == "cpl" else model
        config_file = case.get_value("CONFIG_%s_FILE" % model_str.upper())
        config_dir = os.path.dirname(config_file)
        cmd = os.path.join(config_dir, "buildnml")

        if (logging.getLogger().level == logging.DEBUG):
            run_cmd("PREVIEW_NML=1 %s %s" % (cmd, caseroot))
        else:
            run_cmd("%s %s" % (cmd, caseroot))

    # Save namelists to docdir
    if (not os.path.isdir(docdir)):
        try:
            with open(os.path.join(docdir, "README"), "w") as fd:
                fd.write(" CESM Resolved Namelist Files\n   For documentation only DO NOT MODIFY\n")
        except (OSError, IOError) as e:
            expect(False, "Failed to write %s/README: %s" % (docdir, e))


    for cpglob in ["*_in_[0-9]*", "*modelio*", "*_in",
                   "*streams*txt*", "*stxt", "*maps.rc", "*cism.config*"]:
        for file_to_copy in glob.glob(os.path.join(rundir, cpglob)):
            shutil.copy2(file_to_copy, docdir)

    # Copy over chemistry mechanism docs if they exist
    if (os.path.isdir(os.path.join(casebuild, "camconf"))):
        for file_to_copy in glob.glob(os.path.join(casebuild, "camconf", "*chem_mech*")):
            shutil.copy2(file_to_copy, docdir)
