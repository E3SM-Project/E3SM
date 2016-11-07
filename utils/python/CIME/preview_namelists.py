"""
API for preview namelist
"""

from CIME.XML.standard_module_setup import *

import glob, shutil
logger = logging.getLogger(__name__)

def create_dirs(case):
    """
    Make necessary directories for case
    """
    # Get data from XML
    exeroot = case.get_value("EXEROOT")
    libroot = case.get_value("LIBROOT")
    incroot = case.get_value("INCROOT")
    rundir = case.get_value("RUNDIR")
    caseroot = case.get_value("CASEROOT")

    docdir = os.path.join(caseroot, "CaseDocs")
    dirs_to_make = []
    models = case.get_values("COMP_CLASSES")
    for model in models:
        dirname = "cpl" if model == "DRV" else model.lower()
        dirs_to_make.append(os.path.join(exeroot, dirname, "obj"))

    dirs_to_make.extend([exeroot, libroot, incroot, rundir, docdir])

    for dir_to_make in dirs_to_make:
        if (not os.path.isdir(dir_to_make)):
            try:
                logger.debug("Making dir '%s'" % dir_to_make)
                os.makedirs(dir_to_make)
            except OSError as e:
                expect(False, "Could not make directory '%s', error: %s" % (dir_to_make, e))

    # As a convenience write the location of the case directory in the bld and run directories
    for dir_ in (exeroot, rundir):
        with open(os.path.join(dir_,"CASEROOT"),"w+") as fd:
            fd.write(caseroot+"\n")

def create_namelists(case):

    case.flush()

    casebuild = case.get_value("CASEBUILD")
    caseroot = case.get_value("CASEROOT")
    rundir = case.get_value("RUNDIR")

    docdir = os.path.join(caseroot, "CaseDocs")

    # Load modules
    case.load_env()

    # Create namelists
    models = case.get_values("COMP_CLASSES")
    for model in models:
        model_str = model.lower()
        config_file = case.get_value("CONFIG_%s_FILE" % model_str.upper())
        config_dir = os.path.dirname(config_file)
        cmd = os.path.join(config_dir, "buildnml")
        run_cmd_no_fail("%s %s" % (cmd, caseroot), verbose=True)

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
