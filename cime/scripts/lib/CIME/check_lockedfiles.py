"""
API for checking locked files
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.env_build import EnvBuild
from CIME.XML.env_case import EnvCase
from CIME.XML.env_mach_pes import EnvMachPes
from CIME.XML.env_batch import EnvBatch
from CIME.utils import run_cmd_no_fail

logger = logging.getLogger(__name__)

import glob, shutil

LOCKED_DIR = "LockedFiles"

def lock_file(filename, caseroot=None, newname=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    newname = filename if newname is None else newname
    fulllockdir = os.path.join(caseroot, LOCKED_DIR)
    if not os.path.exists(fulllockdir):
        os.mkdir(fulllockdir)

    logging.debug("Locking file {}".format(filename))

    # JGF: It is extremely dangerous to alter our database (xml files) without
    # going through the standard API. The copy below invalidates all existing
    # GenericXML instances that represent this file and all caching that may
    # have involved this file. We should probably seek a safer way of locking
    # files.
    shutil.copyfile(os.path.join(caseroot, filename), os.path.join(fulllockdir, newname))

def unlock_file(filename, caseroot=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    locked_path = os.path.join(caseroot, LOCKED_DIR, filename)
    if os.path.exists(locked_path):
        os.remove(locked_path)

    logging.debug("Unlocking file {}".format(filename))

def is_locked(filename, caseroot=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    return os.path.exists(os.path.join(caseroot, LOCKED_DIR, filename))

def check_pelayouts_require_rebuild(case, models):
    """
    Create if we require a rebuild, expects cwd is caseroot
    """
    locked_pes = os.path.join(LOCKED_DIR, "env_mach_pes.xml")
    if os.path.exists(locked_pes):
        # Look to see if $comp_PE_CHANGE_REQUIRES_REBUILD is defined
        # for any component
        env_mach_pes_locked = EnvMachPes(infile=locked_pes, components=case.get_values("COMP_CLASSES"))
        for comp in models:
            if case.get_value("{}_PE_CHANGE_REQUIRES_REBUILD".format(comp)):
                # Changing these values in env_mach_pes.xml will force
                # you to clean the corresponding component
                old_tasks   = env_mach_pes_locked.get_value("NTASKS_{}".format(comp))
                old_threads = env_mach_pes_locked.get_value("NTHRDS_{}".format(comp))
                old_inst    = env_mach_pes_locked.get_value("NINST_{}".format(comp))

                new_tasks   = case.get_value("NTASKS_{}".format(comp))
                new_threads = case.get_value("NTHRDS_{}".format(comp))
                new_inst    = case.get_value("NINST_{}".format(comp))

                if old_tasks != new_tasks or old_threads != new_threads or old_inst != new_inst:
                    logging.warning("{} pe change requires clean build {} {}".format(comp, old_tasks, new_tasks))
                    cleanflag = comp.lower()
                    run_cmd_no_fail("./case.build --clean {}".format(cleanflag))

        unlock_file("env_mach_pes.xml", case.get_value("CASEROOT"))

def check_lockedfile(case, filebase):
    caseroot = case.get_value("CASEROOT")

    cfile = os.path.join(caseroot, filebase)
    lfile = os.path.join(caseroot, "LockedFiles", filebase)
    components = case.get_values("COMP_CLASSES")
    if os.path.isfile(cfile):
        objname = filebase.split('.')[0]
        if objname == "env_build":
            f1obj = case.get_env('build')
            f2obj = EnvBuild(caseroot, lfile)
        elif objname == "env_mach_pes":
            f1obj = case.get_env('mach_pes')
            f2obj = EnvMachPes(caseroot, lfile, components=components)
        elif objname == "env_case":
            f1obj = case.get_env('case')
            f2obj = EnvCase(caseroot, lfile)
        elif objname == "env_batch":
            f1obj = case.get_env('batch')
            f2obj = EnvBatch(caseroot, lfile)
        else:
            logging.warning("Locked XML file '{}' is not current being handled".format(filebase))
            return

        diffs = f1obj.compare_xml(f2obj)
        if diffs:

            logging.warning("File {} has been modified".format(lfile))
            toggle_build_status = False
            for key in diffs.keys():
                if key != "BUILD_COMPLETE":
                    print("  found difference in {} : case {} locked {}"
                          .format(key, repr(diffs[key][0]), repr(diffs[key][1])))
                    toggle_build_status = True
            if objname == "env_mach_pes":
                expect(False, "Invoke case.setup --reset ")
            elif objname == "env_case":
                expect(False, "Cannot change file env_case.xml, please"
                       " recover the original copy from LockedFiles")
            elif objname == "env_build":
                if toggle_build_status:
                    logging.warning("Setting build complete to False")
                    case.set_value("BUILD_COMPLETE", False)
                    if "PIO_VERSION" in diffs:
                        case.set_value("BUILD_STATUS", 2)
                        logging.critical("Changing PIO_VERSION requires running "
                                         "case.build --clean-all and rebuilding")
                    else:
                        case.set_value("BUILD_STATUS", 1)
                    f2obj.set_value("BUILD_COMPLETE", False)
            elif objname == "env_batch":
                expect(False, "Batch configuration has changed, please run case.setup --reset")
            else:
                expect(False, "'{}' diff was not handled".format(objname))

def check_lockedfiles(case):
    """
    Check that all lockedfiles match what's in case

    If caseroot is not specified, it is set to the current working directory
    """
    caseroot = case.get_value("CASEROOT")
    lockedfiles = glob.glob(os.path.join(caseroot, "LockedFiles", "*.xml"))
    for lfile in lockedfiles:
        fpart = os.path.basename(lfile)
        # ignore files used for tests such as env_mach_pes.ERP1.xml by looking for extra dots in the name
        if fpart.count('.') > 1:
            continue

        check_lockedfile(case, fpart)
