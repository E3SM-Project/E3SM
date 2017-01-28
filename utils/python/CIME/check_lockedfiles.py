"""
API for checking locked files
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.env_build import EnvBuild
from CIME.XML.env_case import EnvCase
from CIME.XML.env_mach_pes import EnvMachPes
from CIME.XML.env_batch import EnvBatch
from CIME.utils import run_cmd_no_fail

import glob, shutil

LOCKED_DIR = "LockedFiles"

def lock_file(filename, caseroot=None, newname=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    newname = filename if newname is None else newname
    fulllockdir = os.path.join(caseroot, LOCKED_DIR)
    if not os.path.exists(fulllockdir):
        os.mkdir(fulllockdir)
    logging.debug("Locking file %s to %s"%(filename, newname))
    shutil.copyfile(os.path.join(caseroot, filename), os.path.join(fulllockdir, newname))

def unlock_file(filename, caseroot=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    locked_path = os.path.join(caseroot, LOCKED_DIR, filename)
    if os.path.exists(locked_path):
        os.remove(locked_path)

def is_locked(filename, caseroot=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    return os.path.exists(os.path.join(caseroot, LOCKED_DIR, filename))

def restore(filename, caseroot=None, newname=None):
    """
    Restore the locked version of filename into main case dir
    """
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    newname = filename if newname is None else newname
    shutil.copyfile(os.path.join(caseroot, LOCKED_DIR, filename), os.path.join(caseroot, newname))
    # relock the restored file if names diffs
    if newname != filename:
        lock_file(filename, caseroot)

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
            if case.get_value("%s_PE_CHANGE_REQUIRES_REBUILD" % comp):
                # Changing these values in env_mach_pes.xml will force
                # you to clean the corresponding component
                old_tasks   = env_mach_pes_locked.get_value("NTASKS_%s" % comp)
                old_threads = env_mach_pes_locked.get_value("NTHRDS_%s" % comp)
                old_inst    = env_mach_pes_locked.get_value("NINST_%s" % comp)

                new_tasks   = case.get_value("NTASKS_%s" % comp)
                new_threads = case.get_value("NTHRDS_%s" % comp)
                new_inst    = case.get_value("NINST_%s" % comp)

                if old_tasks != new_tasks or old_threads != new_threads or old_inst != new_inst:
                    logging.warn("%s pe change requires clean build %s %s" % (comp, old_tasks, new_tasks))
                    cleanflag = comp.lower()
                    run_cmd_no_fail("./case.build --clean %s" % cleanflag)

        unlock_file("env_mach_pes.xml", case.get_value("CASEROOT"))

def check_lockedfiles(caseroot=None):
    """
    Check that all lockedfiles match what's in case

    If caseroot is not specified, it is set to the current working directory
    """
    caseroot = os.getcwd() if caseroot is None else caseroot
    lockedfiles = glob.glob(os.path.join(caseroot, "LockedFiles", "[^\.]*.xml"))
    for lfile in lockedfiles:
        fpart = os.path.basename(lfile)
        cfile = os.path.join(caseroot, fpart)
        if os.path.isfile(cfile):
            objname = fpart.split('.')[0]
            logging.info("Checking file %s"%objname)
            if objname == "env_build":
                f1obj = EnvBuild(caseroot, cfile)
                f2obj = EnvBuild(caseroot, lfile)
            elif objname == "env_mach_pes":
                f1obj = EnvMachPes(caseroot, cfile)
                f2obj = EnvMachPes(caseroot, lfile)
            elif objname == "env_case":
                f1obj = EnvCase(caseroot, cfile)
                f2obj = EnvCase(caseroot, lfile)
            elif objname == "env_batch":
                f1obj = EnvBatch(caseroot, cfile)
                f2obj = EnvBatch(caseroot, lfile)
            else:
                logging.warn("Locked XML file '%s' is not current being handled" % fpart)
                continue

            diffs = f1obj.compare_xml(f2obj)
            if diffs:
                logging.warn("File %s has been modified"%lfile)
                for key in diffs.keys():
                    print("  found difference in %s : case %s locked %s" %
                          (key, repr(diffs[key][0]), repr(diffs[key][1])))

                if objname == "env_mach_pes":
                    expect(False, "Invoke case.setup --reset ")
                elif objname == "env_case":
                    expect(False, "Cannot change file env_case.xml, please"
                           " recover the original copy from LockedFiles")
                elif objname == "env_build":
                    logging.warn("Setting build complete to False")
                    f1obj.set_value("BUILD_COMPLETE", False)
                    if "PIO_VERSION" in diffs.keys():
                        f1obj.set_value("BUILD_STATUS", 2)
                        f1obj.write()
                        logging.critical("Changing PIO_VERSION requires running "
                                         "case.build --clean-all and rebuilding")
                    else:
                        f1obj.set_value("BUILD_STATUS", 1)
                        f1obj.write()
                elif objname == "env_batch":
                    expect(False, "Batch configuration has changed, please run case.setup --reset")
                else:
                    expect(False, "'%s' diff was not handled" % objname)
