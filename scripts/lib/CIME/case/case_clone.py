"""
create_clone is a member of the Case class from file case.py
"""
import os, glob, shutil
from CIME.XML.standard_module_setup import *
from CIME.utils import expect, check_name, safe_copy
from CIME.user_mod_support import apply_user_mods
from CIME.locked_files         import lock_file
from CIME.simple_compare            import compare_files

logger = logging.getLogger(__name__)


def create_clone(self, newcase, keepexe=False, mach_dir=None, project=None,
                      cime_output_root=None, exeroot=None, rundir=None,
                      user_mods_dir=None):
    """
    Create a case clone

    If exeroot or rundir are provided (not None), sets these directories
    to the given paths; if not provided, uses default values for these
    directories. It is an error to provide exeroot if keepexe is True.
    """
    if cime_output_root is None:
        cime_output_root = self.get_value("CIME_OUTPUT_ROOT")

    newcaseroot = os.path.abspath(newcase)
    expect(not os.path.isdir(newcaseroot),
           "New caseroot directory {} already exists".format(newcaseroot))
    newcasename = os.path.basename(newcaseroot)
    expect(check_name(newcasename),
           "New case name invalid {} ".format(newcasename))
    newcase_cimeroot = os.path.abspath(get_cime_root())

    # create clone from case to case
    clone_cimeroot = self.get_value("CIMEROOT")
    if newcase_cimeroot != clone_cimeroot:
        logger.warning(" case  CIMEROOT is {} ".format(newcase_cimeroot))
        logger.warning(" clone CIMEROOT is {} ".format(clone_cimeroot))
        logger.warning(" It is NOT recommended to clone cases from different versions of CIME.")

    # *** create case object as deepcopy of clone object ***
    srcroot = os.path.join(newcase_cimeroot,"..")
    newcase = self.copy(newcasename, newcaseroot, newsrcroot=srcroot)
    newcase.set_value("CIMEROOT", newcase_cimeroot)

    # if we are cloning to a different user modify the output directory
    olduser = self.get_value("USER")
    newuser = os.environ.get("USER")
    if olduser != newuser:
        cime_output_root = cime_output_root.replace(olduser, newuser)
        newcase.set_value("USER", newuser)
    newcase.set_value("CIME_OUTPUT_ROOT", cime_output_root)

    # try to make the new output directory and raise an exception
    # on any error other than directory already exists.
    if os.path.isdir(cime_output_root):
        expect(os.access(cime_output_root, os.W_OK), "Directory {} is not writable "
               "by this user.  Use the --cime-output-root flag to provide a writable "
               "scratch directory".format(cime_output_root))
    else:
        try:
            os.makedirs(cime_output_root)
        except:
            if not os.path.isdir(cime_output_root):
                raise

    # determine if will use clone executable or not
    if keepexe:
        orig_exeroot = self.get_value("EXEROOT")
        newcase.set_value("EXEROOT", orig_exeroot)
        newcase.set_value("BUILD_COMPLETE","TRUE")
        orig_bld_complete = self.get_value("BUILD_COMPLETE")
        if not orig_bld_complete:
            logger.warning("\nWARNING: Creating a clone with --keepexe before building the original case may cause PIO_TYPENAME to be invalid in the clone")
            logger.warning("Avoid this message by building case one before you clone.\n")
    else:
        newcase.set_value("BUILD_COMPLETE","FALSE")

    # set machdir
    if mach_dir is not None:
        newcase.set_value("MACHDIR", mach_dir)

    # set exeroot and rundir if requested
    if exeroot is not None:
        expect(not keepexe, "create_case_clone: if keepexe is True, "
               "then exeroot cannot be set")
        newcase.set_value("EXEROOT", exeroot)
    if rundir is not None:
        newcase.set_value("RUNDIR", rundir)

    # Set project id
    # Note: we do not just copy this from the clone because it seems likely that
    # users will want to change this sometimes, especially when cloning another
    # user's case. However, note that, if a project is not given, the fallback will
    # be to copy it from the clone, just like other xml variables are copied.
    if project is None:
        project = self.get_value("PROJECT", subgroup=self.get_primary_job())
    if project is not None:
        newcase.set_value("PROJECT", project)

    # create caseroot
    newcase.create_caseroot(clone=True)

    # Many files in the case will be links back to the source tree
    # but users may have broken links to modify files locally.  In this case
    # copy the locally modified file.   We only want to do this for files that
    # already exist in the clone.
    #pylint: disable=protected-access
    self._copy_user_modified_to_clone(self.get_value("CASEROOT"), newcase.get_value("CASEROOT"))
    self._copy_user_modified_to_clone(self.get_value("CASETOOLS"), newcase.get_value("CASETOOLS"))

    newcase.flush(flushall=True)

    # copy user_ files
    cloneroot = self.get_case_root()
    files = glob.glob(cloneroot + '/user_*')

    for item in files:
        safe_copy(item, newcaseroot)

    # copy SourceMod and Buildconf files
    # if symlinks exist, copy rather than follow links
    for casesub in ("SourceMods", "Buildconf"):
        shutil.copytree(os.path.join(cloneroot, casesub),
                        os.path.join(newcaseroot, casesub),
                        symlinks=True)

    # lock env_case.xml in new case
    lock_file("env_case.xml", newcaseroot)

    # apply user_mods if appropriate
    newcase_root = newcase.get_value("CASEROOT")
    if user_mods_dir is not None:
        if keepexe:
            # If keepexe CANNOT change any env_build.xml variables - so make a temporary copy of
            # env_build.xml and verify that it has not been modified
            safe_copy(os.path.join(newcaseroot, "env_build.xml"),
                      os.path.join(newcaseroot, "LockedFiles", "env_build.xml"))

        # Now apply contents of user_mods directory
        apply_user_mods(newcase_root, user_mods_dir, keepexe=keepexe)

        # Determine if env_build.xml has changed
        if keepexe:
            success, comment = compare_files(os.path.join(newcaseroot, "env_build.xml"),
                                             os.path.join(newcaseroot, "LockedFiles", "env_build.xml"))
            if not success:
                logger.warning(comment)
                shutil.rmtree(newcase_root)
                expect(False, "env_build.xml cannot be changed via usermods if keepexe is an option: \n "
                           "Failed to clone case, removed {}\n".format(newcase_root))

    # if keep executable, then remove the new case SourceMods directory and link SourceMods to
    # the clone directory
    if keepexe:
        shutil.rmtree(os.path.join(newcase_root, "SourceMods"))
        os.symlink(os.path.join(cloneroot, "SourceMods"),
                   os.path.join(newcase_root, "SourceMods"))

    # Update README.case
    fclone   = open(cloneroot + "/README.case", "r")
    fnewcase = open(newcaseroot  + "/README.case", "a")
    fnewcase.write("\n    *** original clone README follows ****")
    fnewcase.write("\n " +  fclone.read())

    clonename = self.get_value("CASE")
    logger.info(" Successfully created new case {} from clone case {} ".format(newcasename, clonename))

    newcase.case_setup()

    return newcase
# pylint: disable=unused-argument
def _copy_user_modified_to_clone(self, origpath, newpath):
    """
    If file_ exists and is a link in newpath, and exists but is not a
    link in origpath, copy origpath file to newpath
    """
    for file_ in os.listdir(newpath):
        if (os.path.islink(os.path.join(newpath, file_)) and
            os.path.isfile(os.path.join(origpath, file_)) and
            not os.path.islink(os.path.join(origpath, file_))):
            logger.info("Copying user modified file {} to clone".format(file_))
            os.unlink(os.path.join(newpath, file_))
            safe_copy(os.path.join(origpath, file_), newpath)
