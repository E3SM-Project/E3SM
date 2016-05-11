"""
user_mod_support.py
"""

from CIME.XML.standard_module_setup import *
from CIME.utils             import expect, run_cmd
import shutil, glob, tempfile

logger = logging.getLogger(__name__)

def apply_user_mods(caseroot, user_mods_path, ninst={}):
    '''
    Recursivlely apply user_mods to caseroot
    '''
    include_dirs = build_include_dirs_list(user_mods_path)
    for include_dir in include_dirs:
        for user_nl in glob.iglob(os.path.join(include_dir,"user_nl_*")):
            with open(os.path.join(include_dir, user_nl), "r") as fd:
                contents = fd.read()
            case_user_nl = user_nl.replace(include_dir, caseroot)
            comp = case_user_nl.split('_')[-1]
            if comp in ninst.keys():
                for comp_inst in xrange(1,ninst[comp]):
                    case_user_nl_inst = case_user_nl + "_%4.4d"%comp_inst
                    logger.info("Appending file %s"%case_user_nl_inst)
                    with open(case_user_nl_inst, "a") as fd:
                        fd.write(contents)
            else:
                logger.info("Appending file %s"%case_user_nl)
                with open(case_user_nl, "a") as fd:
                    fd.write(contents)
        for root, dirs, files in os.walk(include_dir,followlinks=True,topdown=False):
            if "src" in os.path.basename(root):
                for sfile in files:
                    source_mods = os.path.join(root,sfile)
                    case_source_mods = source_mods.replace(include_dir, caseroot)
                    if os.path.isfile(case_source_mods):
                        logger.warn("Refusing to overwrite existing SourceMods in %s"%case_source_mods)
                    else:
                        logger.info("Adding SourceMod to case %s"%case_source_mods)
                        try:
                            shutil.copyfile(source_mods, case_source_mods)
                        except:
                            expect(False, "Could not write file %s in caseroot %s"
                                   %(case_source_mods,caseroot))
        for shell_commands_file in glob.iglob(os.path.join(include_dir,"shell_commands")):
            case_shell_commands = shell_commands_file.replace(include_dir, caseroot)
            with open(shell_commands_file,"r") as fd:
                new_shell_commands = fd.read()
            with open(case_shell_commands, "a") as fd:
                fd.write(new_shell_commands)


def build_include_dirs_list(user_mods_path, include_dirs=[]):
    '''
    If user_mods_path has a file "include_user_mods" read that
    file and add directories to the include_dirs, recursively check
    each of those directories for further directories.
    The file may also include comments deleneated with # in the first column
    '''
    expect(os.path.isabs(user_mods_path),
           "Expected full directory path, got '%s'"%user_mods_path)
    logger.info("Adding user mods directory %s"%user_mods_path)
    include_dirs.append(os.path.normpath(user_mods_path))
    include_file = os.path.join(include_dirs[-1],"include_user_mods")
    if os.path.isfile(include_file):
        with open(include_file, "r") as fd:
            for newpath in fd:
                newpath = newpath.rstrip()
                if not newpath.startswith("#"):
                    if not os.path.isabs(newpath):
                        newpath = os.path.join(user_mods_path, newpath)
                    if os.path.isabs(newpath):
                        build_include_dirs_list(newpath, include_dirs)
                    else:
                        logger.warn("Could not resolve path '%s' in file '%s'"%newpath,include_file)
    return include_dirs
