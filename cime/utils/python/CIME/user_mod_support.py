"""
user_mod_support.py
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, run_cmd_no_fail
import shutil, glob

logger = logging.getLogger(__name__)

def apply_user_mods(caseroot, user_mods_path, ninst=None):
    '''
    Recursivlely apply user_mods to caseroot - this includes updating user_nl_xxx,
    updating SourceMods and creating case shell_commands and xmlchange_cmds files

    First remove case shell_commands files if any already exist
    '''
    case_shell_command_files = [os.path.join(caseroot,"shell_commands"),
                           os.path.join(caseroot,"xmlchange_cmnds")]
    for shell_command_file in case_shell_command_files:
        if os.path.isfile(shell_command_file):
            os.remove(shell_command_file)

    include_dirs = build_include_dirs_list(user_mods_path)
    for include_dir in include_dirs:
        # write user_nl_xxx file in caseroot
        for user_nl in glob.iglob(os.path.join(include_dir,"user_nl_*")):
            with open(os.path.join(include_dir, user_nl), "r") as fd:
                newcontents = fd.read()
            if len(newcontents) == 0:
                continue
            case_user_nl = user_nl.replace(include_dir, caseroot)
            comp = case_user_nl.split('_')[-1]
            if ninst is not None and comp in ninst.keys() and ninst[comp] > 1:
                for comp_inst in xrange(1, ninst[comp]+1):
                    contents = newcontents
                    case_user_nl_inst = case_user_nl + "_%4.4d"%comp_inst
                    logger.info("Pre-pending file %s"%case_user_nl_inst)
                    if os.path.isfile(case_user_nl_inst):
                        with open(case_user_nl_inst, "r") as fd:
                            old_contents = fd.read()
                            if old_contents.find(contents) == -1:
                                contents = contents + old_contents
                    with open(case_user_nl_inst, "w") as fd:
                        fd.write(contents)
            else:
                contents = newcontents
                logger.info("Pre-pending file %s"%case_user_nl)
                if os.path.isfile(case_user_nl):
                    with open(case_user_nl, "r") as fd:
                        old_contents = fd.read()
                        if old_contents.find(contents) == -1:
                            contents = contents + old_contents
                with open(case_user_nl, "w") as fd:
                    fd.write(contents)

        # update SourceMods in caseroot
        for root, _, files in os.walk(include_dir,followlinks=True,topdown=False):
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

        # create xmlchange_cmnds and shell_commands in caseroot
        shell_command_files = glob.glob(os.path.join(include_dir,"shell_commands")) +\
                              glob.glob(os.path.join(include_dir,"xmlchange_cmnds"))
        for shell_commands_file in shell_command_files:
            case_shell_commands = shell_commands_file.replace(include_dir, caseroot)
            with open(shell_commands_file,"r") as fd:
                new_shell_commands = fd.read().replace("xmlchange","xmlchange --force")
            with open(case_shell_commands, "a") as fd:
                fd.write(new_shell_commands)

    for shell_command_file in case_shell_command_files:
        if os.path.isfile(shell_command_file):
            os.chmod(shell_command_file, 0777)
            run_cmd_no_fail(shell_command_file)

def build_include_dirs_list(user_mods_path, include_dirs=None):
    '''
    If user_mods_path has a file "include_user_mods" read that
    file and add directories to the include_dirs, recursively check
    each of those directories for further directories.
    The file may also include comments deleneated with # in the first column
    '''
    include_dirs = [] if include_dirs is None else include_dirs
    if user_mods_path is None or user_mods_path == 'UNSET':
        return include_dirs
    expect(os.path.isabs(user_mods_path),
           "Expected full directory path, got '%s'"%user_mods_path)
    expect(os.path.isdir(user_mods_path),
           "Directory not found %s"%user_mods_path)
    logger.info("Adding user mods directory %s"%user_mods_path)
    include_dirs.append(os.path.normpath(user_mods_path))
    include_file = os.path.join(include_dirs[-1],"include_user_mods")
    if os.path.isfile(include_file):
        with open(include_file, "r") as fd:
            for newpath in fd:
                newpath = newpath.rstrip()
                if len(newpath) > 0 and not newpath.startswith("#"):
                    if not os.path.isabs(newpath):
                        newpath = os.path.join(user_mods_path, newpath)
                    if os.path.isabs(newpath):
                        build_include_dirs_list(newpath, include_dirs)
                    else:
                        logger.warn("Could not resolve path '%s' in file '%s'"%newpath,include_file)
    return include_dirs
