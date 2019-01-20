"""
user_mod_support.py
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, run_cmd_no_fail, safe_copy
import glob

logger = logging.getLogger(__name__)

def apply_user_mods(caseroot, user_mods_path, keepexe=None):
    '''
    Recursivlely apply user_mods to caseroot - this includes updating user_nl_xxx,
    updating SourceMods and creating case shell_commands and xmlchange_cmds files

    First remove case shell_commands files if any already exist

    If this function is called multiple times, settings from later calls will
    take precedence over earlier calls, if there are conflicts.

    keepexe is an optional argument that is needed for cases where apply_user_mods is
    called from create_clone
    '''
    case_shell_command_files = [os.path.join(caseroot,"shell_commands"),
                           os.path.join(caseroot,"xmlchange_cmnds")]
    for shell_command_file in case_shell_command_files:
        if os.path.isfile(shell_command_file):
            os.remove(shell_command_file)

    include_dirs = build_include_dirs_list(user_mods_path)
    # If a user_mods dir 'foo' includes 'bar', the include_dirs list returned
    # from build_include_dirs has 'foo' before 'bar'. But with the below code,
    # directories that occur later in the list take precedence over the earlier
    # ones, and we want 'foo' to take precedence over 'bar' in this case (in
    # general: we want a given user_mods directory to take precedence over any
    # mods that it includes). So we reverse include_dirs to accomplish this.
    include_dirs.reverse()
    logger.debug("include_dirs are {}".format(include_dirs))
    for include_dir in include_dirs:
        # write user_nl_xxx file in caseroot
        for user_nl in glob.iglob(os.path.join(include_dir,"user_nl_*")):
            with open(os.path.join(include_dir, user_nl), "r") as fd:
                newcontents = fd.read()
            if len(newcontents) == 0:
                continue
            case_user_nl = user_nl.replace(include_dir, caseroot)
            # If the same variable is set twice in a user_nl file, the later one
            # takes precedence. So by appending the new contents, later entries
            # in the include_dirs list take precedence over earlier entries.
            with open(case_user_nl, "a") as fd:
                fd.write(newcontents)

        # update SourceMods in caseroot
        for root, _, files in os.walk(include_dir,followlinks=True,topdown=False):
            if "src" in os.path.basename(root):
                if keepexe is not None:
                    expect(False,
                           "cannot have any source mods in {} if keepexe is an option".format(user_mods_path))
                for sfile in files:
                    source_mods = os.path.join(root,sfile)
                    case_source_mods = source_mods.replace(include_dir, caseroot)
                    # We overwrite any existing SourceMods file so that later
                    # include_dirs take precedence over earlier ones
                    if os.path.isfile(case_source_mods):
                        logger.warning("WARNING: Overwriting existing SourceMods in {}".format(case_source_mods))
                    else:
                        logger.info("Adding SourceMod to case {}".format(case_source_mods))
                    try:
                        safe_copy(source_mods, case_source_mods)
                    except Exception:
                        expect(False, "Could not write file {} in caseroot {}".format(case_source_mods,caseroot))

        # create xmlchange_cmnds and shell_commands in caseroot
        shell_command_files = glob.glob(os.path.join(include_dir,"shell_commands")) +\
                              glob.glob(os.path.join(include_dir,"xmlchange_cmnds"))
        for shell_commands_file in shell_command_files:
            case_shell_commands = shell_commands_file.replace(include_dir, caseroot)
            # add commands from both shell_commands and xmlchange_cmnds to
            # the same file (caseroot/shell_commands)
            case_shell_commands = case_shell_commands.replace("xmlchange_cmnds","shell_commands")
            # Note that use of xmlchange_cmnds has been deprecated and will soon
            # be removed altogether, so new tests should rely on shell_commands
            if shell_commands_file.endswith("xmlchange_cmnds"):
                logger.warning("xmlchange_cmnds is deprecated and will be removed " +\
                            "in a future release; please rename {} shell_commands".format(shell_commands_file))
            with open(shell_commands_file,"r") as fd:
                new_shell_commands = fd.read().replace("xmlchange","xmlchange --force")
            # By appending the new commands to the end, settings from later
            # include_dirs take precedence over earlier ones
            with open(case_shell_commands, "a") as fd:
                fd.write(new_shell_commands)

    for shell_command_file in case_shell_command_files:
        if os.path.isfile(shell_command_file):
            os.chmod(shell_command_file, 0o777)
            run_cmd_no_fail(shell_command_file,verbose=True)


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
           "Expected full directory path, got '{}'".format(user_mods_path))
    expect(os.path.isdir(user_mods_path),
           "Directory not found {}".format(user_mods_path))
    norm_path = os.path.normpath(user_mods_path)

    for dir_ in include_dirs:
        if norm_path == dir_:
            include_dirs.remove(norm_path)
            break

    logger.info("Adding user mods directory {}".format(norm_path))
    include_dirs.append(norm_path)
    include_file = os.path.join(norm_path,"include_user_mods")
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
                        logger.warning("Could not resolve path '{}' in file '{}'".format(newpath, include_file))

    return include_dirs
