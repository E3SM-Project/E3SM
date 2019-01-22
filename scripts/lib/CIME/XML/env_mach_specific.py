"""
Interface to the env_mach_specific.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase
from CIME.utils import transform_vars, get_cime_root
import string, resource
from collections import OrderedDict

logger = logging.getLogger(__name__)

# Is not of type EntryID but can use functions from EntryID (e.g
# get_type) otherwise need to implement own functions and make GenericXML parent class
class EnvMachSpecific(EnvBase):
    # pylint: disable=unused-argument
    def __init__(self, caseroot=None, infile="env_mach_specific.xml",
                 components=None, unit_testing=False, read_only=False):
        """
        initialize an object interface to file env_mach_specific.xml in the case directory
        """
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_mach_specific.xsd")
        EnvBase.__init__(self, caseroot, infile, schema=schema, read_only=read_only)
        self._allowed_mpi_attributes = ("compiler", "mpilib", "threaded", "unit_testing", "queue")
        self._unit_testing = unit_testing

    def populate(self, machobj):
        """Add entries to the file using information from a Machines object."""
        items = ("module_system", "environment_variables", "resource_limits", "mpirun", "run_exe","run_misc_suffix")
        default_run_suffix = machobj.get_child("default_run_suffix", root=machobj.root)
        default_run_exe_node = machobj.get_child("default_run_exe", root=default_run_suffix)
        default_run_misc_suffix_node = machobj.get_child("default_run_misc_suffix", root=default_run_suffix)

        group_node = self.make_child("group", {"id":"compliant_values"})

        for item in items:
            nodes = machobj.get_first_child_nodes(item)
            if item == "run_exe" or item == "run_misc_suffix":
                if len(nodes) == 0:
                    value = self.text(default_run_exe_node) if item == "run_exe" else self.text(default_run_misc_suffix_node)
                else:
                    value = nodes[0].text

                entity_node = self.make_child("entry", {"id":item, "value":value}, root=group_node)

                self.make_child("type", root=entity_node, text="char")
                self.make_child("desc", root=entity_node, text=("executable name" if item == "run_exe" else "redirect for job output"))

            else:
                for node in nodes:
                    self.add_child(node)

    def _get_modules_for_case(self, case, job=None):
        module_nodes = self.get_children("modules", root=self.get_child("module_system"))
        modules_to_load = None
        if module_nodes is not None:
            modules_to_load = self._compute_module_actions(module_nodes, case, job=job)

        return modules_to_load

    def _get_envs_for_case(self, case, job=None):
        env_nodes = self.get_children("environment_variables")

        envs_to_set = None
        if env_nodes is not None:
            envs_to_set = self._compute_env_actions(env_nodes, case, job=job)

        return envs_to_set

    def load_env(self, case, force_method=None, job=None, verbose=False):
        """
        Should only be called by case.load_env
        """
        # Do the modules so we can refer to env vars set by the modules
        # in the environment_variables block
        modules_to_load = self._get_modules_for_case(case)
        if (modules_to_load is not None):
            self._load_modules(modules_to_load, force_method=force_method, verbose=verbose)

        envs_to_set = self._get_envs_for_case(case, job=job)
        if (envs_to_set is not None):
            self._load_envs(envs_to_set, verbose=verbose)

        self._get_resources_for_case(case)

    def _get_resources_for_case(self, case):
        resource_nodes = self.get_children("resource_limits")
        if resource_nodes is not None:
            nodes = self._compute_resource_actions(resource_nodes, case)
            for name, val in nodes:
                attr = getattr(resource, name)
                limits = resource.getrlimit(attr)
                logger.info("Setting resource.{} to {} from {}".format(name, val, limits))
                limits = (int(val), limits[1])
                resource.setrlimit(attr, limits)

    def _load_modules(self, modules_to_load, force_method=None, verbose=False):
        module_system = self.get_module_system_type() if force_method is None else force_method
        if (module_system == "module"):
            self._load_module_modules(modules_to_load, verbose=verbose)
        elif (module_system == "soft"):
            self._load_modules_generic(modules_to_load, verbose=verbose)
        elif (module_system == "generic"):
            self._load_modules_generic(modules_to_load, verbose=verbose)
        elif (module_system == "none"):
            self._load_none_modules(modules_to_load)
        else:
            expect(False, "Unhandled module system '{}'".format(module_system))

    def list_modules(self):
        module_system = self.get_module_system_type()

        # If the user's login shell is not sh, it's possible that modules
        # won't be configured so we need to be sure to source the module
        # setup script if it exists.
        init_path = self.get_module_system_init_path("sh")
        if init_path:
            source_cmd = "source {} && ".format(init_path)
        else:
            source_cmd = ""

        if (module_system in ["module"]):
            return run_cmd_no_fail("{}module list".format(source_cmd), combine_output=True)
        elif (module_system == "soft"):
            # Does soft really not provide this capability?
            return ""
        elif (module_system == "generic"):
            return run_cmd_no_fail("{}use -lv".format(source_cmd))
        elif (module_system == "none"):
            return ""
        else:
            expect(False, "Unhandled module system '{}'".format(module_system))

    def save_all_env_info(self, filename):
        """
        Get a string representation of all current environment info and
        save it to file.
        """
        with open(filename, "w") as f:
            f.write(self.list_modules())
        run_cmd_no_fail("echo -e '\n' && env", arg_stdout=filename)

    def make_env_mach_specific_file(self, shell, case):
        module_system = self.get_module_system_type()
        sh_init_cmd = self.get_module_system_init_path(shell)
        sh_mod_cmd = self.get_module_system_cmd_path(shell)
        lines = ["# This file is for user convenience only and is not used by the model"]

        lines.append("# Changes to this file will be ignored and overwritten")
        lines.append("# Changes to the environment should be made in env_mach_specific.xml")
        lines.append("# Run ./case.setup --reset to regenerate this file")
        if sh_init_cmd:
            lines.append("source {}".format(sh_init_cmd))

        if "SOFTENV_ALIASES" in os.environ:
            lines.append("source $SOFTENV_ALIASES")
        if "SOFTENV_LOAD" in os.environ:
            lines.append("source $SOFTENV_LOAD")

        modules_to_load = self._get_modules_for_case(case)
        envs_to_set = self._get_envs_for_case(case)
        filename = ".env_mach_specific.{}".format(shell)
        if modules_to_load is not None:
            if module_system == "module":
                lines.extend(self._get_module_commands(modules_to_load, shell))
            else:
                for action, argument in modules_to_load:
                    lines.append("{} {} {}".format(sh_mod_cmd, action, "" if argument is None else argument))

        if envs_to_set is not None:
            for env_name, env_value in envs_to_set:
                if shell == "sh":
                    lines.append("export {}={}".format(env_name, env_value))
                elif shell == "csh":
                    lines.append("setenv {} {}".format(env_name, env_value))
                else:
                    expect(False, "Unknown shell type: '{}'".format(shell))

        with open(filename, "w") as fd:
            fd.write("\n".join(lines))

    # Private API

    def _load_envs(self, envs_to_set, verbose=False):
        for env_name, env_value in envs_to_set:
            logger_func = logger.warning if verbose else logger.debug
            if env_value is None and env_name in os.environ:
                del os.environ[env_name]
                logger_func("Unsetting Environment {}".format(env_name))
            elif env_value is not None:
                os.environ[env_name] = env_value
                logger_func("Setting Environment {}={}".format(env_name, env_value))

    def _compute_module_actions(self, module_nodes, case, job=None):
        return self._compute_actions(module_nodes, "command", case, job=job)

    def _compute_env_actions(self, env_nodes, case, job=None):
        return self._compute_actions(env_nodes, "env", case, job=job)

    def _compute_resource_actions(self, resource_nodes, case, job=None):
        return self._compute_actions(resource_nodes, "resource", case, job=job)

    def _compute_actions(self, nodes, child_tag, case, job=None):
        result = [] # list of tuples ("name", "argument")
        compiler, mpilib = case.get_value("COMPILER"), case.get_value("MPILIB")

        for node in nodes:
            if (self._match_attribs(self.attrib(node), case, job=job)):
                for child in self.get_children(root=node):
                    expect(self.name(child) == child_tag, "Expected {} element".format(child_tag))
                    if (self._match_attribs(self.attrib(child), case, job=job)):
                        val = self.text(child)
                        if val is not None:
                            # We allow a couple special substitutions for these fields
                            for repl_this, repl_with in [("$COMPILER", compiler), ("$MPILIB", mpilib)]:
                                val = val.replace(repl_this, repl_with)

                            val = self.get_resolved_value(val)
                            expect("$" not in val, "Not safe to leave unresolved items in env var value: '{}'".format(val))

                        # intentional unindent, result is appended even if val is None
                        result.append( (self.get(child, "name"), val) )

        return result

    def _match_attribs(self, attribs, case, job=None):
        # check for matches with case-vars
        for attrib in attribs:
            if attrib == "unit_testing": # special case
                if not self._match(self._unit_testing, attribs["unit_testing"].upper()):
                    return False
            elif attrib == "queue":
                if job is not None:
                    val = case.get_value("JOB_QUEUE", subgroup=job)
                    expect(val is not None, "Cannot match attrib '%s', case has no value for it" % attrib.upper())
                    if not self._match(val, attribs[attrib]):
                        return False
            elif attrib == "name":
                pass
            else:
                val = case.get_value(attrib.upper())
                expect(val is not None, "Cannot match attrib '%s', case has no value for it" % attrib.upper())
                if not self._match(val, attribs[attrib]):
                    return False

        return True

    def _match(self, my_value, xml_value):
        if xml_value.startswith("!"):
            result = re.match(xml_value[1:] + "$",str(my_value)) is None
        elif isinstance(my_value, bool):
            if my_value: result = xml_value == "TRUE"
            else: result = xml_value == "FALSE"
        else:
            result = re.match(xml_value + "$",str(my_value)) is not None

        logger.debug("(env_mach_specific) _match {} {} {}".format(my_value, xml_value, result))
        return result

    def _get_module_commands(self, modules_to_load, shell):
        # Note this is independent of module system type
        mod_cmd = self.get_module_system_cmd_path(shell)
        cmds = []
        last_action = None
        last_cmd    = None

        # Normally, we will try to combine or batch module commands together...
        #
        # module load X
        # module load Y
        # module load Z
        #
        # is the same as ...
        #
        # module load X Y Z
        #
        # ... except the latter is significatly faster due to performing 1/3 as
        # many forks.
        #
        # Not all module commands support batching though and we enurmerate those
        # here.
        actions_that_cannot_be_batched = ["swap", "switch"]

        for action, argument in modules_to_load:
            if argument is None:
                argument = ""

            if action == last_action and action not in actions_that_cannot_be_batched:
                last_cmd = "{} {}".format(last_cmd, argument)
            else:
                if last_cmd is not None:
                    cmds.append(last_cmd)

                last_cmd = "{} {} {}".format(mod_cmd, action, "" if argument is None else argument)
                last_action = action

        if last_cmd:
            cmds.append(last_cmd)

        return cmds

    def _load_module_modules(self, modules_to_load, verbose=False):
        logger_func = logger.warning if verbose else logger.debug
        for cmd in self._get_module_commands(modules_to_load, "python"):
            logger_func("module command is {}".format(cmd))
            stat, py_module_code, errout = run_cmd(cmd)
            expect(stat==0 and (len(errout) == 0 or self.allow_error()),
                   "module command {} failed with message:\n{}".format(cmd, errout))
            exec(py_module_code)

    def _load_modules_generic(self, modules_to_load, verbose=False):
        sh_init_cmd = self.get_module_system_init_path("sh")
        sh_mod_cmd = self.get_module_system_cmd_path("sh")
        logger_func = logger.warning if verbose else logger.debug

        # Purpose is for environment management system that does not have
        # a python interface and therefore can only determine what they
        # do by running shell command and looking at the changes
        # in the environment.

        cmd = "source {}".format(sh_init_cmd)

        if "SOFTENV_ALIASES" in os.environ:
            cmd += " && source $SOFTENV_ALIASES"
        if "SOFTENV_LOAD" in os.environ:
            cmd += " && source $SOFTENV_LOAD"

        for action,argument in modules_to_load:
            cmd += " && {} {} {}".format(sh_mod_cmd, action, "" if argument is None else argument)

        # Use null terminated lines to give us something more definitive to split on.
        # Env vars can contain newlines, so splitting on newlines can be ambiguous
        cmd += " && env -0"
        logger_func("cmd: {}".format(cmd))
        output = run_cmd_no_fail(cmd)

        ###################################################
        # Parse the output to set the os.environ dictionary
        ###################################################
        newenv = OrderedDict()
        for line in output.split('\0'):
            if "=" in line:
                key, val = line.split("=", 1)
                newenv[key] = val

        # resolve variables
        for key, val in newenv.items():
            newenv[key] = string.Template(val).safe_substitute(newenv)

        # Set environment with new or updated values
        for key in newenv:
            if key in os.environ and os.environ[key] == newenv[key]:
                pass
            else:
                os.environ[key] = newenv[key]

        for oldkey in list(os.environ.keys()):
            if oldkey not in newenv:
                del os.environ[oldkey]

    def _load_none_modules(self, modules_to_load):
        """
        No Action required
        """
        expect(not modules_to_load,
               "Module system was specified as 'none' yet there are modules that need to be loaded?")

    def _mach_specific_header(self, shell):
        '''
        write a shell module file for this case.
        '''
        header = '''
#!/usr/bin/env {}
#===============================================================================
# Automatically generated module settings for $self->{{machine}}
# DO NOT EDIT THIS FILE DIRECTLY!  Please edit env_mach_specific.xml
# in your CASEROOT. This file is overwritten every time modules are loaded!
#===============================================================================
'''.format(shell)
        header += "source {}".format(self.get_module_system_init_path(shell))
        return header

    def get_module_system_type(self):
        """
        Return the module system used on this machine
        """
        module_system = self.get_child("module_system")
        return self.get(module_system, "type")

    def allow_error(self):
        """
        Return True if stderr output from module commands should be assumed
        to be an error. Default False. This is necessary since implementations
        of environment modules are highlty variable and some systems produce
        stderr output even when things are working fine.
        """
        module_system = self.get_child("module_system")
        value = self.get(module_system, "allow_error")
        return value.upper() == "TRUE" if value is not None else False

    def get_module_system_init_path(self, lang):
        init_nodes = self.get_optional_child("init_path", attributes={"lang":lang}, root=self.get_child("module_system"))
        return self.text(init_nodes) if init_nodes is not None else None

    def get_module_system_cmd_path(self, lang):
        cmd_nodes = self.get_optional_child("cmd_path", attributes={"lang":lang}, root=self.get_child("module_system"))
        return self.text(cmd_nodes) if cmd_nodes is not None else None

    def get_mpirun(self, case, attribs, job, exe_only=False):
        """
        Find best match, return (executable, {arg_name : text})
        """
        mpirun_nodes = self.get_children("mpirun")
        best_match = None
        best_num_matched = -1
        default_match = None
        best_num_matched_default = -1
        args = []
        for mpirun_node in mpirun_nodes:
            xml_attribs = self.attrib(mpirun_node)
            all_match = True
            matches = 0
            is_default = False

            for key, value in attribs.items():
                expect(key in self._allowed_mpi_attributes, "Unexpected key {} in mpirun attributes".format(key))
                if key in xml_attribs:
                    if xml_attribs[key].lower() == "false":
                        xml_attrib = False
                    elif xml_attribs[key].lower() == "true":
                        xml_attrib = True
                    else:
                        xml_attrib = xml_attribs[key]

                    if xml_attrib == value:
                        matches += 1
                    elif key == "mpilib" and value != "mpi-serial" and xml_attrib == "default":
                        is_default = True
                    else:
                        all_match = False
                        break

            if all_match:
                if is_default:
                    if matches > best_num_matched_default:
                        default_match = mpirun_node
                        best_num_matched_default = matches
                else:
                    if matches > best_num_matched:
                        best_match = mpirun_node
                        best_num_matched = matches

        # if there are no special arguments required for mpi-serial it need not have an entry in config_machines.xml
        if "mpilib" in attribs and attribs["mpilib"] == "mpi-serial" and best_match is None:
            return "",[]

        expect(best_match is not None or default_match is not None,
               "Could not find a matching MPI for attributes: {}".format(attribs))

        the_match = best_match if best_match is not None else default_match

        # Now that we know the best match, compute the arguments
        if not exe_only:
            arg_node = self.get_optional_child("arguments", root=the_match)
            if arg_node is not None:
                arg_nodes = self.get_children("arg", root=arg_node)
                for arg_node in arg_nodes:
                    arg_value = transform_vars(self.text(arg_node),
                                               case=case,
                                               subgroup=job,
                                               default=self.get(arg_node, "default"))
                    args.append(arg_value)

        exec_node = self.get_child("executable", root=the_match)
        expect(exec_node is not None,"No executable found")
        executable = self.text(exec_node)

        return executable, args

    def get_type_info(self, vid):
        return "char"
