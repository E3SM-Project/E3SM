"""
Interface to the env_mach_specific.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase
from CIME.utils import transform_vars, get_cime_root
import string

logger = logging.getLogger(__name__)

# Is not of type EntryID but can use functions from EntryID (e.g
# get_type) otherwise need to implement own functions and make GenericXML parent class
class EnvMachSpecific(EnvBase):
    # pylint: disable=unused-argument
    def __init__(self, caseroot=None, infile="env_mach_specific.xml",
                 components=None, unit_testing=False):
        """
        initialize an object interface to file env_mach_specific.xml in the case directory
        """
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_mach_specific.xsd")
        EnvBase.__init__(self, caseroot, infile, schema=schema)
        self._allowed_mpi_attributes = ("compiler", "mpilib", "threaded", "unit_testing")
        self._unit_testing = unit_testing

    def populate(self, machobj):
        """Add entries to the file using information from a Machines object."""
        items = ("module_system", "environment_variables", "mpirun", "run_exe","run_misc_suffix")
        default_run_exe_node = machobj.get_node("default_run_exe")
        default_run_misc_suffix_node = machobj.get_node("default_run_misc_suffix")

        for item in items:
            nodes = machobj.get_first_child_nodes(item)
            if item == "run_exe" or item == "run_misc_suffix":
                newnode = ET.Element(item)
                newnode.tag = "entry"
                newnode.set("id", item)
                if len(nodes) == 0:
                    if item == "run_exe":
                        value = default_run_exe_node.text
                    else:
                        value =  default_run_misc_suffix_node.text
                else:
                    value = nodes[0].text
                newnode.set("value", value)
                newnode = self.add_sub_node(newnode, "type", "char")
                if item == "run_exe":
                    newnode = self.add_sub_node(newnode, "desc", "executable name")
                else:
                    newnode = self.add_sub_node(newnode, "desc", "redirect for job output")

                self.add_child(newnode)
            else:
                for node in nodes:
                    self.add_child(node)

    def _get_modules_for_case(self, case):
        module_nodes = self.get_nodes("modules")
        modules_to_load = None
        if module_nodes is not None:
            modules_to_load = self._compute_module_actions(module_nodes, case)

        return modules_to_load

    def _get_envs_for_case(self, case):
        env_nodes = self.get_nodes("environment_variables")

        envs_to_set = None
        if env_nodes is not None:
            envs_to_set = self._compute_env_actions(env_nodes, case)

        return envs_to_set

    def load_env(self, case):
        """
        Should only be called by case.load_env
        """
        # Do the modules so we can refer to env vars set by the modules
        # in the environment_variables block
        modules_to_load = self._get_modules_for_case(case)
        if (modules_to_load is not None):
            self.load_modules(modules_to_load)

        envs_to_set = self._get_envs_for_case(case)
        if (envs_to_set is not None):
            self.load_envs(envs_to_set)

    def load_modules(self, modules_to_load):
        module_system = self.get_module_system_type()
        if (module_system == "module"):
            self._load_module_modules(modules_to_load)
        elif (module_system == "module_lmod"):
            self._load_modules_generic(modules_to_load)
        elif (module_system == "soft"):
            self._load_modules_generic(modules_to_load)
        elif (module_system == "generic"):
            self._load_modules_generic(modules_to_load)
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

        if (module_system in ["module", "module_lmod"]):
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
        modules_to_load = self._get_modules_for_case(case)
        envs_to_set = self._get_envs_for_case(case)
        filename = ".env_mach_specific.{}".format(shell)
        lines = []
        if modules_to_load is not None:
            lines.extend(self._get_module_commands(modules_to_load, shell))

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

    def load_envs(self, envs_to_set):
        for env_name, env_value in envs_to_set:
            os.environ[env_name] = "" if env_value is None else env_value

    # Private API

    def _compute_module_actions(self, module_nodes, case):
        return self._compute_actions(module_nodes, "command", case)

    def _compute_env_actions(self, env_nodes, case):
        return self._compute_actions(env_nodes, "env", case)

    def _compute_actions(self, nodes, child_tag, case):
        result = [] # list of tuples ("name", "argument")
        compiler, mpilib = case.get_value("COMPILER"), case.get_value("MPILIB")

        for node in nodes:
            if (self._match_attribs(node.attrib, case)):
                for child in node:
                    expect(child.tag == child_tag, "Expected {} element".format(child_tag))
                    if (self._match_attribs(child.attrib, case)):
                        val = child.text
                        if val is not None:
                            # We allow a couple special substitutions for these fields
                            for repl_this, repl_with in [("$COMPILER", compiler), ("$MPILIB", mpilib)]:
                                val = val.replace(repl_this, repl_with)

                            val = self.get_resolved_value(val)
                            expect("$" not in val, "Not safe to leave unresolved items in env var value: '{}'".format(val))

                        # intentional unindent, result is appended even if val is None
                        result.append( (child.get("name"), val) )

        return result

    def _match_attribs(self, attribs, case):
        # check for matches with case-vars
        for attrib in attribs:
            if attrib == "unit_testing": # special case
                if not self._match(self._unit_testing, attribs["unit_testing"].upper()):
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
            result = my_value != xml_value[1:]
        elif isinstance(my_value, bool):
            if my_value: result = xml_value == "TRUE"
            else: result = xml_value == "FALSE"
        else:
            result = my_value == xml_value

        logger.debug("(env_mach_specific) _match {} {} {}".format(my_value, xml_value, result))
        return result

    def _get_module_commands(self, modules_to_load, shell):
        # Note this is independent of module system type
        mod_cmd = self.get_module_system_cmd_path(shell)
        cmds = []
        for action, argument in modules_to_load:
            if argument is None:
                argument = ""
            cmds.append("{} {} {}".format(mod_cmd, action, argument))
        return cmds

    def _load_module_modules(self, modules_to_load):
        for cmd in self._get_module_commands(modules_to_load, "python"):
            logger.debug("module command is {}".format(cmd))
            stat, py_module_code, errout = run_cmd(cmd)
            expect(stat==0 and len(errout) == 0,
                   "module command {} failed with message:\n{}".format(cmd, errout))
            exec(py_module_code)

    def _load_modules_generic(self, modules_to_load):
        sh_init_cmd = self.get_module_system_init_path("sh")
        sh_mod_cmd = self.get_module_system_cmd_path("sh")

        # Purpose is for environment management system that does not have
        # a python interface and therefore can only determine what they
        # do by running shell command and looking at the changes
        # in the environment.

        cmd = "source {}".format(sh_init_cmd)

        if os.environ.has_key("SOFTENV_ALIASES"):
            cmd += " && source $SOFTENV_ALIASES"
        if os.environ.has_key("SOFTENV_LOAD"):
            cmd += " && source $SOFTENV_LOAD"

        for action,argument in modules_to_load:
            cmd += " && {} {} {}".format(sh_mod_cmd, action, argument)

        cmd += " && env"
        output = run_cmd_no_fail(cmd)

        ###################################################
        # Parse the output to set the os.environ dictionary
        ###################################################
        newenv = {}
        dolater = []
        for line in output.splitlines():
            if line.find('$')>0:
                dolater.append(line)
                continue

            m=re.match(r'^(\S+)=(\S+)\s*;*\s*$',line)
            if m:
                key = m.groups()[0]
                val = m.groups()[1]
                newenv[key] = val

        # Now that initial newenv has been set, resolve variables
        for line in dolater:
            m=re.match(r'^(\S+)=(\S+)\s*;*\s*$',line)
            if m:
                key = m.groups()[0]
                valunresolved = m.groups()[1]
                val = string.Template(valunresolved).safe_substitute(newenv)
                expect(val is not None,
                       'string value {} unable to be resolved'.format(valunresolved))
                newenv[key] = val

        # Set environment with new or updated values
        for key in newenv:
            if key in os.environ and key not in newenv:
                del(os.environ[key])
            else:
                os.environ[key] = newenv[key]

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
        module_system = self.get_node("module_system")
        return module_system.get("type")

    def get_module_system_init_path(self, lang):
        init_nodes = self.get_optional_node("init_path", attributes={"lang":lang})
        return init_nodes.text if init_nodes is not None else None

    def get_module_system_cmd_path(self, lang):
        cmd_nodes = self.get_optional_node("cmd_path", attributes={"lang":lang})
        return cmd_nodes.text if cmd_nodes is not None else None

    def get_mpirun(self, case, attribs, check_members=None, job="case.run", exe_only=False):
        """
        Find best match, return (executable, {arg_name : text})
        """
        mpirun_nodes = self.get_nodes("mpirun")
        best_match = None
        best_num_matched = -1
        default_match = None
        best_num_matched_default = -1
        args = []
        for mpirun_node in mpirun_nodes:
            xml_attribs = mpirun_node.attrib
            all_match = True
            matches = 0
            is_default = False

            for key, value in attribs.iteritems():
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
            arg_node = self.get_optional_node("arguments", root=the_match)
            if arg_node is not None:
                arg_nodes = self.get_nodes("arg", root=arg_node)
                for arg_node in arg_nodes:
                    arg_value = transform_vars(arg_node.text,
                                               case=case,
                                               subgroup=job,
                                               check_members=check_members,
                                               default=arg_node.get("default"))
                    args.append(arg_value)

        exec_node = self.get_node("executable", root=the_match)
        expect(exec_node is not None,"No executable found")
        executable = exec_node.text

        return executable, args
