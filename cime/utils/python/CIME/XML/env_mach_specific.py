"""
Interface to the env_mach_specific.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase

import string

logger = logging.getLogger(__name__)

# Is not of type EntryID but can use functions from EntryID (e.g
# get_type) otherwise need to implement own functions and make GenericXML parent class
class EnvMachSpecific(EnvBase):

    def __init__(self, caseroot, infile="env_mach_specific.xml"):
        """
        initialize an object interface to file env_mach_specific.xml in the case directory
        """
        fullpath = infile if os.path.isabs(infile) else os.path.join(caseroot, infile)
        EnvBase.__init__(self, caseroot, fullpath)

    def populate(self, machobj):
        """Add entries to the file using information from a Machines object."""
        items = ("module_system", "environment_variables", "mpirun")
        for item in items:
            nodes = machobj.get_first_child_nodes(item)
            for node in nodes:
                self.add_child(node)

    def get_values(self, item, attribute=None, resolved=True, subgroup=None):
        """Returns the value as a string of the first xml element with item as attribute value.
        <element_name attribute='attribute_value>value</element_name>"""

        logger.debug("(get_values) Input values: %s , %s , %s , %s , %s" , self.__class__.__name__ , item, attribute, resolved, subgroup)

        nodes   = [] # List of identified xml elements
        results = [] # List of identified parameters

        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']
        if item :
            nodes = self.get_nodes("*",{"name" : item})
        else :
            # Return all nodes
            logger.debug("Retrieving all parameter")
            nodes = self.get_nodes("env")

        # Return value for first occurence of node with attribute value = item
        for node in nodes:

            group   = super(EnvMachSpecific, self)._get_group(node)
            val     = node.text
            attr    = node.attrib['name']
            t       = self._get_type(node)
            desc    = self._get_description(node)
            #default = super(EnvBase , self).get_default(node)
            default = self._get_default(node)
            filename    = self.filename

            #t   =  super(EnvBase , self).get_type( node )
            v = { 'group' : group , 'attribute' : attr , 'value' : val , 'type' : t , 'description' : desc , 'default' : default , 'file' : filename}
            logger.debug("Found node with value for %s = %s" , item , v )
            results.append(v)

        return results

    def _get_env_for_case(self, compiler, debug, mpilib):
        module_nodes = self.get_nodes("modules")
        env_nodes    = self.get_nodes("environment_variables")

        modules_to_load = None
        if module_nodes is not None:
            modules_to_load = self._compute_module_actions(module_nodes, compiler, debug, mpilib)

        envs_to_set = None
        if env_nodes is not None:
            envs_to_set = self._compute_env_actions(env_nodes, compiler, debug, mpilib)

        return modules_to_load, envs_to_set

    def load_env_for_case(self, compiler, debug, mpilib):
        modules_to_load, envs_to_set = self._get_env_for_case(compiler, debug, mpilib)

        if (modules_to_load is not None):
            self.load_modules(modules_to_load)
        if (envs_to_set is not None):
            self.load_envs(envs_to_set)

    def load_modules(self, modules_to_load):
        module_system = self.get_module_system_type()
        if (module_system == "module"):
            self._load_module_modules(modules_to_load)
        elif (module_system == "soft"):
            self._load_soft_modules(modules_to_load)
        elif (module_system == "dotkit"):
            self._load_dotkit_modules(modules_to_load)
        elif (module_system == "none"):
            self._load_none_modules(modules_to_load)
        else:
            expect(False, "Unhandled module system '%s'" % module_system)

    def list_modules(self):
        module_system = self.get_module_system_type()

        # If the user's login shell is not sh, it's possible that modules
        # won't be configured so we need to be sure to source the module
        # setup script if it exists.
        init_path = self.get_module_system_init_path("sh")
        if init_path:
            source_cmd = "source %s && " % init_path
        else:
            source_cmd = ""

        if (module_system == "module"):
            return run_cmd_no_fail("%smodule list 2>&1" % source_cmd)
        elif (module_system == "soft"):
            # Does soft really not provide this capability?
            return ""
        elif (module_system == "dotkit"):
            return run_cmd_no_fail("%suse -lv" % source_cmd)
        elif (module_system == "none"):
            return ""
        else:
            expect(False, "Unhandled module system '%s'" % module_system)

    def make_env_mach_specific_file(self, compiler, debug, mpilib, shell):
        modules_to_load, envs_to_set = self._get_env_for_case(compiler, debug, mpilib)

        filename = ".env_mach_specific.%s" % shell
        lines = []
        if modules_to_load is not None:
            lines.extend(self._get_module_commands(modules_to_load, shell))

        if envs_to_set is not None:
            for env_name, env_value in envs_to_set:
                # Let bash do the work on evaluating and resolving env_value
                if shell == "sh":
                    lines.append("export %s=%s" % (env_name, env_value))
                elif shell == "csh":
                    lines.append("setenv %s %s" % (env_name, env_value))
                else:
                    expect(False, "Unknown shell type: '%s'" % shell)

        with open(filename, "w") as fd:
            fd.write("\n".join(lines))

    def load_envs(self, envs_to_set):
        for env_name, env_value in envs_to_set:
            # Let bash do the work on evaluating and resolving env_value
            os.environ[env_name] = run_cmd_no_fail("echo %s" % env_value)

    # Private API

    def _compute_module_actions(self, module_nodes, compiler, debug, mpilib):
        return self._compute_actions(module_nodes, "command", compiler, debug, mpilib)

    def _compute_env_actions(self, env_nodes, compiler, debug, mpilib):
        return self._compute_actions(env_nodes, "env", compiler, debug, mpilib)

    def _compute_actions(self, nodes, child_tag, compiler, debug, mpilib):
        result = [] # list of tuples ("name", "argument")

        for node in nodes:
            if (self._match_attribs(node.attrib, compiler, debug, mpilib)):
                for child in node:
                    expect(child.tag == child_tag, "Expected %s element" % child_tag)
                    if (self._match_attribs(child.attrib, compiler, debug, mpilib)):
                        result.append( (child.get("name"), child.text) )

        return result

    def _match_attribs(self, attribs, compiler, debug, mpilib):
        if ("compiler" in attribs and
            not self._match(compiler, attribs["compiler"])):
            return False
        elif ("mpilib" in attribs and
            not self._match(mpilib, attribs["mpilib"])):
            return False
        elif ("debug" in attribs and
            not self._match("TRUE" if debug else "FALSE", attribs["debug"].upper())):
            return False

        return True

    def _match(self, my_value, xml_value):
        if (xml_value.startswith("!")):
            return my_value != xml_value[1:]
        else:
            return my_value == xml_value

    def _get_module_commands(self, modules_to_load, shell):
        # Note this is independent of module system type
        mod_cmd = self.get_module_system_cmd_path(shell)
        cmds = []
        for action, argument in modules_to_load:
            if argument is None:
                argument = ""
            cmds.append("%s %s %s" % (mod_cmd, action, argument))
        return cmds

    def _load_module_modules(self, modules_to_load):
        for cmd in self._get_module_commands(modules_to_load, "python"):
            py_module_code = run_cmd_no_fail(cmd)
            exec(py_module_code)

    def _load_soft_modules(self, modules_to_load):
        sh_init_cmd = self.get_module_system_init_path("sh")
        sh_mod_cmd = self.get_module_system_cmd_path("sh")

        # Some machines can set the environment
        # variables using a script (such as /etc/profile.d/00softenv.sh
        # on mira or /etc/profile.d/a_softenv.sh on blues)
        # which load the new environment variables using softenv-load.

        # Other machines need to run soft-dec.sh and evaluate the output,
        # which may or may not have unresolved variables such as
        # PATH=/soft/com/packages/intel/16/initial/bin:${PATH}

        cmd = "source %s" % sh_init_cmd

        if os.environ.has_key("SOFTENV_ALIASES"):
            cmd += " && source $SOFTENV_ALIASES"
        if os.environ.has_key("SOFTENV_LOAD"):
            cmd += " && source $SOFTENV_LOAD"

        for action,argument in modules_to_load:
            cmd += " && %s %s %s" % (sh_mod_cmd, action, argument)

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
                       'string value %s unable to be resolved' % valunresolved)
                newenv[key] = val

        # Set environment with new or updated values
        for key in newenv:
            if key in os.environ and key not in newenv:
                del(os.environ[key])
            else:
                os.environ[key] = newenv[key]

    def _load_dotkit_modules(self, _):
        expect(False, "Not yet implemented")

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
#!/usr/bin/env %s
#===============================================================================
# Automatically generated module settings for $self->{machine}
# DO NOT EDIT THIS FILE DIRECTLY!  Please edit env_mach_specific.xml
# in your CASEROOT. This file is overwritten every time modules are loaded!
#===============================================================================
'''%shell
        header += "source %s"%self.get_module_system_init_path(shell)
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

