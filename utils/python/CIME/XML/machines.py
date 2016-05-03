"""
Interface to the config_machines.xml file.  This class inherits from GenericXML.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.utils import expect

import socket

logger = logging.getLogger(__name__)

class Machines(GenericXML):

    def __init__(self, infile=None, files=None, machine=None):
        """
        initialize an object
        if a filename is provided it will be used,
        otherwise if a files object is provided it will be used
        otherwise create a files object from default values
        """

        self.machine_node = None
        self.machine = None
        self.machines_dir = None

        if infile is None:
            if files is None:
                files = Files()
            infile = files.get_value("MACHINES_SPEC_FILE")
            self.machines_dir = os.path.dirname(infile)

        GenericXML.__init__(self, infile)

        """
        Append the contents of $HOME/.cime/config_machines.xml if it exists
        This could cause problems if node matchs are repeated when only one is expected
        """
        infile = os.path.join(os.environ.get("HOME"),".cime","config_machines.xml")
        if os.path.exists(infile):
            GenericXML.read(self, infile)

        if machine is None:
            machine = self.probe_machine_name()

        expect(machine is not None, "Could not initialize machine object")
        self.set_machine(machine)

    def get_machines_dir(self):
        """
        Return the directory of the machines file
        """
        return self.machines_dir

    def get_machine_name(self):
        """
        Return the name of the machine
        """
        return self.machine

    def get_node_names(self):
        """
        Return the names of all the child nodes for the target machine
        """
        nodes = self.machine_node.findall("./")
        node_names = []
        for node in nodes:
            node_names.append(node.tag)
        return node_names

    def get_first_child_nodes(self, nodename):
        """
        Return the names of all the child nodes for the target machine
        """
        nodes = self.machine_node.findall("./" + nodename)
        return nodes

    def list_available_machines(self):
        """
        Return a list of machines defined for a given CIME_MODEL
        """
        machines = []
        nodes  = self.get_nodes("machine")
        for node in nodes:
            mach = node.get("MACH")
            machines.append(mach)
        return machines

    def probe_machine_name(self):
        """
        Find a matching regular expression for hostname
        in the NODENAME_REGEX field in the file.   First match wins.
        """
        machine = None
        nametomatch = socket.gethostname().split(".")[0]
        nodes = self.get_nodes("machine")

        for node in nodes:
            machtocheck = node.get("MACH")
            logger.debug("machine is " + machtocheck)
            regex_str_node = self.get_optional_node("NODENAME_REGEX", root=node)
            regex_str = machtocheck if regex_str_node is None else regex_str_node.text

            if regex_str is not None:
                logger.debug("machine regex string is " + regex_str)
                regex = re.compile(regex_str)
                if regex.match(nametomatch):
                    logger.debug("Found machine: %s matches %s" % (machtocheck, nametomatch))
                    machine = machtocheck
                    break

        if machine is None:
            # Check for a local definition
            logger.warning("Could not probe machine for hostname '%s'" % nametomatch)

        return machine


    def set_machine(self, machine):
        """
        Sets the machine block in the Machines object

        >>> machobj = Machines(machine="melvin")
        >>> machobj.get_machine_name()
        'melvin'
        >>> machobj.set_machine("trump")
        Traceback (most recent call last):
        ...
        SystemExit: ERROR: No machine trump found
        """
        if self.machine != machine or self.machine_node is None:
            self.machine_node = self.get_optional_node("machine", {"MACH" : machine})
            expect(self.machine_node is not None, "No machine %s found" % machine)
            self.machine = machine

        return machine


    def get_value(self, name, resolved=True):
        """
        Get Value of fields in the config_machines.xml file
        """
        expect(self.machine_node is not None, "Machine object has no machine defined")
        value = None

        # COMPILER and MPILIB are special, if called without arguments they get the default value from the
        # COMPILERS and MPILIBS lists in the file.
        if name == "COMPILER":
            value = self.get_default_compiler()
        elif name == "MPILIB":
            value = self.get_default_MPIlib()
        else:
            node = self.get_optional_node(name, root=self.machine_node)
            if node is not None:
                value = node.text

        if value is None:
            # if all else fails
            value = GenericXML.get_value(self, name)

        if resolved:
            if value is not None:
                value = self.get_resolved_value(value)
            elif name in os.environ:
                value = os.environ[name]

        return value

    def get_field_from_list(self, listname, reqval=None):
        """
        Some of the fields have lists of valid values in the xml, parse these
        lists and return the first value if reqval is not provided and reqval
        if it is a valid setting for the machine
        """
        expect(self.machine_node is not None, "Machine object has no machine defined")
        supported_values = self.get_value(listname)
        expect(supported_values is not None,
               "No list found for " + listname + " on machine " + self.machine)
        supported_values = supported_values.split(",")

        if reqval is None or reqval == "UNSET":
            return supported_values[0]
        for val in supported_values:
            if val == reqval:
                return reqval

        expect(False, "%s value %s not supported for machine %s" %
               (listname, reqval, self.machine))

    def get_default_compiler(self):
        """
        Get the compiler to use from the list of COMPILERS
        """
        return self.get_field_from_list("COMPILERS")

    def get_default_MPIlib(self):
        """
        Get the MPILIB to use from the list of MPILIBS
        """
        return self.get_field_from_list("MPILIBS")

    def is_valid_compiler(self,compiler):
        """
        Check the compiler is valid for the current machine

        >>> machobj = Machines(machine="edison")
        >>> machobj.get_default_compiler()
        'intel'
        >>> machobj.is_valid_compiler("cray")
        True
        >>> machobj.is_valid_compiler("nag")
        Traceback (most recent call last):
        ...
        SystemExit: ERROR: COMPILERS value nag not supported for machine edison
        """
        if self.get_field_from_list("COMPILERS", compiler) is not None:
            return True
        return False

    def is_valid_MPIlib(self, mpilib):
        """
        Check the MPILIB is valid for the current machine

        >>> machobj = Machines(machine="edison")
        >>> machobj.is_valid_MPIlib("mpi-serial")
        True
        """
        if self.get_field_from_list("MPILIBS", mpilib) is not None:
            return True
        return False

    def has_batch_system(self):
        """
        Return if this machine has a batch system

        >>> machobj = Machines(machine="edison")
        >>> machobj.has_batch_system()
        True
        >>> machobj.set_machine("melvin")
        'melvin'
        >>> machobj.has_batch_system()
        False
        """
        result = False
        batch_system = self.get_optional_node("batch_system", root=self.machine_node)
        if batch_system is not None:
            result = (batch_system.get("type") != "none")
        logger.debug("Machine %s has batch: %s" % (self.machine, result))
        return result

    def get_batch_system_type(self):
        """
        Return the batch system used on this machine

        >>> machobj = Machines(machine="edison")
        >>> machobj.get_batch_system_type()
        'slurm'
        """
        batch_system = self.get_node("batch_system", root=self.machine_node)
        return batch_system.get("type")

    def get_module_system_type(self):
        """
        Return the module system used on this machine

        >>> machobj = Machines()
        >>> name = machobj.set_machine("edison")
        >>> machobj.get_module_system_type()
        'module'
        """
        module_system = self.get_node("module_system", root=self.machine_node)
        return module_system.get("type")

    def get_module_system_init_path(self, lang):
        init_nodes = self.get_node("init_path", attributes={"lang":lang}, root=self.machine_node)
        return init_nodes.text

    def get_module_system_cmd_path(self, lang):
        cmd_nodes = self.get_node("cmd_path", attributes={"lang":lang}, root=self.machine_node)
        return cmd_nodes.text

    def get_default_queue(self):
        return self.get_optional_node("queue", attributes={"default" : "true"}, root=self.machine_node)

    def get_all_queues(self):
        return self.get_nodes("queue", root=self.machine_node)

    def get_walltimes(self):
        return self.get_nodes("walltime", root=self.machine_node)

    def get_default_walltime(self):
        return self.get_optional_node("walltime", attributes={"default" : "true"}, root=self.machine_node)

    def get_suffix(self, suffix_type):
        node = self.get_optional_node("default_run_suffix")
        if node is not None:
            suffix_node = self.get_optional_node(suffix_type, root=node)
            if suffix_node is not None:
                return suffix_node.text

        return None

    def get_mpirun(self, attribs, batch_maker):
        """
        Find best match, return (executable, {arg_name : text})
        """
        mpirun_nodes = self.get_nodes("mpirun", root=self.machine_node)
        best_match = None
        best_num_matched = -1
        default_match = None
        best_num_matched_default = -1
        args = {}
        for mpirun_node in mpirun_nodes:
            xml_attribs = mpirun_node.attrib
            all_match = True
            matches = 0
            is_default = False
            for key, value in attribs.iteritems():
                if key in xml_attribs:
                    if xml_attribs[key].lower() == "false":
                        xml_attrib = False
                    elif xml_attribs[key].lower() == "true":
                        xml_attrib = True
                    else:
                        xml_attrib = xml_attribs[key]

                    if xml_attrib == value:
                        matches += 1
                    elif key == "mpilib" and xml_attrib == "default":
                        is_default = True
                    else:
                        all_match = False
                        break

            for key in xml_attribs:
                expect(key in attribs, "Unhandled MPI property '%s'" % key)

            if all_match:
                arg_node = self.get_optional_node("arguments", root=mpirun_node)
                if arg_node is not None:
                    arg_nodes = self.get_nodes("arg", root=arg_node)
                    for arg_node in arg_nodes:
                        arg_value = batch_maker.transform_vars(arg_node.text, arg_node.get("default"))
                        args[arg_node.get("name")] = arg_value
                if is_default:
                    if matches > best_num_matched_default:
                        default_match = mpirun_node
                        best_num_matched_default = matches
                else:
                    if matches > best_num_matched:
                        best_match = mpirun_node
                        best_num_matched = matches

        expect(best_match is not None or default_match is not None,
               "Could not find a matching MPI for attributes: %s" % attribs)

        the_match = best_match if best_match is not None else default_match


        executable = self.get_node("executable", root=the_match)

        return executable.text, args

    def print_values(self):
        # write out machines
        machines = self.get_nodes(nodename="machine")
        print  "Machines"
        for machine in machines:
            name = machine.get("MACH")
            desc = machine.find("DESC")
            os_  = machine.find("OS")
            compilers = machine.find("COMPILERS")
            max_tasks_per_node = machine.find("MAX_TASKS_PER_NODE")
            pes_per_node = machine.find("PES_PER_NODE")

            print  "  %s : %s "% (name , desc.text)
            print  "      os             ", os_.text
            print  "      compilers      ",compilers.text
            if pes_per_node is not None:
                print  "      pes/node       ",pes_per_node.text
            if max_tasks_per_node is not None:
                print  "      max_tasks/node ",max_tasks_per_node.text
