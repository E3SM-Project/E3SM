"""
Interface to the config_machines.xml file.  This class inherits from GenericXML.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.utils import convert_to_unknown_type, get_cime_config

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
        schema = None
        if files is None:
            files = Files()
        if infile is None:
            infile = files.get_value("MACHINES_SPEC_FILE")
        schema = files.get_schema("MACHINES_SPEC_FILE")
        logger.debug("Verifying using schema {}".format(schema))

        self.machines_dir = os.path.dirname(infile)

        GenericXML.__init__(self, infile, schema)

        # Append the contents of $HOME/.cime/config_machines.xml if it exists
        # This could cause problems if node matchs are repeated when only one is expected
        local_infile = os.path.join(os.environ.get("HOME"),".cime","config_machines.xml")
        logger.debug("Infile: {}".format(local_infile))
        if os.path.exists(local_infile):
            GenericXML.read(self, local_infile, schema)

        if machine is None:
            if "CIME_MACHINE" in os.environ:
                machine = os.environ["CIME_MACHINE"]
            else:
                cime_config = get_cime_config()
                if cime_config.has_option("main", "machine"):
                    machine = cime_config.get("main", "machine")
                if machine is None:
                    machine = self.probe_machine_name()

        expect(machine is not None, "Could not initialize machine object from {} or {}".format(infile, local_infile))
        self.set_machine(machine)

    def get_child(self, name=None, attributes=None, root=None, err_msg=None):
        if root is None:
            root = self.machine_node
        return super(Machines, self).get_child(name, attributes, root, err_msg)

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
        nodes = self.get_children(root=self.machine_node)
        node_names = []
        for node in nodes:
            node_names.append(self.name(node))
        return node_names

    def get_first_child_nodes(self, nodename):
        """
        Return the names of all the child nodes for the target machine
        """
        nodes = self.get_children(nodename, root=self.machine_node)
        return nodes

    def list_available_machines(self):
        """
        Return a list of machines defined for a given CIME_MODEL
        """
        machines = []
        nodes  = self.get_children("machine")
        for node in nodes:
            mach = self.get(node, "MACH")
            machines.append(mach)
        return machines

    def probe_machine_name(self, warn=True):
        """
        Find a matching regular expression for hostname
        in the NODENAME_REGEX field in the file.   First match wins.
        """

        names_not_found = []

        nametomatch = socket.getfqdn()
        machine = self._probe_machine_name_one_guess(nametomatch)

        if machine is None:
            names_not_found.append(nametomatch)

            nametomatch = socket.gethostname()
            machine = self._probe_machine_name_one_guess(nametomatch)

            if machine is None:
                names_not_found.append(nametomatch)

                names_not_found_quoted = ["'" + name + "'" for name in names_not_found]
                names_not_found_str = ' or '.join(names_not_found_quoted)
                if warn:
                    logger.warning("Could not find machine match for {}".format(names_not_found_str))

        return machine

    def _probe_machine_name_one_guess(self, nametomatch):
        """
        Find a matching regular expression for nametomatch in the NODENAME_REGEX
        field in the file. First match wins. Returns None if no match is found.
        """

        machine = None
        nodes = self.get_children("machine")

        for node in nodes:
            machtocheck = self.get(node, "MACH")
            logger.debug("machine is " + machtocheck)
            regex_str_node = self.get_optional_child("NODENAME_REGEX", root=node)
            regex_str = machtocheck if regex_str_node is None else self.text(regex_str_node)

            if regex_str is not None:
                logger.debug("machine regex string is " + regex_str)
                regex = re.compile(regex_str)
                if regex.match(nametomatch):
                    logger.debug("Found machine: {} matches {}".format(machtocheck, nametomatch))
                    machine = machtocheck
                    break

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
        if machine == "Query":
            self.machine = machine
        elif self.machine != machine or self.machine_node is None:
            self.machine_node = super(Machines,self).get_child("machine", {"MACH" : machine}, err_msg="No machine {} found".format(machine))
            self.machine = machine

        return machine
    #pylint: disable=arguments-differ
    def get_value(self, name, attributes=None, resolved=True, subgroup=None):
        """
        Get Value of fields in the config_machines.xml file
        """
        expect(self.machine_node is not None, "Machine object has no machine defined")
        expect(subgroup is None, "This class does not support subgroups")
        value = None

        # COMPILER and MPILIB are special, if called without arguments they get the default value from the
        # COMPILERS and MPILIBS lists in the file.
        if name == "COMPILER":
            value = self.get_default_compiler()
        elif name == "MPILIB":
            value = self.get_default_MPIlib(attributes)
        else:
            node = self.get_optional_child(name, root=self.machine_node, attributes=attributes)
            if node is not None:
                value = self.text(node)
        if resolved:
            if value is not None:
                value = self.get_resolved_value(value)
            elif name in os.environ:
                value = os.environ[name]

            value = convert_to_unknown_type(value)

        return value

    def get_field_from_list(self, listname, reqval=None, attributes=None):
        """
        Some of the fields have lists of valid values in the xml, parse these
        lists and return the first value if reqval is not provided and reqval
        if it is a valid setting for the machine
        """
        expect(self.machine_node is not None, "Machine object has no machine defined")
        supported_values = self.get_value(listname, attributes=attributes)
        # if no match with attributes, try without
        if supported_values is None:
            supported_values = self.get_value(listname, attributes=None)

        expect(supported_values is not None,
               "No list found for " + listname + " on machine " + self.machine)
        supported_values = supported_values.split(",") #pylint: disable=no-member

        if reqval is None or reqval == "UNSET":
            return supported_values[0]

        for val in supported_values:
            if val == reqval:
                return reqval
        return None

    def get_default_compiler(self):
        """
        Get the compiler to use from the list of COMPILERS
        """
        cime_config = get_cime_config()
        if cime_config.has_option('main','COMPILER'):
            value = cime_config.get('main', 'COMPILER')
            expect(self.is_valid_compiler(value), "User-selected compiler {} is not supported on machine {}".format(value, self.machine))
        else:
            value = self.get_field_from_list("COMPILERS")
        return value

    def get_default_MPIlib(self, attributes=None):
        """
        Get the MPILIB to use from the list of MPILIBS
        """
        return self.get_field_from_list("MPILIBS", attributes=attributes)

    def is_valid_compiler(self,compiler):
        """
        Check the compiler is valid for the current machine

        >>> machobj = Machines(machine="edison")
        >>> machobj.get_default_compiler()
        'intel'
        >>> machobj.is_valid_compiler("gnu")
        True
        >>> machobj.is_valid_compiler("nag")
        False
        """
        return self.get_field_from_list("COMPILERS", reqval=compiler) is not None

    def is_valid_MPIlib(self, mpilib, attributes=None):
        """
        Check the MPILIB is valid for the current machine

        >>> machobj = Machines(machine="edison")
        >>> machobj.is_valid_MPIlib("mpi-serial")
        True
        >>> machobj.is_valid_MPIlib("fake-mpi")
        False
        """
        return mpilib == "mpi-serial" or \
            self.get_field_from_list("MPILIBS", reqval=mpilib, attributes=attributes) is not None

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
        batch_system = self.get_optional_child("BATCH_SYSTEM", root=self.machine_node)
        if batch_system is not None:
            result = (self.text(batch_system) is not None and self.text(batch_system) != "none")
        logger.debug("Machine {} has batch: {}".format(self.machine, result))
        return result

    def get_suffix(self, suffix_type):
        node = self.get_optional_child("default_run_suffix")
        if node is not None:
            suffix_node = self.get_optional_child(suffix_type, root=node)
            if suffix_node is not None:
                return self.text(suffix_node)

        return None

    def set_value(self, vid, value, subgroup=None, ignore_type=True):
        tmproot = self.root
        self.root = self.machine_node
        #pylint: disable=assignment-from-no-return
        result = super(Machines, self).set_value(vid, value, subgroup=subgroup,
                                               ignore_type=ignore_type)
        self.root = tmproot
        return result


    def print_values(self):
        # write out machines
        machines = self.get_children("machine")
        logger.info("Machines")
        for machine in machines:
            name = self.get(machine, "MACH")
            desc = self.get_child("DESC", root=machine)
            os_  = self.get_child("OS", root=machine)
            compilers = self.get_child("COMPILERS", root=machine)
            max_tasks_per_node = self.get_child("MAX_TASKS_PER_NODE", root=machine)
            max_mpitasks_per_node = self.get_child("MAX_MPITASKS_PER_NODE", root=machine)

            print( "  {} : {} ".format(name , self.text(desc)))
            print( "      os             ", self.text(os_))
            print( "      compilers      ",self.text(compilers))
            if max_mpitasks_per_node is not None:
                print("      pes/node       ",self.text(max_mpitasks_per_node))
            if max_tasks_per_node is not None:
                print("      max_tasks/node ",self.text(max_tasks_per_node))

    def return_values(self):
        """ return a dictionary of machine info
        This routine is used by external tools in https://github.com/NCAR/CESM_xml2html
        """
        machines = self.get_children("machine")
        mach_dict = dict()
        logger.debug("Machines return values")
        for machine in machines:
            name = self.get(machine, "MACH")
            desc = self.get_child("DESC", root=machine)
            mach_dict[(name,"description")] = self.text(desc)
            os_  = self.get_child("OS", root=machine)
            mach_dict[(name,"os")] = self.text(os_)
            compilers = self.get_child("COMPILERS", root=machine)
            mach_dict[(name,"compilers")] = self.text(compilers)
            max_tasks_per_node = self.get_child("MAX_TASKS_PER_NODE", root=machine)
            mach_dict[(name,"max_tasks_per_node")] = self.text(max_tasks_per_node)
            max_mpitasks_per_node = self.get_child("MAX_MPITASKS_PER_NODE", root=machine)
            mach_dict[(name,"max_mpitasks_per_node")] = self.text(max_mpitasks_per_node)

        return mach_dict
