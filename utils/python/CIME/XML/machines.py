"""
Interface to the config_machines.xml file.  This class inherits from GenericXML.py
"""
from standard_module_setup import *
import socket
from generic_xml import GenericXML
from files import Files
from CIME.utils import expect

class Machines(GenericXML):

    def __init__(self, infile=None, files=None, machine=None):
        """
        initialize an object
        if a filename is provided it will be used,
        otherwise if a files object is provided it will be used
        otherwise create a files object from default values
        """
        self.machine = None
        self.name    = None

        if (infile is None):
            if (files is None):
                files = Files()
            infile = files.get_value("MACHINES_SPEC_FILE")

        logging.info("Open file " + infile)
        GenericXML.__init__(self, infile)

        if (machine is None):
            machine = self.probe_machine_name()

        self.set_machine(machine)

    def get_machine_name(self):
        """
        Return the name of the machine
        """
        return self.name

    def get_node(self, nodename, attributes=None, root=None):
        if (self.machine is not None and root is None and nodename is not "machine"):
            node = GenericXML.get_node(self, nodename, attributes, root=self.machine)
        else:
            node = GenericXML.get_node(self, nodename, attributes, root)
        return node

    def list_available_machines(self):
        """
        Return a list of machines defined for a given CIME_MODEL
        """
        machines = []
        nodes  = self.get_node("machine")
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
        nodes = self.get_node("machine")

        for node in nodes:
            machtocheck = node.get("MACH")
            logging.debug("machine is " + machtocheck)
            self.set_machine(machtocheck)
            regex_str_nodes = self.get_node("NODENAME_REGEX", root=self.machine)

            if (regex_str_nodes):
                regex_str = regex_str_nodes[0].text
            else:
                regex_str = machtocheck

            if (regex_str is not None):
                logging.debug("machine regex string is " + regex_str)
                regex = re.compile(regex_str)
                if (regex.match(nametomatch)):
                    logging.info("Found machine: %s matches %s" % (machtocheck, nametomatch))
                    machine = machtocheck
                    break

        if (machine is None):
            logging.warning("Could not probe machine for hostname '%s'" % nametomatch)

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
        if (self.machine is not None and self.name is not machine):
            self.machine = None
        mach_nodes = self.get_node("machine", {"MACH":machine})
        expect(mach_nodes, "No machine %s found" % machine)
        self.machine = mach_nodes[0]
        self.name = machine
        return machine

    def get_value(self, name, resolved=True):
        """
        Get Value of fields in the config_machines.xml file
        """
        expect(self.machine is not None, "Machine object has no machine defined")
        value = None

        # COMPILER and MPILIB are special, if called without arguments they get the default value from the
        # COMPILERS and MPILIBS lists in the file.
        if (name == "COMPILER"):
            value = self.get_default_compiler()
        elif (name == "MPILIB"):
            value = self.get_default_MPIlib()
        else:
            nodes = self.get_node(name, root=self.machine)
            if (nodes):
                node = nodes[0]
                expect(node is not None, "No match found for %s in machine %s" % (name, self.name))
                value = node.text

        if (value is None):
            # if all else fails
            value = GenericXML.get_value(self, name)

        if (resolved):
            value = self.get_resolved_value(value)

        return value

    def get_field_from_list(self, listname, reqval=None):
        """
        Some of the fields have lists of valid values in the xml, parse these
        lists and return the first value if reqval is not provided and reqval
        if it is a valid setting for the machine
        """
        expect(self.machine is not None, "Machine object has no machine defined")
        supported_values = self.get_value(listname)
        expect(supported_values is not None,
               "No list found for " + listname + " on machine " + self.name)
        supported_values = supported_values.split(",")

        if (reqval is None or reqval == "UNSET"):
            return supported_values[0]
        for val in supported_values:
            if (val == reqval):
                return reqval

        expect(False, "%s value %s not supported for machine %s" %
               (listname, reqval, self.name))

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
        if (self.get_field_from_list("COMPILERS", compiler) is not None):
            return True
        return False

    def is_valid_MPIlib(self, mpilib):
        """
        Check the MPILIB is valid for the current machine

        >>> machobj = Machines(machine="edison")
        >>> machobj.is_valid_MPIlib("mpi-serial")
        True
        """
        if (self.get_field_from_list("MPILIBS", mpilib) is not None):
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
        batch_system = self.get_node("batch_system")
        if (batch_system):
            result = (batch_system[0].get("type") != "none")
        logging.debug("Machine %s has batch: %s" % (self.name, result))
        return result

    def get_batch_system_type(self):
        """
        Return the batch system using on this machine

        >>> machobj = Machines(machine="edison")
        >>> machobj.get_batch_system_type()
        'slurm'
        """
        batch_system = self.get_node("batch_system")
        return batch_system[0].get("type")
