"""
Interface to the config_machines.xml file.  This class inherits from GenericXML.py
"""
import logging
import re
import socket
from GenericXML import GenericXML
from Files import Files
from CIME.utils import expect

class Machines(GenericXML):
    def __init__(self,infile=None):
        """ initialize an object """
        if(infile is None):
            files = Files()
            infile = files.get_resolved_value(
                files.get_value('MACHINES_SPEC_FILE'))
        logging.info("Open file "+infile)
        GenericXML.__init__(self,infile)
        self.machine = None
        self.name = None

    def get_node(self,nodename,attributes=None,root=None):
        if(self.machine is not None and root is None and nodename is not "machine"):
            node = GenericXML.get_node(self,nodename,attributes,root=self.machine)
        else:
            node = GenericXML.get_node(self,nodename,attributes,root)
        return node

    def list_available_machines(self):
        """
        Return a list of machines defined for a given CIME_MODEL
        """
        machines = []
        nodes  = self.get_node('machine')
        for node in nodes:
            mach = node.get("MACH")
            machines.append(mach)
        return machines

    def probe_machine_name(self):
        """
        Find a matching regular expression for hostname
        in the NODENAME_REGEX field in the file.   First match wins.
        """
        nametomatch = socket.gethostname().split(".")[0]
        nodes = self.get_node('machine')
        for node in nodes:
            machine = node.get('MACH')
            logging.info("machine is "+machine)
            self.set_machine(machine)
            regex_str_nodes =  self.get_node('NODENAME_REGEX',root=self.machine)
            if(len(regex_str_nodes)>0):
                regex_str = regex_str_nodes[0].text
                if (regex_str is not None):
                    logging.info("machine regex string is "+ regex_str)
                    regex = re.compile(regex_str)
                    if (regex.match(nametomatch)):
                        logging.info("Found machine: "+machine)
                        return machine

        return None


    def set_machine(self,machine):
        """
        Sets the machine block in the Machines object
        """
        if(self.machine is not None and self.name is not machine):
            self.machine = None
        mach_nodes = self.get_node('machine',{'MACH':machine})
        expect(mach_nodes, "No machine %s found" % machine)
        self.machine = mach_nodes[0]
        self.name = machine

    def get_value(self,name):
        """
        Get Value of fields in the config_machines.xml file
        """
        expect(self.machine is not None, "Machine object has no machine defined")
        value = None
        """
        COMPILER and MPILIB are special, if called without arguments they get the default value from the
        COMPILERS and MPILIBS lists in the file.
        """
        if(name == "COMPILER"):
            value = self.get_default_compiler()
        elif(name == "MPILIB"):
            value = self.get_default_MPIlib()
        else:
            nodes = self.get_node(name,root=self.machine)
            if(len(nodes)>0):
                node = nodes[0]
                expect(node is not None,"No match found for %s in machine %s" % (name,self.name))
                value = node.text
        if(value is None):
            """ if all else fails """
            value = GenericXML.get_value(self,name)
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
               "No list found for "+listname+" on machine "+self.name)
        supported_values = supported_values.split(',')
        if(reqval is None or reqval == "UNSET"):
            return supported_values[0]
        for val in supported_values:
            if(val == reqval):
                return reqval
        expect(False,"%s value %s not supported for machine %s" %
               (listname, reqval, self.name))

    def get_default_compiler(self):
        """
        Get the compiler to use from the list of COMPILERS
        """
        return self.get_field_from_list('COMPILERS')

    def get_default_MPIlib(self):
        """
        Get the MPILIB to use from the list of MPILIBS
        """
        return self.get_field_from_list('MPILIBS')

    def is_valid_compiler(self,compiler):
        """
        Check the compiler is valid for the current machine
        """
        if(self.get_field_from_list('COMPILERS',compiler) is not None):
            return True
        return False

    def is_valid_MPIlib(self, mpilib):
        """
        Check the MPILIB is valid for the current machine
        """
        if(self.get_field_from_list('MPILIBS',mpilib) is not None):
            return True
        return False

    def has_batch_system(self):
        """
        Return if this machine has a batch system
        """
        batch_system = self.get_node("batch_system")
        return not (batch_system is None or batch_system[0].get('type') == "none")
