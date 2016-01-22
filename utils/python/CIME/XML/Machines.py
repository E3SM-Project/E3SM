"""
Interface to the config_machines.xml file.  This class inherits from GenericXML.py
"""
import xml.etree.ElementTree as ET
import logging
from GenericXML import GenericXML
from CIME.utils import expect

class Machines(GenericXML):
    def __init__(self,infile):
        """ initialize an object """
        logging.info("Open file "+infile)
        self.machine = None
        GenericXML.__init__(self,infile)
    
    def list_available_machines(self):
        machines = []
        nodes  = self.GetNode('machine')
        for node in nodes:
            mach = node.attrib["MACH"]
            machines.append(mach)
        return machines
    
    def set_machine(self,machine):
        self.machine = self.GetNode('machine',{'MACH':machine})
        self.name = machine

    def GetValue(self,name):
        expect(self.machine is not None, "Machine object has no machine defined")
        node = self.GetNode(name,root=self.machine)
        expect(len(node)==0,"Expecting exactly one match for %s, got %s" % (name,len(node)))
        return node.text


    def getMPIlib(self, mpilib=None):
        expect(self.machine is not None, "Machine object has no machine defined")
        supported_mpilibs = self.GetValue("MPILIBS")
        mpilibs = supported_mpilibs.split(',')
        if(mpilib is None or mpilib == "UNSET"):
            return mpilibs[0]
        for lib in mpilibs:
            print lib
            if(lib == mpilib):
                return mpilib
        logging.critical(mpilib + " not defined for machine " +  self.name)




