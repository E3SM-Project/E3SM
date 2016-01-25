"""
Interface to the config_machines.xml file.  This class inherits from GenericXML.py
"""
import logging
from GenericXML import GenericXML
from CIME.utils import expect

class Machines(GenericXML):
    def __init__(self,infile):
        """ initialize an object """
        logging.warn("Open file "+infile)
        self.machine = None
        GenericXML.__init__(self,infile)
    
    def list_available_machines(self):
        machines = []
        nodes  = self.GetNode('machine')
        print "nodes  %d" % len(nodes)

        for node in nodes:
            mach = node.attrib["MACH"]
            machines.append(mach)
        return machines
    
    def set_machine(self,machine):
        self.machine = self.GetNode('machine',{'MACH':machine})[0]
        expect(self.machine is not None,"No machine %s found" % machine)
        
        self.name = machine

    def GetValue(self,name):
        expect(self.machine is not None, "Machine object has no machine defined")
        node = self.GetNode(name,root=self.machine)[0]
        expect(node is not None,"No match found for %s in machine %s" % (name,self.name))
        return node.text

    def getfieldfromlist(self, listname, reqval=None):
        expect(self.machine is not None, "Machine object has no machine defined")
        supported_values = self.GetValue(listname).split(',')

        if(reqval is None or reqval == "UNSET"):
            return supported_values[0]
        for val in supported_values:
            if(val == reqval):
                return reqval
        logging.critical("%s value %s not supported for machine %s" %
                         (listname, reqval, self.name))

    def getCompiler(self, compiler=None):
        return self.getfieldfromlist('Compilers',compiler)

    def getMPIlib(self, mpilib=None):
        return self.getfieldfromlist('MPILIBS',mpilib)




