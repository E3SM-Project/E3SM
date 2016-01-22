"""
Interface to the config_machines.xml file.  This class inherits from GenericXML.py
"""
import xml.etree.ElementTree as ET
import logging
from GenericXML import GenericXML

class Machines(GenericXML):
    def __init__(self,infile):
        """ initialize an object """
        logging.info("Open file "+infile)
        GenericXML.__init__(self,infile)
    
    def list_available_machines(self):
        machines = []
        nodes  = self.GetNode('machine')
        for node in nodes:
            mach = node.attrib["MACH"]
            machines.append(mach)
        return machines
    








