"""
Interface to the env_mach_specific.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from generic_xml import GenericXML

class EnvMachSpecific(GenericXML):

    def __init__(self, caseroot, infile="env_mach_specific.xml"):
        """
        initialize an object interface to file env_mach_specific.xml in the case directory
        """
	if(os.path.isabs(infile)):
	    fullpath = infile
        else:
            fullpath = os.path.join(caseroot,infile)
        GenericXML.__init__(self, fullpath)
