"""
Interface to the env_mach_specific.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

class EnvMachSpecific(EnvBase):

    def __init__(self, caseroot=os.getcwd(), infile="env_mach_specific.xml"):
        """
        initialize an object interface to file env_mach_specific.xml in the case directory
        """
        EnvBase.__init__(self, caseroot, infile)
