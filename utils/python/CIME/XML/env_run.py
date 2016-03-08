"""
Interface to the env_run.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

class EnvRun(EnvBase):
    def __init__(self, case_root=os.getcwd(), infile="env_run.xml"):
        """
        initialize an object interface to file env_run.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)
