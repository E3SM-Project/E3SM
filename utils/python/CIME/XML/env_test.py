"""
Interface to the env_test.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

class EnvTest(EnvBase):
    def __init__(self, case_root=os.getcwd(), infile="env_test.xml"):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)
