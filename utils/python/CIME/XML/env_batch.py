"""
Interface to the env_batch.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

class EnvBatch(EnvBase):
    def __init__(self, case_root=os.getcwd(), infile="env_batch.xml"):
        """
        initialize an object interface to file env_batch.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)
