"""
Interface to the env_build.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from env_base import EnvBase

logger = logging.getLogger(__name__)

class EnvBuild(EnvBase):
    def __init__(self, case_root=os.getcwd(), infile="env_build.xml"):
        """
        initialize an object interface to file env_build.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)
