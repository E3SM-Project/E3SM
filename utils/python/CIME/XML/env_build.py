"""
Interface to the env_build.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase

logger = logging.getLogger(__name__)

class EnvBuild(EnvBase):
    # pylint: disable=unused-argument
    def __init__(self, case_root=None, infile="env_build.xml",components=None):
        """
        initialize an object interface to file env_build.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)






