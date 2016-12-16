"""
Interface to the env_run.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase

logger = logging.getLogger(__name__)

class EnvRun(EnvBase):

    def __init__(self, case_root=None, infile="env_run.xml"):
        """
        initialize an object interface to file env_run.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)
        self._components = []
        self._component_value_list = ["PIO_TYPENAME"]

    def comp_var_split(self, vid):
        parts = vid.split("_", 1)
        if len(parts) == 2 and parts[1] in self._component_value_list:
            return parts[1], parts[0], True
        return vid, None, False



