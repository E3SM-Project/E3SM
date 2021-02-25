"""
Interface to the env_case.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.env_base import EnvBase

logger = logging.getLogger(__name__)

class EnvCase(EnvBase):
    # pylint: disable=unused-argument
    def __init__(self, case_root=None, infile="env_case.xml", components=None, read_only=False):
        """
        initialize an object interface to file env_case.xml in the case directory
        """
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_entry_id.xsd")
        EnvBase.__init__(self, case_root, infile, schema=schema, read_only=read_only)
