"""
Interface to the archive.xml file.  This class inherits from GenericXML.py

"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.utils import expect, get_cime_root, get_model

logger = logging.getLogger(__name__)

class Archive(GenericXML):

    def __init__(self, infile=None):
        """
        initialize an object
        """
        if infile is None:
            infile = os.path.join(get_cime_root(), "cime_config", get_model(), "config_archive.xml")

        GenericXML.__init__(self, infile)

