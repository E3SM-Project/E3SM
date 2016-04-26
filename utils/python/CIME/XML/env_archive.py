"""
Interface to the env_archive.xml file.  This class inherits from EnvBase
"""
from standard_module_setup import *

from CIME.XML.generic_xml import GenericXML
from CIME.XML.archive import Archive
from CIME.XML.headers import Headers

logger = logging.getLogger(__name__)

class EnvArchive(GenericXML):

    def __init__(self, case_root=None, infile="env_archive.xml"):
        """
        initialize an object interface to file env_archive.xml in the case directory
        """
        if case_root is None:
            case_root = os.getcwd()
        if os.path.isabs(infile):
            fullpath = infile
        else:
            fullpath = os.path.join(case_root, infile)
        GenericXML.__init__(self, fullpath)
        if not os.path.isfile(infile):
            headerobj = Headers()
            headernode = headerobj.get_header_node(os.path.basename(fullpath))
            self.root.append(headernode)
            archive = Archive()
            self.root.append(archive.root)

