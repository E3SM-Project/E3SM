"""
Base class for env files .  This class inherits from EntryID.py
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.entry_id import EntryID
from CIME.XML.headers import Headers

logger = logging.getLogger(__name__)

class EnvBase(EntryID):

    def __init__(self, case_root, infile):
        if case_root is None:
            case_root = os.getcwd()

        if os.path.isabs(infile):
            fullpath = infile
        else:
            fullpath = os.path.join(case_root, infile)

        EntryID.__init__(self, fullpath)
        if not os.path.isfile(fullpath):
            headerobj = Headers()
            headernode = headerobj.get_header_node(os.path.basename(fullpath))
            self.root.append(headernode)
