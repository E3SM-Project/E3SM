"""
Interface to the env_test.xml file.  This class inherits from EntryID.py
"""
from standard_module_setup import *

from entry_id import EntryID
from CIME.utils import expect, get_cime_root, get_model
from CIME.XML.headers import Headers

class EnvTest(EntryID):
    def __init__(self, infile=None):
        """
        initialize an object interface to file env_test.xml in the case directory
        """
        if(infile is None):
            infile = "env_test.xml"
        EntryID.__init__(self,infile)
        if(not os.path.isfile(infile)):
            headerobj = Headers()
            headernode = headerobj.get_header_node('env_test.xml')
            self.root.append(headernode)
