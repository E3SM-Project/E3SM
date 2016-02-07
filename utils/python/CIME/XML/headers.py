"""
Interface to the config_headers.xml file.  This class inherits from EntryID.py
"""
from standard_module_setup import *

from entry_id import EntryID
from files import Files
from CIME.utils import expect, get_cime_root, get_model

class Headers(EntryID):
    def __init__(self,infile=None):
        """
        initialize an object

        >>> files = Files()
        >>> files.get_value('CASEFILE_HEADERS',resolved=False)
        '$CIMEROOT/cime_config/config_headers.xml'
        """
        if(infile is None):
            files = Files()
            infile = files.get_value('CASEFILE_HEADERS',resolved=True)
        EntryID.__init__(self,infile)

    def get_header_node(self,fname):
        fnode=self.get_node("file", attributes={"name":fname})
        expect(len(fnode)==1,"Invalid match count in headers for filename '%s'" % fname)
        headernode = self.get_node('header',root=fnode[0])
        expect(len(headernode)==1,"Invalid match count in headers for filename '%s'" % fname)
        return headernode[0]
