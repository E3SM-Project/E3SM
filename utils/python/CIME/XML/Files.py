"""
Interface to the config_files.xml file.  This class inherits from EntryID.py
"""
import logging
import os
from EntryID import EntryID
from CIME.utils import expect, get_cime_root, get_model

class Files(EntryID):
    def __init__(self):
        """
        initialize an object

        >>> files = Files()
        >>> files.get_value('CASEFILE_HEADERS')
        '$CIMEROOT/cime_config/config_headers.xml'
        """
        infile = os.path.join(get_cime_root(),"cime_config",get_model(),"config_files.xml")
        EntryID.__init__(self,infile)
