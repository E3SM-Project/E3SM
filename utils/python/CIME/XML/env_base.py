"""
Base class for env files .  This class inherits from EntryID.py
"""
from standard_module_setup import *

from entry_id import EntryID
from headers import Headers

class EnvBase(EntryID):
    def __init__(self, case_root, infile):
        fullpath = os.path.join(case_root, infile)
        EntryID.__init__(self, fullpath)
        if (not os.path.isfile(fullpath)):
            headerobj = Headers()
            headernode = headerobj.get_header_node(os.path.basename(fullpath))
            self.root.append(headernode)

    def test_reset(self):
        # get nodes who have a child node test_value
        testnodes = self.get_node("test_value/..")
        for node in nodes:
            tvnode = self.get_node("test_value",root=node)
            self.set_value(node, tvnode.text)



