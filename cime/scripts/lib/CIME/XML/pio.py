"""
Class for config_pio files .  This class inherits from EntryID.py
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.entry_id import EntryID
from CIME.XML.files import Files
logger = logging.getLogger(__name__)

class PIO(EntryID):

    def __init__(self, infile=None, files=None):
        if infile is None:
            if files is None:
                files = Files()
            infile = files.get_value("PIO_SPEC_FILE")

        EntryID.__init__(self, infile)

    def get_defaults(self, grid=None, compset=None, mach=None, compiler=None, mpilib=None): # pylint: disable=unused-argument
        # should we have a env_pio file
        defaults = {}

        # Load args into attribute dict
        attributes = {}
        for attrib in ["grid", "compset", "mach", "compiler", "mpilib"]:
            if locals()[attrib] is not None:
                attributes[attrib] = locals()[attrib]

        # Find defauts
        for node in self.get_children("entry"):
            value = self.get_default_value(node, attributes)
            if value:
                defaults[self.get(node, "id")] = value

        return defaults
