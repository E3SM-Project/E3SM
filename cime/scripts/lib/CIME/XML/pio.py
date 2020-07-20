"""
Class for config_pio files .  This class inherits from EntryID.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.entry_id import EntryID
from CIME.XML.files import Files

from collections import OrderedDict

logger = logging.getLogger(__name__)

class PIO(EntryID):

    def __init__(self, comp_classes, infile=None, files=None):
        if infile is None:
            if files is None:
                files = Files()
            infile = files.get_value("PIO_SPEC_FILE")

        EntryID.__init__(self, infile)

        self._components = list(comp_classes)

    def check_if_comp_var(self, vid, attribute=None, node=None):
        comp = None
        new_vid = None
        for comp in self._components:
            if vid.endswith('_'+comp):
                new_vid = vid.replace('_'+comp, '', 1)
            elif vid.startswith(comp+'_'):
                new_vid = vid.replace(comp+'_', '', 1)
            elif '_' + comp + '_' in vid:
                new_vid = vid.replace(comp+'_','', 1)

            if new_vid is not None:
                return new_vid, comp, True

        return vid, None, False

    def get_defaults(self, grid=None, compset=None, mach=None, compiler=None, mpilib=None): # pylint: disable=unused-argument
        # should we have a env_pio file
        defaults = OrderedDict()
        save_for_last = []

        # Load args into attribute dict
        attributes = {}
        for attrib in ["grid", "compset", "mach", "compiler", "mpilib"]:
            if locals()[attrib] is not None:
                attributes[attrib] = locals()[attrib]

        # Find defauts
        for node in self.get_children("entry"):
            value = self.get_default_value(node, attributes)
            if value:
                myid = self.get(node, "id")
                iscompvar = self.check_if_comp_var(myid)[-1]
                if iscompvar:
                    save_for_last.append( (myid, value) )
                else:
                    defaults[myid] = value

        # comp-specific vars must come last so they take precedence over general settings
        for k, v in save_for_last:
            defaults[k] = v

        return defaults
